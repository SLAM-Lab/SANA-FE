import csv
import os
import sys
import numpy as np
import pandas as pd
import dill
# Plotting
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
import matplotlib.patches as patches
import math

import torch
from torchvision import datasets, transforms
from torch.utils.data import DataLoader
import tonic

# Try importing the installed sanafe library. If not installed, require a
#  fall-back to a local build of the sanafe Python library
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
    RUN_PATH =  os.path.abspath((os.path.join(PROJECT_DIR, "runs", "crossbar")))

    print(f"Project dir: {PROJECT_DIR}")
    sys.path.insert(0, PROJECT_DIR)
    sys.path.insert(0, os.path.join(PROJECT_DIR, "etc")) # For helper script
    import sanafe


def ternary_weight(weight):
    threshold = 0.7 * torch.mean(torch.abs(weight))
    output = torch.zeros_like(weight)
    output[weight > threshold] = 1.0
    output[weight < -threshold] = -1.0
    return output.numpy()


# Okabe-Ito color palette (colorblind-friendly)
okabe_ito_colors = [
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermillion
    '#CC79A7',  # Reddish Purple
    '#000000'   # Black
]


def run_spiking_digits(num_inputs, analog_synapses=True):
    print(f"Loading models for binarized spiking digits")
    # Data loading
    transform = tonic.transforms.Compose([
        tonic.transforms.Downsample(time_factor=0.05, spatial_factor=0.1),  # For Jason's circuit-aware trained SNN
        tonic.transforms.ToFrame(sensor_size=(70, 1, 1), time_window=1000)
    ])
    testset = tonic.datasets.SHD(
        save_to="./runs/crossbar/data",
        transform=transform, train=False
    )
    dataloader = torch.utils.data.DataLoader(testset, batch_size=1)

    inputs = []
    labels = []
    count = 0
    for input_frame, label in dataloader:
        input_frame = input_frame.cpu()
        label = label.cpu()
        inputs.append(input_frame.numpy())
        labels.append(int(label))
        count += 1
        if count >= num_inputs:
            break

    if analog_synapses:
        snn_path = "/home/usr1/jboyle/neuro/sana-fe/runs/crossbar/app_models/shd_70_256R_20_bin_crossbar_aware.pt"
        snn = torch.load(snn_path, pickle_module=dill, map_location=torch.device("cpu")).state_dict()
    else:
        # Load a normally trained rate-coded SHD model to deploy on Loihi's
        #  (digital, full-precision) synapses
        snn_path = "/home/usr1/jboyle/neuro/sana-fe/runs/crossbar/app_models/shd_70_256R_20.pt"
        snn = torch.load(snn_path, pickle_module=dill, map_location=torch.device("cpu")).state_dict()

    weights = {}
    if analog_synapses:
        weights["fc1"] = ternary_weight(snn["fc1.weight"])
        weights["fc2"] = ternary_weight(snn["fc2.weight"])
        weights["fcr"] = ternary_weight(snn["fcr.weight"])

        crossbar_biases = {}
        crossbar_biases[1] = ternary_weight(snn["fc1.bias"])
        crossbar_biases[2] = ternary_weight(snn["fc2.bias"])
        crossbar_biases[3] = ternary_weight(snn["fcr.bias"])
    else:
        weights["fc1"] = snn["fc1.weight"]
        weights["fc2"] = snn["fc2.weight"]
        weights["fcr"] = snn["fcr.weight"]

    decays = {}
    decays[1] = float(snn["lif1.beta"])
    decays[2] = float(snn["lif2.beta"])

    thresholds = {}
    thresholds[1] = float(snn["lif1.threshold"])
    thresholds[2] = float(snn["lif2.threshold"])

    # PyTorch stores weights in an array with dims (num out x num in)
    in_neurons = weights["fc1"].shape[1]
    hidden_neurons = weights["fc1"].shape[0]
    out_neurons = weights["fc2"].shape[0]

    print(f"in:{in_neurons}, hidden 1:{hidden_neurons}, our:{out_neurons}")
    #"""
    # Load the LASANA architecture with analog neurons
    if analog_synapses:
        platform = "crossbar"
        arch = sanafe.load_arch("/home/usr1/jboyle/neuro/lasana/crossbar/lasana_crossbar.yaml")
    else:
        platform = "loihi"
        arch = sanafe.load_arch(os.path.join(PROJECT_DIR, "arch", "loihi.yaml"))

    spike_filename = f"spikes_{platform}_shd.csv"
    perf_filename = f"perf_{platform}_shd.csv"
    potential_filename = f"potential_{platform}_shd.csv"


    print("Creating network in SANA-FE")
    network = sanafe.Network()
    in_layer = network.create_neuron_group("in", in_neurons)
    hidden_layer_1 = network.create_neuron_group("hidden", hidden_neurons)
    out_layer = network.create_neuron_group("out", out_neurons)

    # Defaults for Loihi baseline
    synapse_hw = "loihi_dense_synapse"
    dendrite_hw = "loihi_dendrites"

    print("Creating input layer")
    for id, neuron in enumerate(in_layer):
        neuron.set_attributes(dendrite_hw_name=dendrite_hw,
                              log_spikes=True,
                              log_potential=True,
                              soma_hw_name=f"loihi_inputs[{id}]")
        neuron.map_to_core(arch.tiles[0].cores[0])

    print("Creating hidden layer")
    soma_hw = "loihi_lif"

    num_crossbar_rows_per_core = 256
    input_cores = 1
    total_rows_per_neuron = int(math.ceil(in_neurons / 32) + math.ceil(hidden_neurons / 32))
    hidden_neurons_per_core = num_crossbar_rows_per_core // total_rows_per_neuron
    hidden_cores = int(math.ceil(hidden_neurons / hidden_neurons_per_core))
    output_cores = 1

    for id, neuron in enumerate(hidden_layer_1):
        hidden_parameters = {
            "threshold": thresholds[1],
            "leak_decay": decays[1],
            "reset_mode": "hard"
        }
        if analog_synapses:
            synapse_hw = f"synapse_crossbar[{id}]"
            dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=soma_hw,
                            dendrite_hw_name=dendrite_hw,
                            default_synapse_hw_name=synapse_hw,
                            log_spikes=True,
                            log_potential=True,
                            model_attributes=hidden_parameters)
        # We assume that there are 256 crossbar rows available per core
        core_id = input_cores + (id // hidden_neurons_per_core)
        print(f"id:{id} hidden_per_core:{hidden_neurons_per_core} core_id:{core_id}")
        tile_id = core_id // 4
        core_offset = core_id % 4
        assert(tile_id < 4)
        assert(core_offset < 4)
        neuron.map_to_core(arch.tiles[tile_id].cores[core_offset])

    for id, neuron in enumerate(out_layer):
        hidden_parameters = {
            "threshold": thresholds[2],
            "leak_decay": decays[2],
            "reset_mode": "hard"
        }
        if analog_synapses:
            synapse_hw = f"synapse_crossbar[{id}]"
            dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=soma_hw,
                              dendrite_hw_name=dendrite_hw,
                              default_synapse_hw_name=synapse_hw,
                              log_spikes=True,
                              model_attributes=hidden_parameters)
        core_id = input_cores + hidden_cores
        tile_id = core_id // 4
        core_offset = core_id % 4
        assert(tile_id < 4)
        assert(core_offset < 4)
        neuron.map_to_core(arch.tiles[tile_id].cores[core_offset])

    # Connect neurons in both layers
    print("Connecting neurons")
    print("Adding first layer connections")

    total_hidden_1_crossbar_rows = (in_neurons + 31) // 32
    total_hidden_r_crossbar_rows = (hidden_neurons + 31) // 32
    total_out_crossbar_rows = (out_neurons + 31) // 32

    print("Total crossbar rows:")
    print(f"In->Hidden: {total_hidden_1_crossbar_rows}")
    print(f"Hidden->Hidden: {total_hidden_r_crossbar_rows}")
    print(f"Hidden->Out: {total_out_crossbar_rows}")

    for src in range(in_neurons):
        for dst in range(hidden_neurons):
            connection_parameters = {}
            if analog_synapses:
                row = src // 32
                col = src % 32
                connection_parameters["crossbar_position"] = [row, col]
                connection_parameters["weight"] = int(weights["fc1"][dst, src])
                connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
                connection_parameters["crossbar_bias"] = int(crossbar_biases[1][dst])
            else:
                connection_parameters["weight"] = float(weights["fc1"][dst, src])

            network["in"][src].connect_to_neuron(
                network["hidden"][dst], connection_parameters)

    for src in range(hidden_neurons):
        for dst in range(hidden_neurons):
            connection_parameters = {}
            if analog_synapses:
                # Offset from the previously mapped rows for the hidden layer
                row = total_hidden_1_crossbar_rows + (src // 32)
                col = src % 32
                connection_parameters["crossbar_position"] = [row, col]
                connection_parameters["weight"] = int(weights["fcr"][dst, src])
                connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
                connection_parameters["crossbar_bias"] = int(crossbar_biases[3][dst])
            else:
                connection_parameters["weight"] = float(weights["fcr"][dst, src])

            network["hidden"][src].connect_to_neuron(
                network["hidden"][dst], connection_parameters)

    for src in range(hidden_neurons):
        for dst in range(out_neurons):
            connection_parameters = {}
            if analog_synapses:
                row = src // 32
                col = src % 32
                connection_parameters["crossbar_position"] = [row, col]
                connection_parameters["weight"] = int(weights["fc2"][dst, src])
                connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
                connection_parameters["crossbar_bias"] = int(crossbar_biases[2][dst])
            else:
                connection_parameters["weight"] = float(weights["fc2"][dst, src])

            network["hidden"][src].connect_to_neuron(
                network["out"][dst], connection_parameters)

    network.save(os.path.join(RUN_PATH, f"{platform}_shd.yaml"))

    # Run a simulation
    print("Building h/w")
    hw = sanafe.SpikingChip(arch)
    print("Loading SNN")
    hw.load(network)
    print(f"Running simulation")

    timesteps_per_input = []
    mapped_inputs = hw.mapped_neuron_groups["in"]
    print("Setting inputs")
    for input_idx in range(num_inputs):
        print(f"Simulating input: {input_idx}")
        spiking_digit_input = inputs[input_idx].squeeze()
        for id, mapped_neuron in enumerate(mapped_inputs):
            spiketrain = list(spiking_digit_input[:, id])
            mapped_neuron.set_model_attributes(
                model_attributes={"spikes": spiketrain})
        timesteps = len(spiketrain)

        print(f"Simulating for {timesteps} timesteps")
        is_first_input = (input_idx == 0)

        hw.sim(timesteps,
               spike_trace=os.path.join(RUN_PATH, spike_filename),
               perf_trace=os.path.join(RUN_PATH, perf_filename),
               potential_trace=os.path.join(RUN_PATH, potential_filename),
               write_trace_headers=is_first_input,
               processing_threads=4,
               scheduler_threads=8)
        timesteps_per_input.append(timesteps)
        hw.reset()
        #input()

    # The inputs (and therefore timesteps per input) will be the same across
    #  Loihi/analog neuron runs. Only store this once.
    np.savetxt(os.path.join(PROJECT_DIR, "runs", "crossbar", f"crossbar_shd.csv"),
               np.array(timesteps_per_input), fmt="%d")
    print("Simulation finished")


def load_dataset(analog_neurons=True):
    weights = {}

    # Load the Spiking Heidelberg Digits (SHD) network
    if analog_neurons: # Use Indiveri neural circuits
        # Use a SNN trained using circuit-aware methods i.e. was
        #  trained specifically on an IMAC-sim crossbar model.
        spiking_digits_model = torch.load(
                "/home/usr1/jboyle/neuro/sana-fe/runs/crossbar/app_models/shd_70_256R_20_bin_crossbar_aware.pt",
                map_location=torch.device("cpu")).state_dict()

        for attribute_name, param in spiking_digits_model.items():
            weights[attribute_name] = param.detach().numpy()
    else: # Use digital LIF neuron models
        # Use an SNN trained on a normal (linear) synapses
        spiking_digits_model = torch.load("/home/usr1/jboyle/neuro/sana-fe/runs/crossbar/app_models/shd_70_256R_20.pt")

        #for attribute_name, param in spiking_digits.named_parameters():
        #    weights[attribute_name] = param.detach().numpy()

    # Load the spiking digits test inputs from the raw dataset, applying the
    #  same transformations as the training scripts
    frame_transform = tonic.transforms.Compose([
        tonic.transforms.Downsample(time_factor=0.05, spatial_factor=0.1),  # For Jason's circuit-aware trained SNN
        tonic.transforms.ToFrame(sensor_size=(70, 1, 1), time_window=1000)
    ])
    testset = tonic.datasets.SHD(
        save_to="./runs/lasana/data",
        transform=frame_transform, train=False
    )
    dataloader = torch.utils.data.DataLoader(testset, batch_size=1)

    inputs = []
    labels = []
    count = 0
    for input_frame, label in dataloader:
        input_frame = input_frame.cpu()
        label = label.cpu()
        inputs.append(input_frame.numpy())
        labels.append(int(label))
        count += 1
        if count >= num_inputs:
            break

    return inputs, labels, weights



def calculate_accuracy(analog_synapses=True):
    print(f"Calculating accuracy architecture")
    timesteps_per_input = list(np.loadtxt(
        os.path.join(RUN_PATH, "crossbar_shd.csv"), dtype=int, ndmin=1))

    _, labels, weights = load_dataset()
    in_neurons = weights["fc1.weight"].shape[1]
    hidden_neurons = weights["fc1.weight"].shape[0]
    out_neurons = weights["fc2.weight"].shape[0]

    platform = "crossbar" if analog_synapses else "loihi"
    spike_filename = f"spikes_{platform}_shd.csv"
    with open(os.path.join(RUN_PATH, spike_filename)) as spike_csv:
        print("Processing spike data")
        spike_data = csv.DictReader(spike_csv)

        in_spikes = [[] for _ in range(0, in_neurons)]
        hidden_spikes = [[] for _ in range(0, hidden_neurons)]
        out_spikes = [[] for _ in range(0, out_neurons)]
        for spike in spike_data:
            # Spike entry has the format <group_name.neuron_id,timestep>
            timestep, neuron = int(spike["timestep"]), spike["neuron"]
            group_name, neuron_id = neuron.split(".")
            neuron_id = int(neuron_id)
            # Track spikes for all three layers
            if group_name == "in":
                in_spikes[neuron_id].append(timestep)
            elif group_name == "hidden":
                hidden_spikes[neuron_id].append(timestep)
            elif group_name == "out":
                out_spikes[neuron_id].append(timestep)
            else:
                print(f"Warning: Group {group_name} not recognized!")


    counts = np.zeros((out_neurons, num_inputs), dtype=int)
    for digit, spikes in enumerate(out_spikes):
        input_idx = 0
        end_idx = timesteps_per_input[0]
        for spike_timestep in spikes:
            # print(f"digit:{digit} spike timestep:{spike_timestep} end:{end_idx} input_idx:{input_idx}")
            if spike_timestep > end_idx:
                # Go to the next input, assuming spikes are in ascending timestep
                #  order
                input_idx += 1
                end_idx += timesteps_per_input[input_idx]
            assert(input_idx < num_inputs)
            counts[digit, input_idx] += 1

    correct = 0
    for i in range(0, num_inputs):
        # print(f"Spike counts per class for inference of digit:{counts[:, i]} "
        #     f"out:{np.argmax(counts[:, i])} actual:{labels[i]}")
        if np.argmax(counts[:, i]) == labels[i]:
            correct += 1

    accuracy = (correct / num_inputs) * 100
    print(f"Accuracy: {accuracy}% ({correct}/{num_inputs})")
    print(f"Min timesteps: {np.min(timesteps_per_input)}")
    print(f"Max timesteps: {np.max(timesteps_per_input)}")


def plot_spiking_digits():
    print("Plotting experiments")
    timesteps_per_input = list(np.loadtxt(
        os.path.join(RUN_PATH, "crossbar_shd.csv"),
        dtype=int, ndmin=1))
    raster_num_inputs = min(len(timesteps_per_input), 7)
    raster_total_timesteps = sum(timesteps_per_input[0:raster_num_inputs])

    plt.rcParams.update({
        "font.size": 8,
        "font.family": "sans-serif",
        "font.sans-serif": "Arial",
        "pdf.fonttype": 42
    })

    with open(os.path.join(RUN_PATH, "spikes_crossbar_shd.csv")) as spike_csv:
        print("Processing spike data")
        spike_data = csv.DictReader(spike_csv)
        in_neurons = 70
        hidden_neurons = 256
        out_neurons = 20

        in_spikes = [[] for _ in range(0, in_neurons)]
        hidden_spikes = [[] for _ in range(0, hidden_neurons)]
        out_spikes = [[] for _ in range(0, out_neurons)]
        for spike in spike_data:
            # Spike entry has the format <group_name.neuron_id,timestep>
            timestep, neuron = int(spike["timestep"]), spike["neuron"]
            group_name, neuron_id = neuron.split(".")
            neuron_id = int(neuron_id)
            # Track spikes for all three layers
            if group_name == "in":
                in_spikes[neuron_id].append(timestep)
            elif group_name == "hidden":
                hidden_spikes[neuron_id].append(timestep)
            elif group_name == "out":
                out_spikes[neuron_id].append(timestep)
            else:
                print(f"Warning: Group {group_name} not recognized!")

    fig = plt.figure(figsize=(6.5, 2.5))
    gs = GridSpec(2, 1, height_ratios=[1.0, 0.7], hspace=0.1, left=0.08, right=0.92, top=0.99, bottom=0.15)

    ax_spikes = fig.add_subplot(gs[0])
    ax_spikes.set_xlim((0, raster_total_timesteps))
    ax_spikes.set_ylim((0, in_neurons + hidden_neurons + out_neurons + 2))

    total_neurons = 0
    for neuron_id in range(0, in_neurons):
        ax_spikes.scatter(in_spikes[neuron_id],
                          [total_neurons] * len(in_spikes[neuron_id],),
                          rasterized=True,
                          c='k', s=1, marker='.', linewidths=0.2)
        total_neurons += 1

    for neuron_id in range(0, hidden_neurons):
        ax_spikes.scatter(hidden_spikes[neuron_id],
                          [total_neurons] * len(hidden_spikes[neuron_id]),
                          rasterized=True,
                          c='k', s=1, marker='.', linewidths=0.2)
        total_neurons += 1

    for neuron_id in range(0, out_neurons):
        ax_spikes.scatter(out_spikes[neuron_id],
                          [total_neurons] * len(out_spikes[neuron_id]),
                          rasterized=True,
                          c='k', s=1, marker='.', linewidths=0.2)
        total_neurons += 1

    ax_spikes.set_ylabel("Spiking Neuron ID")
    # Add brackets on the right side
    bracket_x = raster_total_timesteps * 1.015  # Slightly beyond the plot edge

    bracket_width = raster_total_timesteps * 0.01
    input_verts = [
        (bracket_x - bracket_width, 3),
        (bracket_x, 3),
        (bracket_x, in_neurons - 3),
        (bracket_x - bracket_width, in_neurons - 3)
    ]
    input_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]
    input_path = Path(input_verts, input_codes)
    input_patch = patches.PathPatch(input_path, facecolor='none', edgecolor='k', linewidth=1.0, clip_on=False)
    ax_spikes.add_patch(input_patch)

    ax_spikes.text(bracket_x + raster_total_timesteps*0.005, in_neurons/2, "Inputs",
                va='center', fontsize=7)

    # Hidden neurons bracket - single path
    hidden_verts = [
        (bracket_x - bracket_width, in_neurons + 3),
        (bracket_x, in_neurons + 3),
        (bracket_x, in_neurons + hidden_neurons - 3),
        (bracket_x - bracket_width, in_neurons + hidden_neurons - 3)
    ]
    hidden_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]
    hidden_path = Path(hidden_verts, hidden_codes)
    hidden_patch = patches.PathPatch(hidden_path, facecolor='none', edgecolor='k', linewidth=1.0, clip_on=False)
    ax_spikes.add_patch(hidden_patch)

    ax_spikes.text(bracket_x + raster_total_timesteps*0.005, in_neurons + hidden_neurons/2, "Hidden\nNeurons",
                va='center', fontsize=7)

    # Output neurons bracket - single path
    hidden_verts = [
        (bracket_x - bracket_width, hidden_neurons + in_neurons + 3),
        (bracket_x, hidden_neurons + in_neurons + 3),
        (bracket_x, total_neurons - 3),
        (bracket_x - bracket_width, total_neurons - 3)
    ]
    hidden_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]
    hidden_path = Path(hidden_verts, hidden_codes)
    hidden_patch = patches.PathPatch(hidden_path, facecolor='none', edgecolor='k', linewidth=1.0, clip_on=False)
    ax_spikes.add_patch(hidden_patch)

    ax_spikes.text(bracket_x + raster_total_timesteps*0.005,
                   in_neurons + hidden_neurons + out_neurons/2, "Outputs",
                va='center', fontsize=7)

    ax_spikes.set_xticks([])

    ax_perf = fig.add_subplot(gs[1])
    ax_perf.set_ylabel("Simulated Energy (µJ)")

    analog_perf_df = pd.read_csv(os.path.join(RUN_PATH, "perf_crossbar_shd.csv"))
    analog_perf_df["total_energy_uj"] = analog_perf_df["total_energy"] * 1.0e6

    loihi_perf_df = pd.read_csv(os.path.join(RUN_PATH, "perf_loihi_shd.csv"))
    analog_perf_df["loihi_energy_uj"] = loihi_perf_df["total_energy"] * 1.0e6

    analog_perf_df.plot(x="timestep", y=["loihi_energy_uj", "total_energy_uj"],
                        ax=ax_perf, color=[okabe_ito_colors[4], okabe_ito_colors[2]],
                        style=["--", "-"], legend=False)
    ax_perf.set_xlim((0, raster_total_timesteps))
    ax_perf.set_xlabel("Time-step")
    #ax_perf.legend(("Loihi", "Loihi-IMAC"), fontsize=7, bbox_to_anchor=(0.5, 1.0), handlelength=1.9)
    ax_perf.minorticks_on()

    # Add colored annotations on the right side for energy subplot
    energy_label_x = raster_total_timesteps * 1.002
    final_loihi_energy = analog_perf_df["loihi_energy_uj"].iloc[-1]
    final_analog_energy = analog_perf_df["total_energy_uj"].iloc[-1]

    ax_perf.text(energy_label_x + raster_total_timesteps*0.005, final_loihi_energy, "Loihi",
                va='center', fontsize=7, color=okabe_ito_colors[4])
    ax_perf.text(energy_label_x + raster_total_timesteps*0.005, final_analog_energy, "Loihi-IMAC",
                va='center', fontsize=7, color=okabe_ito_colors[2])

    plt.savefig(os.path.join(RUN_PATH, "shd_crossbar_raster.png"), dpi=300)
    plt.savefig(os.path.join(RUN_PATH, "shd_crossbar_raster.pdf"), dpi=300)

    analog_perf_df["total_energy_uj"] = analog_perf_df["total_energy"] * 1.0e6
    analog_perf_df["soma_energy_uj"] = analog_perf_df["soma_energy"] * 1.0e6
    print("***")
    mean_energy = analog_perf_df['total_energy'].mean()
    print(f"Mean Total Crossbar energy: {mean_energy} J (100 %)")
    print(f"Mean Soma Crossbar energy: {analog_perf_df['soma_energy'].mean()} J ({100.0 * analog_perf_df['soma_energy'].mean() / mean_energy} %)")
    print(f"Mean Synapse Crossbar energy: {analog_perf_df['synapse_energy'].mean()} J ({100.0 * analog_perf_df['synapse_energy'].mean() / mean_energy} %)")
    print(f"Mean Network Crossbar energy: {analog_perf_df['network_energy'].mean()} J ({100.0 * analog_perf_df['network_energy'].mean() / mean_energy} %)")
    print(f"Mean Crossbar firing neurons: {analog_perf_df['fired'].mean()}")

    print("***")
    mean_energy = loihi_perf_df['total_energy'].mean()
    print(f"Mean Total Loihi energy: {loihi_perf_df['total_energy'].mean()} J (100 %)")
    print(f"Mean Soma Loihi energy: {loihi_perf_df['soma_energy'].mean()} J ({100.0 * loihi_perf_df['soma_energy'].mean() / mean_energy} %)")
    print(f"Mean Synapse Loihi energy: {loihi_perf_df['synapse_energy'].mean()} J ({100.0 * loihi_perf_df['synapse_energy'].mean() / mean_energy} %)")
    print(f"Mean Network Loihi energy: {loihi_perf_df['network_energy'].mean()} J ({100.0 * loihi_perf_df['network_energy'].mean() / mean_energy} %)")
    print(f"Mean Loihi firing neurons: {loihi_perf_df['fired'].mean()}")


if __name__ == "__main__":
    num_inputs = 2264  # Number of inferences
    #num_inputs = 10
    run_experiments = False

    if run_experiments:
        run_spiking_digits(num_inputs=num_inputs, analog_synapses=True)
        run_spiking_digits(num_inputs=num_inputs, analog_synapses=False)

    calculate_accuracy(analog_synapses=False)
    calculate_accuracy(analog_synapses=True)

    plot_spiking_digits()

# Legacy scripts that may be useful at some point to someone, but is not
#  currently needed for any publication
"""

# MNIST Data-set info (circuit-aware trained 23 Sep 2025 for 10 epochs)
# Train set: Average loss: 0.1868, Accuracy: 56680/60000 (94.47%)
# Test set: Accuracy: 9475/10000 (94.75%)

def run_mnist(analog_synapses=True, timesteps=30):
    print(f"Loading models for binarized MNIST")
    # Data loading
    transform = transforms.Compose([
        transforms.ToTensor(),
        transforms.Resize((20, 20)),
        transforms.Normalize((0,), (1,))  # MNIST normalization
    ])
    test_dataset = datasets.MNIST("./runs/lasana/data", train=False,
                                transform=transform, download=True)
    dataloader = DataLoader(test_dataset, batch_size=1, shuffle=False)

    inputs = []
    labels = []
    count = 0
    for input_frame, label in dataloader:
        input_frame = input_frame.cpu()
        label = label.cpu()
        inputs.append(input_frame.numpy())
        labels.append(int(label))
        count += 1
        if count >= num_inputs:
            break

    if analog_synapses:
        snn_path = "/home/usr1/jboyle/neuro/binary_mnist/binarized_mnist_model.pth"
        snn = torch.load(snn_path, weights_only=True)
    else:
        # Load a normally trained rate-coded MNIST model to deploy on Loihi
        #  (normal digital) synapses
        snn = torch.load(
                os.path.join(PROJECT_DIR, "etc", "mnist_400_128_84_10.pt"),
                pickle_module=dill,
                map_location=torch.device("cpu")).state_dict()

    weights = {}
    if analog_synapses:
        weights["fc1"] = ternary_weight(snn["fc1.weight"])
        weights["fc2"] = ternary_weight(snn["fc2.weight"])
        weights["fc3"] = ternary_weight(snn["fc_out.weight"])

        crossbar_biases = {}
        crossbar_biases[1] = ternary_weight(snn["fc1.bias"])
        crossbar_biases[2] = ternary_weight(snn["fc2.bias"])
        crossbar_biases[3] = ternary_weight(snn["fc_out.bias"])
    else:
        weights["fc1"] = snn["fc1.weight"]
        weights["fc2"] = snn["fc2.weight"]
        weights["fc3"] = snn["fc_out.weight"]


    decays = {}
    decays[1] = float(snn["lif1.beta"])
    decays[2] = float(snn["lif2.beta"])
    decays[3] = float(snn["lif_out.beta"])


    thresholds = {}
    thresholds[1] = float(snn["lif1.threshold"])
    thresholds[2] = float(snn["lif2.threshold"])
    thresholds[3] = float(snn["lif_out.threshold"])

    # PyTorch stores weights in an array with dims (num out x num in)
    in_neurons = weights["fc1"].shape[1]
    hidden_1_neurons = weights["fc1"].shape[0]
    hidden_2_neurons = weights["fc2"].shape[0]
    out_neurons = weights["fc3"].shape[0]

    print(f"in:{in_neurons}, hidden 1:{hidden_1_neurons}, hidden 2:{hidden_2_neurons}, out:{out_neurons}")
    # Load the LASANA architecture with analog neurons
    if analog_synapses:
        arch = sanafe.load_arch("/home/usr1/jboyle/neuro/lasana/crossbar/lasana_crossbar.yaml")
    else:
        arch = sanafe.load_arch(os.path.join(PROJECT_DIR, "arch", "loihi.yaml"))

    print("Creating network in SANA-FE")
    network = sanafe.Network()
    in_layer = network.create_neuron_group("in", in_neurons)
    hidden_layer_1 = network.create_neuron_group("hidden_1", hidden_1_neurons)
    hidden_layer_2 = network.create_neuron_group("hidden_2", hidden_2_neurons)
    out_layer = network.create_neuron_group("out", out_neurons)

    # Defaults for Loihi baseline
    synapse_hw = "loihi_dense_synapse"
    dendrite_hw = "loihi_dendrites"

    print("Creating input layer")
    for id, neuron in enumerate(in_layer):
        neuron.set_attributes(dendrite_hw_name=dendrite_hw,
                            soma_hw_name=f"loihi_inputs[{id}]")
        neuron.map_to_core(arch.tiles[0].cores[0])

    print("Creating hidden layer")
    soma_hw = "loihi_lif"
    for id, neuron in enumerate(hidden_layer_1):
        hidden_parameters = {
            "threshold": thresholds[1],
            "leak_decay": decays[1],
            "reset_mode": "hard"
        }
        if analog_synapses:
            hidden_parameters["crossbar_bias"] = int(crossbar_biases[1][id])
            synapse_hw = f"synapse_crossbar[{id}]"
            dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=soma_hw,
                              dendrite_hw_name=dendrite_hw,
                              log_spikes=True,
                              default_synapse_hw_name=synapse_hw,
                              model_attributes=hidden_parameters)
        neuron.map_to_core(arch.tiles[0].cores[1])

    for id, neuron in enumerate(hidden_layer_2):
        hidden_parameters = {
            "threshold": thresholds[2],
            "leak_decay": decays[2],
            "reset_mode": "hard"
        }
        if analog_synapses:
            hidden_parameters["crossbar_bias"] = int(crossbar_biases[2][id])
            synapse_hw = f"synapse_crossbar[{id}]"
            dendrite_hw = synapse_hw


        neuron.set_attributes(soma_hw_name=soma_hw,
                              dendrite_hw_name=dendrite_hw,
                              log_spikes=True,
                              default_synapse_hw_name=synapse_hw,
                              model_attributes=hidden_parameters)
        neuron.map_to_core(arch.tiles[0].cores[2])

    print("Creating output layer")
    for id, neuron in enumerate(out_layer):
        hidden_parameters = {
            "threshold": thresholds[3],
            "leak_decay": decays[3],
            "reset_mode": "hard"
        }
        if analog_synapses:
            hidden_parameters["crossbar_bias"] = int(crossbar_biases[3][id])
            synapse_hw = f"synapse_crossbar[{id}]"
            dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=soma_hw,
                            dendrite_hw_name=dendrite_hw,
                            default_synapse_hw_name=synapse_hw,
                            log_spikes=True,
                            model_attributes=hidden_parameters)
        neuron.map_to_core(arch.tiles[0].cores[3])

    # Connect neurons in both layers

    print("Connecting neurons")
    print("Adding first layer connections")

    for src in range(in_neurons):
        for dst in range(hidden_1_neurons):
            connection_parameters = {}
            if analog_synapses:
                connection_parameters["weight"] = int(weights["fc1"][dst, src])
                connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            else:
                connection_parameters["weight"] = float(weights["fc1"][dst, src])

            network["in"][src].connect_to_neuron(
                network["hidden_1"][dst], connection_parameters)

    for src in range(hidden_1_neurons):
        for dst in range(hidden_2_neurons):
            connection_parameters = {}
            if analog_synapses:
                connection_parameters["weight"] = int(weights["fc2"][dst, src])
                connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            else:
                connection_parameters["weight"] = float(weights["fc2"][dst, src])

            network["hidden_1"][src].connect_to_neuron(
                network["hidden_2"][dst], connection_parameters)

    for src in range(hidden_2_neurons):
        for dst in range(out_neurons):
            connection_parameters = {}
            if analog_synapses:
                connection_parameters["weight"] = int(weights["fc3"][dst, src])
                connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            else:
                connection_parameters["weight"] = float(weights["fc3"][dst, src])

            network["hidden_2"][src].connect_to_neuron(
                network["out"][dst], connection_parameters)

    network.save("snn/crossbar_mnist.yaml")

    # Run a simulation
    print("Building h/w")
    hw = sanafe.SpikingChip(arch)
    print("Loading SNN")
    hw.load(network)
    print(f"Running simulation for {timesteps} timesteps")

    timesteps_per_input = []
    mapped_inputs = hw.mapped_neuron_groups["in"]
    print("Setting inputs")
    for input_idx in range(num_inputs):
        print(f"Simulating input: {input_idx}")
        mnist_input = inputs[input_idx].flatten()
        for id, mapped_neuron in enumerate(mapped_inputs):
            # print(list(inputs[:, id]))
            mapped_neuron.set_model_attributes(#model_attributes={"spikes": list(inputs[:, id])})
                model_attributes={"rate": mnist_input[id]})

        print(f"Simulating for {timesteps} timesteps")
        is_first_input = (input_idx == 0)
        results = hw.sim(timesteps,
                        spike_trace="spikes.csv",
                        perf_trace="perf.csv",
                        write_trace_headers=is_first_input,
                        processing_threads=4,
                        scheduler_threads=8)
        timesteps_per_input.append(timesteps)
        # print(results)
        hw.reset()


def plot_mnist(num_inputs):
    with open("spikes.csv") as spike_csv:
        spike_data = csv.DictReader(spike_csv)

        in_spikes = [[] for _ in range(0, in_neurons)]
        hidden_spikes_1 = [[] for _ in range(0, hidden_neurons)]
        hidden_spikes_2 = [[] for _ in range(0, hidden_2_neurons)]
        out_spikes = [[] for _ in range(0, out_neurons)]
        for spike in spike_data:
            # Spike entry has the format <group_name.neuron_id,timestep>
            timestep, neuron = int(spike["timestep"]), spike["neuron"]
            group_name, neuron_id = neuron.split(".")
            neuron_id = int(neuron_id)
            # Track spikes for all three layers
            if group_name == "in":
                in_spikes[neuron_id].append(timestep)
            elif group_name == "hidden_1":
                hidden_spikes_1[neuron_id].append(timestep)
            elif group_name == "hidden_2":
                hidden_spikes_2[neuron_id].append(timestep)
            elif group_name == "out":
                out_spikes[neuron_id].append(timestep)
            else:
                print(f"Warning: Group {group_name} not recognized!")

    # 1) Count the total output spikes
    counts = np.zeros((out_neurons, num_inputs), dtype=int)
    for digit, spikes in enumerate(out_spikes):
        for spike_timestep in spikes:
            # Get which input we're dealing with assuming we execute for a fixed
            #  number of timesteps on each input
            input_idx = (spike_timestep - 1) // timesteps
            assert(input_idx < num_inputs)
            counts[digit, input_idx] += 1

    correct = 0
    for i in range(0, num_inputs):
        print(f"Spike counts per class for inference of digit:{counts[:, i]} "
            f"out:{np.argmax(counts[:, i])} actual:{labels[i]}")
        if np.argmax(counts[:, i]) == labels[i]:
            correct += 1

    accuracy = (correct / num_inputs) * 100
    print(f"Accuracy: {accuracy}%")
"""