"""
Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Explore the usage of analog neuron circuits [Indiveri, 2003] within a Loihi-
based architecture running rate-encoded MNIST and spiking digits applications.
"""
import csv
import torch
import os
import sys
import dill
import numpy as np
import pandas as pd
# PyTorch and Benchmark data
import torchvision
from torchvision import datasets
import tonic
import tonic.transforms as transforms
# Plotting
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
import matplotlib.patches as patches

# Try importing the installed sanafe library. If not installed, require a
#  fall-back to a local build of the sanafe Python library
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))

    LASANA_DIR = os.path.abspath(os.path.join(
        "/", "home", "usr1", "jboyle", "neuro", "lasana"))
    RUN_PATH =  os.path.abspath((os.path.join(PROJECT_DIR, "runs", "lasana")))
    print(f"Project dir: {PROJECT_DIR}")
    sys.path.insert(0, PROJECT_DIR)
    import sanafe


def load_dataset(dataset, analog_neurons):
    if dataset == "mnist":
        weights = {}
        # Load the MNIST torch network, which is the same for both Loihi and
        #  analog implementations. This was trained with some circuit knowledge,
        #  i.e., by clipping the snntorch "Leaky" neuron to be nonnegative. This
        #  clipping behavior is supported in Loihi, so we have a common SNN that
        #  works reasonably well for both Indiveri and Loihi neurons. We found
        #  that MNIST is generally not super sensitive to small circuit
        #  variations (unlike spiking digits)
        #
        # This network was trained using the script: TODO
        mnist_model = torch.load(
                os.path.join(PROJECT_DIR, "etc", "mnist_784_128_10.pt"),
                pickle_module=dill,
                map_location=torch.device("cpu"))
        for attribute_name, param in mnist_model.named_parameters():
            weights[attribute_name] = param.detach().numpy()

        transform = torchvision.transforms.Compose([
            torchvision.transforms.Resize((28, 28)),
            torchvision.transforms.Grayscale(),
            torchvision.transforms.ToTensor(),
            torchvision.transforms.Normalize((0,), (1,))])
        test_dataset = datasets.MNIST(root="./runs/lasana/data", train=False,
                                      download=True, transform=transform)
        labels = test_dataset.targets

        dataloader = torch.utils.data.DataLoader(
            test_dataset, batch_size=1, shuffle=False)
        count = 0
        inputs = np.zeros((num_inputs, 28*28))
        for data, targets in dataloader:
            data = data.cpu()
            inputs[count, :] = data.numpy().flatten()
            count += 1
            if count >= num_inputs:
                break

    elif dataset == "shd":
        weights = {}

        # Load the Spiking Heidelberg Digits (SHD) network
        if analog_neurons: # Use Indiveri neural circuits
            # Use a SNN trained using circuit-aware methods i.e. was
            #  trained specifically on an Indiveri neuron circuit model.
            #  This training script (spiking_digits_indiveri.py) was provided by
            #  Jason Ho and uses LASANA surrogate model layers in between
            #  regular Linear synapses.
            spiking_digits_model = torch.load(
                    os.path.join(PROJECT_DIR, "etc", "spiking_digits_indiveri.pt"),
                    weights_only=True,
                    map_location=torch.device("cpu"))

            for attribute_name, param in spiking_digits_model.items():
                weights[attribute_name] = param.detach().numpy()
        else: # Use digital LIF neuron models
            # Use an SNN trained on a generic/abstract LIF neuron, applying clipping
            #  to mimic the Indiveri behavior as close as is possible with an LIF
            #  (Leaky in SNNTorch) neuron
            spiking_digits = torch.load(
                    os.path.join(PROJECT_DIR, "etc", "spiking_digits_lif.pt"),
                    pickle_module=dill,
                    map_location=torch.device("cpu"))

            for attribute_name, param in spiking_digits.named_parameters():
                weights[attribute_name] = param.detach().numpy()

        # Load the spiking digits test inputs from the raw dataset, applying the
        #  same transformations as the training scripts
        frame_transform = transforms.Compose([
            transforms.Downsample(time_factor=0.05, spatial_factor=0.1),  # For Jason's circuit-aware trained SNN
            transforms.ToFrame(sensor_size=(70, 1, 1), time_window=1000)
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
    else:
        print(f"Dataset not recognized: {dataset}")

    return inputs, labels, weights


def create_net(arch, dataset, weights, analog_neurons):
    # PyTorch stores weights in an array with dims (num out x num in)
    in_neurons = weights["fc1.weight"].shape[1]
    hidden_neurons = weights["fc1.weight"].shape[0]
    out_neurons = weights["fc2.weight"].shape[0]
    print(f"in:{in_neurons}, hidden:{hidden_neurons}, out:{out_neurons}")

    print("Creating network in SANA-FE")
    network = sanafe.Network()
    print(f"in:{in_neurons} hidden:{hidden_neurons} out:{out_neurons}")
    in_layer = network.create_neuron_group("in", in_neurons)
    hidden_layer = network.create_neuron_group("hidden", hidden_neurons)
    out_layer = network.create_neuron_group("out", out_neurons)

    print("Creating input layer")
    for id, neuron in enumerate(in_layer):
        neuron.set_attributes(soma_hw_name=f"input[{id}]",
                                log_spikes=True)
        neuron.map_to_core(arch.tiles[0].cores[0])

    print("Creating hidden layer")
    for id, neuron in enumerate(hidden_layer):
        hidden_parameters = {
            "threshold": 1.0,
            "leak_decay": 0.85,
            "reset_mode": "hard",
        }
        if dataset == "shd":
            # In the previously trained network lif1.alpha was available in the
            #  trained SNN
            if analog_neurons:
                hidden_parameters["input_decay"] = 0.0
            else:
                hidden_parameters["input_decay"] = weights["lif1.alpha"][id]

        if analog_neurons:
            neuron.set_attributes(
                soma_hw_name=f"analog_lif[{id + hidden_neurons}]",
                log_spikes=True,
                log_potential=True,
                model_attributes=hidden_parameters
            )
        else:
            if dataset == "shd":
                # Clip anything below 0
                hidden_parameters["reverse_threshold"] = 0.0
                hidden_parameters["reverse_reset_mode"] = "saturate"
            neuron.set_attributes(
                soma_hw_name=f"loihi_lif",
                log_spikes=True,
                model_attributes=hidden_parameters
            )
        neuron.map_to_core(arch.tiles[0].cores[0])

    print("Creating output layer")
    for id, neuron in enumerate(out_layer):
        if dataset == "shd":
            neuron.set_attributes(
                soma_hw_name=f"loihi_lif",
                log_potential=True,
                model_attributes={
                    "leak_decay": 0.85,
                    "reset_mode": "none",
                    "input_decay": weights["lif2.alpha"][id]
                }
            )
        elif dataset == "mnist":
            if analog_neurons:
                neuron.set_attributes(
                    soma_hw_name=f"analog_lif[{id}]",
                    log_spikes=True
                )
            else:
                neuron.set_attributes(
                    soma_hw_name=f"loihi_lif",
                    log_spikes=True,
                    model_attributes={
                        "threshold": 1.0,
                        "leak_decay": 0.85,
                        "reset_mode": "hard"
                    }
                )
        else:
            raise Exception(f"Dataset: {dataset} not recognized!")
        neuron.map_to_core(arch.tiles[0].cores[0])

    # Connect neurons in both layers
    print("Connecting neurons")
    min_weight = 1.0e-10
    print("Adding first layer connections")
    for src in range(in_neurons):
        for dst in range(hidden_neurons):
            weight = weights["fc1.weight"][dst, src]
            if abs(weight) > min_weight:
                network["in"][src].connect_to_neuron(
                    network["hidden"][dst], {"weight": weight})

    if dataset == "shd":  # Add recurrent connections only for spiking digits
        print("Adding recurrent connections")
        for src in range(hidden_neurons):
            for dst in range(hidden_neurons):
                weight = weights["fcr.weight"][dst, src]
                if abs(weight) > min_weight:
                    network["hidden"][src].connect_to_neuron(
                        network["hidden"][dst], {"weight": weight}
                    )

    print("Adding second layer connections")
    for src in range(hidden_neurons):
        for dst in range(out_neurons):
            weight = weights["fc2.weight"][dst, src]
            if abs(weight) > min_weight:
                network["hidden"][src].connect_to_neuron(
                    network["out"][dst], {"weight": weight})

    return network


def run_experiment(num_inputs, dataset="shd", analog_neurons=True):
    print(f"Loading models for {dataset}")

    inputs, labels, weights = load_dataset(dataset, analog_neurons)
    # Load the LASAGNA architecture with analog neurons
    arch = sanafe.load_arch(
        os.path.abspath(os.path.join(LASANA_DIR, "indiveri", "lasana.yaml")))
    snn = create_net(arch, dataset, weights, analog_neurons)

    platform = "analog" if analog_neurons else "loihi"
    spike_filename = f"spikes_{platform}_{dataset}.csv"
    perf_filename = f"perf_{platform}_{dataset}.csv"
    potential_filename = f"potential_{platform}_{dataset}.csv"

    # Run a simulation
    print("Building h/w")
    hw = sanafe.SpikingChip(arch)
    print("Loading SNN")
    hw.load(snn)

    timesteps_per_input = []
    mapped_inputs = hw.mapped_neuron_groups["in"]
    print("Setting inputs")
    for input in range(num_inputs):
        print(f"Simulating input: {input}")

        if dataset == "mnist":
            timesteps = 100
            mnist_input = inputs[input, :]
            for id, mapped_neuron in enumerate(mapped_inputs):
                mapped_neuron.set_model_attributes(
                    model_attributes={"poisson": mnist_input[id]})
        elif dataset == "shd":
            spiking_digit_input = inputs[input].squeeze()
            for id, mapped_neuron in enumerate(mapped_inputs):
                spiketrain = list(spiking_digit_input[:, id])
                mapped_neuron.set_model_attributes(
                    model_attributes={"spikes": spiketrain})
            timesteps = len(spiketrain)
        else:
            raise Exception(f"Dataset: {dataset} not recognized!")

        print(f"Simulating for {timesteps} timesteps")
        is_first_input = (input == 0)
        hw.sim(timesteps, spike_trace=os.path.join(RUN_PATH, spike_filename),
               potential_trace=os.path.join(RUN_PATH, potential_filename),
               perf_trace=os.path.join(RUN_PATH, perf_filename),
               write_trace_headers=is_first_input, processing_threads=16)
        timesteps_per_input.append(timesteps)
        hw.reset()

    platform = "analog" if analog_neurons else "loihi"
    snn.save(os.path.join(PROJECT_DIR, "runs", "lasana", f"indiveri_{platform}_{dataset}.yaml"))
    # The inputs (and therefore timesteps per input) will be the same across
    #  Loihi/analog neuron runs. Only store this once.
    np.savetxt(os.path.join(PROJECT_DIR, "runs", "lasana", f"indiveri_{dataset}.csv"),
               np.array(timesteps_per_input), fmt="%d")
    return


def calculate_cumulative_spikes(out_spikes, out_neurons, total_timesteps,
                                timesteps_per_input):
    # Initialize array to store cumulative counts for each output neuron at each
    #  timestep
    cumulative_counts = np.zeros((out_neurons, total_timesteps))

    # For each output neuron
    for neuron_idx, spikes in enumerate(out_spikes):
        # Sort spikes by timestep for cumulative counting
        sorted_spikes = sorted([s for s in spikes if s <= total_timesteps])
        # Process each timestep
        for t in range(total_timesteps):
            # Calculate which input period we're in
            input_image = t // timesteps_per_input
            period_start = input_image * timesteps_per_input

            # Count spikes from start of current input period up to current
            #  timestep
            count = sum(1 for spike in sorted_spikes
                    if spike <= t + 1 and spike > period_start)

            cumulative_counts[neuron_idx, t] = count

    return cumulative_counts


# we need code to
#  1) work out and count the total spikes for each output for all runs
#  2) use that to calculate accuracy
#  3) separately, use other data to plot energy per timestep (just take the first 10)

# We can store the spike and perf data somewhere in the run directory, and plot from this

def calculate_accuracy(dataset, analog_neurons, timesteps_per_input):
    arch = "analog" if analog_neurons else "loihi"
    print(f"Calculating accuracy for {arch} architecture")

    _, labels, weights = load_dataset(dataset, analog_neurons)
    in_neurons = weights["fc1.weight"].shape[1]
    hidden_neurons = weights["fc1.weight"].shape[0]
    out_neurons = weights["fc2.weight"].shape[0]

    spike_filename = f"spikes_{arch}_{dataset}.csv"
    perf_filename = f"perf_{arch}_{dataset}.csv"
    potential_filename = f"potential_{arch}_{dataset}.csv"

    if dataset == "mnist":
        # Read in simulation results
        with open(os.path.join(RUN_PATH, spike_filename)) as spike_csv:
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
            for spike_timestep in spikes:
                input_idx = (spike_timestep - 1) // timesteps_per_input[0]
                assert(input_idx < num_inputs)
                counts[digit, input_idx] += 1

        correct = 0
        for i in range(0, num_inputs):
            print(f"Spike counts per class for inference of digit:{counts[:, i]} "
                f"out:{np.argmax(counts[:, i])} actual:{labels[i]}")
            if np.argmax(counts[:, i]) == labels[i]:
                correct += 1

    elif dataset == "shd":
        # This code is to plot a raster not to calculate accuracy
        columns_to_extract = [
            "neuron out.0", "neuron out.1", "neuron out.2", "neuron out.3",
            "neuron out.4", "neuron out.5", "neuron out.6", "neuron out.7",
            "neuron out.8", "neuron out.9", "neuron out.10", "neuron out.11",
            "neuron out.12", "neuron out.13", "neuron out.14", "neuron out.15",
            "neuron out.16", "neuron out.17", "neuron out.18", "neuron out.19"]
        potentials_df = pd.read_csv(os.path.join(RUN_PATH, potential_filename))
        potentials = potentials_df[columns_to_extract].to_numpy()

        correct = 0
        start_timestep = 0
        for i, timesteps in enumerate(timesteps_per_input):
            timestep_potentials = potentials[start_timestep:start_timestep+timesteps, :]
            max_potentials = np.max(timestep_potentials, axis=0)
            category = np.argmax(max_potentials, axis=0)

            # print(f"Max potentials per class for inference:{max_potentials}\n"
            #     f"out:{category} actual:{labels[i]}")
            if category == labels[i]:
                correct += 1
            start_timestep += timesteps

    accuracy = (correct / num_inputs) * 100
    print(f"Accuracy: {accuracy}% ({correct}/{num_inputs})")


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


def plot_experiments(dataset):
    # Plot both sets of experiments:
    #  Raster plot of a small number of inputs e.g. 10 inputs
    #  Time series plot of a small number of inputs
    #  Statistics for a longer run showing accuracy
    print("Plotting experiments")
    timesteps_per_input = list(np.loadtxt(
        os.path.join(PROJECT_DIR, "runs", "lasana", f"indiveri_{dataset}.csv"),
        dtype=int, ndmin=1))

    inputs, labels, weights = load_dataset(dataset, True)
    if dataset == "shd":
        plot_shd(timesteps_per_input, labels, weights)
    elif dataset == "mnist":
        plot_mnist(timesteps_per_input[0], inputs, weights)

    # Calculate accuracy for both sets of runs
    calculate_accuracy(dataset, True, timesteps_per_input)
    calculate_accuracy(dataset, False, timesteps_per_input)


def plot_shd(timesteps_per_input, labels, weights):
    analog_spike_filename = f"spikes_analog_shd.csv"
    analog_perf_filename = f"perf_analog_shd.csv"
    analog_potential_filename = f"potential_analog_shd.csv"

    loihi_perf_filename = f"perf_loihi_shd.csv"

    in_neurons = weights["fc1.weight"].shape[1]
    hidden_neurons = weights["fc1.weight"].shape[0]
    out_neurons = weights["fc2.weight"].shape[0]

    plt.rcParams.update({
        "font.size": 6,
        "font.family": "sans-serif",
        "font.sans-serif": "Arial",
        "pdf.fonttype": 42
    })

    # Create an array tracking the input corresponding to each timestep
    columns_to_extract = [
        "neuron out.0", "neuron out.1", "neuron out.2", "neuron out.3",
        "neuron out.4", "neuron out.5", "neuron out.6", "neuron out.7",
        "neuron out.8", "neuron out.9", "neuron out.10", "neuron out.11",
        "neuron out.12", "neuron out.13", "neuron out.14", "neuron out.15",
        "neuron out.16", "neuron out.17", "neuron out.18", "neuron out.19"]
    potentials_df = pd.read_csv(os.path.join(RUN_PATH, analog_potential_filename))
    potentials = potentials_df[columns_to_extract].to_numpy()

    analog_perf_df = pd.read_csv(os.path.join(RUN_PATH, analog_perf_filename))
    analog_perf_df["total_energy_uj"] = analog_perf_df["total_energy"] * 1.0e6
    loihi_perf_df = pd.read_csv(os.path.join(RUN_PATH, loihi_perf_filename))
    analog_perf_df["total_energy_loihi_uj"] = (loihi_perf_df["total_energy"] - loihi_perf_df["soma_energy"]) * 1.0e6

    raster_num_inputs = min(len(timesteps_per_input), 10)
    raster_total_timesteps = sum(timesteps_per_input[0:raster_num_inputs])

    with open(os.path.join(RUN_PATH, analog_spike_filename)) as spike_csv:
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

    fig = plt.figure(figsize=(7.2, 2.5))
    gs = GridSpec(2, 1, height_ratios=[1.0, 1.0], hspace=0.1, left=0.06, right=0.94, top=0.99, bottom=0.12)

    ax_spikes = fig.add_subplot(gs[0])
    ax_spikes.set_xlim((0, raster_total_timesteps))
    ax_spikes.set_ylim((0, in_neurons + hidden_neurons))

    total_neurons = 0
    for neuron_id in range(0, in_neurons):
        ax_spikes.scatter(in_spikes[neuron_id],
                          [total_neurons] * len(in_spikes[neuron_id]),
                          c='k', s=1, marker='.', linewidths=0.3)
        total_neurons += 1

    for neuron_id in range(0, hidden_neurons):
        ax_spikes.scatter(hidden_spikes[neuron_id],
                          [total_neurons] * len(hidden_spikes[neuron_id]),
                          c='k', s=1, marker='.', linewidths=0.3)
        total_neurons += 1

    ax_spikes.set_ylabel("Spiking Neuron ID")
    # Add brackets on the right side
    bracket_x = raster_total_timesteps * 1.015  # Slightly beyond the plot edge

    bracket_width = raster_total_timesteps * 0.01
    input_verts = [
        (bracket_x - bracket_width, 2),
        (bracket_x, 2),
        (bracket_x, in_neurons-2),
        (bracket_x - bracket_width, in_neurons-2)
    ]
    input_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]
    input_path = Path(input_verts, input_codes)
    input_patch = patches.PathPatch(input_path, facecolor='none', edgecolor='k', linewidth=1.0, clip_on=False)
    ax_spikes.add_patch(input_patch)

    ax_spikes.text(bracket_x + raster_total_timesteps*0.005, in_neurons/2, 'Input',
                va='center', fontsize=5)

    # Hidden neurons bracket - single path
    hidden_verts = [
        (bracket_x - bracket_width, in_neurons+2),
        (bracket_x, in_neurons+2),
        (bracket_x, total_neurons-2),
        (bracket_x - bracket_width, total_neurons-2)
    ]
    hidden_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]
    hidden_path = Path(hidden_verts, hidden_codes)
    hidden_patch = patches.PathPatch(hidden_path, facecolor='none', edgecolor='k', linewidth=1.0, clip_on=False)
    ax_spikes.add_patch(hidden_patch)

    ax_spikes.text(bracket_x + raster_total_timesteps*0.005, in_neurons + hidden_neurons/2, 'Hidden',
                va='center', fontsize=5)

    ax_spikes.set_xticks([])

    ax_perf = fig.add_subplot(gs[1])
    ax_perf.set_ylabel("Simulated Energy (µJ)")
    analog_perf_df.plot(x="timestep", y=["total_energy_loihi_uj", "total_energy_uj"],
                        ax=ax_perf, color=[okabe_ito_colors[0], okabe_ito_colors[4]],
                        style=["--", "-"])
    ax_perf.set_xlim((0, raster_total_timesteps))
    ax_perf.set_xlabel("Time-step")
    ax_perf.legend(("Loihi", "Loihi-Ind"), fontsize=5, bbox_to_anchor=(0.5, 1.0), handlelength=1.9)
    ax_perf.minorticks_on()

    # Currently, I'm not planning to put this figure anywhere. Other data shows
    #  that the vast majority of the energy consumption is due to everything
    #  apart from the soma. The majority of energy savings by using the Indiveri
    #  model is due to reduced spiking activity (likely due to the SNN being
    #  trained slightly differently). The energy savings seen in this plot are
    #  due to reduced synaptic energy usage...
    plt.savefig(os.path.join("runs", "lasana", "shd_indiveri_raster.png"), dpi=300)
    plt.savefig(os.path.join("runs", "lasana", "shd_indiveri_raster.pdf"))

    # Aggregate the metrics per inference
    per_input_metrics = []
    start_timestep = 0
    for input_idx in range(num_inputs):
        print(".", end="")
        end_timestep = start_timestep + timesteps_per_input[input_idx]

        # Filter data for this input's timestep range
        analog_mask = (analog_perf_df['timestep'] >= start_timestep) & (analog_perf_df['timestep'] < end_timestep)
        loihi_mask = (loihi_perf_df['timestep'] >= start_timestep) & (loihi_perf_df['timestep'] < end_timestep)

        # Calculate energy totals for this input
        indiveri_total = analog_perf_df.loc[analog_mask, 'total_energy'].sum()
        indiveri_synapse = analog_perf_df.loc[analog_mask, 'synapse_energy'].sum()
        indiveri_soma = analog_perf_df.loc[analog_mask, 'soma_energy'].sum()
        indiveri_network = analog_perf_df.loc[analog_mask, 'network_energy'].sum()

        loihi_total = loihi_perf_df.loc[loihi_mask, 'total_energy'].sum()
        loihi_synapse = loihi_perf_df.loc[loihi_mask, 'synapse_energy'].sum()
        loihi_soma = loihi_perf_df.loc[loihi_mask, 'soma_energy'].sum()
        loihi_network = loihi_perf_df.loc[loihi_mask, 'network_energy'].sum()

        # Count total spikes for this input
        indiveri_spikes = analog_perf_df.loc[analog_mask, 'fired'].sum()
        loihi_spikes = loihi_perf_df.loc[loihi_mask, 'fired'].sum()

        per_input_metrics.append({
            'input_idx': input_idx,
            'timesteps': timesteps_per_input[input_idx],
            'indiveri_total_energy': indiveri_total,
            'indiveri_synapse_energy': indiveri_synapse,
            'indiveri_soma_energy': indiveri_soma,
            'indiveri_network_energy': indiveri_network,
            'indiveri_spikes': indiveri_spikes,
            'loihi_total_energy': loihi_total,
            'loihi_synapse_energy': loihi_synapse,
            'loihi_soma_energy': loihi_soma,
            'loihi_network_energy': loihi_network,
            'loihi_spikes': loihi_spikes,
        })
        start_timestep = end_timestep

    # Print per-inference statistics summaries
    spiking_neurons = 70 + hidden_neurons
    print("")
    print(f"Per-Inference Averages (across {num_inputs} inputs):")
    avg_indiveri_total = sum(m['indiveri_total_energy'] for m in per_input_metrics) / num_inputs
    avg_indiveri_synapse = sum(m['indiveri_synapse_energy'] for m in per_input_metrics) / num_inputs
    avg_indiveri_soma = sum(m['indiveri_soma_energy'] for m in per_input_metrics) / num_inputs
    avg_indiveri_network = sum(m['indiveri_network_energy'] for m in per_input_metrics) / num_inputs
    avg_indiveri_spikes = (sum(m['indiveri_spikes'] / spiking_neurons for m in per_input_metrics)) / num_inputs

    print(f"Indiveri - Avg total energy per inference: {avg_indiveri_total:.6e} J (100.0 %)")
    print(f"Indiveri - Avg synapse energy per inference: {avg_indiveri_synapse:.6e} J ({100.0 * avg_indiveri_synapse / avg_indiveri_total:.1f} %)")
    print(f"Indiveri - Avg soma energy per inference: {avg_indiveri_soma:.6e} J ({100.0 * avg_indiveri_soma / avg_indiveri_total:.1f} %)")
    print(f"Indiveri - Avg network energy per inference: {avg_indiveri_network:.6e} J ({100.0 * avg_indiveri_network / avg_indiveri_total:.1f} %)")
    print(f"Indiveri - Avg spikes out per neuron per inference: {avg_indiveri_spikes:.1f}")
    print("***")

    avg_loihi_total = sum(m['loihi_total_energy'] for m in per_input_metrics) / num_inputs
    avg_loihi_synapse = sum(m['loihi_synapse_energy'] for m in per_input_metrics) / num_inputs
    avg_loihi_soma = sum(m['loihi_soma_energy'] for m in per_input_metrics) / num_inputs
    avg_loihi_network = sum(m['loihi_network_energy'] for m in per_input_metrics) / num_inputs
    avg_loihi_spikes = (sum(m['loihi_spikes'] / spiking_neurons for m in per_input_metrics)) / num_inputs

    print(f"Loihi - Avg total energy per inference: {avg_loihi_total:.6e} J (100.0 %)")
    print(f"Loihi - Avg synapse energy per inference: {avg_loihi_synapse:.6e} J ({100.0 * avg_loihi_synapse / avg_loihi_total:.1f} %)")
    print(f"Loihi - Avg soma energy per inference: {avg_loihi_soma:.6e} J ({100.0 * avg_loihi_soma / avg_loihi_total:.1f} %)")
    print(f"Loihi - Avg network energy per inference: {avg_loihi_network:.6e} J ({100.0 * avg_loihi_network / avg_loihi_total:.1f} %)")
    print(f"Loihi - Avg spikes out per neuron per inference: {avg_loihi_spikes:.1f}")
    print("***")

    # Print per-timestep statistics. I will probably not include this in the paper
    #  either but its otherwise interesting/useful to me
    total_indiveri_energy = analog_perf_df['total_energy'].mean()
    print("Per-timestep Averages")
    print(f"Indiveri - Mean per-timestep total energy: {total_indiveri_energy} J (100 %)")
    print(f"Indiveri - Mean per-timestep synapse energy: {analog_perf_df['synapse_energy'].mean()} J ({100.0 * analog_perf_df['synapse_energy'].mean() / total_indiveri_energy} %)")
    print(f"Indiveri - Mean per-timestep soma energy: {analog_perf_df['soma_energy'].mean()} J ({100.0 * analog_perf_df['soma_energy'].mean() / total_indiveri_energy} %)")
    print(f"Indiveri - Mean per-timestep network energy: {analog_perf_df['network_energy'].mean()} J ({100.0 * analog_perf_df['network_energy'].mean() / total_indiveri_energy} %)")
    print(f"Indiveri - Mean per-timestep fired: {analog_perf_df['fired'].mean()}")
    print("***")

    total_loihi_energy = loihi_perf_df['total_energy'].mean()
    print(f"Loihi - Mean per-timestep Loihi total energy: {total_loihi_energy} J (100 %)")
    print(f"Loihi - Mean per-timestep Loihi synapse energy: {loihi_perf_df['synapse_energy'].mean()} J ({100.0 * loihi_perf_df['synapse_energy'].mean() / total_loihi_energy} %)")
    print(f"Loihi - Mean per-timestep Loihi soma energy: {loihi_perf_df['soma_energy'].mean()} J ({100.0 * loihi_perf_df['soma_energy'].mean() / total_loihi_energy} %)")
    print(f"Loihi - Mean per-timestep Loihi network energy: {loihi_perf_df['network_energy'].mean()} J ({100.0 * loihi_perf_df['network_energy'].mean() / total_loihi_energy} %)")
    print(f"Loihi - Mean per-timestep Loihi fired: {loihi_perf_df['fired'].mean()}")
    print("***")

    # *** Other misc plots ***
    fig = plt.figure(figsize=(12.0, 8.0))
    gs = GridSpec(3, 1, height_ratios=[4, 4, 4], hspace=0.1)
    # Plot a detailed set of figures with raster, potentials and energies
    #  This isn't for publication but is useful/interesting either way
    ax_spikes = fig.add_subplot(gs[0])
    ax_spikes.set_title('Mixed-Signal Architecture Classifying Spiking Digits')
    ax_spikes.set_xlim((0, raster_total_timesteps))
    ax_spikes.set_ylim((0, in_neurons + hidden_neurons))

    total_neurons = 0
    for neuron_id in range(0, in_neurons):
        ax_spikes.scatter(in_spikes[neuron_id],
                    [total_neurons]*len(in_spikes[neuron_id]), c='r', s=2,
                    marker='.', linewidths=1)
        total_neurons += 1

    for neuron_id in range(0, hidden_neurons):
        ax_spikes.scatter(hidden_spikes[neuron_id],
                    [total_neurons]*len(hidden_spikes[neuron_id]), c='b', s=2,
                    marker='.', linewidths=1)
        total_neurons += 1

    # The output layer is a leaky-integrator layer with *no* spiking (the output
    #  is taken from the maximum potential on the integrators). Therefore,
    #  don't add output spikes (if any) to the raster
    ax_spikes.set_ylabel("Neuron")
    ax_spikes.set_xticks([])

    # Plot the output neuron spike counts
    ax_potentials = fig.add_subplot(gs[1])
    ax_potentials.set_xlim((0, raster_total_timesteps))
    ax_potentials.set_ylabel("Output Potentials")
    ax_potentials.set_xticks([])
    ax_potentials.set_xlabel("")  # Remove label
    ax_potentials.plot(potentials)

    ax_perf = fig.add_subplot(gs[2])
    ax_perf.set_ylabel("Simulated Energy (µJ)")
    analog_perf_df.plot(x="timestep", y=["total_energy_uj", "total_energy_loihi_uj"], ax=ax_perf)
    ax_perf.set_xlim((0, raster_total_timesteps))
    ax_perf.set_xlabel("Time-step")

    # Save a larger plot with more information for debug. This figure isn't
    #  intended for publication or presentation.
    plt.savefig(os.path.join("runs", "lasana", "shd.png"), dpi=300)


def plot_mnist(timesteps_per_input, inputs, weights):
    spike_filename = f"spikes_analog_mnist.csv"
    perf_filename = f"perf_analog_mnist.csv"

    in_neurons = weights["fc1.weight"].shape[1]
    hidden_neurons = weights["fc1.weight"].shape[0]
    out_neurons = weights["fc2.weight"].shape[0]

    # Read in simulation results
    with open(os.path.join(RUN_PATH, spike_filename)) as spike_csv:
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


    # 2) Create a raster plot of all spikes
    cumulative_counts = calculate_cumulative_spikes(
        out_spikes, out_neurons, timesteps_per_input * num_inputs,
        timesteps_per_input)

    fig = plt.figure(figsize=(12.0, 8.0))
    gs = GridSpec(4, 1, height_ratios=[2, 4, 4, 4], hspace=0.1)

    # Top subplot for MNIST digits
    ax_digits = fig.add_subplot(gs[0])
    ax_digits.set_xticks([])
    ax_digits.set_yticks([])

    # Calculate positions for each digit
    digit_width = 28  # MNIST digits are 28x28
    time_per_digit = timesteps_per_input

    # Create a blank canvas for all digits
    total_timesteps = timesteps_per_input * num_inputs
    display_height = 20  # Adjust this factor to change digit height

    # Place each digit at its corresponding time position
    for i, digit in enumerate(inputs[0:num_inputs, :]):
        start_time = i * time_per_digit
        # Calculate extent for each digit: [left, right, bottom, top]
        # Width of each digit display is set to match the time window
        digit_extent = [start_time, start_time + time_per_digit, 0, display_height]
        ax_digits.imshow(digit.reshape(digit_width, digit_width),
                        cmap='gray', aspect='auto', extent=digit_extent)

    # Display the digits
    ax_digits.set_xlim(0, total_timesteps)
    ax_digits.set_ylim(0, display_height)
    ax_digits.set_title('Mixed-Signal Architecture Classifying MNIST')

    ax_spikes = fig.add_subplot(gs[1])
    ax_spikes.set_xlim((0, total_timesteps))
    ax_spikes.set_ylim((0, in_neurons + hidden_neurons + out_neurons))

    total_neurons = 0
    for neuron_id in range(0, in_neurons):
        ax_spikes.scatter(in_spikes[neuron_id],
                    [total_neurons]*len(in_spikes[neuron_id]), c='r', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1

    for neuron_id in range(0, hidden_neurons):
        ax_spikes.scatter(hidden_spikes[neuron_id],
                    [total_neurons]*len(hidden_spikes[neuron_id]), c='b', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1

    for neuron_id in range(0, out_neurons):
        ax_spikes.scatter(out_spikes[neuron_id],
                    [total_neurons]*len(out_spikes[neuron_id]), c='k', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1

    # Add vertical lines to show digit presentation boundaries
    for i in range(num_inputs + 1):
        ax_spikes.axvline(x=i * time_per_digit, color='gray', linestyle='--', alpha=0.3)
    ax_spikes.set_ylabel("Neuron")
    ax_spikes.set_xticks([])

    # Plot the output neuron spike counts
    ax_potentials = fig.add_subplot(gs[2])
    ax_potentials.set_xlim((0, total_timesteps))
    ax_potentials.plot(cumulative_counts.transpose())
    # Add vertical lines to show digit presentation boundaries
    for i in range(num_inputs + 1):
        ax_potentials.axvline(x=i*time_per_digit + 1, color='gray', linestyle='--', alpha=0.3)
    ax_potentials.set_ylabel("Spike Counts")
    ax_potentials.set_xticks([])
    ax_potentials.set_xlabel("")

    ax_perf = fig.add_subplot(gs[3])
    ax_perf.set_xlim((0, total_timesteps))
    perf_df = pd.read_csv(os.path.join(RUN_PATH, perf_filename))
    perf_df["total_energy_uj"] = perf_df["total_energy"] * 1.0e6
    perf_df["soma_energy_uj"] = perf_df["soma_energy"] * 1.0e6
    perf_df.plot(x="timestep", y=["total_energy_uj", "soma_energy_uj"], ax=ax_perf)
    #perf_df.plot(x="timestep", y=["soma_energy_uj",], ax=ax_perf)
    for i in range(num_inputs + 1):
        ax_perf.axvline(x=i*time_per_digit + 1, color='gray', linestyle='--', alpha=0.3)
    ax_perf.set_ylabel("Simulated Energy (uJ)")
    ax_perf.get_legend().remove()

    ax_perf.set_xlabel("Time-step")
    plt.savefig("runs/lasana/mnist_raster.png", dpi=300)
    plt.savefig("runs/lasana/mnist_raster.pdf")

    print("***")
    mean_energy = perf_df['total_energy'].mean()
    print(f"Mean Total Indiveri energy: {mean_energy} J (100 %)")
    print(f"Mean Soma Indiveri energy: {perf_df['soma_energy'].mean()} J ({100.0 * perf_df['soma_energy'].mean() / mean_energy} %)")
    print(f"Mean Synapse Indiveri energy: {perf_df['synapse_energy'].mean()} J ({100.0 * perf_df['synapse_energy'].mean() / mean_energy} %)")
    print(f"Mean Network Indiveri energy: {perf_df['network_energy'].mean()} J ({100.0 * perf_df['network_energy'].mean() / mean_energy} %)")
    print(f"Mean Indiveri firing neurons: {perf_df['fired'].mean()}")

    print("***")
    perf_filename = f"perf_loihi_mnist.csv"
    perf_df = pd.read_csv(os.path.join(RUN_PATH, perf_filename))
    mean_energy = perf_df['total_energy'].mean()
    print(f"Mean Total Loihi energy: {perf_df['total_energy'].mean()} J (100 %)")
    print(f"Mean Soma Loihi energy: {perf_df['soma_energy'].mean()} J ({100.0 * perf_df['soma_energy'].mean() / mean_energy} %)")
    print(f"Mean Synapse Loihi energy: {perf_df['synapse_energy'].mean()} J ({100.0 * perf_df['synapse_energy'].mean() / mean_energy} %)")
    print(f"Mean Network Loihi energy: {perf_df['network_energy'].mean()} J ({100.0 * perf_df['network_energy'].mean() / mean_energy} %)")
    print(f"Mean Loihi firing neurons: {perf_df['fired'].mean()}")


# Rate-encoded MNIST as a toy demo
#dataset = "mnist"
# Spiking Heidelberg Digits (SHD) is the second neuromorphic application.
#  I have trained two SNNs for SHD for this script:
#  The first SNN was trained using a modified LIF (Leaky) behavioral model (that
#  captured the effects of positive and negative clipping) and then converting
#  to the analog Indiveri neurons. The abstract SNN was trained to ~70% accuracy
#  but after conversion and deploying onto analog neurons, accuracy drops to
#  under ~50%.
#
#  The second SNN was trained using the LASANA model directly (by Jason). This
#   circuit-aware training has no accuracy degredation and so achieves ~70%
#   accuracy, which is comparable to the SHD benchmark paper.
dataset = "shd"

run_experiments = False
create_plots = True

if dataset == "mnist":
    num_inputs = 10
    #num_inputs = 1
    #num_inputs = 10000 # Entire test set
elif dataset == "shd":
    timesteps = None
    #num_inputs = 1
    #num_inputs = 10
    #num_inputs = 100
    num_inputs = 2264  # Entire test set

if run_experiments:
    run_experiment(num_inputs, dataset, analog_neurons=True)
    run_experiment(num_inputs, dataset, analog_neurons=False)

if create_plots:
    plot_experiments(dataset)

print("Finished.")
#plt.show()  # Enable if running on the host and we want to display