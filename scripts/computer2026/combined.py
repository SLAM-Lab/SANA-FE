import csv
import os
import sys
import numpy as np
import pandas as pd
import dill
import argparse
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

# from indiveri_model import Indiveri

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))

parser = argparse.ArgumentParser(
        prog="Loihi-Crossbar-Indiveri",
        description="Explore mixed-signal architecture with analog crossbar synapses and Indiveri neurons.")
parser.add_argument("data_path", nargs="?",
                    default=os.path.abspath((os.path.join(PROJECT_DIR, "runs", "combined"))))
parser.add_argument("run_path", nargs="?",
                    default=os.path.abspath((os.path.join(PROJECT_DIR, "runs", "combined"))))
parser.add_argument("lasana_dir", nargs="?",
                    default=os.path.abspath(os.path.join("/", "home", "usr1", "jboyle", "neuro", "lasana", "build")))
parser.add_argument("--quick", action="store_true")
parser.add_argument("--run", action="store_true")
parser.add_argument("--plot", action="store_true")
args = parser.parse_args()

DATA_PATH = args.data_path
RUN_PATH = args.run_path
LASANA_DIR = args.lasana_dir
QUICK_RUN = args.quick
RUN_EXPERIMENTS = args.run
PLOT_EXPERIMENTS = args.plot


# Try importing the installed sanafe library. If not installed, require a
#  fall-back to a local build of the sanafe Python library
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
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
        save_to=os.path.join(RUN_PATH, "data"),
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

    # snn_path = os.path.join(RUN_PATH, "app_models", "spiking_digits_indiveri_crossbar_latest.pt")
    # snn_path = os.path.join(RUN_PATH, "app_models", "shd_70_256R_20_bin_crossbar_indiveri_aware_60.pt.state")
    # snn_path = os.path.join(RUN_PATH, "app_models", "shd_70_256R_20_bin_crossbar_indiveri_aware_55.pt.state")
    snn_path = os.path.join(DATA_PATH, "app_models", "shd_70_256R_20_crossbar_indiveri_aware.pt")

    # After loading
    # snn = torch.load(snn_path,
    #                  pickle_module=dill, map_location=torch.device("cpu"))
    # torch.save(snn.state_dict(), "shd_70_256R_20_bin_crossbar_indiveri_aware_60.pt.state")
    # print("Saved!")
    # exit()

    snn = torch.load(snn_path, weights_only=True,
                     map_location=torch.device("cpu"))

    weights = {}
    weights["fc1"] = ternary_weight(snn["fc1.weight"])
    weights["fc2"] = ternary_weight(snn["fc2.weight"])
    weights["fcr"] = ternary_weight(snn["fcr.weight"])

    crossbar_biases = {}
    crossbar_biases[1] = ternary_weight(snn["fc1.bias"])
    crossbar_biases[2] = ternary_weight(snn["fc2.bias"])
    crossbar_biases[3] = ternary_weight(snn["fcr.bias"])

    # print(weights)
    # print(crossbar_biases)
    # exit()

    # PyTorch stores weights in an array with dims (num out x num in)
    in_neurons = weights["fc1"].shape[1]
    hidden_neurons = weights["fc1"].shape[0]
    out_neurons = weights["fc2"].shape[0]

    print(f"in:{in_neurons}, hidden 1:{hidden_neurons}, our:{out_neurons}")
    #"""
    # Load the LASANA architecture with analog neurons
    arch = sanafe.load_arch(os.path.join(LASANA_DIR, "crossbar", "lasana_crossbar_indiveri.yaml"))

    spike_filename = f"spikes_crossbar_indiveri_shd.csv"
    perf_filename = f"perf_crossbar_indiveri_shd.csv"
    potential_filename = f"potential_crossbar_indiveri_shd.csv"


    print("Creating network in SANA-FE")
    network = sanafe.Network()
    in_layer = network.create_neuron_group("in", in_neurons)
    hidden_layer_1 = network.create_neuron_group("hidden", hidden_neurons)
    out_layer = network.create_neuron_group("out", out_neurons)

    synapse_hw = "loihi_dense_synapse"
    dendrite_hw = "loihi_dendrites"
    print("Creating input layer")
    for id, neuron in enumerate(in_layer):
        neuron.set_attributes(soma_hw_name=f"loihi_inputs[{id}]",
                              dendrite_hw_name=dendrite_hw,
                              log_spikes=True,
                              log_potential=True)
        neuron.map_to_core(arch.tiles[0].cores[0])

    print("Creating hidden layer")

    num_crossbar_rows_per_core = 256
    input_cores = 1
    total_rows_per_neuron = int(math.ceil(in_neurons / 32) + math.ceil(hidden_neurons / 32))
    hidden_neurons_per_core = num_crossbar_rows_per_core // total_rows_per_neuron
    hidden_cores = int(math.ceil(hidden_neurons / hidden_neurons_per_core))
    output_cores = 1

    for id, neuron in enumerate(hidden_layer_1):
        hidden_parameters = {
            "threshold": 1.0,
            "leak_decay": 0.85,
            "input_decay": 0.0,
            "reset_mode": "hard",
            "scale_input": 0.5,
            "input_offset": 1.0
        }
        synapse_hw = f"synapse_crossbar[{id}]"
        dendrite_hw = synapse_hw


        neuron.set_attributes(soma_hw_name=f"analog_lif[{id}]",
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
            "threshold": 1.0,
            "leak_decay": 0.85,
            "reset_mode": "hard",
            "input_offset": 1.0
        }
        synapse_hw = f"synapse_crossbar[{id}]"
        dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=f"analog_lif[{id + hidden_neurons}]",
                            dendrite_hw_name=dendrite_hw,
                            default_synapse_hw_name=synapse_hw,
                            log_spikes=True,
                            log_potential=True,
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
            # Offset from the previously mapped rows for the hidden layer
            row = total_hidden_1_crossbar_rows + (src // 32)
            col = src % 32
            connection_parameters["crossbar_position"] = [row, col]
            connection_parameters["weight"] = int(weights["fcr"][dst, src])
            connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            connection_parameters["crossbar_bias"] = int(crossbar_biases[3][dst])

            network["hidden"][src].connect_to_neuron(
                network["hidden"][dst], connection_parameters)

    for src in range(hidden_neurons):
        for dst in range(out_neurons):
            connection_parameters = {}
            row = src // 32
            col = src % 32
            connection_parameters["crossbar_position"] = [row, col]
            connection_parameters["weight"] = int(weights["fc2"][dst, src])
            connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            connection_parameters["crossbar_bias"] = int(crossbar_biases[2][dst])

            network["hidden"][src].connect_to_neuron(
                network["out"][dst], connection_parameters)

    network.save(os.path.join(RUN_PATH, "crossbar_indiveri_shd.yaml"))

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
        # Make sure we add an extra timestep at the end. This turned out to
        #  be really important for accuracy when aligning against SNNTorch
        #  trained models! Perhaps the margins are quite fine...
        timesteps = len(spiketrain) + 1

        print(f"Simulating for {timesteps} timesteps")
        is_first_input = (input_idx == 0)

        hw.sim(timesteps,
               spike_trace=os.path.join(RUN_PATH, spike_filename),
               perf_trace=os.path.join(RUN_PATH, perf_filename),
               potential_trace=os.path.join(RUN_PATH, potential_filename),
               write_trace_headers=is_first_input,
               processing_threads=8,
               scheduler_threads=4)
        timesteps_per_input.append(timesteps)
        hw.reset()
        #input()

    # The inputs (and therefore timesteps per input) will be the same across
    #  Loihi/analog neuron runs. Only store this once.
    np.savetxt(os.path.join(RUN_PATH, f"crossbar_indiveri_shd.csv"),
               np.array(timesteps_per_input), fmt="%d")
    print("Simulation finished")


def load_dataset(num_inputs):
    weights = {}

    # Load the Spiking Heidelberg Digits (SHD) network
    # Use a SNN trained using circuit-aware methods i.e. was
    #  trained specifically on an IMAC-sim crossbar model.
    spiking_digits_model = torch.load(
            os.path.join(DATA_PATH, "app_models", "shd_70_256R_20_crossbar_indiveri_aware.pt"),
            map_location=torch.device("cpu"), weights_only=True)

    for attribute_name, param in spiking_digits_model.items():
        weights[attribute_name] = param.detach().numpy()

    # Load the spiking digits test inputs from the raw dataset, applying the
    #  same transformations as the training scripts
    frame_transform = tonic.transforms.Compose([
        tonic.transforms.Downsample(time_factor=0.05, spatial_factor=0.1),  # For Jason's circuit-aware trained SNN
        tonic.transforms.ToFrame(sensor_size=(70, 1, 1), time_window=1000)
    ])
    testset = tonic.datasets.SHD(
        save_to=os.path.join(RUN_PATH, "data"),
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


def calculate_shd_accuracy(num_inputs):
    print(f"Calculating accuracy for SHD")
    timesteps_per_input = list(np.loadtxt(
        os.path.join(RUN_PATH, "crossbar_indiveri_shd.csv"), dtype=int, ndmin=1))

    _, labels, weights = load_dataset(num_inputs)
    in_neurons = weights["fc1.weight"].shape[1]
    hidden_neurons = weights["fc1.weight"].shape[0]
    out_neurons = weights["fc2.weight"].shape[0]

    platform = "crossbar_indiveri"
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

    # print(out_spikes)

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
    print(f"Accuracy ({platform}): {accuracy}% ({correct}/{num_inputs})")
    print(f"Min timesteps: {np.min(timesteps_per_input)}")
    print(f"Max timesteps: {np.max(timesteps_per_input)}")
    return accuracy


def plot_spiking_digits(num_inputs):
    print("Plotting spiking digit experiments")
    timesteps_per_input = list(np.loadtxt(
        os.path.join(RUN_PATH, "crossbar_indiveri_shd.csv"),
        dtype=int, ndmin=1))

    imac_indiveri_accuracy = calculate_shd_accuracy(num_inputs=num_shd_inputs)
    analog_perf_df = pd.read_csv(os.path.join(RUN_PATH, "perf_crossbar_indiveri_shd.csv"))
    analog_perf_df["total_energy_uj"] = analog_perf_df["total_energy"] * 1.0e6

    # Per time-step
    imac_mean_energy = analog_perf_df['total_energy'].mean()
    def per_timestep_row(df, label, mean_energy):
        return {
            "Platform": label,
            "Total Energy (uJ)":  df['total_energy'].mean()   * 1e6,
            "Total Energy (%)": 100.0,
            "Soma Energy (uJ)": df['soma_energy'].mean()    * 1e6,
            "Soma Energy (%)": 100.0 * df['soma_energy'].mean()    / mean_energy,
            "Synapse Energy (uJ)": df['synapse_energy'].mean() * 1e6,
            "Synapse Energy (%)": 100.0 * df['synapse_energy'].mean() / mean_energy,
            "Network Energy (uJ)": df['network_energy'].mean() * 1e6,
            "Network Energy (%)": 100.0 * df['network_energy'].mean() / mean_energy,
            "Mean Fired": df['fired'].mean(),
            "Mean Latency (ms)": df['sim_time'].mean() * 1e3,
        }

    per_timestep_results = pd.DataFrame([
        per_timestep_row(analog_perf_df, "Loihi-IMAC-Indiveri", imac_mean_energy),
    ]).set_index("Platform")

    # print("=" * 80)
    # print("Per Time-step Results")
    # print("=" * 80)
    # print(per_timestep_results.to_string())
    # print()

    # Per-inference
    imac_mean_energy = analog_perf_df['total_energy'].sum() / num_inputs

    def per_inference_row(df, label, mean_energy, accuracy):
        return {
            "Platform":            label,
            "Total Energy (uJ)":   df['total_energy'].sum()   / num_inputs * 1e6,
            "Total Energy (%)":    100.0,
            "Soma Energy (uJ)":    df['soma_energy'].sum()    / num_inputs * 1e6,
            "Soma Energy (%)":     100.0 * df['soma_energy'].sum()    / (num_inputs * mean_energy),
            "Synapse Energy (uJ)": df['synapse_energy'].sum() / num_inputs * 1e6,
            "Synapse Energy (%)":  100.0 * df['synapse_energy'].sum() / (num_inputs * mean_energy),
            "Network Energy (uJ)": df['network_energy'].sum() / num_inputs * 1e6,
            "Network Energy (%)":  100.0 * df['network_energy'].sum() / (num_inputs * mean_energy),
            "Firing Neurons":      df['fired'].sum()   / num_inputs,
            "Latency (ms)":             df['sim_time'].sum() / num_inputs * 1e3,
            "Accuracy (%)": accuracy
        }

    per_inference_results = pd.DataFrame([
        per_inference_row(analog_perf_df, "Loihi-IMAC-Indiveri", imac_mean_energy, imac_indiveri_accuracy),
    ]).set_index("Platform")

    print("=" * 80)
    print("SHD Per-Inference Results")
    print("=" * 80)
    print(per_inference_results.to_string())
    print()

    # Save to CSV
    combined = pd.concat(
        [per_inference_results, per_timestep_results],
        keys=["Per-Inference", "Per-Timestep"]
    )
    combined.to_csv(os.path.join(RUN_PATH, "imac_shd_results.csv"))
    print("Results saved to imac_shd_results.csv")

#"""

# MNIST Data-set info (circuit-aware trained 23 Sep 2025 for 10 epochs)
# Train set: Average loss: 0.1868, Accuracy: 56680/60000 (94.47%)
# Test set: Accuracy: 9475/10000 (94.75%)

def run_mnist(num_inputs, timesteps=100):
    print(f"Loading models for binarized MNIST")
    platform = "crossbar_indiveri"
    spike_filename = f"spikes_{platform}_mnist.csv"
    perf_filename = f"perf_{platform}_mnist.csv"
    potential_filename = f"potential_{platform}_mnist.csv"

    # Data loading
    transform = transforms.Compose([
        transforms.ToTensor(),
        transforms.Resize((20, 20)),
        transforms.Normalize((0,), (1,))  # MNIST normalization
    ])
    test_dataset = datasets.MNIST(os.path.join(RUN_PATH, "data"), train=False,
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

    snn_path = os.path.join(DATA_PATH, "app_models", "mnist_10_128_84_10_crossbar_indiveri_aware.pt")
    snn = torch.load(snn_path, map_location=torch.device("cpu"))

    weights = {}
    weights["fc1"] = ternary_weight(snn["fc1.weight"])
    weights["fc2"] = ternary_weight(snn["fc2.weight"])
    weights["fc3"] = ternary_weight(snn["fc_out.weight"])

    crossbar_biases = {}
    crossbar_biases[1] = ternary_weight(snn["fc1.bias"])
    crossbar_biases[2] = ternary_weight(snn["fc2.bias"])
    crossbar_biases[3] = ternary_weight(snn["fc_out.bias"])

    # PyTorch stores weights in an array with dims (num out x num in)
    in_neurons = weights["fc1"].shape[1]
    hidden_1_neurons = weights["fc1"].shape[0]
    hidden_2_neurons = weights["fc2"].shape[0]
    out_neurons = weights["fc3"].shape[0]

    print(f"in:{in_neurons}, hidden 1:{hidden_1_neurons}, hidden 2:{hidden_2_neurons}, out:{out_neurons}")
    # Load the LASANA architecture with analog neurons
    arch = sanafe.load_arch(os.path.join(LASANA_DIR, "crossbar", "lasana_crossbar_indiveri.yaml"))

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
            "reset_mode": "hard",
            "input_offset": 1.0,
            "digital_to_analog_converter": False,
            "analog_to_digital_converter": False
        }
        #hidden_parameters["crossbar_bias"] = int(crossbar_biases[1][id])
        synapse_hw = f"synapse_crossbar[{id}]"
        dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=f"analog_lif[{id}]",
                              dendrite_hw_name=dendrite_hw,
                              log_spikes=True,
                              default_synapse_hw_name=synapse_hw,
                              model_attributes=hidden_parameters)
        neuron.map_to_core(arch.tiles[0].cores[1])

    for id, neuron in enumerate(hidden_layer_2):
        hidden_parameters = {
            "reset_mode": "hard",
            "input_offset": 1.0,
            "digital_to_analog_converter": False,
            "analog_to_digital_converter": False
        }
        #hidden_parameters["crossbar_bias"] = int(crossbar_biases[2][id])
        synapse_hw = f"synapse_crossbar[{id}]"
        dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=f"analog_lif[{id}]",
                              dendrite_hw_name=dendrite_hw,
                              log_spikes=True,
                              default_synapse_hw_name=synapse_hw,
                              model_attributes=hidden_parameters)
        neuron.map_to_core(arch.tiles[0].cores[2])

    print("Creating output layer")
    for id, neuron in enumerate(out_layer):
        hidden_parameters = {
            "reset_mode": "hard",
            "input_offset": 1.0,
            "digital_to_analog_converter": False,
            "analog_to_digital_converter": False
        }
        #hidden_parameters["crossbar_bias"] = int(crossbar_biases[3][id])
        synapse_hw = f"synapse_crossbar[{id}]"
        dendrite_hw = synapse_hw

        neuron.set_attributes(soma_hw_name=f"analog_lif[{id}]",
                            dendrite_hw_name=dendrite_hw,
                            default_synapse_hw_name=synapse_hw,
                            log_spikes=True,
                            model_attributes=hidden_parameters)
        neuron.map_to_core(arch.tiles[0].cores[3])

    # Connect neurons in both layers

    print("Connecting neurons")
    print("Adding first layer connections")
    total_hidden_1_crossbar_rows = (in_neurons + 31) // 32
    total_hidden_2_crossbar_rows = (hidden_1_neurons + 31) // 32
    total_out_crossbar_rows = (out_neurons + 31) // 32

    print("Total crossbar rows:")
    print(f"In->Hidden: {total_hidden_1_crossbar_rows}")
    print(f"Hidden->Hidden: {total_hidden_2_crossbar_rows}")
    print(f"Hidden->Out: {total_out_crossbar_rows}")


    for src in range(in_neurons):
        for dst in range(hidden_1_neurons):
            connection_parameters = {}
            row = src // 32
            col = src % 32
            connection_parameters["crossbar_position"] = [row, col]
            connection_parameters["weight"] = int(weights["fc1"][dst, src])
            connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            connection_parameters["crossbar_bias"] = int(crossbar_biases[1][dst])

            network["in"][src].connect_to_neuron(
                network["hidden_1"][dst], connection_parameters)

    for src in range(hidden_1_neurons):
        for dst in range(hidden_2_neurons):
            connection_parameters = {}
            row = src // 32
            col = src % 32
            connection_parameters["crossbar_position"] = [row, col]
            connection_parameters["weight"] = int(weights["fc2"][dst, src])
            connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            connection_parameters["crossbar_bias"] = int(crossbar_biases[2][dst])

            network["hidden_1"][src].connect_to_neuron(
                network["hidden_2"][dst], connection_parameters)

    for src in range(hidden_2_neurons):
        for dst in range(out_neurons):
            connection_parameters = {}
            row = src // 32
            col = src % 32
            connection_parameters["crossbar_position"] = [row, col]
            connection_parameters["weight"] = int(weights["fc3"][dst, src])
            connection_parameters["synapse_hw_name"] = f"analog_crossbar[{dst}]"
            connection_parameters["crossbar_bias"] = int(crossbar_biases[3][dst])

            network["hidden_2"][src].connect_to_neuron(
                network["out"][dst], connection_parameters)

    network.save(os.path.join(RUN_PATH, "crossbar_mnist.yaml"))

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
                        spike_trace=os.path.join(RUN_PATH, spike_filename),
                        perf_trace=os.path.join(RUN_PATH, perf_filename),
                        potential_trace=os.path.join(RUN_PATH, potential_filename),
                        write_trace_headers=is_first_input,
                        processing_threads=8,
                        scheduler_threads=4)
        timesteps_per_input.append(timesteps)
        # print(results)
        hw.reset()


def calculate_mnist_accuracy(num_inputs, analog_synapses=True):
    print("Calculating MNIST accuracy for "
          f"{'Loihi-IMAC' if analog_synapses else 'Loihi'}")
    filename = "spikes_crossbar_indiveri_mnist.csv"

    # Data loading
    transform = transforms.Compose([
        transforms.ToTensor(),
        transforms.Resize((20, 20)),
        transforms.Normalize((0,), (1,))  # MNIST normalization
    ])
    test_dataset = datasets.MNIST(os.path.join(RUN_PATH, "data"), train=False,
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

    with open(os.path.join(RUN_PATH, filename)) as spike_csv:
        spike_data = csv.DictReader(spike_csv)

        # PyTorch stores weights in an array with dims (num out x num in)
        in_neurons = 400
        hidden_1_neurons = 128
        hidden_2_neurons = 84
        out_neurons = 10
        timesteps = 100

        in_spikes = [[] for _ in range(0, in_neurons)]
        hidden_spikes_1 = [[] for _ in range(0, hidden_1_neurons)]
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
        # print(f"Spike counts per class for inference of digit:{counts[:, i]} "
        #     f"out:{np.argmax(counts[:, i])} actual:{labels[i]}")
        if np.argmax(counts[:, i]) == labels[i]:
            correct += 1

    accuracy = (correct / num_inputs) * 100
    print(f"Loihi-IMAC-Indiveri Accuracy: {accuracy}%")
    return accuracy


def plot_mnist(num_inputs):
    print("Outputting MNIST experiment results")
    analog_perf_df = pd.read_csv(
        os.path.join(RUN_PATH, "perf_crossbar_indiveri_mnist.csv"))

    imac_indiveri_accuracy = calculate_mnist_accuracy(
        num_inputs=num_mnist_inputs, analog_synapses=True)

    # Per time-step
    imac_mean_energy = analog_perf_df['total_energy'].mean()

    def per_timestep_row(df, label, mean_energy):
        return {
            "Platform": label,
            "Total Energy (uJ)": df['total_energy'].mean() * 1e6,
            "Total Energy (%)": 100.0,
            "Soma Energy (uJ)": df['soma_energy'].mean() * 1e6,
            "Soma Energy (%)": 100.0 * df['soma_energy'].mean() / mean_energy,
            "Synapse Energy (uJ)": df['synapse_energy'].mean() * 1e6,
            "Synapse Energy (%)": 100.0 * df['synapse_energy'].mean() / mean_energy,
            "Network Energy (uJ)": df['network_energy'].mean() * 1e6,
            "Network Energy (%)": 100.0 * df['network_energy'].mean() / mean_energy,
            "Mean Fired": df['fired'].mean(),
            "Total Spikes": df['spikes'].sum(),
        }

    per_timestep_results = pd.DataFrame([
        per_timestep_row(analog_perf_df, "Loihi-IMAC-Indiveri", imac_mean_energy),
    ]).set_index("Platform")

    # print("=" * 80)
    # print("Per Time-step Results")
    # print("=" * 80)
    # print(per_timestep_results.to_string())
    # print()

    # Per-inference
    imac_mean_energy = analog_perf_df['total_energy'].sum() / num_inputs

    def per_inference_row(df, label, mean_energy, accuracy):
        return {
            "Platform": label,
            "Total Energy (uJ)": df['total_energy'].sum() / num_inputs * 1e6,
            "Total Energy (%)": 100.0,
            "Soma Energy (uJ)": df['soma_energy'].sum() / num_inputs * 1e6,
            "Soma Energy (%)": 100.0 * df['soma_energy'].sum() / (num_inputs * mean_energy),
            "Synapse Energy (uJ)": df['synapse_energy'].sum() / num_inputs * 1e6,
            "Synapse Energy (%)":  100.0 * df['synapse_energy'].sum() / (num_inputs * mean_energy),
            "Network Energy (uJ)": df['network_energy'].sum() / num_inputs * 1e6,
            "Network Energy (%)": 100.0 * df['network_energy'].sum() / (num_inputs * mean_energy),
            "Firing Neurons": df['fired'].sum()  / num_inputs,
            "Latency (ms)": df['sim_time'].sum() / num_inputs * 1e3,
            "Accuracy (%)": accuracy
        }

    per_inference_results = pd.DataFrame([
        per_inference_row(analog_perf_df, "Loihi-IMAC-Indiveri", imac_mean_energy, imac_indiveri_accuracy),
    ]).set_index("Platform")

    print("=" * 80)
    print("MNIST Per-Inference Results")
    print("=" * 80)
    print(per_inference_results.to_string())
    print()

    # Save to CSV
    combined = pd.concat(
        [per_inference_results, per_timestep_results],
        keys=["Per-Inference", "Per-Timestep"]
    )
    combined.to_csv(os.path.join(RUN_PATH, "imac_indiveri_mnist_results.csv"))
    print("Results saved to imac_indiveri_mnist_results.csv")


if __name__ == "__main__":
    print(f"Launching Loihi-IMAC-Indiveri, run:{RUN_EXPERIMENTS} plot:{PLOT_EXPERIMENTS}")
    if QUICK_RUN:
        num_shd_inputs = 100
        num_mnist_inputs = 100
    else:
        num_shd_inputs = 2264  # Number of inferences
        num_mnist_inputs = 10000

    if RUN_EXPERIMENTS:
        run_spiking_digits(num_inputs=num_shd_inputs)
        run_mnist(num_inputs=num_mnist_inputs, timesteps=100)

    if PLOT_EXPERIMENTS:
        plot_spiking_digits(num_inputs=num_shd_inputs)
        plot_mnist(num_inputs=num_mnist_inputs)

    print("Loihi-IMAC-Indiveri finished!")