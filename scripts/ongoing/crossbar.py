import csv
import os
import sys
import numpy as np
import pandas as pd
# Plotting
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import torch
from torchvision import datasets, transforms
from torch.utils.data import DataLoader


def ternary_weight(weight):
    threshold = 0.7 * torch.mean(torch.abs(weight))
    output = torch.zeros_like(weight)
    output[weight > threshold] = 1.0
    output[weight < -threshold] = -1.0
    return output.numpy()


# Try importing the installed sanafe library. If not installed, require a
#  fall-back to a local build of the sanafe Python library
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
    print(f"Project dir: {PROJECT_DIR}")
    sys.path.insert(0, PROJECT_DIR)
    import sanafe

analog_synpases = True

#timesteps = 100
timesteps = 30
num_inputs = 100

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

# TODO: hacky just to debug
# import pandas as pd
# df = pd.read_csv("/home/usr1/jboyle/neuro/binary_mnist/inputs.csv")
# inputs = df.to_numpy()
# print(inputs)

snn_path = "/home/usr1/jboyle/neuro/binary_mnist/binarized_mnist_model.pth"
snn = torch.load(snn_path, weights_only=True)
weights = {}
weights["fc1"] = ternary_weight(snn["fc1.weight"])
weights["fc2"] = ternary_weight(snn["fc2.weight"])
weights["fc3"] = ternary_weight(snn["fc_out.weight"])

biases = {}
biases[1] = ternary_weight(snn["fc1.bias"])
biases[2] = ternary_weight(snn["fc2.bias"])
biases[3] = ternary_weight(snn["fc_out.bias"])

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
#"""
# Load the LASANA architecture with analog neurons
arch = sanafe.load_arch("/home/usr1/jboyle/neuro/lasana/lasana_crossbar.yaml")

print("Creating network in SANA-FE")
network = sanafe.Network()
in_layer = network.create_neuron_group("in", in_neurons)
hidden_layer_1 = network.create_neuron_group("hidden_1", hidden_1_neurons)
hidden_layer_2 = network.create_neuron_group("hidden_2", hidden_2_neurons)
out_layer = network.create_neuron_group("out", out_neurons)

print("Creating input layer")
for id, neuron in enumerate(in_layer):
    neuron.set_attributes(dendrite_hw_name="default_dendrite",
                          soma_hw_name=f"input[{id}]", log_spikes=True)
    neuron.map_to_core(arch.tiles[0].cores[0])

print("Creating hidden layer")
for id, neuron in enumerate(hidden_layer_1):
    hidden_parameters = {
        "threshold": thresholds[1],
        "leak_decay": decays[1],
        "reset_mode": "hard",
        "crossbar_bias": int(biases[1][id])
    }

    neuron.set_attributes(soma_hw_name="loihi_lif",
                          dendrite_hw_name=f"synapse_crossbar[{id}]",
                          default_synapse_hw_name=f"synapse_crossbar[{id}]",
                          log_spikes=True,
                          log_potential=True,
                          model_attributes=hidden_parameters)
    neuron.map_to_core(arch.tiles[0].cores[1])

for id, neuron in enumerate(hidden_layer_2):
    hidden_parameters = {
        "threshold": thresholds[2],
        "leak_decay": decays[2],
        "reset_mode": "hard",
        "crossbar_bias": int(biases[2][id])
    }
    neuron.set_attributes(soma_hw_name=f"loihi_lif",
                          dendrite_hw_name=f"synapse_crossbar[{id}]",
                          default_synapse_hw_name=f"synapse_crossbar[{id}]",
                          log_spikes=True, log_potential=True,
                          model_attributes=hidden_parameters)
    neuron.map_to_core(arch.tiles[0].cores[2])

print("Creating output layer")
for id, neuron in enumerate(out_layer):
    neuron.set_attributes(soma_hw_name=f"loihi_lif",
                          dendrite_hw_name=f"synapse_crossbar[{id}]",
                          default_synapse_hw_name=f"synapse_crossbar[{id}]",
                          log_spikes=True,
                          model_attributes={
                            "threshold": thresholds[3],
                            "leak_decay": decays[3],
                            "reset_mode": "hard",
                            "crossbar_bias": int(biases[3][id])
                        })
    neuron.map_to_core(arch.tiles[0].cores[3])

# Connect neurons in both layers

print("Connecting neurons")
print("Adding first layer connections")

for src in range(in_neurons):
    for dst in range(hidden_1_neurons):
        weight = int(weights["fc1"][dst, src])
        network["in"][src].connect_to_neuron(
            network["hidden_1"][dst], {"weight": weight, "synapse_hw_name": f"analog_crossbar[{dst}]"})

for src in range(hidden_1_neurons):
    for dst in range(hidden_2_neurons):
        weight = int(weights["fc2"][dst, src])
        network["hidden_1"][src].connect_to_neuron(
            network["hidden_2"][dst], {"weight": weight, "synapse_hw_name": f"analog_crossbar[{dst}]"})

for src in range(hidden_2_neurons):
    for dst in range(out_neurons):
        weight = int(weights["fc3"][dst, src])
        network["hidden_2"][src].connect_to_neuron(
            network["out"][dst], {"weight": weight, "synapse_hw_name": f"analog_crossbar[{dst}]"})

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
            model_attributes={"poisson": mnist_input[id]})

    print(f"Simulating for {timesteps} timesteps")
    is_first_input = (input_idx == 0)
    results = hw.sim(timesteps, spike_trace="spikes.csv",
                     potential_trace="potentials.csv", perf_trace="perf.csv",
                     write_trace_headers=is_first_input, processing_threads=1,
                     scheduler_threads=8)
    timesteps_per_input.append(timesteps)
    # print(results)
    hw.reset()
#"""

with open("spikes.csv") as spike_csv:
    spike_data = csv.DictReader(spike_csv)

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
    print(f"Spike counts per class for inference of digit:{counts[:, i]} "
        f"out:{np.argmax(counts[:, i])} actual:{labels[i]}")
    if np.argmax(counts[:, i]) == labels[i]:
        correct += 1

accuracy = (correct / num_inputs) * 100
print(f"Accuracy: {accuracy}%")

# def calculate_cumulative_spikes(out_spikes, out_neurons, total_timesteps, timesteps_per_input):
#     # Initialize array to store cumulative counts for each output neuron at each timestep
#     cumulative_counts = np.zeros((out_neurons, total_timesteps))

#     # For each output neuron
#     for neuron_idx, spikes in enumerate(out_spikes):
#         # Sort spikes by timestep for cumulative counting
#         sorted_spikes = sorted([s for s in spikes if s <= total_timesteps])
#         # Process each timestep
#         for t in range(total_timesteps):
#             # Calculate which input period we're in
#             input_image = t // timesteps_per_input
#             period_start = input_image * timesteps_per_input

#             # Count spikes from start of current input period up to current timestep
#             count = sum(1 for spike in sorted_spikes
#                     if spike <= t + 1 and spike > period_start)

#             cumulative_counts[neuron_idx, t] = count

#     return cumulative_counts


# cumulative_counts = calculate_cumulative_spikes(
#     out_spikes, out_neurons, timesteps * num_inputs, timesteps)
