import csv
import os
import sys
import numpy as np
import pandas as pd
import dill
# Plotting
import matplotlib.pyplot as plt

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

analog_synapses = True
#analog_synapses = False

timesteps = 30
num_inputs = 100
#num_inputs = 1

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

print(weights)
print(decays)
print(thresholds)

# PyTorch stores weights in an array with dims (num out x num in)
in_neurons = weights["fc1"].shape[1]
hidden_1_neurons = weights["fc1"].shape[0]
hidden_2_neurons = weights["fc2"].shape[0]
out_neurons = weights["fc3"].shape[0]

print(f"in:{in_neurons}, hidden 1:{hidden_1_neurons}, hidden 2:{hidden_2_neurons}, out:{out_neurons}")
#"""
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



# MNIST Data-set info (circuit-aware trained 23 Sep 2025 for 10 epochs)
# Train set: Average loss: 0.1868, Accuracy: 56680/60000 (94.47%)
# Test set: Accuracy: 9475/10000 (94.75%)

# TODO: run both with and without analog synapses
# TODO: plot and compare performance of both crossbar/synaptic memory
