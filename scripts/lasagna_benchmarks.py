import torch
import os
import sys
import dill
import numpy as np
import matplotlib.pyplot as plt
import csv

# Try importing the installed sanafe library. If not installed, require a
#  fall-back to a local build of the sanafe Python library
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
    sys.path.insert(0, PROJECT_DIR)
    import sanafe

# Setup constants and read in the spiking networks from file
timesteps = 100
# Load models previously trained on Google Collab
mnist_model = torch.load(
        os.path.join(PROJECT_DIR, "etc", "mnist.pt"),
        pickle_module=dill,
        map_location=torch.device("cpu"))
spiking_digits_model = torch.load(
        os.path.join(PROJECT_DIR, "etc", "spiking_digits.pt"),
        pickle_module=dill,
        map_location=torch.device("cpu"))
mnist_weights = {}
for param_name, param in mnist_model.named_parameters():
    mnist_weights[param_name] = param.detach().numpy()

mnist_inputs = np.loadtxt(os.path.join(PROJECT_DIR, "etc", "mnist.csv"),
                          delimiter=',')
#spiking_digits_inputs = np.loadtxt(os.path.join(PROJECT_DIR, "etc", "spiking_digits.npz"))

# PyTorch stores weights in an array with dims (num out x num in)
in_neurons = mnist_weights["fc1.weight"].shape[1]
hidden_neurons = mnist_weights["fc1.weight"].shape[0]
out_neurons = mnist_weights["fc2.weight"].shape[0]

#"""
# Load the LASAGNA architecture with analog neurons
arch = sanafe.load_arch("/home/james/code/lasagne/lasagne.yaml")

# Load the network into SANA-FE
network = sanafe.Network()
# TODO: mapping order still seems messed up..
print(f"in:{in_neurons} hidden:{hidden_neurons} out:{out_neurons}")
in_layer = network.create_neuron_group("in", in_neurons)
hidden_layer = network.create_neuron_group("hidden", hidden_neurons,
                                           soma_hw_name="soma")
out_layer = network.create_neuron_group("out", out_neurons, soma_hw_name="soma")

# Select only the first MNIST digit for now
mnist_input = mnist_inputs[0, :]
#plt.imshow(mnist_input.reshape(28, 28), cmap="gray")
#plt.colorbar()
#plt.show()

# Map neurons to one core for now
for id, neuron in enumerate(in_layer.neurons):
    # Setup input neurons
    neuron.set_attributes(soma_hw_name=f"input[{id}]",
                          model_parameters={"poisson": mnist_input[id]},
                          log_spikes=True)
    neuron.map_to_core(arch.tiles[0].cores[0])

for id, neuron in enumerate(hidden_layer.neurons):
    neuron.set_attributes(soma_hw_name=f"loihi", log_spikes=True,
                          model_parameters={"threshold": 1.0, "leak_decay": 0.85})
    #neuron.set_attributes(soma_hw_name=f"analog_lif[{id}]", log_spikes=True)
    neuron.map_to_core(arch.tiles[0].cores[0])

for id, neuron in enumerate(out_layer.neurons):
    #neuron.set_attributes(soma_hw_name=f"analog_lif[{id + hidden_neurons}]", log_spikes=True, log_potential=True)
    neuron.set_attributes(soma_hw_name=f"loihi", log_spikes=True,
                          model_parameters={"threshold": 1.0, "leak_decay": 0.85})
    neuron.map_to_core(arch.tiles[0].cores[0])


# Connect neurons in both layers
for src in range(0, in_neurons):
    for dst in range(0, hidden_neurons):
        weight = mnist_weights["fc1.weight"][dst, src]
        if weight > 0.0:
            network.groups["in"].neurons[src].connect_to_neuron(
                network.groups["hidden"].neurons[dst], {"weight": weight})

for src in range(0, hidden_neurons):
    for dst in range(0, out_neurons):
        weight = mnist_weights["fc2.weight"][dst, src]
        if weight > 0.0:
            print(weight)
            network.groups["hidden"].neurons[src].connect_to_neuron(
                network.groups["out"].neurons[dst], {"weight": weight})

# Run a simulation
hw = sanafe.SpikingHardware(arch, record_spikes=True, record_potentials=True)
hw.load(network)
results = hw.sim(timesteps)
print(results)
#"""


# Plot simulation results
with open("spikes.csv") as spike_csv:
    spike_data = csv.DictReader(spike_csv)
    plt.xlim((0, timesteps))
    plt.ylim((0, in_neurons + hidden_neurons + out_neurons))

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

# 1) Raster plot of spikes
total_neurons = 0
for neuron_id in range(0, in_neurons):
    plt.scatter(in_spikes[neuron_id],
                [total_neurons]*len(in_spikes[neuron_id]), c='r', s=2,
                marker='.', linewidths=0.5)
    total_neurons += 1

for neuron_id in range(0, hidden_neurons):
    plt.scatter(hidden_spikes[neuron_id],
                [total_neurons]*len(hidden_spikes[neuron_id]), c='b', s=2,
                marker='.', linewidths=0.5)
    total_neurons += 1

for neuron_id in range(0, out_neurons):
    plt.scatter(out_spikes[neuron_id],
                [total_neurons]*len(out_spikes[neuron_id]), c='k', s=2,
                marker='.', linewidths=0.5)
    total_neurons += 1

#print("timesteps: {0}".format(timesteps))
plt.xlabel("Time-step")
plt.ylabel("Neuron")
plt.show()

# TODO: load the spike data and plot raster
# TODO: break raster plot into script that can ship with sana-fe python library

print("Finished")
