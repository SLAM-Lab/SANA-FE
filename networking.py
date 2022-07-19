"""First attempt at simulation run script.

Run a few basic experiments, show how we might interface a python
script with the simulator kernel.
"""
import matplotlib
matplotlib.use('Agg')

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt

MAX_COMPARTMENTS = 1024
MAX_CORES = 128
NETWORK_FILENAME = "connected_layer.csv"
TECH_FILENAME = "loihi.tech"

def run_sim(network, core_count):
    fields = ["Neuron ID", "Core ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
    with open(NETWORK_FILENAME, "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(fields)
        writer.writerows(network)

    timesteps = 10
    command = ("./sim", "{0}".format(TECH_FILENAME), "{0}".format(timesteps),
               NETWORK_FILENAME)
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    with open("results.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)

    return results

import random
def fully_connected(layer_neurons, spiking=True, max_connections=100):
    # Two layers, fully connected
    network = []

    reset = 0
    # Set up input into network
    for i in range(0, layer_neurons):
        # TODO: this is a bug in the csv handling code (in network.c)
        #  I can't seem to have fields that have 0 characters in, they all get
        #  treated as one field
        # neuron = ['i', None, None, None, None, None]
        neuron = ['i', ' ', ' ', ' ', ' ', ' ', i, 2.0]
        #for dest in range(i*4096, min((i+1)*4096, layer_neurons)):
        #    neuron.extend((i, random.random()))  # Same weight for all connections
        network.append(neuron)

    for n in range(0, layer_neurons):
        core_id = n / MAX_COMPARTMENTS

        target_core = random.choice(list(range(64,128)))
        min_neuron = target_core*1024
        neuron_list = list(range(min_neuron,min_neuron+1024))
        #neuron_list = list(range(layer_neurons, 2*layer_neurons))
        dest_neurons = random.sample(neuron_list, max_connections)
        #print(dest_neurons)

        neuron = [n, core_id, 1.0, reset, 0, 0]
        for dest in dest_neurons:
            neuron.extend((dest, random.random()))  # Same weight for all connections
        network.append(neuron)
        #print("neuron has {0} fields".format(len(neuron)))

        if (n % 10000) == 0:
            print(n)

    if spiking:
        threshold = 3.0
    else:  # never spike
        threshold = 2*layer_neurons

    for n in range(layer_neurons, 2*layer_neurons):
        core_id = n / MAX_COMPARTMENTS
        neuron = [n, core_id, threshold, reset, 0, 0]
        network.append(neuron)

    core_count = core_id + 1

    #print(network)

    return network, core_count


def empty(neurons, max_compartments=MAX_COMPARTMENTS):
    network = []
    threshold = -1.0  # never spike
    reset = 0.0

    core_id = 0
    compartment = 0
    # Map a number of neurons onto cores, but we may not necessarily use all
    #  compartments of each core
    # TODO: how to do 0 compartments yet still activate the core?
    # Maybe the simulator can take number of cores to simulate as an arg
    #  then we assume all simulated cores are powered
    for n in range(0, neurons):
        if compartment == max_compartments:
            core_id += 1
            compartment = 0

        neuron = [n, core_id, threshold, reset, 0, 0]
        network.append(neuron)
        compartment += 1

    core_count = core_id + 1

    return network, core_count


if __name__ == "__main__":
    neurons = []
    spiking_times = []
    spiking_update_energy = []
    spiking_spike_gen_energy = []
    spiking_synapse_energy = []
    spiking_network_energy = []

    layer_neurons = 65536
    #layer_neurons = 4096

    network, core_count = fully_connected(layer_neurons, True)
    print("Testing network with {0} neurons and inputs".format(len(network)))
    results = run_sim(network, core_count)

    neurons.append(len(network))
    spiking_times.append(results["time"])
    spiking_update_energy.append(results["update_energy"])
    spiking_spike_gen_energy.append(results["spike_gen_energy"])
    spiking_synapse_energy.append(results["synapse_energy"])
    spiking_network_energy.append(results["network_energy"])

    # Write all the simulation data to csv
    with open("sim_network.csv", "w") as spiking_csv:
        spiking_writer = csv.DictWriter(spiking_csv,
                                        ("neurons", "update_energy", "spike_gen_energy",
                                         "synapse_energy", "network_energy"))
        spiking_writer.writeheader()
        for neuron_count, update_energy, spike_gen_energy, synapse_energy, network_energy in zip(neurons, spiking_update_energy,
                                                  spiking_spike_gen_energy, spiking_synapse_energy,
                                                  spiking_network_energy):
            spiking_writer.writerow({"neurons": neuron_count,
                                     "update_energy": update_energy,
                                     "spike_gen_energy": spike_gen_energy,
                                     "synapse_energy": synapse_energy,
                                     "network_energy": network_energy})

    plt.figure(figsize=(5.5, 5.5))
    plt.bar(1, spiking_update_energy)
    plt.bar(2, spiking_spike_gen_energy, color="orange")
    plt.bar(3, spiking_synapse_energy, color="red")
    plt.bar(4, spiking_network_energy, color="green")
    plt.xticks([1,2,3,4], ["Update", "Spike Gen", "Synapse", "Network"])
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_energy.png")
    #plt.show()
