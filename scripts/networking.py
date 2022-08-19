"""First attempt at simulation run script.

Run a few basic experiments, show how we might interface a python
script with the simulator kernel.
"""
#import matplotlib
#matplotlib.use('Agg')

import csv
import subprocess
import yaml
import pandas as pd
from matplotlib import pyplot as plt

MAX_COMPARTMENTS = 1024
MAX_CORES = 128
ARCH_FILENAME = "loihi.list"
NETWORK_FILENAME = "connected_layer.csv"
TECH_FILENAME = "loihi.tech"

def run_sim(network):
    # Create the network to run
    fields = ["Neuron ID", "Core ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
    with open(NETWORK_FILENAME, "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(fields)
        writer.writerows(network)

    timesteps = 10
    command = ("./sim", "{0}".format(TECH_FILENAME), ARCH_FILENAME,
               NETWORK_FILENAME, "{0}".format(timesteps),)
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    # Parse the detailed perf statistics
    print("Reading performance data")
    stats = pd.read_csv("perf.csv")
    analysis = parse_stats(stats)

    return analysis


def parse_stats(stats):
    print("Parsing statistics")
    total = stats.sum()
    print(total)
    update_keys, synapse_keys, spike_gen_energy, network_energy = [], [], [], []
    for k in total.keys():
        if "update_energy" in k: update_keys.append(k)
        elif "synapse_energy" in k: synapse_keys.append(k)
        elif "o[" in k: spike_gen_energy.append(k)
        elif "r[" in k: network_energy.append(k)

    analysis = {} 
    analysis["update_energy"] = total[update_keys].sum()
    analysis["synapse_energy"] = total[synapse_keys].sum()
    analysis["spike_gen_energy"] = total[spike_gen_energy].sum()
    analysis["network_energy"] = total[network_energy].sum()

    return analysis


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
        neuron = ['i', ' ', ' ', ' ', ' ', ' ', ' ', i, 2.0]
        #for dest in range(i*4096, min((i+1)*4096, layer_neurons)):
        #    neuron.extend((i, random.random()))  # Same weight for all connections
        network.append(neuron)

    for n in range(0, layer_neurons):
        #target_core = random.choice(list(range(64,128)))
        #min_neuron = target_core*1024
        #neuron_list = list(range(min_neuron,min_neuron+1024))
        neuron_list = list(range(layer_neurons, 2*layer_neurons))
        dest_neurons = random.sample(neuron_list, max_connections)
        #print(dest_neurons)

        neuron = [n, n, 1.0, reset, 0, 0, 0]
        for dest in dest_neurons:
            neuron.extend((dest, random.random()))  # Same weight for all connections
        network.append(neuron)
        #print("neuron has {0} fields".format(len(neuron)))

        if (n % 10000) == 0:
            print("Created {0} neurons".format(n))

    if spiking:
        threshold = 3.0
    else:  # never spike
        threshold = 2*layer_neurons

    for n in range(layer_neurons, 2*layer_neurons):
        neuron = [n, n, threshold, reset, 0, 0, 0]
        network.append(neuron)

    #print(network)
    return network


if __name__ == "__main__":
    neurons = []
    spiking_times = []
    spiking_update_energy = []
    spiking_spike_gen_energy = []
    spiking_synapse_energy = []
    spiking_network_energy = []

    layer_neurons = 65536
    #layer_neurons = 4096
    #layer_neurons = 256
    total_neurons = 2 * layer_neurons

    network = fully_connected(layer_neurons, True)
    print("Testing network with {0} neurons and inputs".format(len(network)))
    analysis = run_sim(network)

    neurons.append(len(network))
    spiking_update_energy.append(analysis["update_energy"])
    spiking_spike_gen_energy.append(analysis["spike_gen_energy"])
    spiking_synapse_energy.append(analysis["synapse_energy"])
    spiking_network_energy.append(analysis["network_energy"])

    # Write all the simulation data to csv
    #with open("sim_network.csv", "w") as spiking_csv:
    #    spiking_writer = csv.DictWriter(spiking_csv,
    #                                    ("neurons", "update_energy", "spike_gen_energy",
    #                                     "synapse_energy", "network_energy"))
    #    spiking_writer.writeheader()
    #    for neuron_count, update_energy, spike_gen_energy, synapse_energy, network_energy in zip(neurons, spiking_update_energy,
    #                                              spiking_spike_gen_energy, spiking_synapse_energy,
    #                                              spiking_network_energy):
    #        spiking_writer.writerow({"neurons": neuron_count,
    #                                 "update_energy": update_energy,
    #                                 "spike_gen_energy": spike_gen_energy,
    #                                 "synapse_energy": synapse_energy,
    #                                 "network_energy": network_energy})

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
    plt.show()
