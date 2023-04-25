"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Calibrating SANA-FE against real-world hardware

Use several partitions of a small benchmark to calibrate
the simulator.
"""

import matplotlib
matplotlib.use('Agg')

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt
import pickle
import math
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '/home/usr1/jboyle/neuro/sana-fe')
import utils

np.random.seed(1)

MAX_TILES = 32
MAX_CORES = 4
MAX_COMPARTMENTS = 1024
NETWORK_FILENAME = "runs/connected_layers.net"
ARCH_FILENAME = "loihi.arch"

def run_sim(network, timesteps):
    network.save(NETWORK_FILENAME)
    run_command = ("./sim", ARCH_FILENAME, NETWORK_FILENAME,
               "{0}".format(timesteps))
    print("Command: {0}".format(" ".join(run_command)))
    subprocess.call(run_command)

    with open("stats.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)

    return results


import random
random.seed(1)
def connected_layers(weights):
    network = utils.Network()
    loihi_compartments = utils.init_compartments(32, 4, 1024)

    layer_neuron_count = len(weights)
    threshold = 2.0

    reset = 0
    force_update = True
    log_spikes = False
    log_voltage = False
    leak = 1.0

    neurons_per_core = [0, 0, 0, 0]
    neurons_per_core[0] = layer_neurons

    layer_mapping = [(0, 0) for _ in range(0, neurons_per_core[0])]
    for _ in range(0, neurons_per_core[1]):
        layer_mapping.append((0, 1))
    for _ in range(0, neurons_per_core[2]):
        layer_mapping.append((0, 2))
    for _ in range(0, neurons_per_core[3]):
        layer_mapping.append((0, 3))

    layer_1 = utils.create_layer(network, layer_neuron_count,
                                     loihi_compartments, log_spikes,
                                     log_voltage, force_update, threshold,
                                     reset, leak, mappings=layer_mapping)
    for neuron in layer_1.neurons:
        neuron.add_bias((random.random()) - 0.75)

    force_update = False
    neurons_per_core = [0, 0, 0, 0, 0]
    neurons_per_core[0] = min(layer_neurons, 1024 - layer_neurons)
    neurons_per_core[1] = layer_neurons - neurons_per_core[0]

    layer_mapping = [(0, 0) for _ in range(0, neurons_per_core[0])]
    for _ in range(0, neurons_per_core[1]):
        layer_mapping.append((0, 1))
    for _ in range(0, neurons_per_core[2]):
        layer_mapping.append((0, 2))
    for _ in range(0, neurons_per_core[3]):
        layer_mapping.append((0, 3))
    for _ in range(0, neurons_per_core[4]):
        layer_mapping.append((1, 0))

    layer_2 = utils.create_layer(network, layer_neuron_count,
                                     loihi_compartments, log_spikes,
                                     log_voltage, force_update, threshold,
                                     reset, leak, mappings=layer_mapping)

    for src in layer_1.neurons:
        for dest in layer_2.neurons:
            # Take the ID of the neuron in the 2nd layer
            weight = float(weights[src.id][dest.id]) / 256
            if abs(weight) >= (1.0 / 256):
                # Zero weights are pruned i.e. removed
                src.add_connection(dest, weight)

    return network


if __name__ == "__main__":
    # This experiment looks at two fully connected layers, spiking

    with open("sandia_data/weights_loihi.pkl", "rb") as weights_file:
        weights = pickle.load(weights_file)

    neuron_counts = []
    spiking_times = []
    spiking_energy = []

    timesteps = 64
    n = 30
    layer_neurons = n*n

    #network = fully_connected(layer_neurons, spiking=True, probability=connection_probabilities[i-1])
    commands = connected_layers(weights[n-1].transpose())
    print("Testing network with {0} neurons".format(2*layer_neurons))
    results = run_sim(commands, timesteps)

    neuron_counts.append(layer_neurons*2)
    spiking_times.append(results["time"])
    spiking_energy.append(results["energy"])

    # Write all the simulation data to csv
    with open("runs/interception.csv", "w") as spiking_csv:
        spiking_writer = csv.DictWriter(spiking_csv,
                                        ("neuron_counts", "energy", "time"))
        spiking_writer.writeheader()
        for count, time, energy_val in zip(neuron_counts, spiking_times,
                                                  spiking_energy):
            spiking_writer.writerow({"neuron_counts": count,
                                     "energy": energy_val,
                                     "time": time})

    plt.rcParams.update({'font.size': 9, 'lines.markersize': 3})

    # Plot the latency
    run_data = pd.read_csv("perf.csv")
    # There is a weird effect, that the first sample of all inputs > 1 is
    #  a 0 value. Just ignore the entries for both arrays (so we have
    #  timestep-1)
    firing_neurons = run_data.loc[:, "total_neuron_updates"]
    synapse_reads = run_data.loc[:, "total_synapse_reads"]
    energy = run_data.loc[:, "total_energy"]
    latency = run_data.loc[:, "time"]

    print(firing_neurons)

    plt.figure(figsize=(8.4, 2.2))
    plt.xlabel("Time-step")
    ax1 = plt.subplot()
    line1, = ax1.plot(np.arange(1, (timesteps+1)), firing_neurons/1e3, 'x', color="#1f77b4")
    line2, = ax1.plot(np.arange(1, (timesteps+1)), synapse_reads/1e3, 'o', color="darkorange")
    ax1.set_ylabel("Activity Counts (x$10^3$)")
    #ax1.ticklabel_format(style="sci", axis="y", scilimits=(0,0))

    ax2 = ax1.twinx()
    line3, = ax2.plot(np.arange(1, (timesteps+1)), energy*1e6, 'k-')
    #ax2.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    ax2.set_ylabel("Energy (uJ)")
    ax2.set_ylim((0.0, 0.3))

    ax3 = ax1.twinx()
    line4, = ax3.plot(np.arange(1, (timesteps+1)), latency * 1e6, 'r--')
    ax3.spines["right"].set_position(('outward', 50))
    ax3.set_ylabel("Latency (us)")
    ax3.set_ylim((0.0, 20.0))

    plt.legend((line1, line2, line3, line4),
               ("Neuron Updates", "Synapse Reads", "Dynamic Energy", "Time-step Latency"),
               loc="upper left", prop={"size": 8})
    plt.tight_layout(pad=0)
    plt.savefig("runs/time_series.pdf")
    plt.savefig("runs/time_series.png")
