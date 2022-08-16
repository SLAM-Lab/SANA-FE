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
import pickle
import math

MAX_COMPARTMENTS = 1024
MAX_CORES = 128
NETWORK_FILENAME = "connected_layer.net"
ARCH_FILENAME = "loihi.list"

def run_sim(network):
    fields = ["Neuron ID", "Compartment ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
    #with open(NETWORK_FILENAME, "w") as csv_file:
    #    writer = csv.writer(csv_file)
    #    writer.writerow(fields)
    #    writer.writerows(network)

    with open(NETWORK_FILENAME, "w") as out_file:
        for line in network["groups"]:
            els = [str(s) for s in line]
            print(" ".join(els), file=out_file)
        for line in network["mappings"]:
            els = [str(s) for s in line]
            print(" ".join(els), file=out_file)
        for line in network["neurons"]:
            els = [str(s) for s in line]
            print(" ".join(els), file=out_file)

    timesteps = 2
    command = ("./sim", ARCH_FILENAME,
               NETWORK_FILENAME, "{0}".format(timesteps))
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    with open("stats.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)

    return results


def create_layer(network, compartments_per_core, layer_neurons, threshold, reset):
    # Create a layer of neurons, which may be one or more groups
    layer_groups = []
    neurons_left = layer_neurons
    compartments_left = 1024

    # TODO: thinking about this, forcing each neuron group to map to the same
    #  hardware makes handling the group kind of nasty, need to think about this
    while neurons_left:
        if ((len(compartments_per_core) == 0) or
                           (compartments_per_core[-1] == MAX_COMPARTMENTS)):
            compartments_per_core.append(0)

        compartments_left = MAX_COMPARTMENTS - compartments_per_core[-1]

        # Figure out the total neurons TODO: this is a mess
        first_neuron_id = 0
        for g in network["groups"]:
            count = g[1]
            first_neuron_id += count 
        group_neurons = min(neurons_left, compartments_left)
        group = ['g', group_neurons, threshold, reset]
        network["groups"].append(group)
    
        group_id = len(network["groups"]) - 1

        layer_groups.append((group_id, first_neuron_id, group_neurons))

        # Map this group to hardware
        core = len(compartments_per_core) - 1
        core_id = core % 4
        tile_id = core // 4  # Cores per tile

        # TODO: we need to know which group this is
        mapping = ['&', group_id, tile_id, core_id, 0, 0, 0, 0, 0]
        network["mappings"].append(mapping)

        neurons_left -= group_neurons

    return layer_groups


# TODO: this is horribly convoluted, I should probably simplify and just have
#  a mapping of every neuron in the simulation / layer to a neuron group
def get_neuron_id(groups, neuron):
    for info in groups:
        gid, first_neuron_id, group_neurons = info
        if (neuron >= first_neuron_id) and (neuron < (first_neuron_id + group_neurons)):
            nid = neuron % group_neurons
            return gid, nid

    print("Error couldn't find this group")
    exit()
    return None, None


def connected_layer(weights, spiking=True):
    network = {"neurons": [], "groups": [], "mappings": []}
    compartments_per_core = []

    layer_neurons = len(weights)
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neurons

    reset = 0
    force_update = True
    log_spikes = False
    log_voltage = False

    layer_1 = create_layer(network, compartments_per_core, layer_neurons, threshold, reset)
    layer_2 = create_layer(network, compartments_per_core, layer_neurons, threshold, reset)

    # *** layer 1 ***
    for n in range(0, layer_neurons):
        gid, nid = get_neuron_id(layer_1, n)
        neuron = ['n', gid, nid, int(log_spikes), int(log_voltage), int(force_update)]
        for dest in range(0, layer_neurons):
            # Take the ID of the neuron in the 2nd layer
            weight = float(weights[n][dest]) / 255
            if weight != 0:
                # Zero weights are pruned i.e. removed
                dest_gid, dest_nid = get_neuron_id(layer_2, (layer_neurons+dest))
                neuron.extend((dest_gid, dest_nid, weight))
        network["neurons"].append(neuron)

    # *** layer 2 ***
    for n in range(layer_neurons, 2*layer_neurons):
        gid, nid = get_neuron_id(layer_2, n)
        neuron = ['n', gid, nid, int(log_spikes), int(force_update), int(force_update)]
        network["neurons"].append(neuron)

    #print(network)
    return network


import random
def fully_connected(layer_neurons, spiking=True, probability=1.0):
    # Two layers, fully connected
    network = {"neurons": [], "groups": [], "mappings": []}
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neurons

    reset = 0
    force_update = True
    log_spikes = False
    log_voltage = False

    layer_1 = create_layer(network, layer_neurons, threshold, reset)
    layer_2 = create_layer(network, layer_neurons, threshold, reset)

    weight = 1.0
    reset = 0
    force_update = True
    for n in range(0, layer_neurons):
        gid, nid = get_neuron_id(layer_1, n)
        neuron = ['n', gid, nid, int(log_spikes), int(force_update), int(force_update)]

        for dest in range(layer_neurons, 2*layer_neurons):
            if random.random() < probability:
                neuron.extend((dest, weight))  # Same weight for all connections
        network["neurons"].append(neuron)

    for n in range(0, layer_neurons):
        gid, nid = get_neuron_id(layer_2, n)
        neuron = ['n', gid, nid, int(log_spikes), int(force_update), int(force_update)]
        network["neurons"].append(neuron)

    return network


def empty(neurons, max_compartments=MAX_COMPARTMENTS):
    network = {"neurons": [], "groups": [], "mappings": []}
    threshold = -1.0  # never spike
    reset = 0.0

    compartment = 0
    # Map a number of neurons onto cores, but we may not necessarily use all
    #  compartments of each core
    layer = create_layer(network, neurons, threshold, reset)
    for n in range(0, neurons):
        gid, nid = get_neuron_id(layer, n)
        neuron = ['n', gid, nid, int(log_spikes), int(force_update), int(force_update)]
        network["neurons"].append(neuron)
        compartment += 1

    return network


if __name__ == "__main__":
    #core_count = [1, 2, 4, 8, 16, 32, 64, 128]
    times = {0: [], 256: [], 512: [], 768: [], 1024: []}
    energy = {0: [], 256: [], 512: [], 768: [], 1024: []}
    """
    for cores in core_count:
        for compartments in range(0, MAX_COMPARTMENTS+1, 256):
            n = compartments * cores
            network = empty(n, compartments)
            results = run_sim(network, cores)

            times[compartments].append(results["time"])
            energy[compartments].append(results["energy"])

    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(5.5, 5.5))
    for compartments in range(0, MAX_COMPARTMENTS+1, 256):
        plt.plot(core_count, times[compartments], "-o")
    plt.ylabel("Time (s)")
    plt.xlabel("Cores Used")
    plt.legend(("0", "256", "512", "768", "1024"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("empty_times.png")

    plt.figure(figsize=(5.5, 5.5))
    for compartments in range(0, MAX_COMPARTMENTS+1, 256):
        plt.plot(core_count, energy[compartments], "-o")
    plt.ylabel("Energy (J)")
    plt.xlabel("Cores Used")
    plt.legend(("0", "256", "512", "768", "1024"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("empty_energy.png")
    """
    # This experiment looks at two fully connected layers, spiking

    with open("sandia_data/weights.pkl", "rb") as weights_file:
        weights = pickle.load(weights_file)

    neurons = []
    spiking_times = []
    spiking_energy = []

    #for i in range(1, 2):
    for i in range(1, 31):
        layer_neurons = i*i

        #network = fully_connected(layer_neurons, spiking=True, probability=connection_probabilities[i-1])
        network = connected_layer(weights[i-1], spiking=True)
        print("Testing network with {0} neurons".format(2*layer_neurons))
        results = run_sim(network)

        neurons.append(layer_neurons*2)
        spiking_times.append(results["time"])
        spiking_energy.append(results["energy"])

    # Write all the simulation data to csv
    with open("sim_spiking.csv", "w") as spiking_csv:
        spiking_writer = csv.DictWriter(spiking_csv,
                                        ("neurons", "energy", "time"))
        spiking_writer.writeheader()
        for neuron_count, time, energy_val in zip(neurons, spiking_times,
                                                  spiking_energy):
            spiking_writer.writerow({"neurons": neuron_count,
                                     "energy": energy_val,
                                     "time": time})

    """
    neurons = []
    nonspiking_times = []
    nonspiking_energy = []

    # The second experiment looks at two fully connected layers, not spiking
    for i in range(1, 31):
        layer_neurons = i*i

        network = fully_connected(layer_neurons, spiking=False)
        print("Testing network with {0} neurons".format(len(network)))
        results = run_sim(network)

        neurons.append(layer_neurons*2)
        nonspiking_times.append(results["time"])
        nonspiking_energy.append(results["energy"])

    with open("sim_nonspiking.csv", "w") as nonspiking_csv:
        nonspiking_writer = csv.DictWriter(nonspiking_csv,
                                           ("neurons", "energy", "time"))
        nonspiking_writer.writeheader()
        for neuron_count, time, energy_val in zip(neurons, nonspiking_times,
                                                  nonspiking_energy):
            nonspiking_writer.writerow({"neurons": neuron_count,
                                        "energy": energy_val,
                                        "time": time})
    """
    # **************************************************************************
    # Read Loihi measurement data from csv, this is only available to me locally
    #  since this is restricted data!
    neurons = []
    loihi_times_spikes = []
    loihi_energy_spikes = []

    spiking_energy = []
    spiking_times = []
    with open("sim_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            spiking_times.append(float(row["time"]))
            spiking_energy.append(float(row["energy"]))
            neurons.append(int(row["neurons"]))

    nonspiking_times = []
    nonspiking_energy = []
    with open("sim_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            nonspiking_times.append(float(row["time"]))
            nonspiking_energy.append(float(row["energy"]))

    with open("loihi_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            loihi_times_spikes.append(float(row["time"]))
            loihi_energy_spikes.append(float(row["energy"]))

    loihi_times_no_spikes = []
    loihi_energy_no_spikes = []
    with open("loihi_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            loihi_times_no_spikes.append(float(row["time"]))
            loihi_energy_no_spikes.append(float(row["energy"]))

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, spiking_times, "-o")
    plt.plot(neurons, loihi_times_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, spiking_energy, "-o")
    plt.plot(neurons, loihi_energy_spikes, "--x", color="orange")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_energy.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, nonspiking_energy, "-o")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_energy_sim_only.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, nonspiking_times, "-o")
    plt.plot(neurons, loihi_times_no_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_not_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, nonspiking_energy, "-o")
    plt.plot(neurons, loihi_energy_no_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_not_spiking_energy.png")

    # Some additional plots to highlight trends
    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, loihi_times_spikes, "--x", color="orange")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Measured",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_time_loihi_only.png")

    #plt.show()
