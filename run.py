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

MAX_TILES = 32
MAX_CORES = 4
MAX_COMPARTMENTS = 1024
NETWORK_FILENAME = "runs/connected_layers.net"
ARCH_FILENAME = "loihi.arch"


class NeuronGroup:
    def __init__(self, tile_id, core_id):
        self.tile_id = tile_id
        self.core_id = core_id
        self.neuron_count = 0


def run_sim(commands, timesteps):
    fields = ["Neuron ID", "Compartment ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]

    with open(NETWORK_FILENAME, "w") as out_file:
        for command in commands:
            field_strs = [str(field) for field in command]
            print(" ".join(field_strs), file=out_file)

    command = ("./sim", ARCH_FILENAME, NETWORK_FILENAME,
               "{0}".format(timesteps))
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    with open("stats.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)

    return results


# TODO: generalize this to map networks to any hardware. Probably using the
#  same architecture description
def loihi_init_compartments():
    compartments = []
    for tile in range(0, MAX_TILES):
        c = []
        for core in range(0, MAX_CORES):
            c.append(MAX_COMPARTMENTS)
        compartments.append(c)

    return compartments


def loihi_map_neuron_to_compartment(loihi_compartments):
    for tile in range(0, MAX_TILES):
        for core in range(0, MAX_CORES):
            if loihi_compartments[tile][core] > 0:
                loihi_compartments[tile][core] -= 1
                return tile, core

    # No free compartments left
    return None, None


def create_layer(layer_neuron_count, groups, loihi_compartments):
    # Create a layer of neurons, which may be one or more groups
    print("Creating layer with {0} neurons", layer_neuron_count)
    print("Compartments free: {0}".format(loihi_compartments))
    prev_core = (None, None)  # (tile id, core id)
    layer_neurons = []
    for _ in range(0, layer_neuron_count):
        core = loihi_map_neuron_to_compartment(loihi_compartments)
        # If we use up all compartments in the core we need to spill over to the
        #  next one. *At the moment* we need a different group for each core.
        #  I.e. Neuron groups need to share common hardware.
        if core != prev_core:
            tile_id, core_id = core
            # We always add neurons to the last group
            groups.append(NeuronGroup(tile_id, core_id))
            # Now track the current core being mapped to
            prev_core = core

        groups[-1].neuron_count += 1
        # The neuron group is tracked out of all groups created for the spiking
        #  network. One layer might require multiple groups if that layer
        #  spreads over multiple cores. TODO: is to decide whether we allow
        #  one neuron group to be spread over multiple cores, then each layer
        #  can just have its own group
        group_id = len(groups) - 1
        # The neuron ID is just the neuron # in the group. If the layer
        #  only uses one group, then that ID corresponds to the neuron #
        #  in that layer. This is something only needed internally. When
        #  connecting neurons in scripts we index them without thinking about
        #  which group they are in.
        neuron_id = groups[-1].neuron_count - 1

        layer_neurons.append((group_id, neuron_id))

    return layer_neurons


def connected_layers(weights, spiking=True):
    commands = []
    loihi_compartments = loihi_init_compartments()

    layer_neurons = len(weights)
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neurons

    reset = 0
    force_update = True
    log_spikes = False
    log_voltage = False

    neuron_groups = []
    layer_1 = create_layer(layer_neurons, neuron_groups, loihi_compartments)
    layer_2 = create_layer(layer_neurons, neuron_groups, loihi_compartments)

    for group_id, group in enumerate(neuron_groups):
        # Create the neuron group with the right number of neurons
        commands.append(['g', group.neuron_count, threshold, reset])
        # Map the group to the right tile and core
        # TODO: support architectures that can have multiple hardware units
        #  in the same core e.g. it might have multiple different soma
        #  processors
        commands.append(['&', group_id, group.tile_id, group.core_id,
                         0, 0, 0, 0, 0])

    # *** layer 1 ***
    for n in range(0, layer_neurons):
        gid, nid = layer_1[n]
        neuron = ['n', gid, nid, int(log_spikes), int(log_voltage),
                  int(force_update)]
        for dest in range(0, layer_neurons):
            # Take the ID of the neuron in the 2nd layer
            weight = float(weights[n][dest]) / 255
            if weight != 0:
                # Zero weights are pruned i.e. removed
                dest_gid, dest_nid = layer_2[dest]
                neuron.extend((dest_gid, dest_nid, weight))
        commands.append(neuron)

    # *** layer 2 ***
    for n in range(0, layer_neurons):
        gid, nid = layer_2[n]
        neuron = ['n', gid, nid, int(log_spikes), int(force_update),
                  int(force_update)]
        commands.append(neuron)

    #print(network)
    return commands


import random
def fully_connected(layer_neurons, spiking=True, force_update=False,
                    probability=1.0):
    # Two layers, fully connected
    commands = []
    loihi_compartments = loihi_init_compartments()

    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neurons

    reset = 0
    log_spikes = False
    log_voltage = False

    neuron_groups = []
    layer_1 = create_layer(layer_neurons, neuron_groups, loihi_compartments)
    layer_2 = create_layer(layer_neurons, neuron_groups, loihi_compartments)

    for group_id, group in enumerate(neuron_groups):
        # Create the neuron group with the right number of neurons
        commands.append(['g', group.neuron_count, threshold, reset])
        # Map the group to the right tile and core
        # TODO: support architectures that can have multiple hardware units
        #  in the same core e.g. it might have multiple different soma
        #  processors
        commands.append(['&', group_id, group.tile_id, group.core_id,
                         0, 0, 0, 0, 0])

    weight = 1.0
    reset = 0
    for n in range(0, layer_neurons):
        gid, nid = layer_1[n]
        neuron = ['n', gid, nid, int(log_spikes), int(log_voltage),
                  int(force_update)]

        for dest in range(layer_neurons, 2*layer_neurons):
            if random.random() < probability:
                neuron.extend((dest, weight))  # Same weight for all connections
        commands.append(neuron)

    for n in range(0, layer_neurons):
        gid, nid = layer_2[n]
        neuron = ['n', gid, nid, int(log_spikes), int(log_voltage),
                  int(force_update)]
        commands.append(neuron)

    return commands


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

    neuron_counts = []
    spiking_times = []
    spiking_energy = []

    #for i in range(1, 4):
    timesteps = 2
    for i in range(1, 31):
        layer_neurons = i*i

        #network = fully_connected(layer_neurons, spiking=True, probability=connection_probabilities[i-1])
        commands = connected_layers(weights[i-1], spiking=True)
        print("Testing network with {0} neurons".format(2*layer_neurons))
        results = run_sim(commands, timesteps)

        neuron_counts.append(layer_neurons*2)
        spiking_times.append(results["time"])
        spiking_energy.append(results["energy"])

    # Write all the simulation data to csv
    with open("runs/sim_spiking.csv", "w") as spiking_csv:
        spiking_writer = csv.DictWriter(spiking_csv,
                                        ("neuron_counts", "energy", "time"))
        spiking_writer.writeheader()
        for count, time, energy_val in zip(neuron_counts, spiking_times,
                                                  spiking_energy):
            spiking_writer.writerow({"neuron_counts": count,
                                     "energy": energy_val,
                                     "time": time})

    neuron_counts = []
    nonspiking_times = []
    nonspiking_energy = []

    # The second experiment looks at two fully connected layers, not spiking
    for i in range(1, 31):
        layer_neurons = i*i

        commands = fully_connected(layer_neurons, spiking=False,
                                   force_update=True)
        print("Testing network with {0} neurons".format(layer_neurons*2))
        results = run_sim(commands, timesteps)

        neuron_counts.append(layer_neurons*2)
        nonspiking_times.append(results["time"] / timesteps)
        nonspiking_energy.append(results["energy"] / timesteps)

    with open("runs/sim_nonspiking.csv", "w") as nonspiking_csv:
        nonspiking_writer = csv.DictWriter(nonspiking_csv,
                                           ("neurons", "energy", "time"))
        nonspiking_writer.writeheader()
        for neuron_count, time, energy_val in zip(neuron_counts, nonspiking_times,
                                                  nonspiking_energy):
            nonspiking_writer.writerow({"neurons": neuron_count,
                                        "energy": energy_val,
                                        "time": time})

    # **************************************************************************
    # Read Loihi measurement data from csv, this is only available to me locally
    #  since this is restricted data!
    neuron_counts = []
    loihi_times_spikes = []
    loihi_energy_spikes = []

    spiking_energy = []
    spiking_times = []
    with open("runs/sim_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            spiking_times.append(float(row["time"]))
            spiking_energy.append(float(row["energy"]))
            neuron_counts.append(int(row["neuron_counts"]))

    nonspiking_times = []
    nonspiking_energy = []
    with open("runs/sim_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            nonspiking_times.append(float(row["time"]))
            nonspiking_energy.append(float(row["energy"]))

    with open("sandia_data/loihi_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            loihi_times_spikes.append(float(row["time"]))
            loihi_energy_spikes.append(float(row["energy"]))

    loihi_times_no_spikes = []
    loihi_energy_no_spikes = []
    with open("sandia_data/loihi_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            loihi_times_no_spikes.append(float(row["time"]))
            loihi_energy_no_spikes.append(float(row["energy"]))

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, spiking_times, "-o")
    plt.plot(neuron_counts, loihi_times_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/connected_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, spiking_energy, "-o")
    plt.plot(neuron_counts, loihi_energy_spikes, "--x", color="orange")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/connected_spiking_energy.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, nonspiking_energy, "-o")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/connected_spiking_energy_sim_only.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, nonspiking_times, "-o")
    plt.plot(neuron_counts, loihi_times_no_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/connected_not_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, nonspiking_energy, "-o")
    plt.plot(neuron_counts, loihi_energy_no_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/connected_not_spiking_energy.png")

    # Some additional plots to highlight trends
    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, loihi_times_spikes, "--x", color="orange")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Measured",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/connected_spiking_time_loihi_only.png")

    #plt.show()
