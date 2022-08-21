"""First attempt at simulation run script.

Run a few basic experiments, show how we might interface a python
script with the simulator kernel.
"""

import matplotlib
matplotlib.use('Agg')

import csv
import subprocess
import yaml
import pickle
import math

from matplotlib import pyplot as plt

MAX_TILES = 32
MAX_CORES = 4
MAX_COMPARTMENTS = 1024
NETWORK_FILENAME = "../runs/connected_layers.net"
ARCH_FILENAME = "../loihi.arch"


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

    command = ("../sim", ARCH_FILENAME, NETWORK_FILENAME,
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
    #  Returns a list of tuples (gid, nid) for each neuron in the layer
    print("Creating layer with {0} neurons".format(layer_neuron_count))
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



def create_groups(neuron_groups, threshold, reset):
    commands = []
    for group_id, group in enumerate(neuron_groups):
        # Create the neuron group with the right number of neurons
        commands.append(['g', group.neuron_count, threshold, reset])
        # Map the group to the right tile and core
        # TODO: support architectures that can have multiple hardware units
        #  in the same core e.g. it might have multiple different soma
        #  processors
        commands.append(['&', group_id, group.tile_id, group.core_id,
                         0, 0, 0, 0, 0])

    return commands



def connected_layers(weights, spiking=True):
    commands = []
    loihi_compartments = loihi_init_compartments()

    layer_neuron_count = len(weights)
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neuron_count

    reset = 0
    force_update = True
    log_spikes = False
    log_voltage = False

    neuron_groups = []
    layer_1 = create_layer(layer_neuron_count, neuron_groups,
                           loihi_compartments)
    layer_2 = create_layer(layer_neuron_count, neuron_groups,
                           loihi_compartments)

    commands = create_groups(neuron_groups, reset, threshold)
    # *** layer 1 ***
    for n in range(0, layer_neuron_count):
        gid, nid = layer_1[n]
        neuron = ['n', gid, nid, int(log_spikes), int(log_voltage),
                  int(force_update)]
        for dest in range(0, layer_neuron_count):
            # Take the ID of the neuron in the 2nd layer
            weight = float(weights[n][dest]) / 255
            if weight != 0:
                # Zero weights are pruned i.e. removed
                dest_gid, dest_nid = layer_2[dest]
                neuron.extend((dest_gid, dest_nid, weight))
        commands.append(neuron)

    # *** layer 2 ***
    for n in range(0, layer_neuron_count):
        gid, nid = layer_2[n]
        neuron = ['n', gid, nid, int(log_spikes), int(force_update),
                  int(force_update)]
        commands.append(neuron)

    #print(network)
    return commands

