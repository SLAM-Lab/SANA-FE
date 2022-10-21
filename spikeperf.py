"""Python API

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
NETWORK_FILENAME = "runs/connected_layers.net"
ARCH_FILENAME = "loihi.arch"


class Network:
    def __init__(self, external_inputs=0):
        self.external_inputs = external_inputs
        self.inputs = []
        self.groups = []

    def create_group(self, threshold, reset):
        group_id = len(self.groups)
        group = NeuronGroup(group_id, threshold, reset)
        self.groups.append(group)
        return group

    def create_input(self):
        input_id = len(self.inputs)
        input_node = Input(input_id)
        self.inputs.append(input_node)
        return input_node

    def save(self, filename):
        with open(filename, 'w') as network_file:
            if self.external_inputs > 0:
                network_file.write("e {0} rate\n".format(self.external_inputs))
            for group in self.groups:
                network_file.write(group.to_command())

            for group in self.groups:
                for neuron in group.neurons:
                    network_file.write(neuron.to_command())

            for input_node in self.inputs:
                network_file.write(input_node.to_command())


class NeuronGroup:
    def __init__(self, group_id, threshold, reset):
        self.id = group_id
        self.neurons = []
        self.threshold = threshold
        self.reset = reset

    def to_command(self):
        return "g {0} {1} {2}\n".format(len(self.neurons), self.threshold,
                                        self.reset)

    def create_neuron(self, log_spikes, log_voltage, force_update):
        neuron_id = len(self.neurons)
        neuron = Neuron(self.id, neuron_id, log_spikes, log_voltage,
                        force_update)
        self.neurons.append(neuron)
        #print(self.neurons)
        return neuron


class Input:
    def __init__(self, input_id):
        self.id = input_id
        self.connections = []

    def add_connection(self, dest, weight):
        self.connections.append((dest, weight))

    def to_command(self):
        command = "< {0}".format(self.id)
        for connection in self.connections:
            dest_neuron, weight = connection
            command += " {0} {1} {2}".format(
                dest_neuron.group_id, dest_neuron.id, weight)
        command += '\n'
        return command

class Neuron:
    def __init__(self, group_id, neuron_id, log_spikes=False, log_voltage=False,
               force_update=False):
        self.group_id = group_id
        self.id = neuron_id
        self.log_spikes = log_spikes
        self.log_voltage = log_voltage
        self.force_update = force_update
        self.connections = []
        self.tile = None
        self.core = None

    def add_connection(self, dest, weight):
        self.connections.append((dest, weight))

    def to_command(self):
        command = "n {0} {1} {2} {3} {4}".format(
            self.group_id, self.id, int(self.log_spikes), int(self.log_voltage),
            int(self.force_update))

        for connection in self.connections:
            dest_neuron, weight = connection
            command += " {0} {1} {2}".format(
                dest_neuron.group_id, dest_neuron.id, weight)

        command += '\n'
        command += "& {0} {1} {2} {3}\n".format(self.group_id, self.id,
                                                self.tile, self.core)
        return command


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


def create_layer(network, layer_neuron_count, loihi_compartments,
                 log_spikes, log_voltage, force_update, threshold, reset,
                 mappings=None):
    print("Creating layer with {0} neurons".format(layer_neuron_count))
    print("Compartments free: {0}".format(loihi_compartments))
    layer_group = network.create_group(threshold, reset)

    if mappings is not None:
        assert(len(mappings) == layer_neuron_count)

    for i in range(0, layer_neuron_count):
        neuron = layer_group.create_neuron(log_spikes, log_voltage,
                                           force_update)

        if mappings is not None:
            tile, core = mappings[i]
        else:
            tile, core = loihi_map_neuron_to_compartment(loihi_compartments)
        neuron.tile, neuron.core = tile, core

    return layer_group
