"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Python API helper functions, mostly to aid with network generation
"""

import matplotlib
matplotlib.use('Agg')

import csv
import subprocess
import yaml
import pickle
import math

from matplotlib import pyplot as plt

NETWORK_FILENAME = "runs/connected_layers.net"
ARCH_FILENAME = "loihi.arch"


class Network:
    def __init__(self, external_inputs=0):
        self.external_inputs = external_inputs
        self.inputs = []
        self.groups = []
        self._max_tiles = None
        self._max_cores = None
        self._max_compartments = None

    def create_group(self, threshold, reset, leak):
        group_id = len(self.groups)
        group = NeuronGroup(group_id, threshold, reset, leak)
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
    def __init__(self, group_id, threshold, reset, leak):
        self.id = group_id
        self.neurons = []
        self.threshold = threshold
        self.reset = reset
        # TODO: remove defaults
        self.reverse_reset = 0.0
        self.reverse_threshold = 0.0
        self.leak = leak

    def to_command(self):
        neuron_count = len(self.neurons)
        return (f"g {neuron_count} {self.threshold} {self.reset} "
                f"{self.reverse_threshold} {self.reverse_reset} {self.leak} "
                f"hard none\n")

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
        self.bias = 0.0

    def add_connection(self, dest, weight):
        self.connections.append((dest, weight))

    def add_bias(self, bias):
        self.bias = bias

    def to_command(self):
        log_spikes = int(self.log_spikes)
        log_voltage = int(self.log_voltage)
        force_update = int(self.force_update)
        self.leak = 1.0
        command = (f"n {self.group_id} {self.id} {self.bias} {log_spikes} "
                   f"{log_voltage} {force_update}")

        for connection in self.connections:
            dest_neuron, weight = connection
            command += " {0} {1} {2}".format(
                dest_neuron.group_id, dest_neuron.id, weight)

        command += '\n'
        command += "& {0} {1} {2} {3}\n".format(self.group_id, self.id,
                                                self.tile, self.core)
        return command


def init_compartments(max_tiles, max_cores, max_compartments):
    compartments = []
    for tile in range(0, max_tiles):
        c = []
        for core in range(0, max_cores):
            c.append(max_compartments)
        compartments.append(c)

    return compartments


def map_neuron_to_compartment(compartments):
    for tile, cores in enumerate(compartments):
        for core, _ in enumerate(cores):
            if compartments[tile][core] > 0:
                compartments[tile][core] -= 1
                return tile, core

    # No free compartments left
    return None, None


def create_layer(network, layer_neuron_count, compartments,
                 log_spikes, log_voltage, force_update, threshold, reset,
                 leak, mappings=None):
    print("Creating layer with {0} neurons".format(layer_neuron_count))
    #print("Compartments free: {0}".format(compartments))
    layer_group = network.create_group(threshold, reset, leak)

    if mappings is not None:
        #print("mappings: {0} len({1}) layer_neuron_count: {2}".format(
        #    mappings, len(mappings), layer_neuron_count))
        assert(len(mappings) == layer_neuron_count)

    for i in range(0, layer_neuron_count):
        if (i % 10000) == 0:
            print(f"Creating neuron {i}")
        neuron = layer_group.create_neuron(log_spikes, log_voltage,
                                           force_update)

        if mappings is not None:
            tile, core = mappings[i]
        else:
            tile, core = map_neuron_to_compartment(compartments)
        neuron.tile, neuron.core = tile, core

    return layer_group
