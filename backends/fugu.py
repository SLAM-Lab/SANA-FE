#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
isort:skip_file
"""
import os
import sys

# fmt: off
from collections import deque
from warnings import warn

import fugu.simulators.SpikingNeuralNetwork as snn
import fugu
from fugu import Scaffold, Brick
from fugu.bricks import Vector_Input
from fugu.backends import snn_Backend

from .backend import Backend
from ..utils.export_utils import results_df_from_dict
from ..utils.misc import CalculateSpikeTimes

"""
keep this file structure-
Fugu
├── backends
│   ├── sanafe_backend.py
│   │...
│...
SANAFE
├── sim.py
│...
"""

BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
FUGU_DIR = os.path.dirname(BACKEND_DIR)
FUGU_PROJECT_DIR = os.path.dirname(FUGU_DIR)
PROJECT_DIR = os.path.dirname(FUGU_PROJECT_DIR)
sys.path.append(PROJECT_DIR)
import sim as sanafe
# Set some flags for the dynamic linking library
# Important to do before importing the simcpp .so library!
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)
import sanafecpp

ARCH_FILENAME = "arch/loihi.yaml"

class sanafe_Backend(Backend):

    def map_neurons(self):
        return 0

    def _build_network(self):
        self.nn = snn.NeuralNetwork()
        neuron_dict = {}
        """
        Add Neurons
        """
        # Add input neurons, as identified by circuit information.
        record_all = (self.record == 'all')
        for node, vals in self.fugu_circuit.nodes.data():
            if vals.get('layer') != 'input': continue
            for list in vals['output_lists']:
                for neuron in list:
                    n = snn.InputNeuron(neuron, record=record_all)
                    neuron_dict[neuron] = n
                    self.nn.add_neuron(n)

        # Add all other neurons.
        for neuron, props in self.fugu_graph.nodes.data():
            if neuron in neuron_dict: continue
            voltage = props.get('voltage', 0.0)
            threshold = props.get('threshold', 1.0)
            reset_voltage = props.get('reset_voltage', 0.0)
            leakage = 1.0 - props.get('decay', 0.0)
            bias = props.get('bias', 0.0)
            p = props.get('p', 1.0)
            if 'potential' in props: Vinit   = props['potential']
            if 'leakage_constant' in props: Vretain = props['leakage_constant']
            n = snn.LIFNeuron(neuron, voltage=voltage, threshold=threshold,
                              reset_voltage=reset_voltage,
                              leakage_constant=leakage,
                              bias=bias, p=p, record=record_all)
            neuron_dict[neuron] = n
            self.nn.add_neuron(n)

        # Tag output neurons based on circuit information.
        if not record_all:
            for node, vals in self.fugu_circuit.nodes.data():
                if vals.get("layer") != "output": continue
                for list in vals["output_lists"]:
                    for neuron in list:
                        neuron_dict[neuron].record = True

        for n1, n2, props in self.fugu_graph.edges.data():
            delay  = int(props.get("delay",  1))
            weight = float(props.get("weight", 1.0))
            syn = snn.Synapse(neuron_dict[n1], neuron_dict[n2], delay=delay,
                              weight=weight)
            self.nn.add_synapse(syn)

        del neuron_dict
        '''
        Set initial input values
        '''
        self.set_input_spikes()

        if self.debug_mode:
            self.nn.list_neurons()

    def compile(self, scaffold, compile_args={}):
        # creates neuron populations and synapses
        self.scaffold = scaffold
        self.fugu_circuit = scaffold.circuit
        self.fugu_graph = scaffold.graph
        self.brick_to_number = scaffold.brick_to_number
        if "record" in compile_args:
            self.record = compile_args['record']
        else:
            self.record = False
        if "ds_format" in compile_args:
            self.ds_format = compile_args['ds_format']
        else:
            self.ds_format = True
        if "debug_mode" in compile_args:
            self.debug_mode = compile_args['debug_mode']
        else:
            self.debug_mode = False

        self._build_network()

    def _find_node(name):
        return None

    def run(self, n_steps=10, return_potential=False, debug_mode=False):
        "Convert to SANA-FE network and save to file"
        arch = sim.load_arch_yaml(ARCH_FILENAME)
        network = sanafecpp.Network(arch)

        id_dict = {} # Group name -> gid
        nid_dict = [] # [Group index] Neuron name -> SANA-FE nid
        nid_dict_idx = {} # Neuron Name -> gid
        neurons = list(enumerate(sorted(self.scaffold.graph.nodes)))

        # Add the neurons that belong to the current group
        #TODO: optimize move out of this for loop
        group = network.create_neuron_group(len(neurons))
        for j, neuron in neurons:
            neuron_properties = self.scaffold.graph.nodes[neuron]
            threshold = neuron_properties["threshold"]
            decay = neuron_properties["decay"]
            reset = neuron_properties.get("reset", 0.0)

            i = j
            #TODO: threshold, reset, leak -> supposed to be neuron specific,
            #  make new group for each setting type?
            sanafe_neuron = group.define_neuron(j,
                {"threshold": threshold, "reset": reset, "decay": decay})
            index = sanafe_neuron.id
            name = neuron
            toadd = {name: index}
            nid_dict_idx[name] = i
            #neurons.remove((j, neuron))

            # Map neuron
            sim.map_neuron_to_compartment(arch, neuron)
            nid_dict.append(toadd)

        # Add synapses
        for i, synapse in enumerate(self.scaffold.graph.edges):
                #print(str(synapse) + " | " + str(self.scaffold.graph.edges[synapse]))
                group_idx_src = nid_dict_idx[synapse[0]]
                neuron_idx_src = nid_dict[group_idx_src][synapse[0]]

                group_idx_dest = nid_dict_idx[synapse[1]]
                neuron_idx_dest = nid_dict[group_idx_dest][synapse[1]]

                weight = self.scaffold.graph.edges[synapse]['weight']
                network.groups[group_idx_src].neurons[neuron_idx_src].add_connection(
                    network.groups[group_idx_dest].neurons[neuron_idx_dest],
                    {"weight": weight})

        # print(id_dict)
        # print(nid_dict)
        # print(nid_dict_idx)
        # self.scaffold.summary(verbose = 1)

        #save
        sim = sanafecpp.Simulation(arch, network)
        sim.run(n_steps)
        return

    def cleanup(self):
        """
        Deletes/frees neurons and synapses
        """
        pass

    def reset(self):
        """
        Resets time-step to 0 and resets neuron/synapse properties
        """
        self._build_network()
        pass

    def set_properties(self, properties={}):
        """
        Set properties for specific neurons and synapses
        Args:
            properties: dictionary of parameter for bricks

        Example:
           for brick in properties:
               neuron_props, synapse_props = self.circuit[brick].get_changes(properties[brick])
               for neuron in neuron_props:
                   set neuron properties
               for synapse in synapse_props:
                   set synapse properties

        @NOTE: Currently, this function behaves differently for Input Bricks
           * Instead of returning the changes, they change internally and reset the iterator
           * This is because of how initial spike times are calculated using said bricks
           * I have not yet found a way of incorporating my proposed method (above) into these bricks yet
        """
        pass

    def set_input_spikes(self):
        pass
