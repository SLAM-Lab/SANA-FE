# Copyright (c) 2025 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
SANA-FE back-end for the Fugu framework

A SANA-FE back-end for executing applications created in Fugu, using scaffolds
and bricks. SANA-FE supports a simple mapping of bricks to hardware, with
basic spike and voltage probes.
"""

# TODO: support recordInGraph
# TODO: read from arch object in the future
from fugu.backends import Backend
from collections import defaultdict
import sanafe
import pandas as pd

class sanafe_Backend(Backend):
    _net = None
    _arch = None

    def _map_to_cores(self):
        """
        Assigns neurons to cores based on capacity and hardware compatibility.
        """
        MAX_NEURONS_PER_CORE = 1024
        cores = self.arch.cores()
        total_cores = len(cores)
        neurons_per_core = {i: 0 for i in range(total_cores)}

        current_core_id = 0
        # We assume self.node_map contains {fugu_id: sanafe_neuron_object}
        for fugu_node_id, neuron in self.node_map.items():
            if neurons_per_core[current_core_id] >= MAX_NEURONS_PER_CORE:
                # There is space in the core
                current_core_id += 1
                assert(current_core_id < total_cores)
            neurons_per_core[current_core_id] += 1

            # Assign spike inputs to custom hardware
            if fugu_node_id in self.input_map:
                input_id = neurons_per_core[current_core_id] - 1
                neuron.set_attributes(soma_hw_name=f"loihi_inputs[{input_id}]")

            # Perform the actual mapping in SANA-FE
            target_core = cores[current_core_id]
            neuron.map_to_core(target_core)


    def _convert_props(self, fugu_props):
        """
        Convert Fugu properties to equivalent SANA-FE attributes
        """
        param_map = {
            "decay": "leak_decay",
            "reset_voltage": "reset",
        }
        sanafe_props = {param_map.get(k, k): v for k, v in fugu_props.items()}
        del sanafe_props["index"]
        del sanafe_props["brick"]
        del sanafe_props["neuron_number"]
        del sanafe_props["p"]

        return sanafe_props

    def _build_network(self):
        """
        Build a Fugu SNN in SANA-FE
        """
        self.net = sanafe.Network()
        self.node_map = {}
        self.fugu_name_to_neuron_number = {}
        self.input_map = set()

        # --- STEP 1: Sort Neurons by Brick ---
        # We create a dictionary where Key = Brick Tag, Value = List of Neuron IDs
        self.brick_groups = defaultdict(list)
        neurons_to_record = set()
        input_neurons = set()
        record_all = True if self.record == "all" else False

        # --- STEP 2: Create input spike trains and outputs probes ---
        for brick_id, props in self.fugu_circuit.nodes.data():
            if props.get("layer") == "input":
                for timestep, neurons in enumerate(props['brick']):
                    for n in neurons:
                        node = self.fugu_graph.nodes[n]
                        if "spikes" not in node:
                            node["spikes"] = []
                        spikes = node["spikes"]
                        spikes.append(timestep + 1)
                        input_neurons.add(n)
            elif props.get("layer") == "output":
                if self.debug_mode:
                    print(f"Found Output Brick: {props.get('name')}")
                # Dig into the ports to find the specific data neurons
                # 'ports' maps port_name -> PortData object
                if "ports" in props:
                    for port in props["ports"].values():
                        # 'channels' maps channel_name ('data', 'complete') -> ChannelData
                        if "data" in port.channels:
                            # Get the list of low-level Fugu names (e.g., 'AND_0')
                            data_neurons = port.channels["data"].neurons
                            neurons_to_record.update(data_neurons)

        # --- STEP 2: Build up list of Fugu nodes and their unique names ---
        for n, props in self.fugu_graph.nodes.data():
            # Get the unique tag (e.g., 'Brick-2', 'Input0-0')
            # Default to 'Misc' if a neuron somehow has no brick parent
            brick_tag = props.get("brick", "Misc")
            self.brick_groups[brick_tag].append(n)
            self.fugu_name_to_neuron_number[n] = props.get("neuron_number")
            if record_all:
                # If recording all neurons, add this to this list in addition
                #  to any specific output nodes
                neurons_to_record.add(n)

        # --- STEP 3: Create SANA-FE groups ---
        for brick_tag, neuron_list in self.brick_groups.items():
            # Lookup the "Pretty Name" (e.g., 'AND') instead of 'Brick-2'
            # The high-level info is stored in fugu_circuit using the tag as the key
            group_name = brick_tag
            if brick_tag in self.fugu_circuit.nodes:
                # Get the human-readable name you gave it in the script
                group_name = self.fugu_circuit.nodes[brick_tag].get('name', brick_tag)

            # Create the neuron group in SANA-FE
            if self.debug_mode:
                print(f"Creating Group: {group_name} with {len(neuron_list)} neurons")
            sanafe_group = self.net.create_neuron_group(group_name, len(neuron_list), {})

            # Configure neurons & map IDs
            for i, fugu_node_id in enumerate(neuron_list):
                fugu_props = self.fugu_graph.nodes[fugu_node_id]

                # Remap properties which have different names in Fugu vs SANA-FE
                sanafe_props = self._convert_props(fugu_props)

                if fugu_node_id in input_neurons:
                    self.input_map.add(fugu_node_id)
                sanafe_group[i].set_attributes(model_attributes=sanafe_props)

                if fugu_node_id in neurons_to_record:
                    # Enable logging for this neuron
                    #if self.debug_mode:
                    print(f"Enabling logging for neuron {fugu_node_id}")
                    sanafe_group[i].set_attributes(log_spikes=True,
                                                   log_potential=True)

                # Map the Fugu node ID to the corresponding SANA-FE neuron
                #  This ensures edges can find the correct destination next
                self.node_map[fugu_node_id] = sanafe_group[i]

        for n1, n2, props in self.fugu_graph.edges.data():
            if n1 in self.node_map and n2 in self.node_map:
                src = self.node_map[n1]
                dst = self.node_map[n2]
                src.connect_to_neuron(dst, props)


    def compile(self, scaffold, compile_args={}):
        """
        Compile Fugu scaffold as SANA-FE network
        Args:
            scaffold: The Fugu scaffold containing one or more bricks
            compile_args (Dict, optional): Arguments controlling compilation. Defaults to {}.
        """
        # creates neuron populations and synapses
        self.scaffold = scaffold
        self.fugu_circuit = scaffold.circuit
        self.fugu_graph = scaffold.graph
        self.brick_to_number = scaffold.brick_to_number
        self.recordInGraph = "recordInGraph" in compile_args
        if "record" in compile_args:
            self.record = compile_args["record"]
        else:
            self.record = False
        if "ds_format" in compile_args:
            self.ds_format = compile_args["ds_format"]
        else:
            self.ds_format = True
        if "debug_mode" in compile_args:
            self.debug_mode = compile_args["debug_mode"]
        else:
            self.debug_mode = False
        if "arch" in compile_args:
            self.arch_name = compile_args["arch"]

        self._build_network()

    def _find_node(name):
        return None

    def run(self, n_steps, return_potentials=False, debug_mode=False):
        """
        Convert to SANA-FE network and save to file
        Args:
            n_steps (int): Number of timesteps to simulate.
            return_potentials (bool, optional): Log neuron potentials. Defaults to False.
            debug_mode (bool, optional): Enable for verbose prints. Defaults to False.
        """

        self.return_potentials = return_potentials
        if debug_mode:
            self.scaffold.summary(verbose=1)
        self.arch = sanafe.load_arch(self.arch_name)
        self._map_to_cores()

        # Run simulation
        chip = sanafe.SpikingChip(self.arch)
        chip.load(self.net)

        with open("spikes.csv", "w") as spike_trace, open("potentials.csv", "w") as potential_trace:
            chip.sim(n_steps, spike_trace=spike_trace, potential_trace=potential_trace)

        # Parse output data
        spikes_out_df = pd.read_csv("spikes.csv")
        potentials_out_df = pd.read_csv("potentials.csv")

        if self.recordInGraph:
            for n, node in self.fugu_graph.nodes.data():
                if not "outputs" in node:
                    continue
                # TODO: read the spike trace and attach spiking info to each
                #  node in the Fugu graph
            return True

        # Otherwise, reformat the spike dataframe to match what Fugu expects
        spikes_out_df["time"] = spikes_out_df["timestep"] - 1.0
        # Remap SANA-FE neuron IDs to Fugu neuron numbers
        neuron_numbers = []
        for n in spikes_out_df["neuron"]:
            brick_tag, offset_str = n.split('.')
            offset = int(offset_str)
            # Get the Fugu node name
            brick = self.brick_groups[brick_tag]
            fugu_name = brick[offset]
            # Map neuron ID to number
            neuron_numbers.append(self.fugu_name_to_neuron_number[fugu_name])
        spikes_out_df["neuron_number"] = neuron_numbers

        spikes_out_df = spikes_out_df.drop(columns=["neuron", "timestep"])
        if not self.return_potentials:
            return spikes_out_df
        return spikes_out_df, potentials_out_df

    def cleanup(self):
        """
        Deletes/frees neurons and synapses
        """
        del self.brick_groups
        del self.fugu_name_to_neuron_number
        del self.net
        del self.arch

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
            properties (dict, optional): dictionary of properties for bricks. Defaults to {}.
        """
        for brick in properties:
            if brick != 'compile_args':
                brick_id = self.brick_to_number[brick]
                self.fugu_circuit.nodes[brick_id]['brick'].set_properties(properties[brick])
        # must call run() for changes to take effect

    def set_input_spikes(self):
        """
        Reset input spike trains
        """
        # Clean out old spike structures.
        for n, node in self.fugu_graph.nodes.data():
            if 'spikes' in node:
                del node['spikes']  # Allow list to be built from scratch.
        # When run() is called, network will be rebuilt with new spike times.
