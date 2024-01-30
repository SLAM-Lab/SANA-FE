"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

sim.py - Simulator script and utility functionality
"""
import subprocess
import os
import yaml

NETWORK_FILENAME = "runs/connected_layers.net"
ARCH_FILENAME = "loihi.arch"

### SNN utility functions ###

class Network:
    def __init__(self, external_inputs=0, save_mappings=False):
        self.external_inputs = external_inputs
        self.inputs = []
        self.groups = []
        self._max_tiles = None
        self._max_cores = None
        self._max_compartments = None
        self._save_mappings = save_mappings

    def create_group(self, threshold, reset, leak, log_spikes=False,
                     log_potential=False, force_update=False,
                     connections_out=None, reverse_threshold=None,
                     reverse_reset_mode=None, neuron_model=None,
                     default_synapse_model=None):
        group_id = len(self.groups)
        group = NeuronGroup(group_id, threshold, reset, leak, log_spikes,
                            log_potential, force_update, connections_out,
                            reverse_threshold, reverse_reset_mode,
                            neuron_model=neuron_model,
                            default_synapse_model=default_synapse_model)
        self.groups.append(group)
        return group

    def create_input(self):
        input_id = len(self.inputs)
        input_node = Input(input_id)
        self.inputs.append(input_node)
        return input_node

    def save(self, filename, group_idx=None):
        if group_idx is None:
            group_idx = slice(0, len(self.groups))
        with open(filename, 'w') as network_file:
            if self.external_inputs > 0:
                network_file.write("x {0} rate\n".format(self.external_inputs))
            for group in self.groups[group_idx]:
                network_file.write(str(group))

            for group in self.groups[group_idx]:
                for neuron in group.neurons:
                    neuron._save_mappings = self._save_mappings
                    network_file.write(str(neuron))

            for input_node in self.inputs:
                network_file.write(str(input_node))

    def load(self, filename):
        with open(filename, 'r') as network_file:
            for line in network_file:
                fields = line.split()
                if fields and fields[0] == 'g':
                    # TODO: support other fields to be loaded
                    neuron_count = int(fields[1])
                    group = self.create_group(0.0, 0.0, 0)
                    for _ in range(0, neuron_count):
                        group.create_neuron()

                elif fields and fields[0] == 'n':
                    pass
                    #neuron_address = fields[1]
                    #gid = int(neuron_address.split('.')[0])
                    #group = self.groups[gid]
                    #group.create_neuron()
                elif fields and fields[0] == 'e':
                    edge_info = fields[1]
                    src_address = edge_info.split("->")[0]
                    dest_address = edge_info.split("->")[1]

                    src_gid = int(src_address.split(".")[0])
                    src_nid = int(src_address.split(".")[1])
                    src = self.groups[src_gid].neurons[src_nid]

                    dest_gid = int(dest_address.split(".")[0])
                    dest_nid = int(dest_address.split(".")[1])
                    dest = self.groups[dest_gid].neurons[dest_nid]

                    weight = None
                    for f in fields:
                        if "w=" in f or "weight=" in f:
                            weight = float(f.split("=")[1])
                    src.add_connection(dest, weight)
                elif fields and fields[0] == '&':
                    # TODO: for now ignore the mappings, the whole reason I'm
                    #  trying this code is to explore different mappings
                    pass

        """
        for g in self.groups:
            for n in g.neurons:
                print(n)
        """
        return


class NeuronGroup:
    def __init__(self, group_id, threshold, reset, leak, log_spikes=None,
                 log_potential=None, force_update=None, connections_out=None,
                 reverse_reset=None, reverse_reset_mode=None,
                 neuron_model=None, default_synapse_model=None):
        # TODO: support all features here
        self.id = group_id
        self.neurons = []
        self.threshold = threshold
        self.reset = reset
        self.reset_mode = None
        self.reverse_reset = reverse_reset
        self.reverse_reset_mode = reverse_reset_mode
        self.reverse_threshold = None
        self.leak_decay = leak
        self.connections_out = connections_out
        self.log_spikes = log_spikes
        self.log_potential = log_potential
        self.force_update = force_update
        self.neuron_model = neuron_model
        self.default_synapse_model = default_synapse_model

    def __str__(self):
        neuron_count = len(self.neurons)

        group_str = f"g {neuron_count}"
        if self.neuron_model is not None:
            group_str += f" soma_model={self.neuron_model}"
        if self.threshold is not None:
            group_str += f" threshold={self.threshold}"
        if self.reset is not None:
            group_str += f" reset={self.reset}"
        if self.reverse_threshold is not None:
            group_str += f" reverse_threshold={self.reverse_threshold}"
        if self.reverse_reset is not None:
            group_str += f" reverse_reset={self.reverse_reset}"
        if self.leak_decay is not None:
            group_str += f" leak_decay={self.leak_decay}"
        if self.reset_mode is not None:
            group_str += f" reset_mode={self.reset_mode}"
        if self.reverse_reset_mode is not None:
            group_str += f" reverse_reset_mode={self.reverse_reset_mode}"
        if self.log_spikes is not None:
            group_str += f" log_spikes={int(self.log_spikes)}"
        if self.log_potential is not None:
            group_str += f" log_v={int(self.log_potential)}"
        if self.force_update is not None:
            group_str += f" force_update={int(self.force_update)}"
        if self.connections_out is not None:
            group_str += f" connections_out={self.connections_out}"
        if self.default_synapse_model is not None:
            group_str += f" synapse_model={self.default_synapse_model}"

        group_str += "\n"
        return group_str

    def create_neuron(self, log_spikes=None, log_potential=None,
                      force_update=None):
        neuron_id = len(self.neurons)
        neuron = Neuron(self, neuron_id, log_spikes=log_spikes,
                        log_potential=log_potential, force_update=force_update)
        self.neurons.append(neuron)
        return neuron


class Input:
    def __init__(self, input_id):
        self.id = input_id
        self.connections = []

    def add_connection(self, dest, weight):
        self.connections.append((dest, weight))

    def __str__(self):
        line = "< {0}".format(self.id)
        for connection in self.connections:
            dest_neuron, weight = connection
            line += " {0} {1} {2}".format(
                dest_neuron.group_id, dest_neuron.id, weight)
        line += '\n'
        return line


class Neuron:
    def __init__(self, group, neuron_id, log_spikes=None,
                 log_potential=None, force_update=None):
        self.group = group
        self.id = neuron_id
        self.log_spikes = log_spikes
        self.log_potential = log_potential
        self.force_update = force_update
        self.connections = []
        self.tile = None
        self.core = None
        self.bias = None
        self._save_mappings = False

    def add_connection(self, dest, weight):
        self.connections.append((dest, weight))

    def add_bias(self, bias):
        self.bias = bias

    def __str__(self, map_neuron=True):
        neuron_str = f"n {self.group.id}.{self.id}"
        if self.bias is not None:
            neuron_str += f" bias={self.bias}"
        if self.log_spikes is not None:
            self.log_spikes = int(self.log_spikes)
            neuron_str += f" log_spikes={int(self.log_spikes)}"
        if self.log_potential is not None:
            neuron_str += f" log_v={int(self.log_potential)}"
        if self.force_update is not None:
            neuron_str += f" force_update={int(self.force_update)}"
        if (self.group.connections_out is None or (self.connections and
            (len(self.connections) > self.group.connections_out))):
           neuron_str += f" connections_out={len(self.connections)}"
        neuron_str += "\n"

        for connection in self.connections:
            dest_neuron, weight = connection
            neuron_str += f"e {self.group.id}.{self.id}->"
            neuron_str += f"{dest_neuron.group.id}.{dest_neuron.id}"
            if isinstance(weight, float):
                neuron_str += f" w={weight:.5e}"
            else:
                neuron_str += f" w={weight}"
            neuron_str += "\n"

        if self._save_mappings:
            neuron_str += "& {0}.{1}@{2}.{3}\n".format(self.group.id, self.id,
                                                    self.tile, self.core)
        return neuron_str


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
                 log_spikes=False, log_potential=False, force_update=False,
                 threshold=1.0, reset=0.0, leak=1.0, mappings=None,
                 connections_out=None, reverse_threshold=None,
                 reverse_reset_mode=None, neuron_model=None,
                 synapse_model=None):
    print("Creating layer with {0} neurons".format(layer_neuron_count))
    layer_group = network.create_group(threshold, reset, leak, log_spikes,
                                       log_potential, force_update,
                                       connections_out=connections_out,
                                       reverse_threshold=reverse_threshold,
                                       reverse_reset_mode=reverse_reset_mode,
                                       neuron_model=neuron_model,
                                       default_synapse_model=synapse_model)

    if mappings is not None:
        assert(len(mappings) == layer_neuron_count)

    for i in range(0, layer_neuron_count):
        if (i % 10000) == 0:
            print(f"Creating neuron {i}")
        neuron = layer_group.create_neuron()

        if mappings is not None:
            tile, core = mappings[i]
        else:
            tile, core = map_neuron_to_compartment(compartments)
        neuron.tile, neuron.core = tile, core

    return layer_group


### Architecture description parsing ###
def parse_file(input_filename, output_filename):
    with open(input_filename, "r") as arch_file:
        arch_dict = yaml.safe_load(arch_file)

    if "architecture" not in arch_dict:
        raise Exception("Error: no architecture defined")

    parse_arch(arch_dict["architecture"])
    arch_elements = _entry_list

    with open(output_filename, "w") as list_file:
        for line in arch_elements:
            list_file.write(line + '\n')
    return


def parse_arch(arch):
    global _tiles
    global _cores_in_tile
    global _entry_list

    _tiles = 0
    _cores_in_tile = []
    _entry_list = []

    arch_name = arch["name"]
    if "[" in arch_name:
        raise Exception("Error: multiple architectures not supported")

    if "tile" not in arch:
        raise Exception("Error: No tiles defined, must be at least one tile")

    tiles = arch["tile"]
    for tile in tiles:
        parse_tile(tile)

    create_noc(arch)
    return


def parse_tile(tile_dict):
    tile_name = tile_dict["name"]
    # Work out how many instances of this tile to create
    if "[" in tile_name:
        # Can use notation [min..max] to indicate range of elements
        range_min, range_max = parse_range(tile_name)
    else:
        range_min, range_max = 0, 0

    for instance in range(range_min, range_max+1):
        tile_name = tile_name.split("[")[0] + "[{0}]".format(instance)
        tile_id = create_tile(tile_dict, tile_name)

        # Add any elements local to this h/w structure. They have access to any
        #  elements in the parent structures
        if "core" not in tile_dict:
            raise Exception("Error: No cores defined, "
                            "must be at least one core")
        cores = tile_dict["core"]
        for _, core_dict in enumerate(cores):
            parse_core(core_dict, tile_id)

    return

def parse_core(core_dict, tile_id):
    core_name = core_dict["name"]

    # Work out how many instances of this tile to create
    if "[" in core_name:
        # Can use notation [min..max] to indicate range of elements
        range_min, range_max = parse_range(core_name)
    else:
        range_min, range_max = 0, 0

    elements = ("axon_in", "synapse", "dendrite", "soma", "axon_out")
    for el in elements:
        if el not in core_dict:
            raise Exception("Error: {0} not defined in core {1}".format(
                el, core_name))

    for instance in range(range_min, range_max+1):
        core_name = core_name.split("[")[0] + "[{0}]".format(instance)
        core_id = create_core(tile_id, core_name, core_dict)
        create_axon_in(tile_id, core_id, core_dict["axon_in"][0])

        for synapse in core_dict["synapse"]:
            create_synapse(tile_id, core_id, synapse)
        create_dendrite(tile_id, core_id, core_dict["dendrite"][0])
        for soma in core_dict["soma"]:
            create_soma(tile_id, core_id, soma)
        create_axon_out(tile_id, core_id, core_dict["axon_out"][0])


def parse_range(range_str):
    range_str = range_str.replace("]", "")
    range_str = range_str.split("[")[1]
    range_min = int(range_str.split("..")[0])
    range_max = int(range_str.split("..")[1])

    return range_min, range_max


def get_instances(element_dict):
    element_name = element_dict["name"]

    if "[" in element_name:
        range_min, range_max = parse_range(element_name)
        instances = (range_max - range_min) + 1
    else:
        instances = 1

    return instances


_tiles = 0
_cores_in_tile = []
_entry_list = []


def format_attributes(attributes):
    line = ""
    if attributes is None:
        attributes = {}

    for key in attributes:
        line += (f" {key}={attributes[key]}")
    return line


def create_tile(tile, name):
    global _tiles
    tile_id = _tiles
    _tiles += 1

    tile = f"t {name}" + format_attributes(tile["attributes"])
    _entry_list.append(tile)
    # Track how many cores are in this tile
    _cores_in_tile.append(0)

    return tile_id


def create_core(tile_id, name, core_dict):
    core_id = _cores_in_tile[tile_id]
    core = f"c {name} {tile_id}" + format_attributes(core_dict["attributes"])
    _entry_list.append(core)
    _cores_in_tile[tile_id] += 1

    return core_id


def create_synapse(tile_id, core_id, synapse_dict):
    name = synapse_dict["name"]
    synapse = f"s {name} {tile_id} {core_id}" + format_attributes(synapse_dict["attributes"])
    _entry_list.append(synapse)
    return


def create_dendrite(tile_id, core_id, dendrite_dict):
    name = dendrite_dict["name"]
    dendrite = f"d {name} {tile_id} {core_id}" + format_attributes(dendrite_dict["attributes"])
    _entry_list.append(dendrite)
    return


def create_soma(tile_id, core_id, soma_dict):
    name = soma_dict["name"]
    soma = (f"+ {name} {tile_id} {core_id}" + format_attributes(soma_dict["attributes"]))
    _entry_list.append(soma)
    return


def create_axon_in(tile_id, core_id, axon_dict):
    name = axon_dict["name"]
    axon = f"i {name} {tile_id} {core_id}" + format_attributes(axon_dict["attributes"])
    _entry_list.append(axon)
    return


def create_axon_out(tile_id, core_id, axon_dict):
    name = axon_dict["name"]
    axon = f"o {name} {tile_id} {core_id}" + format_attributes(axon_dict["attributes"])
    _entry_list.append(axon)
    return


def create_noc(noc_dict):
    if "attributes" not in noc_dict:
        raise Exception("Error: NoC not defined for architecture (please add this under attributes)")
    _entry_list.append("@" + format_attributes(noc_dict["attributes"]))
    return


project_dir = os.path.dirname(os.path.abspath(__file__))
def run(arch_path, network_path, timesteps,
        run_dir=os.path.join(project_dir, "runs"),
        perf_trace=True, spike_trace=False, potential_trace=False,
        message_trace=False):
    parsed_filename = os.path.join(run_dir,
                                   os.path.basename(arch_path) + ".parsed")
    parse_file(arch_path, parsed_filename)
    # Parse inputs and run simulation
    args = []
    if spike_trace:
        args.append("-s",)
    if potential_trace:
        args.append("-v")
    if perf_trace:
        args.append("-p")
    if message_trace:
        args.append("-m")
    command = [os.path.join(project_dir, "sim"),] + args + [parsed_filename,
               network_path, f"{timesteps}"]

    print("Command: {0}".format(" ".join(command)))
    ret = subprocess.call(command)
    if ret != 0:
        raise RuntimeError(f"Error: Simulator kernel failed (code={ret}).")

    with open("run_summary.yaml", "r") as run_summary:
        results = yaml.safe_load(run_summary)

    return results


if __name__ == "__main__":
    # Run SANA-FE from the command-line
    import argparse

    parser = argparse.ArgumentParser(prog="python sim.py",
                                    description="Simulating Advanced Neuromorphic Architectures for Fast Exploration")
    parser.add_argument("architecture", help="Architecture description (YAML) file path", type=str)
    parser.add_argument("snn", help="Spiking Neural Network description file path", type=str)
    parser.add_argument("timesteps", help="Number of timesteps to simulate", type=int)
    parser.add_argument("-s", "--spikes", help="Trace spikes", action="store_true")
    parser.add_argument("-v", "--voltages", help="Trace membrane voltages", action="store_true")

    args = parser.parse_args()
    print(args)

    run(args.architecture, args.snn, args.timesteps,
        spike_trace=args.spikes, potential_trace=args.voltages)
    print("sim finished")
