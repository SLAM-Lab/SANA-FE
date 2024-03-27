"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

sim.py - Simulator script and utility functionality
"""
import sys
import os
import yaml

NETWORK_FILENAME = "runs/connected_layers.net"
ARCH_FILENAME = "loihi.arch"

### SNN utility functions ###
class Architecture:
    def __init__(self, arch=None):
        self._arch_dict = arch
        # TODO: parameterize
        self.max_compartments = 1024
        total_cores = 128

        self.core_to_address = []
        for core in range(0, total_cores):
            tile, offset = (core // 4), (core % 4)
            self.core_to_address.append((tile, offset))
            self.compartments = [self.max_compartments] * total_cores


class Network:
    def __init__(self, save_mappings=False):
        self.inputs = []
        self.groups = []
        # TODO: parameterize
        self._save_mappings = save_mappings

    def create_group(self, threshold, reset, leak, log_spikes=False,
                     log_potential=False, force_update=False,
                     connections_out=None, reverse_threshold=None,
                     reverse_reset_mode=None, soma_hw_name=None,
                     synapse_hw_name=None):
        group_id = len(self.groups)
        group = NeuronGroup(group_id, threshold, reset, leak, log_spikes,
                            log_potential, force_update, connections_out,
                            reverse_threshold, reverse_reset_mode,
                            soma_hw_name=soma_hw_name,
                            synapse_hw_name=synapse_hw_name)
        self.groups.append(group)
        return group

    def create_input(self):
        input_id = len(self.inputs)
        input_node = Input(input_id)
        self.inputs.append(input_node)
        return input_node

    def save(self, filename, group_idx=None, save_mappings=None):
        if group_idx is None:
            group_idx = slice(0, len(self.groups))
        if save_mappings is None:
            save_mappings = self._save_mappings

        with open(filename, 'w') as network_file:
            for group in self.groups[group_idx]:
                network_file.write(str(group))

            for group in self.groups[group_idx]:
                for neuron in group.neurons:
                    neuron._save_mappings = save_mappings
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
                 soma_hw_name=None, synapse_hw_name=None):
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
        self.soma_hw_name = soma_hw_name
        self.synapse_hw_name = synapse_hw_name

    def __str__(self):
        neuron_count = len(self.neurons)

        group_str = f"g {neuron_count}"
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
        if self.soma_hw_name is not None:
            group_str += f" soma_hw_name={self.soma_hw_name}"
        if self.synapse_hw_name is not None:
            group_str += f" synapse_hw_name={self.synapse_hw_name}"

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


def map_neuron_to_compartment(arch, core=None):
    if core is not None:
        assert(core < len(arch.compartments) and arch.compartments[core] > 0)
        arch.compartments[core] -= 1
        return core

    for curr, _ in enumerate(arch.compartments):
        if arch.compartments[curr] > 0:
            arch.compartments[curr] -= 1
            return core

    # No free compartments left
    return None, None


def map_neuron_group_to_cores(neuron_group, arch, core_count):
    return map_neurons_to_cores(neuron_group.neurons, arch, core_count)


def map_neurons_to_cores(neurons, arch, core_count):
    """Map neurons to one or more cores"""
    neuron_count = len(neurons)
    neurons_per_core = neuron_count // core_count
    print(f"Mapping {neuron_count} neurons to {core_count} cores")

    # Map N sequential cores. Find the first empty core
    first_unused_core = arch.compartments.index(arch.max_compartments)

    # Create a list of the cores to map to (absolute ID)
    cores = list(range(first_unused_core, (first_unused_core+core_count)))

    # Go through list and assign compartments of each entry
    map_neurons = [neurons_per_core] * (core_count)
    map_neurons[-1] += (neuron_count - (neurons_per_core*core_count))
    assert(sum(map_neurons) == neuron_count)

    mapping = dict(zip(cores, map_neurons))
    #print(f"Mapping dictionary: {mapping}")

    curr_neuron = 0
    for core in mapping.keys():
        # For the last core, map all remaining neurons
        for _ in range(0, mapping[core]):
            n = neurons[curr_neuron]
            mapped_core = map_neuron_to_compartment(arch,
                                                    core=core)
            n.tile, n.core = arch.core_to_address[mapped_core]
            curr_neuron += 1

    assert(curr_neuron == neuron_count)
    return mapping


def create_connected_layer(network, prev_layer, input_shape, weights,
                           log_spikes=False, log_potential=False,
                           force_update=False, threshold=1.0, reset=0.0,
                           leak=1.0, reverse_threshold=None,
                           reverse_reset_mode=None, soma_hw_name=None,
                           synapse_hw_name=None):
    layer_neuron_count = weights.shape[1]
    layer = create_layer(network, layer_neuron_count,
                         log_spikes=log_spikes, log_potential=log_potential,
                         force_update=force_update,
                         threshold=threshold, reset=reset, leak=leak,
                         connections_out=None,
                         reverse_threshold=reverse_threshold,
                         reverse_reset_mode=reverse_reset_mode,
                         soma_hw_name=soma_hw_name,
                         synapse_hw_name=synapse_hw_name)

    for src in prev_layer.neurons:
        for dest in layer.neurons:
            # Take the ID of the neuron in the 2nd layer
            weight = weights[src.id, dest.id]
            src.add_connection(dest, weight)

    return layer


def create_conv_layer(network, input_layer, input_shape, filters,
                      stride=1, pad=0, log_spikes=False, log_potential=False,
                      force_update=False, threshold=1.0, reset=0.0,
                      leak=1.0, reverse_threshold=None,
                      reverse_reset_mode=None, soma_hw_name=None,
                      synapse_hw_name=None):

    # Calculate dimensions of this layer
    input_width = input_shape[0]
    input_height = input_shape[1]
    input_channels = input_shape[2]

    kernel_width = filters.shape[0]
    kernel_height = filters.shape[1]
    kernel_channels = filters.shape[2]
    kernel_count = filters.shape[3]

    output_width = ((input_width + (2*pad) - kernel_width) // stride) + 1
    output_height = ((input_height + (2*pad) - kernel_height) // stride) + 1

    output_channels = kernel_count
    assert(kernel_channels == input_channels)

    layer_neuron_count = output_channels * output_width * output_height
    output_layer = create_layer(network, layer_neuron_count,
                                log_spikes=log_spikes,
                                log_potential=log_potential,
                                force_update=force_update,
                                threshold=threshold, reset=reset, leak=leak,
                                connections_out=None,
                                reverse_threshold=reverse_threshold,
                                reverse_reset_mode=reverse_reset_mode,
                                soma_hw_name=soma_hw_name,
                                synapse_hw_name=synapse_hw_name)

    # Create the convolutional connections
    for c_out in range(0, output_channels):
        for y_out in range(0, output_height):
            for x_out in range(0, output_width):
                dest_idx = c_out * output_width * output_height
                dest_idx += y_out * output_width
                dest_idx += x_out

                dest = output_layer.neurons[dest_idx]
                for c_in in range(0, input_channels):
                    for y_kernel in range(0, kernel_height):
                       if not 0 <= y_out*stride + y_kernel < input_height:
                            continue
                       for x_kernel in range(0, kernel_width):
                            if not 0 <= x_out*stride + x_kernel < input_width:
                                continue

                            # Calculate the index of the src neuron
                            src_idx = c_in * input_width * input_height
                            src_idx += ((y_out*stride) + y_kernel)*input_width
                            src_idx += ((x_out*stride) + x_kernel)
                            src = input_layer.neurons[src_idx]

                            weight = filters[y_kernel, x_kernel, c_in, c_out]
                            src.add_connection(dest, weight)

    return output_layer


def create_layer(network, layer_neuron_count,
                 log_spikes=False, log_potential=False, force_update=False,
                 threshold=1.0, reset=0.0, leak=1.0,
                 connections_out=None, reverse_threshold=None,
                 reverse_reset_mode=None, soma_hw_name=None,
                 synapse_hw_name=None, biases=None):
    print("Creating layer with {0} neurons".format(layer_neuron_count))
    layer_group = network.create_group(threshold, reset, leak, log_spikes,
                                       log_potential, force_update,
                                       connections_out=connections_out,
                                       reverse_threshold=reverse_threshold,
                                       reverse_reset_mode=reverse_reset_mode,
                                       soma_hw_name=soma_hw_name,
                                       synapse_hw_name=synapse_hw_name)

    for i in range(0, layer_neuron_count):
        if (i % 1000) == 0:
            print(f"Creating neuron {i}")
        layer_group.create_neuron()

    if biases is not None:
        assert(len(biases) == layer_neuron_count)
        for bias, neuron in zip(biases, layer_group.neurons):
            neuron.add_bias(bias)

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
        run_dir=os.path.join(project_dir, "runs"), perf_trace=False,
        spike_trace=False, potential_trace=False, message_trace=False,
        run_alive=False):
    parsed_filename = os.path.join(run_dir,
                                   os.path.basename(arch_path) + ".parsed")
    parse_file(arch_path, parsed_filename)

    # Set some flags for the dynamic linking library
    # Important to do before importing the simcpp .so library!
    sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)
    import simcpp as sim

    # Parse inputs and run simulation
    sana_fe = sim.SANA_FE()
    if spike_trace:
        sana_fe.set_spike_flag()
    if potential_trace:
        sana_fe.set_pot_flag()
    if perf_trace:
        sana_fe.set_perf_flag()
    if message_trace:
        sana_fe.set_mess_flag()
    
    sana_fe.set_arch(parsed_filename)
    sana_fe.set_net(network_path)

    if run_alive:
        while True:
            if timesteps > 0:
                print('-----------Inter Run Summary-----------')
                sana_fe.run_timesteps(timesteps)
                sana_fe.run_summary()
                print('---------------------------------------')
            print("Enter timesteps to run: ", end="")
            user_in = input()

            if user_in == "q" or user_in == "quit":
                break
            if user_in.startswith("u"):
                try:
                    group_id = int(user_in[2])
                except ValueError:
                    print(f"Error: Expected int. Got \"{user_in[2]}\".")
                    exit(1)

                try:
                    n_id = int(user_in[4])
                except ValueError:
                    print(f"Error: Expected int. Got \"{user_in[4]}\".")
                    exit(1)

                user_in = user_in[6:]
                kwargs = user_in.split(" ")
                
                print(group_id, n_id, kwargs, len(kwargs))
                sana_fe.update_neuron(group_id, n_id, kwargs, len(kwargs))

                timesteps = 0
                continue
            try:
                timesteps = int(user_in)
            except ValueError:
                print(f"Error: Expected int. Got {user_in}.")
                exit(1)
    else:
        if timesteps < 1:
            print(f"Error: Given {timesteps} timesteps, require int > 1.")
            exit(1)
        sana_fe.run_timesteps(timesteps)
    
    print('-----------Total Run Summary-----------')
    sana_fe.sim_summary()

    with open("run_summary.yaml", "r") as run_summary:
        results = yaml.safe_load(run_summary)

    sana_fe.clean_up()

    return results


if __name__ == "__main__":
    # Run SANA-FE from the command-line
    import argparse

    parser = argparse.ArgumentParser(
        prog="python sim.py",
        description="Simulating Advanced Neuromorphic Architectures for Fast Exploration"
    )
    parser.add_argument("architecture", help="Architecture description (YAML) file path", type=str)
    parser.add_argument("snn", help="Spiking Neural Network description file path", type=str)
    parser.add_argument("timesteps", help="Number of timesteps to simulate", type=int)
    parser.add_argument("-s", "--spikes", help="Trace spikes", action="store_true")
    parser.add_argument("-v", "--voltages", help="Trace neuron voltages", action="store_true")
    parser.add_argument("-p", "--perf", help="Trace perf", action="store_true")
    parser.add_argument("-m", "--messages", help="Trace messages", action="store_true")
    parser.add_argument("-v", "--voltages", help="Trace membrane voltages", action="store_true")
    parser.add_argument("-r", "--run", help="Keep simulation alive", action="store_true")

    args = parser.parse_args()
    print(args)

    run(args.architecture, args.snn, args.timesteps, spike_trace=args.spikes,
        potential_trace=args.voltages, perf_trace=args.perf,
        message_trace=args.messages, run_alive=args.run)
    print("sim finished")
