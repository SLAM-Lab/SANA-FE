"""
Copyright (c) 2024 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

sim.py - Simulator script and utility functionality
"""
import sys
import os
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, os.path.join(PROJECT_DIR))
# Set some flags for the dynamic linking library
# Important to do before importing the simcpp .so library!
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)
import sanafecpp

### SNN utility functions ###

# END OF DUPLICATED CODE

"""
class Architecture:
    def __init__(self, arch=None):
        self._arch = sanafecpp.Architecture()


class Network:
    def __init__(self, save_mappings=False):
        self._net = sanafecpp.Network()
        self._save_mappings = save_mappings

    def create_neuron_group(self, neuron_count, attributes):
        attributes =_attribute_values_to_str(attributes)
        group = NeuronGroup(self._net.create_neuron_group(
            neuron_count, attributes))
        return group

    def save(self, filename, group_idx=None, save_mappings=None):
        if group_idx is None:
            group_idx = slice(0, len(self.groups))
        if save_mappings is None:
            save_mappings = self._save_mappings

class NeuronGroup:
    def __init__(self, neuron_count, attributes):
        attributes =_attribute_values_to_str(attributes)
        self._neuron_group = sanafecpp.NeuronGroup(neuron_count, attributes)
        return

        group_str += "\n"
        return ""

    def define_neuron(self, neuron_id, attributes):
        attributes =_attribute_values_to_str(attributes)
        return self._neuron_group.define_neuron(neuron_id, attributes)


class Neuron:
    def __init__(self, neuron_id):
        self._neuron = sanafecpp.Neuron(neuron_id)
        self._save_mappings = False
        return

        if self._save_mappings:
            neuron_str += "& {0}.{1}@{2}.{3}\n".format(self.group.id, self.id,
                                                    self.tile, self.core)
        return ""
"""

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
            mapped_core = core.map_neuron(n)
            n.tile, n.core = arch.core_to_address[mapped_core]
            curr_neuron += 1

    assert(curr_neuron == neuron_count)
    return mapping


def create_connected_layer(network, prev_layer, input_shape, weights, attributes):
    layer_neuron_count = weights.shape[1]
    layer = create_layer(network, layer_neuron_count, attributes)

    for src in prev_layer.neurons:
        for dest in layer.neurons:
            # Take the ID of the neuron in the 2nd layer
            weight = weights[src.id, dest.id]
            src.connect_to_neuron(dest, {"w": weight})

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


def create_layer(net, layer_neuron_count, layer_attributes):
    print("Creating layer with {0} neurons".format(layer_neuron_count))
    layer_group = net.create_group(layer_neuron_count, layer_attributes)

    for i in range(0, layer_neuron_count):
        if (i % 10000) == 0:
            print(f"Defining neuron {i}")
        layer_group.define_neuron(i)

    return layer_group


### Architecture description parsing ###
def load_arch_from_yaml(yaml_filename):
    arch = sanafecpp.Architecture()
    with open(yaml_filename, "r") as arch_file:
        arch_dict = yaml.safe_load(arch_file)

    if "architecture" not in arch_dict:
        raise Exception("Error: no architecture section defined")

    parse_arch(arch, arch_dict["architecture"])

    return arch


def parse_arch(arch, arch_dict):
    arch_name = arch_dict["name"]
    if "[" in arch_name:
        raise Exception("Error: multiple architectures not supported")

    if "tile" not in arch_dict:
        raise Exception("Error: No tiles defined, must be at least one tile")

    tiles = arch_dict["tile"]
    for tile_dict in tiles:
        parse_tile(arch, tile_dict)

    print(arch_dict)
    arch.set_noc_attributes(arch_dict["attributes"])
    return


def parse_tile(arch, tile_dict):
    tile_name = tile_dict["name"]
    tile_name = tile_name.replace(" ", "_")
    tile_name = tile_name.replace("\t", "_")
    # Work out how many instances of this tile to create
    if "[" in tile_name:
        # Can use notation [min..max] to indicate range of elements
        range_min, range_max = parse_range(tile_name)
    else:
        range_min, range_max = 0, 0

    for instance in range(range_min, range_max+1):
        tile_name = tile_name.split("[")[0] + "[{0}]".format(instance)
        tile = arch.create_tile(tile_dict["attributes"])

        # Add any elements local to this h/w structure. They have access to any
        #  elements in the parent structures
        if "core" not in tile_dict:
            raise Exception("Error: No cores defined, "
                            "must be at least one core")
        cores = tile_dict["core"]
        for _, core_dict in enumerate(cores):
            parse_core(arch, tile.get_id(), core_dict)

    return


def parse_core(arch, tile_id, core_dict):
    core_name = core_dict["name"]
    core_name = core_name.replace(" ", "_")
    core_name = core_name.replace("\t", "_")

    # Work out how many instances of this tile to create
    if "[" in core_name:
        # Can use notation [min..max] to indicate range of elements
        range_min, range_max = parse_range(core_name)
    else:
        range_min, range_max = 0, 0

    elements = ("axon_in", "synapse", "dendrite", "soma", "axon_out")
    for el in elements:
        if el not in core_dict:
            raise Exception(f"Error: {el} not defined in core {core_name}")

    tile = arch.tiles[tile_id]
    for instance in range(range_min, range_max+1):
        core_name = core_name.split("[")[0] + "[{0}]".format(instance)
        core_attributes = core_dict.get("attributes", {})
        if core_attributes is None:
            core_attributes = {}
        core = arch.create_core(tile.get_id(), core_attributes)
        for axon_in in core_dict["axon_in"]:
            axon_in_attributes = axon_in.get("attributes", {})
            if axon_in_attributes is None:
                axon_in_attributes = {}
            core.create_axon_in(axon_in["name"], axon_in_attributes)
        for synapse in core_dict["synapse"]:
            synapse_attributes = synapse.get("attributes", {})
            if synapse_attributes is None:
                synapse_attributes = {}
            core.create_synapse(synapse["name"], synapse_attributes)

        #create_dendrite(core, core_dict["dendrite"][0])

        for soma in core_dict["soma"]:
            soma_attributes = soma.get("attributes", {})
            if soma_attributes is None:
                soma_attributes = {}
            print(f"soma:{soma_attributes}")
            core.create_soma(soma["name"], soma_attributes)
            #print(f"tile_id:{tile_id}")
            #print(f"core_id:{core.get_offset()}")
            #print(arch.tiles[tile_id].cores[core.get_offset()].soma)
            #print(core.soma)
            #print(id(core))
            #print(id(arch.tiles[tile_id].cores[core.get_offset()]))
            #exit(1)
        for axon_out in core_dict["axon_out"]:
            axon_out_attributes = soma.get("attributes", {})
            if axon_out_attributes is None:
                axon_out_attributes = {}
            core.create_axon_out(axon_out["name"], axon_out_attributes)
    return


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


## TODO: move these functions into the C++ kernel, this should be responsible
##  for loading and saving this raw format. The Python script should build
##  objects using the PyBind11 interface
def format_attributes(attributes):
    line = ""
    if attributes is None:
        attributes = {}

    for key in attributes:
        line += (f" {key}={attributes[key]}")
    return line

"""
def create_tile(tile, name):
    global _tiles
    tile_id = _tiles
    _tiles += 1

    if "attributes" not in tile:
        raise Exception(f"attributes section not defined for tile {tile_id}. "
                        "Did you mispell attributes?")

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
"""

project_dir = os.path.dirname(os.path.abspath(__file__))
def run_kernel(arch_path, network_path, timesteps,
        perf_trace=False, spike_trace=False,
        potential_trace=False, message_trace=False, out_dir=None):
    print("Loading architecture\n")
    try:
        arch = load_arch_from_yaml(arch_path)
        #arch = sanafecpp.Architecture()
        #arch.load_arch_file("runs/example.yaml.parsed")
    except yaml.parser.ParserError as yaml_parse_exc:
        print("Error: Invalid YAML file given.")
        if "expected '<document start>'" in str(yaml_parse_exc):
            print("Did you mix up the order of simulator files?")
            print("Run 'python3 sim.py' to check usage.")
    except yaml.scanner.ScannerError as yaml_parse_exc:
        print("Error: Problem parsing YAML file.")
        print(f"Full Error: '{yaml_parse_exc}'")
        exit()

    print("Loading network\n")
    #net = load_from_net_file(network_path, arch)
    net = sanafecpp.Network()
    net.load_net_file(network_path, arch)
    # Parse inputs and run simulation

    if out_dir is None:
        out_dir = "."

    print("Creating simulation")
    sim = sanafecpp.Simulation(arch, net, output_dir=out_dir,
                               record_spikes=spike_trace,
                               record_potentials=potential_trace,
                               record_perf=perf_trace,
                               record_messages=message_trace)

    sim.run(timesteps=timesteps, heartbeat=100)
    results = sim.get_run_summary()
    print(f"Results: {results}")

    return results


if __name__ == "__main__":
    # Run SANA-FE from the command-line
    import argparse

    parser = argparse.ArgumentParser(
        prog="python sim.py",
        description="Simulating Advanced Neuromorphic Architectures for Fast Exploration"
    )
    parser = argparse.ArgumentParser(
        prog="python sim.py",
        description="Simulating Advanced Neuromorphic Architectures for Fast Exploration"
    )
    parser.add_argument("architecture", help="Architecture description (YAML) file path", type=str)
    parser.add_argument("snn", help="Spiking Neural Network description file path", type=str)
    parser.add_argument("timesteps", help="Number of timesteps to simulate", type=int)
    parser.add_argument("-o", "--out_dir", help="Output directory", type=str, required=False)
    parser.add_argument("-s", "--spikes", help="Trace spikes", action="store_true")
    parser.add_argument("-v", "--voltages", help="Trace neuron voltages", action="store_true")
    parser.add_argument("-p", "--perf", help="Trace perf", action="store_true")
    parser.add_argument("-m", "--messages", help="Trace messages", action="store_true")

    args = parser.parse_args()
    try:
        run_kernel(args.architecture, args.snn, args.timesteps,
                   spike_trace=args.spikes, potential_trace=args.voltages,
                   perf_trace=args.perf, message_trace=args.messages,
                   out_dir=args.out_dir)
    except RuntimeError as exc:
        print(exc)
    else:
        print("sim finished")


    # Old code: replace this with the code below soon
    """
    import subprocess
    args = []
    if spike_trace:
        args.append("-s",)
    if potential_trace:
        args.append("-v")
    if perf_trace:
        args.append("-p")
    if message_trace:
        args.append("-m")
    if out_dir is not None:
        args.append("-o")
        args.append(f"{out_dir}")
    command = [os.path.join(project_dir, "sim"),] + args + [parsed_filename,
               network_path, f"{timesteps}"]

    print("Command: {0}".format(" ".join(command)))
    ret = subprocess.call(command)

    if out_dir is None:
        out_dir = ""
    summary_file = os.path.join(out_dir, "run_summary.yaml")
    with open(summary_file, "r") as run_summary:
        results = yaml.safe_load(run_summary)

    return results
    """
    # Capstone code for live demo - figure out how to integrate or split into
    #  a new script
    """
    # Code implemented during capstone project. Remove for now
    if run_alive:
        while True:
            if timesteps > 0:
                sim.run(timesteps, heartbeat=100)
                sim.get_run_summary()
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
                # print(group_id, n_id, kwargs, len(kwargs))
                #sim.update_neuron(group_id, n_id, kwargs, len(kwargs))

                timesteps = 0
                continue
            if user_in.startswith("s"):
                try:
                    group_id = int(user_in[2])
                except ValueError:
                    print(f"Error: Expected int. Got \"{user_in[2]}\".")
                    exit(1)

                #print(sim.get_status(group_id))

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
        sim.run_timesteps(timesteps)
    """