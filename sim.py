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
# Important to do before importing the sanafecpp .so library!
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)
import sanafecpp


### SNN utility functions ###
def map_neuron_group(neuron_group, arch, mappings):
    assert(len(neuron_group.neurons) == len(mappings))
    for n, m in zip(neuron_group.neurons, mappings):
        tile_id, core_offset = m
        # Look and find the tile and core objects
        core = arch.tiles[tile_id].cores[core_offset]
        core.map_neuron(n)
        #print(f"Mapping n:{n.get_id()} to "
        #      f"core {tile_id}.{core_offset}")


def map_neurons_to_cores(neurons, arch, core_count):
    """Map neurons uniformly to a given number of cores"""
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


def create_connected_layer(net, prev_layer, weights, attributes):
    import numpy as np
    layer_neuron_count = weights.shape[1]
    print(f"Creating new layer with {layer_neuron_count} neurons")
    layer = net.create_neuron_group(layer_neuron_count, attributes)

    connections = np.empty((layer_neuron_count * len(prev_layer.neurons)),
                           dtype=object)
    flattened_weights = weights.flatten()
    i = 0
    for src in range(0, len(prev_layer.neurons)):
        for dest in range(0, len(layer.neurons)):
            connections[i] = (src, dest)
            i += 1

    print("Connecting new layer with prev layer")
    prev_layer.connect_neurons(layer, connections, {"w": flattened_weights})
    return layer


def create_conv_layer(net, input_layer, input_shape, filters,
                      stride=1, pad=0, layer_attributes={}):
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
    output_layer = net.create_neuron_group(layer_neuron_count,
                                           layer_attributes)

    # Create the convolutional connections
    connections = []
    weights = []
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
                            connections.append((src.get_id(), dest.get_id()),)
                            weights.append(weight)

    input_layer.connect_neurons(output_layer, connections, {"w": weights})

    return output_layer


def load_net(net_filename, arch):
    net = sanafecpp.Network()
    net.load_net_description(net_filename, arch)
    return net


def load_arch_yaml(yaml_filename):
    arch = sanafecpp.Architecture()
    with open(yaml_filename, "r") as arch_file:
        arch_dict = yaml.safe_load(arch_file)

    if "architecture" not in arch_dict:
        raise Exception("Error: no architecture section defined")

    parse_arch_yaml(arch, arch_dict["architecture"])

    return arch


### Architecture description parsing ###
def parse_arch_yaml(arch, arch_dict):
    arch_name = arch_dict["name"]
    if "[" in arch_name:
        raise Exception("Error: multiple architectures not supported")

    if "tile" not in arch_dict:
        raise Exception("Error: No tiles defined, must be at least one tile")

    tiles = arch_dict["tile"]
    for tile_dict in tiles:
        parse_tile_yaml(arch, tile_dict)

    arch.set_noc_attributes(arch_dict["attributes"])
    return


def parse_tile_yaml(arch, tile_dict):
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
        tile = arch.create_tile(tile_name, tile_dict["attributes"])

        # Add any elements local to this h/w structure. They have access to any
        #  elements in the parent structures
        if "core" not in tile_dict:
            raise Exception("Error: No cores defined, "
                            "must be at least one core")
        cores = tile_dict["core"]
        for _, core_dict in enumerate(cores):
            parse_core_yaml(arch, tile.get_id(), core_dict)

    return


def parse_core_yaml(arch, tile_id, core_dict):
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
        core = arch.create_core(core_name, tile.get_id(), core_attributes)
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
        for dendrite in core_dict["dendrite"]:
            dendrite_attributes = dendrite.get("attributes", {})
            if dendrite_attributes is None:
                dendrite_attributes = {}
            core.create_synapse(dendrite["name"], dendrite_attributes)
        for soma in core_dict["soma"]:
            soma_attributes = soma.get("attributes", {})
            if soma_attributes is None:
                soma_attributes = {}
            core.create_soma(soma["name"], soma_attributes)
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
#def format_attributes(attributes):
#    line = ""
#    if attributes is None:
#        attributes = {}

#    for key in attributes:
#        line += (f" {key}={attributes[key]}")
#    return line


project_dir = os.path.dirname(os.path.abspath(__file__))
def run_kernel(arch_path, network_path, timesteps,
        perf_trace=False, spike_trace=False,
        potential_trace=False, message_trace=False, out_dir=None):
    print("Loading architecture\n")
    try:
        arch = load_arch_yaml(arch_path)
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
