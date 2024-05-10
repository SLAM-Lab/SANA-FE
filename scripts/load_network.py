"""
Copyright (c) 2024 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

load_network.py - Simulator script and utility functionality
"""
import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, os.path.join(PROJECT_DIR))
# Set some flags for the dynamic linking library
# Important to do before importing the simcpp .so library!
sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)
import sanafecpp

def load_from_net_file(filename, arch, heartbeat=100):
    neurons_loaded = 0
    edges_loaded = 0
    mappings_loaded = 0

    net = sanafecpp.Network()
    with open(filename, 'r') as network_file:
        for line in network_file:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == 'g':
                neuron_count = int(fields[1])
                group_attributes = dict([f.split('=') for f in fields[2:]])
                group = net.create_neuron_group(neuron_count,
                                                group_attributes)
                print("Loaded group")
            elif fields[0] == 'n':
                group_id = int(fields[1].split('.')[0])
                neuron_id = int(fields[1].split('.')[1])
                neuron_attributes = dict([f.split('=') for f in fields[2:]])
                group = net.groups[group_id]
                group.define_neuron(neuron_id, neuron_attributes)
                neurons_loaded += 1
                if (neurons_loaded % heartbeat) == 0:
                    print(f"Loaded {neurons_loaded} neurons")

            elif fields[0] == 'e':
                edge_info = fields[1]
                src_address = edge_info.split("->")[0]
                dest_address = edge_info.split("->")[1]

                src_gid = int(src_address.split(".")[0])
                src_nid = int(src_address.split(".")[1])
                src = net.groups[src_gid].neurons[src_nid]

                dest_gid = int(dest_address.split(".")[0])
                dest_nid = int(dest_address.split(".")[1])
                dest = net.groups[dest_gid].neurons[dest_nid]

                edge_attributes = dict([f.split('=') for f in fields[2:]])
                src.connect_to_neuron(dest, edge_attributes)
                edges_loaded += 1
                if (edges_loaded % heartbeat) == 0:
                    print(f"Loaded {edges_loaded} edges")

            elif fields[0] == '&':
                mapping_info = fields[1]
                neuron_address = mapping_info.split("@")[0]
                core_address = mapping_info.split("@")[1]

                group_id = int(neuron_address.split(".")[0])
                neuron_id = int(neuron_address.split(".")[1])
                neuron = net.groups[group_id].neurons[neuron_id]

                tile_id = int(core_address.split(".")[0])
                core_offset = int(core_address.split(".")[1])
                core = arch.tiles[tile_id].cores[core_offset]

                core.map_neuron(neuron)
                mappings_loaded += 1
                if (mappings_loaded % heartbeat) == 0:
                    print(f"Loaded {mappings_loaded} mappings")

    return net
