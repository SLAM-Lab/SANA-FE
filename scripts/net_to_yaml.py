"""
Copyright (c) 2024 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

net_to_yaml.py - Convert from v1 .net file to v2 .yaml description
"""
import sys
import os
import yaml

# Hacks to get pyyaml to print with a mix of flow and block styles
class block_dict(dict):
    pass

def block_dict_rep(dumper, data):
    return dumper.represent_mapping(u"tag:yaml.org,2002:map", data,
                                    flow_style=False)

class flow_dict(dict):
    pass

def flow_dict_rep(dumper, data):
    return dumper.represent_mapping(u"tag:yaml.org,2002:map", data,
                                    flow_style=True)

class block_list(list):
    pass

def block_list_rep(dumper, data):
    return dumper.represent_sequence(u"tag:yaml.org,2002:seq", data,
                                    flow_style=False)

class flow_list(list):
    pass

def flow_list_rep(dumper, data):
    return dumper.represent_sequence(u"tag:yaml.org,2002:seq", data,
                                    flow_style=True)

yaml.add_representer(block_dict, block_dict_rep)
yaml.add_representer(flow_dict, flow_dict_rep)
yaml.add_representer(block_list, block_list_rep)
yaml.add_representer(flow_list, flow_list_rep)
# End of hacks

def parse_attributes(attributes):
    float_attributes = ("w", "weight", "bias", "threshold", "reset", "leak_decay")
    bool_attributes = ("log_spikes", "log_v", "force_update")
    for key, val in attributes.items():
        if key in float_attributes:
            val = float(val)
            if val.is_integer():
                val = int(val)
            attributes[key] = val
        elif key in bool_attributes:
            if val.is_integer():
                val = int(val)
                if val != 0:
                    val = True
                else:
                    val = False
            elif val == "true" or val == "True":
                val = True
            else:
                val = False
            attributes[key] = val

    if "connections_out" in attributes:
        del(attributes["connections_out"])

    return attributes

def load_from_net_file(filename, heartbeat=100000):
    # Use short description formats for neurons and edges
    neurons_loaded = 0
    edges_loaded = 0
    mappings_loaded = 0

    net = {"network": {}, "mappings": []}
    neurons_in_group = []
    with open(filename, 'r') as net_file:
        net_name = filename.split('.')[0]
        net["network"] = {"name": net_name, "groups": [], "edges": {}}
        for line in net_file:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == 'g':
                group_attributes = dict([f.split('=') for f in fields[2:]])
                group_attributes = flow_dict(parse_attributes(group_attributes))
                group_id = len(net["network"]["groups"])
                net["network"]["groups"].append({"name": group_id,
                                                 "attributes": group_attributes,
                                                 "neurons": {}})
                neurons_in_group.append([])
                print(f"Loaded group: {group_id}")
            elif fields[0] == 'n':
                group_id = int(fields[1].split('.')[0])
                neuron_id = int(fields[1].split('.')[1])
                neuron_attributes = dict([f.split('=') for f in fields[2:]])
                neuron_attributes = flow_dict(
                    parse_attributes(neuron_attributes))

                # Add some temporary logic to compress
                if ((len(neurons_in_group[group_id]) > 0) and
                    (neuron_attributes == neurons_in_group[group_id][-1]["attributes"])):
                    # Neuron is identical, don't redefine
                    neurons_in_group[group_id][-1]["_last"] += 1
                else:
                    print("Creating new neuron definition")
                    neurons_in_group[group_id].append(
                        {"_first": neuron_id, "_last": neuron_id,
                         "attributes": neuron_attributes})
                neurons_loaded += 1
                if heartbeat is not None and (neurons_loaded % heartbeat) == 0:
                    print(f"Loaded {neurons_loaded} neurons")

            elif fields[0] == 'e':
                edge_info = fields[1]
                src_address = edge_info.split("->")[0]
                dest_address = edge_info.split("->")[1]

                edge_description = f"{src_address} -> {dest_address}"
                edge_attributes = dict([f.split('=') for f in fields[2:]])
                edge_attributes = flow_dict(parse_attributes(edge_attributes))
                net["network"]["edges"][edge_description] = edge_attributes

                edges_loaded += 1
                if heartbeat is not None and (edges_loaded % heartbeat) == 0:
                    print(f"Loaded {edges_loaded} edges")

            # Mapping now goes in separate file, figure out
            elif fields[0] == '&':
                mapping_info = fields[1]
                neuron_address = mapping_info.split("@")[0]
                core_address = mapping_info.split("@")[1]

                mapping = flow_dict({neuron_address: {"core": core_address}})
                net["mappings"].append(mapping)
                mappings_loaded += 1
                if (mappings_loaded % heartbeat) == 0:
                    print(f"Loaded {mappings_loaded} mappings")

    # Go through structs and convert to YAML dictionary
    for group_id, neurons in enumerate(neurons_in_group):
        for n in neurons:
            first, last = n["_first"], n["_last"]
            name = first
            if last > first:
                name = f"{first}..{last}"

            attributes = n["attributes"]
            net["network"]["groups"][group_id]["neurons"][name] = attributes

    return net

if len(sys.argv) < 2:
    print("Usage: python net_to_yaml.py [file.net]")
    exit(1)
filename = sys.argv[1]
print(f"Old filename: {filename}")
net = load_from_net_file(filename)

if len(sys.argv) > 2:
    new_filename = sys.argv[2]
else:
    new_filename = filename.split('.')[0] + ".yaml"
with open(new_filename, "w") as yaml_file:
    print("Writing to YAML file.")
    yaml.dump(net, yaml_file, sort_keys=False)
print(f"YAML file {new_filename} written.")
