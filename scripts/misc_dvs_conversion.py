"""misc_dvs_conversion.py

A temporary script to help with moving from the old network format to the new
YAML-based format. I've already created a script to convert netlists to YAML
files, so this script just outputs the convolutional filters in the new,
more compact format, and can update inputs.

Going forward, I need to make sure the SNNToolbox outputs networks in the new
compact format, using conv2d and dense layers, and efficiently outputs the
mappings.
"""
import yaml
import numpy as np
import sys

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

if (len(sys.argv) != 4):
    print("Usage: convert_dvs_edges.py <snn> <in> <out>")
    exit()
script_name, yaml_filename, np_filename, out_filename = sys.argv

#with open(yaml_filename, "r") as yaml_file:
#    snn = yaml.safe_load(yaml_file)
info = np.load(np_filename)

edges = [
    {"0 -> 1": {"type": "conv2d", "weight": flow_list(info["conv1"].astype(int).flatten().tolist())}},
    {"1 -> 2": {"type": "conv2d", "weight": flow_list(info["conv2"].astype(int).flatten().tolist())}},
    {"2 -> 3": {"type": "conv2d", "weight": flow_list(info["conv3"].astype(int).flatten().tolist())}},
    {"3 -> 4": {"type": "conv2d", "weight": flow_list(info["conv4"].astype(int).flatten().tolist())}},
    {"4 -> 5": {"type": "conv2d", "weight": flow_list(info["dense1"].astype(int).flatten().tolist())}},
]

input_neurons = []
for id, bias in enumerate(info["inputs"].astype(int).tolist()):
    input_neurons.append(flow_dict({id: {"bias": bias}}))

print(yaml.dump(edges))
with open(out_filename, "w") as description_file:
    yaml.dump({"network": {"edges": edges}, "neurons": input_neurons},
              description_file, default_flow_style=False)
