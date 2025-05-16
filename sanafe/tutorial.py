"""
Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

tutorial.py - Tutorial helper scripts, mostly for checking answers
"""
import yaml
import re

green_text = "\033[92m"
red_text = "\033[31m"
default_text = "\033[0m"


def check_arch(arch_filename):
    with open("arch.yaml", "r") as arch_file:
        arch_details = yaml.safe_load(arch_file)
    check_arch_exercise_1(arch_details)
    check_arch_exercise_2(arch_details)
    check_arch_exercise_3(arch_details)


def check_arch_exercise_1(arch_details):
    tiles = arch_details["architecture"]["tile"]
    cores = tiles[0]["core"]
    somas = cores[0]["soma"]
    soma_energy = somas[0]["attributes"]["energy_update_neuron"]
    soma_latency = somas[0]["attributes"]["latency_update_neuron"]
    if soma_energy == 2.0e-12 and soma_latency == 2.0e-9:
        print(f"{green_text}Exercise 1: PASS{default_text}")
    else:
        print(f"{red_text}Exercise 1: FAIL - Soma energy ({soma_energy}J) "
              f"and/or latency ({soma_latency}s) not set correctly{default_text}")


def parse_name_range(s):
    match = re.match(r"(\w+)\[(\d+)(?:\.\.(\d+))?\]", s)
    if match is None:
        return None, None
    else:
        return int(match.group(2)), int(match.group(3) or match.group(2))


def check_arch_exercise_2(arch_details):
    tiles = arch_details["architecture"]["tile"]
    tile_name = tiles[0]["name"]
    range_start, range_end = parse_name_range(tile_name)
    passed = True
    total_tiles = 2
    if range_start is None:
        print(f"{red_text}Exercise 2: FAIL - Tile not duplicated{default_text}")
        return
    elif (range_end - range_start) + 1 != total_tiles:
        print(f"{red_text}Exercise 2: FAIL - Tile duplicated {1+(range_end-range_start)} "
              f"times, should be {total_tiles} times{default_text}")
        return

    cores = tiles[0]["core"]
    core_name = cores[0]["name"]
    range_start, range_end = parse_name_range(core_name)
    total_cores = 4
    if range_start is None:
        print(f"{red_text}Exercise 2: FAIL - Cores not duplicated{default_text}")
        passed = False
    elif (range_end - range_start) + 1 != total_cores:
        print(f"{red_text}Exercise 2: FAIL - Cores duplicated {1+(range_end-range_start)} "
              f"times, should be {total_cores} times{default_text}")
        passed = False

    if passed:
        print(f"{green_text}Exercise 2: PASS{default_text}")


def check_arch_exercise_3(arch_details):
    tiles = arch_details["architecture"]["tile"]
    cores = tiles[0]["core"]
    synapses = cores[0]["synapse"]
    if len(synapses) != 2:
        print(f"{red_text}Exercise 3: FAIL - Expected to see 2 synapse units, "
              f"only found {len(synapses)}")
    else:
        # Get the new soma unit
        synapse = synapses[0]
        if synapse["name"] == "tutorial_synapse_uncompressed":
            synapse = synapses[1]
        synapse_energy = synapse["attributes"]["energy_process_spike"]
        synapse_latency = synapse["attributes"]["latency_process_spike"]

        if synapse_energy == 0.5e-12 and synapse_latency == 2.0e-9:
            print(f"{green_text}Exercise 3: PASS{default_text}")
        else:
            print(f"{red_text}Exercise 3: FAIL - New synapse energy ({synapse_energy}J) "
                  f"and/or latency ({synapse_latency}s) not set correctly{default_text}")


def check_snn(snn_filename):
    with open(snn_filename, "r") as snn_file:
        snn = yaml.safe_load(snn_file)
    check_exercise_snns_1(snn)
    check_exercise_snns_2(snn)
    check_exercise_snns_3(snn)
    check_exercise_snns_4(snn)


def check_exercise_snns_1(snn):
    net = snn["network"]
    group = net["groups"][1]
    neurons_found = len(group["neurons"])
    if len(group["neurons"]) != 2:
        print(f"{red_text}Exercise 1: FAIL - Should be 2 neurons in group 1, found {neurons_found}{default_text}")
        return

    # SANA-FE will check other aspects of the mapping, if it runs, it should be
    #  fine
    print(f"{green_text}Exercise 1: PASS{default_text}")


def check_exercise_snns_2(snn):
    net = snn["network"]
    edges = net["edges"]
    if len(edges) < 3:
        print(f"{red_text}Exercise 2: FAIL - Expected 3 edges but got {len(edges)}{default_text}")
        return

    # TODO: check weights are correct

    print(f"{green_text}Exercise 2: PASS{default_text}")


def check_exercise_snns_3(snn):
    net = snn["network"]
    group = net["groups"][0]
    neuron = group["neurons"][1]
    attributes = list(neuron.values())[0]
    if ("bias" not in attributes or attributes["bias"] != 0.5):
        print(f"{red_text}Exercise 3: FAIL - Neuron 0.1 bias not set to "
              f"0.5{default_text}")
    else:
        print(f"{green_text}Exercise 3: PASS{default_text}")


def check_exercise_snns_4(snn):
    net = snn["network"]
    group = net["groups"][1]
    from functools import reduce
    attributes = reduce(lambda a, b: {**a, **b}, group["attributes"])

    if attributes["synapse_hw_name"] == "tutorial_synapse_uncompressed":
        print(f"{red_text}Exercise 4: FAIL - Set group 1 synapse h/w to your "
              f"new synapse H/W unit{default_text}")
    else:
        print(f"{green_text}Exercise 4: PASS{default_text}")


def check_api(snn):
    check_exercise_api_1(snn)
    check_exercise_api_2(snn)


def check_exercise_api_1(snn):
    group = snn["out"]
    neurons_found = len(group)
    if len(group) != 2:
        print(f"{red_text}Exercise 1: FAIL - Should be 2 neurons in output "
              f"layer, found {neurons_found}{default_text}")
        return

    print(f"{green_text}Exercise 1: PASS{default_text}")


def check_exercise_api_2(snn):
    in_layer = snn["in"]
    neuron = in_layer[0]
    connections_out = in_layer[0].edges_out

    if (len(connections_out) != 2 or
        connections_out[0].post_neuron.group_name != "out" or
        connections_out[1].post_neuron.group_name != "out"):
        print(f"{red_text}Exercise 2: FAIL - Should be 2 edges out of in.0, "
              f"to the output layer, found {len(connections_out)}{default_text}")
        return

    synapse_attributes = connections_out[1].synapse_attributes
    if (synapse_attributes.get("w") not in (-2, -2.0) and
        synapse_attributes.get("weight") not in (-2, -2.0)):
        print(f"{red_text}Exercise 2: FAIL - in.0 weight should be -2{default_text}")
        return

    connections_out = in_layer[1].edges_out
    if (len(connections_out) != 1 or
        connections_out[0].post_neuron.group_name != "out" or
        connections_out[0].post_neuron.neuron_offset != 1):
        print(f"{red_text}Exercise 2: FAIL - Should be 1 edge out of in.1 to "
              f"out.1, found {len(connections_out)}{default_text}")
        return

    synapse_attributes = connections_out[0].synapse_attributes
    if (synapse_attributes.get("w") not in (3, 3.0) and
        synapse_attributes.get("weight") not in (3, 3.0)):
        print(f"{red_text}Exercise 2: FAIL - in.0 weight should be 3{default_text}")
        return

    print(f"{green_text}Exercise 2: PASS{default_text}")
