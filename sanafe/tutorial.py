"""
Copyright (c) 2026 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

tutorial.py - Tutorial helper scripts, mostly for checking answers
"""
import re
import shutil
from importlib.resources import files
from functools import reduce

import yaml
from sanafecpp import load_arch, load_net

GREEN_TEXT = "\033[92m"
RED_TEXT = "\033[31m"
DEFAULT_TEXT = "\033[0m"


def copy_arch():
    """Copy the tutorial architecture into the working dir"""
    resource_path = files("sanafe.examples").joinpath("tutorial_arch.yaml")
    with resource_path.open("rb") as src, open("tutorial_arch.yaml", "wb") as dst:
        shutil.copyfileobj(src, dst)


def copy_snn():
    """Copy the tutorial SNN application into the working dir"""
    resource_path = files("sanafe.examples").joinpath("tutorial_snn.yaml")
    with resource_path.open("rb") as src, open("tutorial_snn.yaml", "wb") as dst:
        shutil.copyfileobj(src, dst)


def load():
    """Load the tutorial architecture and SNN"""
    base = files("sanafe.examples")
    arch = load_arch(base / "tutorial_arch.yaml")
    net  = load_net(base / "tutorial_snn.yaml", arch)
    return arch, net


def get_dvs_data():
    """Get raw (binary) weight data for the DVS gesture challenge."""
    return files("sanafe.examples").joinpath("dvs_challenge.npz").open("rb")


def check_arch():
    """Check architecture modification exercises."""
    with open("tutorial_arch.yaml", "r", encoding="utf-8") as arch_file:
        arch_details = yaml.safe_load(arch_file)
    check_arch_exercise_1(arch_details)
    check_arch_exercise_2(arch_details)
    check_arch_exercise_3(arch_details)


def check_arch_exercise_1(arch_details):
    """Check architecture against expected modifications (1)."""
    tiles = arch_details["architecture"]["tile"]
    cores = tiles[0]["core"]
    somas = cores[0]["soma"]
    soma_energy = somas[0]["attributes"]["energy_update_neuron"]
    soma_latency = somas[0]["attributes"]["latency_update_neuron"]
    if soma_energy == 2.0e-12 and soma_latency == 2.0e-9:
        print(f"{GREEN_TEXT}Exercise 1: PASS{DEFAULT_TEXT}")
    else:
        print(f"{RED_TEXT}Exercise 1: FAIL - Soma energy ({soma_energy}J) "
              f"and/or latency ({soma_latency}s) not set correctly{DEFAULT_TEXT}")


def _parse_name_range(s):
    """Extract range from tile or core name."""
    match = re.match(r"(\w+)\[(\d+)(?:\.\.(\d+))?\]", s)
    if match is None:
        return None, None
    return int(match.group(2)), int(match.group(3) or match.group(2))


def check_arch_exercise_2(arch_details):
    """Check architecture against expected modifications (2)."""
    tiles = arch_details["architecture"]["tile"]
    tile_name = tiles[0]["name"]
    range_start, range_end = _parse_name_range(tile_name)
    total_tiles = 2
    if range_start is None:
        print(f"{RED_TEXT}Exercise 2: FAIL - Tile not duplicated{DEFAULT_TEXT}")
        return

    duplicated_tiles = 0 if range_end is None else 1 + (range_end - range_start)
    if duplicated_tiles != total_tiles:
        print(f"{RED_TEXT}Exercise 2: FAIL - Tile duplicated {duplicated_tiles}"
              f" times, should be {total_tiles} times{DEFAULT_TEXT}")
        return

    cores = tiles[0]["core"]
    core_name = cores[0]["name"]
    range_start, range_end = _parse_name_range(core_name)
    total_cores = 4
    if range_start is None:
        print(f"{RED_TEXT}Exercise 2: FAIL - Cores not duplicated{DEFAULT_TEXT}")
        return

    duplicated_cores = 0 if range_end is None else 1 + (range_end - range_start)
    if duplicated_cores != total_cores:
        print(f"{RED_TEXT}Exercise 2: FAIL - Cores duplicated {duplicated_cores} "
              f"times, should be {total_cores} times{DEFAULT_TEXT}")
        return

    print(f"{GREEN_TEXT}Exercise 2: PASS{DEFAULT_TEXT}")


def check_arch_exercise_3(arch_details):
    """Check architecture against expected modifications (3)."""
    tiles = arch_details["architecture"]["tile"]
    cores = tiles[0]["core"]
    synapses = cores[0]["synapse"]
    if len(synapses) != 2:
        print(f"{RED_TEXT}Exercise 3: FAIL - Expected to see 2 synapse units, "
              f"only found {len(synapses)}")
    else:
        # Get the new soma unit
        synapse = synapses[0]
        if synapse["name"] == "tutorial_synapse_uncompressed":
            synapse = synapses[1]
        synapse_energy = synapse["attributes"]["energy_process_spike"]
        synapse_latency = synapse["attributes"]["latency_process_spike"]

        if synapse_energy == 0.5e-12 and synapse_latency == 2.0e-9:
            print(f"{GREEN_TEXT}Exercise 3: PASS{DEFAULT_TEXT}")
        else:
            print(f"{RED_TEXT}Exercise 3: FAIL - New synapse energy ({synapse_energy}J) "
                  f"and/or latency ({synapse_latency}s) not set correctly{DEFAULT_TEXT}")


def check_snn():
    """Check tutorial SNN against set of modifications."""
    with open("tutorial_snn.yaml", "r", encoding="utf-8") as snn_file:
        snn = yaml.safe_load(snn_file)
    check_exercise_snns_1(snn)
    check_exercise_snns_2(snn)
    check_exercise_snns_3(snn)
    check_exercise_snns_4(snn)


def check_exercise_snns_1(snn):
    """Check SNN against expected change to neuron group size."""
    net = snn["network"]
    group = net["groups"][1]
    neurons_found = len(group["neurons"])
    if len(group["neurons"]) != 2:
        print(f"{RED_TEXT}Exercise 1: FAIL - Should be 2 neurons in group 1, "
              f"found {neurons_found}{DEFAULT_TEXT}")
        return

    # SANA-FE checks other aspects of the mapping, if it runs, it should be fine
    print(f"{GREEN_TEXT}Exercise 1: PASS{DEFAULT_TEXT}")


def check_exercise_snns_2(snn):
    """Check SNN against expected added connection."""
    net = snn["network"]
    edges = net["edges"]
    if len(edges) < 3:
        print(f"{RED_TEXT}Exercise 2: FAIL - Expected 3 edges but got {len(edges)}{DEFAULT_TEXT}")
        return

    # TODO: check weight values are correct

    print(f"{GREEN_TEXT}Exercise 2: PASS{DEFAULT_TEXT}")


def check_exercise_snns_3(snn):
    """Check SNN against expected attribute modifications."""
    net = snn["network"]
    group = net["groups"][0]
    neuron = group["neurons"][1]
    attributes = list(neuron.values())[0]
    if ("bias" not in attributes or attributes["bias"] != 0.5):
        print(f"{RED_TEXT}Exercise 3: FAIL - Neuron 0.1 bias not set to "
              f"0.5{DEFAULT_TEXT}")
    else:
        print(f"{GREEN_TEXT}Exercise 3: PASS{DEFAULT_TEXT}")


def check_exercise_snns_4(snn):
    """Check SNN against expected mapping modifications."""
    net = snn["network"]
    group = net["groups"][1]
    attributes = reduce(lambda a, b: {**a, **b}, group["attributes"])

    if attributes["synapse_hw_name"] == "tutorial_synapse_uncompressed":
        print(f"{RED_TEXT}Exercise 4: FAIL - Set group 1 synapse h/w to your "
              f"new synapse H/W unit{DEFAULT_TEXT}")
    else:
        print(f"{GREEN_TEXT}Exercise 4: PASS{DEFAULT_TEXT}")


def check_api(snn):
    """Check SNN against expected API modifications."""
    check_exercise_api_1(snn)
    check_exercise_api_2(snn)


def check_exercise_api_1(snn):
    """Check SNN against expected API modifications (1)."""
    group = snn["out"]
    neurons_found = len(group)
    if len(group) != 2:
        print(f"{RED_TEXT}Exercise 1: FAIL - Should be 2 neurons in output "
              f"layer, found {neurons_found}{DEFAULT_TEXT}")
        return

    print(f"{GREEN_TEXT}Exercise 1: PASS{DEFAULT_TEXT}")


def check_exercise_api_2(snn):
    """Check SNN against expected API modifications (2)."""
    in_layer = snn["in"]
    connections_out = in_layer[0].edges_out

    if (len(connections_out) != 2 or
        connections_out[0].post_neuron.group_name != "out" or
        connections_out[1].post_neuron.group_name != "out"):
        print(f"{RED_TEXT}Exercise 2: FAIL - Should be 2 edges out of in.0, "
              f"to the output layer, found {len(connections_out)}{DEFAULT_TEXT}")
        return

    synapse_attributes = connections_out[1].synapse_attributes
    if (synapse_attributes.get("w") not in (-2, -2.0) and
        synapse_attributes.get("weight") not in (-2, -2.0)):
        print(f"{RED_TEXT}Exercise 2: FAIL - in.0 weight should be -2{DEFAULT_TEXT}")
        return

    connections_out = in_layer[1].edges_out
    if (len(connections_out) != 1 or
        connections_out[0].post_neuron.group_name != "out" or
        connections_out[0].post_neuron.neuron_offset != 1):
        print(f"{RED_TEXT}Exercise 2: FAIL - Should be 1 edge out of in.1 to "
              f"out.1, found {len(connections_out)}{DEFAULT_TEXT}")
        return

    synapse_attributes = connections_out[0].synapse_attributes
    if (synapse_attributes.get("w") not in (3, 3.0) and
        synapse_attributes.get("weight") not in (3, 3.0)):
        print(f"{RED_TEXT}Exercise 2: FAIL - in.0 weight should be 3{DEFAULT_TEXT}")
        return

    print(f"{GREEN_TEXT}Exercise 2: PASS{DEFAULT_TEXT}")
