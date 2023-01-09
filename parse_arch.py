"""Parse architecture description YAML file, output machine-readable list"""
import yaml

def parse(arch_dict):
    print(arch_dict)

    # On the first pass, parse anything that isn't "sim"
    #  The user can define custom blocks that can be reused
    #  A globally defined block can't contain other structures? lets see

    if "architecture" not in arch_dict:
        print("Error: no architecture defined")
    else:
        parse_arch(arch_dict["architecture"])

    return


def parse_range(range_str):
    range_str = range_str.replace("]", "")
    range_str = range_str.split("[")[1]
    range_min = int(range_str.split("..")[0])
    range_max = int(range_str.split("..")[1])

    return range_min, range_max


def parse_arch(arch_dict):
    arch_name = arch_dict["name"]
    if "[" in arch_name:
        raise Exception("Error: multiple architectures not supported")

    if "attributes" not in arch_dict:
        raise Exception("Error: NoC not defined for architecture "
                        "(please add this under attributes)")
    attributes = arch_dict["attributes"]
    if "topology" not in attributes:
        raise Exception("Error: topology not defined for NoC")

    topology = attributes["topology"]
    if topology == "mesh":
        dimensions = int(attributes["dimensions"])
        width = int(attributes["width"])
    cost = attributes["cost"]

    network_attributes = {"topology": topology, "dimensions": dimensions,
                          "width": width, "cost": cost}

    if "tile" not in arch_dict:
        raise Exception("Error: No tiles defined, must be at least one tile")

    tiles = arch_dict["tile"]
    for tile in tiles:
        parse_tile(tile, network_attributes)

    create_noc(8, 4)
    return


def parse_tile(tile_dict, network_attributes):
    tile_name = tile_dict["name"]
    # Work out how many instances of this tile to create
    if "[" in tile_name:
        # Can use notation [min..max] to indicate range of elements
        range_min, range_max = parse_range(tile_name)
    else:
        range_min, range_max = 0, 0

    tile_elements = list()
    for instance in range(range_min, range_max+1):
        tile_name = tile_name.split("[")[0] + "[{0}]".format(instance)
        tile_id = create_tile(network_attributes)
        print("Parsing struct {0}".format(tile_name))

        # Add any elements local to this h/w structure. They have access to any
        #  elements in the parent structures
        if "core" not in tile_dict:
            raise Exception("Error: No cores defined, "
                            "must be at least one core")
        cores = tile_dict["core"]
        for core_id, core_dict in enumerate(cores):
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
        core_id = create_core(tile_id)
        print("Parsing struct {0}".format(core_name))

        axon_inputs = [parse_axon_in(el, tile_id, core_id)
                       for el in core_dict["axon_in"]]
        synapse_processors = [parse_synapse(el, tile_id, core_id) for
                              el in core_dict["synapse"]]
        dendrite_processors = [parse_dendrite(el, tile_id, core_id) for
                               el in core_dict["dendrite"]]
        soma_processors = [parse_soma(el, tile_id, core_id) for
                           el in core_dict["soma"]]
        axon_outputs = [parse_axon_out(el, tile_id, core_id) for
                        el in core_dict["axon_out"]]


def parse_synapse(element_dict, tile_id, core_id):
    # TODO: parse the number of elements to create
    attributes = element_dict["attributes"]

    print(attributes)
    model = attributes["model"]
    weight_bits = attributes["weight_bits"]

    synaptic_op_energy, synaptic_op_time = 0.0, 0.0
    if "cost" in attributes:
        costs = attributes["cost"]
        synaptic_op_energy, synaptic_op_time = costs["energy"], costs["time"]

    return create_synapse(tile_id, core_id, model, weight_bits,
                          synaptic_op_energy, synaptic_op_time)


def parse_axon_out(element_dict, tile_id, core_id):
    attributes = element_dict["attributes"]
    spike_energy, spike_time = 0.0, 0.0
    if "cost" in attributes:
        costs = attributes["cost"]
        spike_energy, spike_time = costs["energy"], costs["time"]

    return create_axon_out(tile_id, core_id, spike_energy, spike_time)


def parse_axon_in(element_dict, tile_id, core_id):
    return create_axon_in(tile_id, core_id)


def parse_soma(element_dict, tile_id, core_id):
    active_energy, active_time = 0.0, 0.0
    inactive_energy, inactive_time = 0.0, 0.0

    attributes = element_dict["attributes"]
    if "cost" in attributes:
        costs = attributes["cost"]
        active_energy, active_time = costs["active"]["energy"], costs["active"]["time"]
        inactive_energy, inactive_time = costs["inactive"]["energy"], costs["inactive"]["time"]
        spiking_energy, spiking_time = costs["spiking"]["energy"], costs["spiking"]["time"]

    return create_soma(tile_id, core_id, active_energy, active_time,
                        inactive_energy, inactive_time, spiking_energy,
                        spiking_time)


def parse_dendrite(element_dict, tile_id, core_id):
    attributes = element_dict["attributes"]
    update_energy, update_time = 0.0, 0.0
    if "cost" in attributes:
        costs = attributes["cost"]
        update_energy, update_time = costs["energy"], costs["time"]
    return create_dendrite(tile_id, core_id, update_energy, update_time)


def get_instances(element_dict):
    element_name = element_dict["name"]

    if "[" in element_name:
        range_min, range_max = parse_range(element_name)
        instances = (range_max - range_min) + 1
    else:
        instances = 1

    print("Parsing element {0}".format(element_name))
    return instances

"""
Functions to add elements. These can be swapped out for python API calls in
the future.
"""
_tiles = 0
_cores_in_tile = []
_axon_inputs_in_core = [[]]
_synapses_in_core = [[]]
_dendrites_in_core = [[]]
_somas_in_core = [[]]
_axon_outputs_in_core = [[]]
_command_list = []

def create_tile(network_attributes):
    global _tiles
    tile_id = _tiles
    _tiles += 1

    assert(network_attributes["topology"] == "mesh") # TODO: for now
    #if network_attributes["topology"] == "mesh":
    assert(network_attributes["dimensions"] == 2) # TODO: for now
    x = tile_id % network_attributes["width"]
    y = int(tile_id / network_attributes["width"])

    east_west_energy = 0.0
    east_west_time = 0.0
    north_south_energy = 0.0
    north_south_time = 0.0
    barrier_time = 0.0
    if "cost" in network_attributes:
        east_west_energy = network_attributes["cost"]["east_west"]["energy"]
        east_west_time = network_attributes["cost"]["east_west"]["time"]
        north_south_energy = network_attributes["cost"]["north_south"]["energy"]
        north_south_time = network_attributes["cost"]["north_south"]["time"]
        barrier_time = network_attributes["cost"]["global_barrier"]["time"]
        # TODO: this will have to be generalized for other topologies
    costs = "{0} {1} {2} {3} {4}".format(east_west_energy, east_west_time,
            north_south_energy, north_south_time, barrier_time)

    tile = "t {0} {1} {2}".format(x, y, costs)
    _command_list.append(tile)
    # Track how many cores are in this tile
    _cores_in_tile.append(0)
    _axon_inputs_in_core.append(list())
    _synapses_in_core.append(list())
    _dendrites_in_core.append(list())
    _somas_in_core.append(list())
    _axon_outputs_in_core.append(list())

    return tile_id


def create_core(tile_id):
    core_id = _cores_in_tile[tile_id]
    core = "c {0}".format(tile_id)
    _command_list.append(core)
    _cores_in_tile[tile_id] += 1

    _axon_inputs_in_core[tile_id].append(0)
    _synapses_in_core[tile_id].append(0)
    _dendrites_in_core[tile_id].append(0)
    _somas_in_core[tile_id].append(0)
    _axon_outputs_in_core[tile_id].append(0)

    return core_id


def create_synapse(tile_id, core_id, synapse_model, bits, synaptic_op_energy,
                   synaptic_op_time):
    synapse_id = _synapses_in_core[tile_id][core_id]
    _synapses_in_core[tile_id][core_id] += 1

    # TODO: more features
    #models = { "cuba": 0 }
    #model_id = models[synapse_model]

    synapse = "s {0} {1} {2} {3}".format(tile_id, core_id, synaptic_op_energy,
                                         synaptic_op_time)
    _command_list.append(synapse)

    return synapse_id


def create_dendrite(tile_id, core_id, update_energy, update_time):
    dendrite_id = _dendrites_in_core[tile_id][core_id]
    _dendrites_in_core[tile_id][core_id] += 1
    dendrite = "d {0} {1} {2} {3}".format(tile_id, core_id,
                                                update_energy, update_time)
    _command_list.append(dendrite)

    return dendrite_id

_somas = []
def create_soma(tile_id, core_id, active_energy, active_time, inactive_energy,
                    inactive_time, spiking_energy, spiking_time):
    soma_id = _somas_in_core[tile_id][core_id]
    _somas_in_core[tile_id][core_id] += 1
    soma = "+ {0} {1} {2} {3} {4} {5} {6} {7}".format(tile_id, core_id,
            active_energy, active_time, inactive_energy, inactive_time,
            spiking_energy, spiking_time)
    _command_list.append(soma)


def create_axon_in(tile_id, core_id):
    axon_id = _axon_inputs_in_core[tile_id][core_id]
    axon = "i {0} {1}".format(tile_id, core_id)
    _axon_inputs_in_core[tile_id][core_id] += 1
    _command_list.append(axon)

    return axon_id


def create_axon_out(tile_id, core_id, spike_energy, spike_time):
    axon_id = _axon_outputs_in_core[tile_id][core_id]
    axon = "o {0} {1} {2} {3}".format(tile_id, core_id, spike_energy, spike_time)
    _axon_outputs_in_core[tile_id][core_id] += 1
    _command_list.append(axon)

    return axon_id


def create_noc(width, height):
    _command_list.append("@ mesh 2 {0} {1}".format(width, height))
    return width*height


if __name__ == "__main__":
    arch_list = None
    with open("loihi.yaml", "r") as arch_file:
        arch_dict = yaml.safe_load(arch_file)
        parse(arch_dict)

    arch_elements = _command_list
    print(arch_elements)

    with open("out", "w") as list_file:
        for line in arch_elements:
            print(line)
            list_file.write(line + '\n')
