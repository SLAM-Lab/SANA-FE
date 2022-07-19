"""Parse architecture description YAML file, output machine-readable list"""
# TODO: It seems obviously better to send a number of instances to the API
#  rather than just calling the API a number of times with the exact same
#  element. But does it make sense for every element. Or just compartments?
import yaml

MAX_RECURSION = 32

def parse(arch_dict):
    print(arch_dict)

    # On the first pass, parse anything that isn't "sim"
    #  The user can define custom blocks that can be reused
    #  A globally defined block can't contain other structures? lets see

    if "architecture" not in arch_dict:
        print("Error: no architecture under sim defined")
    else:
        parse_struct(arch_dict["architecture"], [], 1)

    return


def parse_range(range_str): 
    range_str = range_str.replace("]", "")
    range_str = range_str.split("[")[1]
    range_min = int(range_str.split("..")[0])
    range_max = int(range_str.split("..")[1])

    return range_min, range_max


def parse_struct(struct_dict, parent_elements, recursion_depth):
    if recursion_depth > MAX_RECURSION:
        raise Exception("Error: Exceeded max recursion depth ({0}), "
                        "stopping!".format(MAX_RECURSION))

    struct_name = struct_dict["name"]
    # Work out how many instances of this structure to create
    if "[" in struct_name:
        # Can use notation [min..max] to indicate range of elements
        range_min, range_max = parse_range(struct_name)
    else:
        range_min, range_max = 0, 0

    struct_elements = list()
    for instance in range(range_min, range_max+1):
        struct_name = struct_name.split("[")[0] + "[{0}]".format(instance)
        print("Parsing struct {0}".format(struct_name))
 
        # Add any elements local to this h/w structure. They have access to any
        #  elements in the parent structures
        if "local" in struct_dict:
            local_elements = parse_local(struct_dict["local"], parent_elements)
        else:
            local_elements = dict()
        
        # When looking at the subtree, elements may use information about this
        #  structure and any parent structures 
        subtree_elements = list()
        elements = combine_elements(local_elements, parent_elements)
        if "subtree" in struct_dict:
            parse_subtree(struct_dict["subtree"], elements, recursion_depth+1)

    return


def combine_elements(local, parent):
    print(local)
    combined = dict(local)
    for key in parent:
        if key not in combined:
            combined[key] = []

        combined[key].append(parent[key])
    return combined


def parse_subtree(subtree, parent_elements, recursion_depth):
    print("Parsing subtree, parent: {0}".format(parent_elements))
    if not isinstance(subtree, list):
        raise Exception("Subtree must have list of branches")

    for branch in subtree:
        parse_struct(branch, parent_elements, recursion_depth)


def parse_local(local, parent_ids):
    if not isinstance(local, list):
        raise Exception("Local vars must be list of elements")

    # Now create all the local elements in the simulation or arch description
    # Create all non-compartment elements first, then pass these as a reference
    #  to the compartment 
    elements = {"memory": [], "synapse": [], "compartment": [], "soma": [],
                "dendrite": [], "axon_in": [], "axon_out": [], "router": [],}
    for el in local:
        # Group elements of the same class together   
        # TODO: rename to get_instances ? 
        el["instances"] = get_instances(el)
        class_name = el["class"]
        elements[class_name].append(el)

    routers = [parse_router(el) for el in elements["router"]]
    somas = [parse_soma(el) for el in elements["soma"]]
    memories = [parse_memory(el) for el in elements["memory"]]
    dendrites = [parse_dendrite(el) for el in elements["dendrite"]]
    synapses = [parse_synapse(el) for el in elements["synapse"]]

    local_ids = { "routers": routers, "somas": somas, "memories": memories,
                  "dendrites": dendrites, "synapses": synapses }

    # These elements need to know about other elements that are local or higher
    #  in the hierarchy
    if "routers" in parent_ids:
        routers = routers + parent_ids["routers"]
    local_ids["axon_inputs"] = []
    for el in elements["axon_in"]:
        el["routers"] = routers
        local_ids["axon_inputs"].append(parse_axon_in(el))

    local_ids["axon_outputs"] = []
    for el in elements["axon_out"]:
        el["routers"] = routers
        local_ids["axon_outputs"].append(parse_axon_out(el))

    # Need to combine both parent and local elements and pass these to the
    #  compartment
    dependencies = combine_elements(local_ids, parent_ids)
    compartments = []
    for el in elements["compartment"]:
        el["synapses"] = dependencies["synapses"]
        el["somas"] = dependencies["somas"]
        el["dendrites"] = dependencies["dendrites"]
        el["synapses"] = dependencies["synapses"]
        el["axon_inputs"] = dependencies["axon_inputs"]
        el["axon_outputs"] = dependencies["axon_outputs"]
        compartments = parse_compartment(el)

    return local_ids


def get_instances(element_dict):
    element_name = element_dict["name"]

    if "[" in element_name:
        range_min, range_max = parse_range(element_name)
        instances = (range_max - range_min) + 1
    else:
        instances = 1

    print("Parsing element {0}".format(element_name)) 
    return instances


def parse_compartment(element_dict):
    global _compartment_count

    instances = element_dict["instances"]
    attributes = element_dict["attributes"]


    compartment_ids = create_compartment(instances, element_dict["synapses"],
                                          element_dict["dendrites"],
                                          element_dict["somas"],
                                          element_dict["axon_inputs"],
                                          element_dict["axon_outputs"])

    return [compartment_ids]


def parse_router(element_dict):
    router_ids = []

    assert(element_dict["instances"] == 1)

    attributes = element_dict["attributes"]
    if "connection" in attributes:
        connection_type = attributes["connection"]
        dimensions = int(attributes["dimensions"])
        width = int(attributes["width"])

    return create_router(dimensions, width, connection_type)


def parse_synapse(element_dict):
    attributes = element_dict["attributes"]

    print(attributes)
    model = attributes["model"]
    bits = attributes["bits"]

    return create_synapse(model, bits)


def parse_soma(element_dict):
    return create_soma()


def parse_axon_out(element_dict):
    routers = element_dict["routers"]
    assert(len(routers) == 1)
    return create_axon_out(routers[0])


def parse_axon_in(element_dict):
    routers = element_dict["routers"]
    assert(len(routers) == 1)
    return create_axon_in(routers[0])


def parse_memory(element_dict):
    attributes = element_dict["attributes"]
    mem_size = attributes["size"]
    return create_memory(mem_size)
"""
Functions to add elements. These can be swapped out for python API calls in
the future.
"""
_compartments = []
_compartment_count = 0
def create_compartment(instances, synapses, dendrites, somas, axon_inputs,
                       axon_outputs):
    global _compartment_count

    compartment_id = _compartment_count
    _compartment_count += instances
    if instances > 1:
        id_str= "{0}..{1}".format(compartment_id,
                                  (compartment_id+instances) - 1)
    else:
        id_str = str(compartment_id)

    compartment = "c {0} {1} {2} {3}".format(id_str, synapses[0],
                                             axon_inputs[0], axon_outputs[0])
    _compartments.append(compartment)

    # TODO: return list of ids
    return compartment_id


_routers = []
def create_router(dimensions, width, connection_type):
    router_id = len(_routers)

    if connection_type == "mesh":
        assert(dimensions == 2)
        x = router_id % width
        y = int(router_id / width)

    router = "r {0} {1} {2}".format(router_id, x, y) 
    _routers.append(router)
    return router_id


_synapses = []
def create_synapse(synapse_model, bits):
    synapse_id = len(_synapses)

    models = { "cuba": 0 }
    model_id = models[synapse_model]

    synapse = "s {0} {1} {2}".format(synapse_id, model_id, bits)
    _synapses.append(synapse)

    return synapse_id


_somas = []
def create_soma():
    soma_id = len(_somas)
    soma = "+ {0}".format(soma_id)
    _somas.append(soma)

    return soma_id


_dendrites = []
def parse_dendrite(element_dict):
    return create_dendrite()

def create_dendrite():
    dendrite_id = len(_dendrites)
    dendrite = "d {0}".format(dendrite_id)
    _dendrites.append(dendrite)

    return dendrite_id


_memories = []
def create_memory(mem_size):
    memory_id = len(_memories)
    memory = "m {0} {1}".format(memory_id, mem_size, )
    _memories.append(memory)

    return memory_id


_axon_inputs = []
def create_axon_in(router_id):
    axon_id = len(_axon_inputs)
    axon = "i {0} {1}".format(axon_id, router_id)
    _axon_inputs.append(axon)

    return axon_id


_axon_outputs = []
def create_axon_out(router_id):
    axon_id = len(_axon_outputs)
    axon = "o {0} {1}".format(axon_id, router_id)
    _axon_outputs.append(axon)

    return axon_id


if __name__ == "__main__":
    arch_list = None
    with open("loihi.yaml", "r") as arch_file:
        arch_dict = yaml.safe_load(arch_file)
        parse(arch_dict)

    #arch_elements = _compartments + _routers
    arch_elements = _compartments + _routers + _memories + _synapses + \
                    _axon_inputs + _axon_outputs + _somas + _dendrites
    #print(arch_elements)

    with open("out", "w") as list_file:
        for line in arch_elements:
            print(line) 
            list_file.write(line + '\n')
