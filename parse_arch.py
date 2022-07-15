"""Parse architecture description YAML file, output machine-readable list"""
import yaml

classes = ("synapse", "soma", "dendrite", "router", "memory",
           "axon_in", "axon_out")
MAX_RECURSION = 32


def parse(arch_dict):
    print(arch_dict)

    # On the first pass, parse anything that isn't "sim"
    #  The user can define custom blocks that can be reused
    #  A globally defined block can't contain other structures? lets see

    if "architecture" not in arch_dict:
        print("Error: no architecture under sim defined")
    else:
        arch_elements = parse_struct(arch_dict["architecture"], [], 1)

    return arch_elements


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
        range_min, range_max = 0, 1

    struct_elements = list()
    for instance in range(range_min, range_max+1):
        struct_name = struct_name.split("[")[0] + "[{0}]".format(instance)
        print("Parsing struct {0}".format(struct_name))
 
        # Add any elements local to this h/w structure. They have access to any
        #  elements in the parent structures
        if "local" in struct_dict:
            local_elements = parse_local(struct_dict["local"], parent_elements)
        else:
            local_elements = list()
        
        # When looking at the subtree, elements may use information about this
        #  structure and any parent structures 
        subtree_elements = list()
        elements = parent_elements + local_elements
        if "subtree" in struct_dict:
            parse_subtree(struct_dict["subtree"], elements, recursion_depth+1)

    return


def parse_subtree(subtree, parent_elements, recursion_depth):
    print("Parsing subtree, parent: {0}".format(parent_elements))
    if not isinstance(subtree, list):
        raise Exception("Subtree must have list of branches")

    for branch in subtree:
        parse_struct(branch, parent_elements, recursion_depth)


def parse_local(local, parent_elements):
    if not isinstance(local, list):
        raise Exception("Local vars must be list of elements")

    elements = list()
    for el in local:
        elements = elements + parse_element(el) 

    return elements


def parse_element(element_dict):
    element_name = element_dict["name"]
    if "[" in element_name:
        range_min, range_max = parse_range(element_name)
        instances = (range_max+1) - range_min 
    else:
        instances = 1

    element_class = element_dict["class"]
    if element_class == "compartment":
        new_elements = create_compartment(element_dict, instances)
    elif element_class == "synapse":
        new_elements = create_synapse(element_dict, instances)
    elif element_class == "soma":
        new_elements = create_soma(element_dict, instances)
    elif element_class == "dendrite":
        new_elements = create_dendrite(element_dict, instances)
    elif element_class == "router":
        new_elements = create_router(element_dict, instances)
    elif element_class == "memory":
        new_elements = create_memory(element_dict, instances)
    elif element_class == "axon_in":
        new_elements = create_axon_in(element_dict, instances)
    elif element_class == "axon_out":
        new_elements = create_axon_out(element_dict, instances)
    else:
        raise Exception("Class {0} not recognized".format(element_class))

    return new_elements


def create_compartment(element_dict, instances):
    return [element_dict["name"]] * instances


def create_router(element_dict, instances):
    return [element_dict["name"]] * instances


def create_synapse(element_dict, instances):
    return [element_dict["name"]] * instances


def create_soma(element_dict, instances):
    return [element_dict["name"]] * instances


def create_dendrite(element_dict, instances):
    return [element_dict["name"]] * instances


def create_memory(element_dict, instances):
    return [element_dict["name"]] * instances


def create_axon_in(element_dict, instances):
    return [element_dict["name"]] * instances


def create_axon_out(element_dict, instances):
    return [element_dict["name"]] * instances


# TODO: this is only needed until we use the Python api instead
def format_element(element):
    """Print element as line to file"""
    return element


if __name__ == "__main__":
    arch_list = None
    with open("loihi.yaml", "r") as arch_file:
        arch_dict = yaml.safe_load(arch_file)
        arch_elements = parse(arch_dict)

    with open("out", "w") as list_file:
        for element in arch_elements:
            line = format_element(element)
            print(line)
            list_file.write(line + '\n')

