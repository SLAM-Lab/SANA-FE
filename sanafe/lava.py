#LAVA to SANAFE UTILS

# This module converts LAVA processes or serialization
# object to the SNN representation runnable on SANA-FE

# Implemented by Lance Lui as part of the capstone senior design project

# TODO: adapt to work with updated SANA-FE
# TODO: update to actually pull parameters from the Lava process
import lava.utils.serialization
from lava.magma.core.process.process import AbstractProcess
from lava.magma.compiler.executable import Executable

from lava.utils.serialization import load

import os
import sys
import importlib

"""
keep this file structure-
lava
├── src
│   ├── lava
│   │   ├── utils
│   │   │   ├── sanafe.py
│   │   │   │...
│   │   │...
│   │...
│...
SANAFE
├── sim.py
│...
"""

UTILS_DIR = os.path.dirname(os.path.abspath(__file__))
LAVA_DIR = os.path.dirname(UTILS_DIR)
SRC_DIR = os.path.dirname(LAVA_DIR)
LAVA_PROJECT_DIR = os.path.dirname(SRC_DIR)
PROJECT_DIR = os.path.dirname(LAVA_PROJECT_DIR)
sys.path.append(PROJECT_DIR)
sim = importlib.import_module('SANA-FE.snn')

NETWORK_FILENAME = "runs/random/random.net"
ARCH_FILENAME = "arch/loihi.yaml"
LOIHI_NEURONS_PER_CORE = 1024
LOIHI_CORES = 128
LOIHI_CORES_PER_TILE = 4
LOIHI_TILES = int(LOIHI_CORES / LOIHI_CORES_PER_TILE)

import os
import sys
import importlib


def serial_to_sanafe(filename: str, name:str = "converted_abstract_process.net"):
    p = load(filename)

    network = sim.Network(save_mappings=True)
    neurons_per_core = LOIHI_NEURONS_PER_CORE
    compartments = sim.init_compartments(LOIHI_TILES, LOIHI_CORES_PER_TILE,
                                        neurons_per_core)

    num_proc = len(p[0])
    print(num_proc)

    if type(p[0]) == AbstractProcess:
        process_to_sanafe(p[0])
        return
    else:
        for i in range(len(p[0])):
            dict = p[0][i].proc_params._parameters
            dim = dict["shape"]
            dim1 = dim[0]
            dim2 = 1
            if(len(dim) > 1): dim2 = dim[1]
            group = sim.create_layer_type(network=network,
                                            layer_neuron_count=dim1*dim2,
                                            compartments=compartments,
                                            neuron_parameters=dict)
            if group.id > 0: sim.connect_layers(network, group.id-1, group.id)
            prev = group.id

    network.save(filename=name)
    return

def process_to_sanafe(process: AbstractProcess, name:str = "converted_abstract_process.net"):

    network = sim.Network(save_mappings=True)
    neurons_per_core = LOIHI_NEURONS_PER_CORE
    compartments = sim.init_compartments(LOIHI_TILES, LOIHI_CORES_PER_TILE,
                                        neurons_per_core)

    dict = process.proc_params._parameters
    dim = dict['shape']
    dim1 = dim[0]
    dim2 = 1
    for i in range(dim[1]):
        group = sim.create_layer(network=network, layer_neuron_count=dim[0], compartments=compartments)
        if group.id > 0: sim.connect_layers(network, group.id-1, group.id)
        prev = group.id

    network.save(filename=name)
    return
