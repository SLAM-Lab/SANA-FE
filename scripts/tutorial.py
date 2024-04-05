"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

tutorial.py: NICE 2024 Tutorial
"""
import sys
import os
#from tensorflow.keras import models
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
NETWORK_PATH = os.path.join(PROJECT_DIR, "snn", "dvs_challenge.net")
ARCH_PATH = os.path.join(PROJECT_DIR, "arch", "loihi.yaml")
TIMESTEPS = 100

sys.path.insert(0, PROJECT_DIR)
import sim

def check_mapping(network):
    # Check that all neurons are mapped to a core
    neurons_mapped_to_cores = [0 for _ in range(0, 128)]
    for n in network.neurons:
        if n.core is None:
            print("Error: Neuron n{n.id} is not mapped!")
            return False
        else:
            neurons_mapped_to_cores[n.core] += 1

    for core in range(0, 128):
        if neurons_mapped_to_cores[core] > 1024:
            print("Error: Core {core} has > 1024 neurons mapped")
            return False

    return True

network = sim.Network()
arch = sim.Architecture()

kernels = np.load("dvs_challenge.npz")
# TODO: save the updated model so we just have to load the model and the weights
# Convert the DVS gesture categorization model to SANA-FE's format

biases = np.zeros((1024))
layer0 = sim.create_layer(network, 1024, threshold=255, biases=biases)
# TODO: set the biases for one input frame
layer1 = sim.create_conv_layer(network, layer0, (32, 32, 1),  # 3600 neurons
                               kernels["conv1"], stride=2, threshold=420,
                               )
layer2 = sim.create_conv_layer(network, layer1, (15, 15, 16),  # 5408 neurons
                               kernels["conv2"], stride=1, threshold=351)
layer3 = sim.create_conv_layer(network, layer2, (13, 13, 32),  # 7744 neurons
                               kernels["conv3"], stride=1, threshold=276)
layer4 = sim.create_conv_layer(network, layer3, (11, 11, 64),  # 891 neurons
                               kernels["conv4"], stride=1, threshold=371)
layer5 = sim.create_connected_layer(network, layer4, (9, 9, 11),  # 11 neurons
                                    kernels["dense1"], threshold=341)

sim.map_neuron_group_to_cores(layer0, arch, 1)
sim.map_neuron_group_to_cores(layer1, arch, 4)
sim.map_neuron_group_to_cores(layer2, arch, 16)
sim.map_neuron_group_to_cores(layer3, arch, 16)
sim.map_neuron_group_to_cores(layer4, arch, 4)
sim.map_neuron_group_to_cores(layer5, arch, 1)

# Save the network to file. Note that this file may require a large amount of
#  disk space
network.save(os.path.join(PROJECT_DIR, "runs", "dvs_challenge.net"),
             save_mappings=True)


sim.run(ARCH_PATH, NETWORK_PATH, TIMESTEPS)