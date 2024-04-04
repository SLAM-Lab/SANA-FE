"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

tutorial.py: NICE 2024 Tutorial
"""
import sys
import os
from tensorflow.keras import models
import numpy as np
import logging

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, PROJECT_DIR)
import sim

def check_mapping(network):
    # Check that all neurons are mapped to a core
    for n in network.neurons:
        if n.core is None:
            logging.error("Neuron n{n.id} is not mapped!")
            return False
        else:
            neurons_mapped_to_cores[n.core] += 1

    for core in range(0, 128):
        if neurons_mapped_to_cores[core] > 1024:
            logging.error("Core {core} has > 1024 neurons mapped")
            return False

    return True

scale = (420.05236577257483, 351.1046444780251, 276.6147837631879, 371.60317670987195, 341.41679600239286)

# Load the Keras model for the DVS gesture categorization benchmark
model = models.load_model(os.path.join(PROJECT_DIR, "etc", "dvs_challenge.h5"))
network = sim.Network()
compartments = sim.init_compartments(32, 128, 1024)


# TODO: save the updated model so we just have to load the model and the weights
# Convert the DVS gesture categorization model to SANA-FE's format
print(model.summary())
layer0 = sim.create_layer(network, 1024, compartments, threshold=scale[0])

filters = np.rint(model.layers[0].get_weights()[0] * scale[0])
layer1 = sim.create_conv_layer(network, layer0, (32, 32, 1),
                               filters, compartments, stride=2,
                               threshold=255)
filters = np.rint(model.layers[1].get_weights()[0] * scale[1])
layer2 = sim.create_conv_layer(network, layer1, (15, 15, 16),
                               filters, compartments, stride=1,
                               threshold=420)
filters = np.rint(model.layers[2].get_weights()[0] * scale[2])
layer3 = sim.create_conv_layer(network, layer2, (13, 13, 32),
                               filters, compartments, stride=1,
                               threshold=351)
filters = np.rint(model.layers[3].get_weights()[0] * scale[3])
layer4 = sim.create_conv_layer(network, layer3, (11, 11, 64),
                               filters, compartments, stride=1,
                               threshold=277)

weights = np.rint(model.layers[5].get_weights()[0] * scale[4])
# Rearrange the Keras weights to be channel second
weights = np.reshape(weights, (9, 9, 11, 11))
weights = np.swapaxes((weights), 0, 1)
weights = np.swapaxes(weights, 2, 1)
weights = np.reshape(weights, (891, 11), order='C')
layer5 = sim.create_connected_layer(network, layer4, (9, 9, 11), weights,
                                    compartments)

# Save the network to file. Note that this file may require a large amount of
#  disk space
network.save(os.path.join("etc", "dvs_challenge.net"))
