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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, PROJECT_DIR)
import sim

def check_mapping(network):
    return True

# Load the Keras model for the DVS gesture categorization benchmark
model = models.load_model(os.path.join(PROJECT_DIR, "etc", "dvs_challenge.h5"))
network = sim.Network()
compartments = sim.init_compartments(32, 128, 1024)

# Convert the DVS gesture categorization model to SANA-FE's format
print(model.summary())
layer0 = sim.create_layer(network, 1024, compartments)

filters = model.layers[0].get_weights()[0]
layer1 = sim.create_conv_layer(network, layer0, (32, 32, 1),
                               filters, compartments, stride=2)
filters = model.layers[1].get_weights()[0]
layer2 = sim.create_conv_layer(network, layer1, (15, 15, 16),
                               filters, compartments, stride=1)
filters = model.layers[2].get_weights()[0]
layer3 = sim.create_conv_layer(network, layer2, (13, 13, 32),
                               filters, compartments, stride=1)
filters = model.layers[3].get_weights()[0]
layer4 = sim.create_conv_layer(network, layer3, (11, 11, 64),
                               filters, compartments, stride=1)

weights = model.layers[5].get_weights()[0]
# Rearrange the Keras weights to be channel second
weights = np.reshape(weights, (9, 9, 11, 11))
weights = np.swapaxes(weights, 2, 1)
weights = np.reshape(weights, (891, 11), order='C')
layer5 = sim.create_connected_layer(network, layer4, (9, 9, 11), weights,
                                    compartments)

# Save the network to file. Note that this file may require a large amount of
#  disk space
network.save(os.path.join("dvstmp.net"))
