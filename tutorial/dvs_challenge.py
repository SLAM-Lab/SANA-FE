"""
Copyright (c) 2024 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

dvs_challenge.py: NICE 2024 Tutorial Challenge
"""
import sys
import os

# NumPy is the only external Python dependency for this script
import numpy as np

# Useful directories, relative to the current path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
NETWORK_PATH = os.path.join(PROJECT_DIR, "tutorial", "dvs_challenge.net")
ARCH_PATH = os.path.join(PROJECT_DIR, "tutorial", "loihi.yaml")
TIMESTEPS = 1000

sys.path.insert(0, PROJECT_DIR)
import sim

network = sim.Network()
arch = sim.Architecture()

# Load the convolutional kernel weights, thresholds and input biases from file.
#  If using the Docker container, this file is included in the image.
#  Otherwise, this file is also hosted on Google Drive and can be downloaded
#  prior to running this script using the command:
#     wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1WkbJZFasTe-v8vTYXrUaz_1e-p_xHMEj'
snn = np.load(os.path.join(PROJECT_DIR, "tutorial", "dvs_challenge.npz"))

# Convert the DVS gesture categorization model to SANA-FE's SNN format
biases = snn["inputs"]
thresholds = snn["thresholds"]
layer0 = sim.create_layer(network, 1024, threshold=thresholds[0], biases=biases)
layer1 = sim.create_conv_layer(network, layer0, (32, 32, 1),  # 3600 neurons
                               snn["conv1"], stride=2, threshold=thresholds[1],
                               )
layer2 = sim.create_conv_layer(network, layer1, (15, 15, 16),  # 5408 neurons
                               snn["conv2"], stride=1, threshold=thresholds[2])
layer3 = sim.create_conv_layer(network, layer2, (13, 13, 32),  # 7744 neurons
                               snn["conv3"], stride=1, threshold=thresholds[3])
layer4 = sim.create_conv_layer(network, layer3, (11, 11, 64),  # 891 neurons
                               snn["conv4"], stride=1, threshold=thresholds[4])
layer5 = sim.create_connected_layer(network, layer4, (9, 9, 11),  # 11 neurons
                                    snn["dense1"], threshold=thresholds[5])

# Map the SNN to Loihi cores. Each layer is represented by its own neuron group
### CHANGE THESE LINES TO TEST DIFFERENT MAPPINGS ###
sim.map_neuron_group_to_cores(layer0, arch, 1)
sim.map_neuron_group_to_cores(layer1, arch, 4)
sim.map_neuron_group_to_cores(layer2, arch, 16)
sim.map_neuron_group_to_cores(layer3, arch, 16)
sim.map_neuron_group_to_cores(layer4, arch, 4)
sim.map_neuron_group_to_cores(layer5, arch, 1)
### END OF LINES TO CHANGE ###

# Save the network to file. Note that this file may require a large amount of
#  disk space (~100 MB)
network.save(NETWORK_PATH, save_mappings=True)

# Run the network you just generated on Loihi
# Comment out this line to stop simulation runs
sim.run(ARCH_PATH, NETWORK_PATH, TIMESTEPS)
