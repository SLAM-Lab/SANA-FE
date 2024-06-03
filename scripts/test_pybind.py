import numpy as np
import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, os.path.join(PROJECT_DIR))
sys.path.insert(0, os.path.join(PROJECT_DIR))
import sanafecpp
import utils

net = sanafecpp.Network()
weights = np.ones((5000, 5000))
layer1 = net.create_neuron_group(5000, {})
layer2 = sim.create_connected_layer(net, layer1, weights, {})

print(layer1)
print(layer2)