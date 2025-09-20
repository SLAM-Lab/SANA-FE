import numpy as np
import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, os.path.join(PROJECT_DIR))
import sanafe

net = sanafe.Network()
weights = np.ones((9, 9))
layer1 = net.create_neuron_group("in", 2, {})
layer2 = net.create_neuron_group("out", 2, {})

print(layer1)
print(layer2)

layer1.connect_neurons_sparse(layer2, {}, [(0, 0), (0, 1)])

groups = net.groups
names = list(net.groups.keys())
print(groups)
print(names)
print(groups["in"].neurons)
print(groups["in"].neurons[0])
