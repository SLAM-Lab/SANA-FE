"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Run Latin Square solver benchmark (CSP solver)
"""
# TODO: First try to get the application working roughly, as proof of concept.
#  Next I will create an equivalent model for Loihi and compare the two.
import matplotlib
matplotlib.use('Agg')

import csv
import subprocess
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
import networkx as nx

# SANA-FE libraries
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, os.path.join(PROJECT_DIR))
import utils

NETWORK_FILENAME = "runs/latin/latin_square.net"
ARCH_FILENAME = "arch/loihi.arch"
LOIHI_CORES = 128
LOIHI_CORES_PER_TILE = 4
LOIHI_TILES = int(LOIHI_CORES / LOIHI_CORES_PER_TILE)
N = 4
TIMESTEPS = 100

def calculate_graph_index(row, col, digit):
    return ((row*N + col)*N) + digit

network = utils.Network()
compartments = utils.init_compartments(LOIHI_TILES, LOIHI_CORES_PER_TILE,
                                        1024)
log_spikes = 1
log_potential = 1
force_update = 1
threshold = 1.0
reset = 0.0
leak = 0.0

print(f"Creating WTA networks for {N} digits")

G = nx.DiGraph()
G.add_nodes_from(range(0, N**3))
# For every position in the square, create a WTA layer representing all possible
#  digit choices
square = []
for i in range(0, N):
    row = []
    for j in range(0, N):
        row.append(utils.create_layer(network, N, compartments,
                          log_spikes, log_potential, force_update,
                          threshold, reset, leak))
    square.append(row)


# TODO: take into account known positions. This should basically remove some
#  neuron connections since we know that *has* to be the answer. That neuron
#  should fire with a bias and not be inhibited. We can immediately prune away
#  any options that would conflict. This neuron can spike on every timestep as
#  well.

# Connect such that every digit in one position inhibits all other digits in
#  that position
connections = 0
for row in range(0, N):
    for col in range(0, N):
        pos = square[row][col]
        for digit in range(0, N):
            pre_neuron = pos.neurons[digit]
            for d in range(0, N):
                if d != digit:
                    # Add inhibiting connection for all other digits at this
                    #  position
                    post_neuron = pos.neurons[d]
                    pre_neuron.add_connection(post_neuron, -1)
                    i = calculate_graph_index(row, col, digit)
                    j = calculate_graph_index(row, col, d)
                    G.add_edge(i, j, weights=-1)
                    connections += 1

            for r in range(0, N):
                if r != row:
                    # Add inhibiting connection for this digit at all other rows
                    dest = square[r][col]
                    post_neuron = dest.neurons[digit]
                    pre_neuron.add_connection(post_neuron, -1)
                    j = calculate_graph_index(r, col, digit)
                    G.add_edge(i, j, weights=-1)
                    connections += 1

            for c in range(0, N):
                if c != col:
                    # Add inhibiting connection for this digit at other cols
                    dest = square[row][c]
                    post_neuron = dest.neurons[digit]
                    pre_neuron.add_connection(post_neuron, -1)
                    j = calculate_graph_index(row, c, digit)
                    G.add_edge(i, j, weights=-1)
                    connections += 1

print(f"Latin square network has {connections} connections")
network_path = os.path.join(PROJECT_DIR, NETWORK_FILENAME)
network.save(network_path)

if N < 4:
    pos = nx.nx_agraph.graphviz_layout(G)
    nx.draw_networkx(G, pos)
    plt.savefig(os.path.join(PROJECT_DIR, "runs/latin/latin_net.png"))

# Now execute the network using SANA-FE and extract the spike timings
arch_path = os.path.join(PROJECT_DIR, ARCH_FILENAME)
run_command = ("./sim", "-s", "-v", arch_path, network_path, f"{TIMESTEPS}")
print("sana-fe command: {0}".format(" ".join(run_command)))
subprocess.call(run_command)

# Use spiking data to create the grid solution produced by the Loihi run
#with open
with open(os.path.join(PROJECT_DIR, "spikes.trace")) as spikes:
    reader = csv.reader(spikes)
    header = next(reader)

    spike_counts = np.zeros((N, N, N))
    for spike in reader:
        gid_str, nid_str = spike[0].split(".")
        gid, nid = int(gid_str), int(nid_str)
        timestep = int(spike[1])

        digit = nid
        col = gid % N
        row = gid // N
        assert(digit < N)
        assert(col < N)
        assert(r < N)
        spike_counts[row][col][digit] += 1

print(spike_counts)
chosen_digits = np.argmax(spike_counts, axis=2)

# Verify the solution is correct

# Plot a grid and fill in the numbers based on the largest number of
#  spikes collected after a fixed point
plt.figure(figsize=(1, 1))
fig, ax = plt.subplots(1, 1)
ax.axis('tight')
ax.axis('off')
ax.table(cellText=chosen_digits, colWidths=[0.1] * N*N, cellLoc="center",
         loc="center")
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig(os.path.join(PROJECT_DIR, "runs/latin/latin_square.png"))

df = pd.read_csv("potential.trace")
plt.figure()
df.plot()
plt.savefig(os.path.join(PROJECT_DIR, "runs/latin/latin_potentials.png"))