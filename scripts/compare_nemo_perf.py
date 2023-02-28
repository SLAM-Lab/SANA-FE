import matplotlib
matplotlib.use('Agg')

import subprocess
from matplotlib import pyplot as plt
import random

import sys
sys.path.insert(0, '/home/usr1/jboyle/neuro/sana-fe')
import utils

random.seed(1)

TRUENORTH_COMPARTMENTS = 256
TRUENORTH_AXONS = TRUENORTH_COMPARTMENTS
#TRUENORTH_TILES = 4096
TRUENORTH_TILES = 32
SPIKE_INTRA_CORE_PROB = 0.8
NETWORK_FILENAME = "runs/nemo_randomized.net"
ARCH_FILENAME = "truenorth.arch"

# Create a random truenorth network, 80% connected to neuron within same
# core, 20% connected to neurons outside
def create_nemo_network():
    network = utils.Network()
    compartments = utils.init_compartments(TRUENORTH_TILES, 1,
                                           TRUENORTH_COMPARTMENTS)
    print("Creating neuron population")
    
    mappings = []
    for i in range(0, TRUENORTH_TILES):
        m = (i, 0)
        mappings.extend((m,) * TRUENORTH_COMPARTMENTS)
    print(len(mappings))
    # Create neurons to fill every TrueNorth compartment, with a negative
    #  threshold and forced updates i.e., spikes every timestep
    population = utils.create_layer(network,
                                    TRUENORTH_TILES*TRUENORTH_COMPARTMENTS,
                                    compartments, 0, 0, 1, 0.0, -1.0, 0.0,
                                    mappings=mappings)


    print("Generating randomized network connections")
    weight = 1.0
    for c in range(0, TRUENORTH_TILES):
        if (c % 32) == 0:
            print(f"Generating synaptic connections for core {c}")
        for n in range(0, TRUENORTH_AXONS):
            if random.random() < SPIKE_INTRA_CORE_PROB:
                possible_cores = list(range(0, TRUENORTH_TILES))
                del(possible_cores[c])
                dest_core = random.choice(possible_cores)
            else:  # 20% chance of picking the same core
                dest_core = c
            dest_axon = random.randrange(0, TRUENORTH_AXONS)
            src = population.neurons[(c*TRUENORTH_AXONS) + n]
            dest = population.neurons[(dest_core*TRUENORTH_AXONS) + dest_axon]
            src.add_connection(dest, weight)

    network.save(NETWORK_FILENAME)


# Run the simulation on SANA-FE
def run_sim_sanafe(timesteps):
    run_command = ("./sim", ARCH_FILENAME, NETWORK_FILENAME,
                   "{0}".format(timesteps))
    print("Command: {0}".format(" ".join(run_command)))
    subprocess.call(run_command)

    return


def run_sim_nemo(timesteps):
    pass


# 3) plot runtime results against runtime results for NeMo, varying the core
#    count


if __name__ == "__main__":
    create_nemo_network()
    run_sim_sanafe(100)

