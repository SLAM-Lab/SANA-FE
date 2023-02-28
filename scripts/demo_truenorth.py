"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Validate and demo TrueNorth neuron models similar to the method from:
"NeMo: A Massively Parallel Discrete-Event SimulationModel for Neuromorphic
Architectures" M. Plagge (2016)

Two different experiments showing 1) Izhikevich phasic spiking and
2) Izhikevich tonic bursting. The main purpose of this experiment is to
cross-validate my TrueNorth implementation.
"""
import matplotlib
matplotlib.use('Agg')

import sys
import subprocess
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
sys.path.insert(0, '/home/usr1/jboyle/neuro/sana-fe')
import utils

ARCH_FILENAME = "truenorth.arch"
def run_sim(network_filename, timesteps, plot_filename):
    run_command = ("./sim", ARCH_FILENAME, network_filename,
               "{0}".format(timesteps))
    print("Running command: {0}".format(" ".join(run_command)))
    subprocess.call(run_command)

    potential_data = pd.read_csv("probe_potential.csv")
    spike_data = pd.read_csv("probe_spikes.csv")

    potentials = potential_data.loc[:, "1.0"]
    # TODO: remove this kind of hacky change! We make the smallest possible
    #  value 0 I guess?
    #norm_potentials = potentials + 15.0
    norm_potentials = potentials
    spike_idx = spike_data.loc[:, "1.0"].to_numpy().nonzero()[0]
    #spike_vals = 20.0 * np.ones(spike_idx.shape)
    spike_vals = np.ones(spike_idx.shape)

    plt.figure()
    plt.plot(norm_potentials)
    plt.scatter(spike_idx, spike_vals, marker='^', color='red')
    plt.xlabel("Simulation Ticks")
    plt.ylabel("Normalized Membrane Potential")
    #plt.ylim((0, 22))
    plt.tight_layout()
    # Need to save
    plt.savefig(plot_filename)

    # Probe potential and probe spikes to get the data
    return (potential_data, spike_data)


def run():
    #run_sim("examples/truenorth_phasic.net", 1000, "runs/phasic.png")
    run_sim("examples/truenorth_bursting.net", 1000, "runs/bursting.png")


if __name__ == "__main__":
    run()
