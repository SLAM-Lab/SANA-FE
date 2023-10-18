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
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import os

import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
sys.path.insert(0, os.path.join(PROJECT_DIR))
import sim

ARCH_PATH = os.path.join(PROJECT_DIR, "arch", "truenorth.yaml")


def run_sim(network_path, timesteps, plot_filename):
    print(ARCH_PATH)
    sim.run(ARCH_PATH, network_path, timesteps, potential_trace=True,
            spike_trace=True)

    potential_data = pd.read_csv("potential.trace")
    spike_data = pd.read_csv("spikes.trace")

    offset=200
    potentials = potential_data.loc[offset:timesteps, "1.0"]

    spikes_in = spike_data.loc[spike_data["gid.nid"] == 0.0]
    spikes_out = spike_data.loc[spike_data["gid.nid"] == 1.0]

    plt.rcParams.update({'font.size': 7, "lines.markersize": 2})
    plt.figure(figsize=(3.2, 1.6))
    plt.plot(np.arange(0, timesteps-offset), potentials)
    spike_idx = spikes_out.loc[:, "timestep"]
    height = max(potentials) + 2
    spike_vals = height * np.ones(spike_idx.shape)
    plt.scatter(spike_idx-offset, spike_vals, marker='^', color='red')
    spike_idx = spikes_in.loc[:, "timestep"]
    spike_vals = (min(potentials) - 1.0) * np.ones(spike_idx.shape)
    plt.scatter(spike_idx-offset, spike_vals, marker='^', color='black')
    plt.xlabel("Simulation Ticks")
    plt.ylabel("Membrane Potential")
    #plt.ylim((0, 22))
    plt.tight_layout()
    # Need to save
    plt.savefig(plot_filename)

    # Probe potential and probe spikes to get the data
    return (potential_data, spike_data)


def run_experiment():
    run_sim("snn/nemo/truenorth_phasic.net", 1200, "runs/nemo/phasic.png")
    run_sim("snn/nemo/truenorth_bursting.net", 1200, "runs/nemo/bursting.pdf")

if __name__ == "__main__":
    run_experiment()
