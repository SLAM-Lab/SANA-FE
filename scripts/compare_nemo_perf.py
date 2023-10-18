"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.
"""
# TODO: I seem to be getting an issue where NeMo is freezing after a number of
#  timestep / ticks

# External libraries, plotting
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# Other external libraries
import numpy as np
import pandas as pd

# Python built-in libraries
import subprocess
import random
import time
import csv
import sys
import os

# SANA-FE libraries
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))

sys.path.insert(0, PROJECT_DIR)
import sim

# Use a dumb seed to get consistent results
random.seed(1)

# Global experiment parameters
TRUENORTH_COMPARTMENTS = 256
TRUENORTH_AXONS = TRUENORTH_COMPARTMENTS
SPIKE_INTRA_CORE_PROB = 0.8
NETWORK_FILENAME = os.path.join(PROJECT_DIR, "runs", "nemo",
                                "nemo_randomized.net")
ARCH_FILENAME = os.path.join(PROJECT_DIR, "arch", "truenorth.yaml")
TIMESTEPS = 10  # i.e., ticks, where each tick is 1 ms of wall-time on the chip
NEMO_BIN_PATH = "/home/usr1/jboyle/neuro/nemo/NeMo/bin/NeMo"
CSV_RESULTS_FILENAME = os.path.join(PROJECT_DIR, "runs", "nemo",
                                    "compare_sanafe_nemo.csv")

# Create a random truenorth network, 80% connected to neuron within same
# core, 20% connected to neurons outside
def create_nemo_network(cores):
    network = sim.Network()
    compartments = sim.init_compartments(cores, 1,
                                           TRUENORTH_COMPARTMENTS)
    print("Creating neuron population")

    mappings = []
    for i in range(0, cores):
        m = (i, 0)
        mappings.extend((m,) * TRUENORTH_COMPARTMENTS)
    # Create neurons to fill every TrueNorth compartment, with a negative
    #  threshold and forced updates i.e., spikes every timestep
    population = sim.create_layer(network, cores*TRUENORTH_COMPARTMENTS,
                                    compartments, 0, 0, 1, 0.0, -1.0, 0.0,
                                    mappings=mappings, connections_out=1)

    print("Generating randomized network connections")
    weight = 1
    for c in range(0, cores):
        if (c % 32) == 0:
            print(f"Generating synaptic connections for core {c}")
        for n in range(0, TRUENORTH_AXONS):
            if random.random() < SPIKE_INTRA_CORE_PROB:
                possible_cores = list(range(0, cores))
                del(possible_cores[c])
                dest_core = random.choice(possible_cores)
            else:  # 20% chance of picking the same core
                dest_core = c
            dest_axon = random.randrange(0, TRUENORTH_AXONS)
            src = population.neurons[(c*TRUENORTH_AXONS) + n]
            dest = population.neurons[(dest_core*TRUENORTH_AXONS) + dest_axon]
            src.add_connection(dest, weight)

    network.save(NETWORK_FILENAME)


# Run the simulation on SANA-FE, generating the network and immediately using it
#  Return the total runtime measured by Python, including setup and processing
#  time.
def run_sim_sanafe(cores, timesteps):
    create_nemo_network(cores)
    start = time.time()
    sim.run(ARCH_FILENAME, NETWORK_FILENAME, timesteps)
    end = time.time()
    run_time = end - start
    print(f"sanafe runtime for {cores} cores was {run_time} s")
    return run_time


# Run the same simulation on NeMo, for a given number of cores and timesteps
#  Return the runtime measured by Python.
# TODO: should we add changes to measure the runtime of simulation and ignore
#  setup time? Is this even possible
def run_sim_nemo(cores, timesteps, debug=True):
    run_command = ["mpirun", NEMO_BIN_PATH,
                   f"--cores={cores}", f"--end={timesteps}",  "--sync=3",
                   "--rand"]
    if debug:
        run_command.append("--svouts")
    print("NeMo command: {0}".format(" ".join(run_command)))
    start = time.time()
    subprocess.call(run_command)
    end = time.time()
    run_time = end - start
    print(f"nemo runtime for {cores} cores was {run_time} s")

    if debug:
        # See how many spikes were sent, i.e., how many line entries are in the
        #  firing logs (1 spike per line). We just want to check we generate
        #  roughly the same number of spikes as sana-fe
        num_spikes = 0
        for i in range(0, 4):
            num_spikes += sum(1 for line in open(f"fire_record_rank_{i}.csv"))
        print(f"NeMo {num_spikes} total spikes")

    return run_time


def plot_results():
    df = pd.read_csv(CSV_RESULTS_FILENAME, index_col="cores")
    plt.rcParams.update({'font.size': 7, 'lines.markersize': 1})
    df.plot.bar(rot=0, figsize=(3.5, 1.6), color=("#ff7f0e", "#1f77b4"))
    plt.xlabel("TrueNorth Core Count")
    plt.ylabel("Run-time (s)")
    plt.minorticks_on()
    plt.tight_layout(pad=0.3)
    plt.savefig(os.path.join(PROJECT_DIR, "runs", "nemo", "compare_sanafe_nemo.png"))
    plt.savefig(os.path.join(PROJECT_DIR, "runs", "nemo", "compare_sanafe_nemo.pdf"))
    return


if __name__ == "__main__":
    run_experiments = False
    plot = True
    if run_experiments:
        core_counts = (32, 64, 128, 256, 512, 1024)
        print(f"Running experiments with following core counts: {core_counts}")

        experimental_runs = len(core_counts)
        sanafe_runtimes = np.zeros(experimental_runs)
        nemo_runtimes = np.zeros(experimental_runs)

        for i, cores in enumerate(core_counts):
            print(f"Running simulation of {cores} cores")
            sanafe_runtimes[i] = run_sim_sanafe(cores, TIMESTEPS)
            #nemo_runtimes[i] = run_sim_nemo(cores, TIMESTEPS)

        print(sanafe_runtimes)
        #with open(CSV_RESULTS_FILENAME, "w") as csv_file:
        #    writer = csv.writer(csv_file)
        #    writer.writerow(("cores", "SANA-FE", "NeMo"))
        #    for i, cores in enumerate(core_counts):
        #        writer.writerow((cores, sanafe_runtimes[i], nemo_runtimes[i]))

        print("Saved results to file")

    if plot:
        plot_results()
