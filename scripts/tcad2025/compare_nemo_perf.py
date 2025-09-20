"""
Copyright (c) 2024 - The University of Texas at Austin
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
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))

sys.path.insert(0, PROJECT_DIR)
import utils

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
    network = sim.Network(save_mappings=True)
    compartments = sim.init_compartments(cores, 1, TRUENORTH_COMPARTMENTS)
    print("Creating neuron population")

    mappings = []
    for i in range(0, cores):
        m = (i, 0)
        mappings.extend((m,) * TRUENORTH_COMPARTMENTS)
    # Create neurons to fill every TrueNorth compartment, with a negative
    #  threshold and forced updates i.e., spikes every timestep
    population = sim.create_layer(network, cores*TRUENORTH_COMPARTMENTS,
                                    compartments, 0, 0, 1, 0.0, -1.0, 0.0,
                                    mappings=mappings, connections_out=1,
                                    soma_hw_name="core_soma",
                                    synapse_hw_name="core_synapses")

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
    run_command = ["mpirun", "-np", "12", NEMO_BIN_PATH,
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
    plt.rcParams.update({'font.size': 6, 'lines.markersize': 1})
    times = np.array(df.values)
    cores = np.array(df.index)
    entries = len(cores)
    #df.plot.bar(rot=0, figsize=(3.5, 1.4), color=("#ff7f0e", "#1f77b4"))

    #fig = plt.figure(figsize=(3.5, 1.4))
    fig = plt.figure(figsize=(3.7, 1.4))
    nemo_throughput = TIMESTEPS / times[:, 1]
    sanafe_throughput = TIMESTEPS / times[:, 0]
    bar1 = plt.bar(np.arange(entries) - 0.15, nemo_throughput, width=0.3)
    bar2 = plt.bar(np.arange(entries) + 0.15, sanafe_throughput, width=0.3,
            alpha=.99)
    plt.legend(("NeMo", "SANA-FE"))

    for i, rect in enumerate(bar1):
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, height, f"{nemo_throughput[i]:.1f}", ha="center", va="bottom")

    for i, rect in enumerate(bar2):
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2.0, height, f"{sanafe_throughput[i]:.1f}", ha="center", va="bottom")

    ax = plt.gca()
    plt.xlabel("TrueNorth Core Count / Total Neurons")
    ax.set_xticks(np.arange(entries))
    neuron_counts = ("8k", "16k", "32k", "64k", "128k", "256k")
    core_labels = []
    for core, neuron in zip(cores, neuron_counts):
        core_labels.append(f"{core}/{neuron}")
    ax.set_xticklabels(core_labels)
    #plt.xticks(rotation=30)
    #plt.ylabel("Run-time (s)")
    plt.ylabel("Throughput (steps per s)")
    #plt.yscale("log")
    #plt.minorticks_on()
    plt.ylim((0, 25))
    ax.tick_params(axis='y', which='minor', labelbottom=False)
    plt.tight_layout(pad=0.3)
    plt.savefig(os.path.join(PROJECT_DIR, "runs", "nemo", "compare_sanafe_nemo.png"))
    plt.savefig(os.path.join(PROJECT_DIR, "runs", "nemo", "compare_sanafe_nemo.pdf"))
    return


if __name__ == "__main__":
    run_experiments = False
    plot_experiments = True
    if run_experiments:
        core_counts = (32, 64, 128, 256, 512, 1024)
        print(f"Running experiments with following core counts: {core_counts}")

        experimental_runs = len(core_counts)
        sanafe_runtimes = np.zeros(experimental_runs)
        nemo_runtimes = np.zeros(experimental_runs)

        for i, cores in enumerate(core_counts):
            print(f"Running simulation of {cores} cores")
            sanafe_runtimes[i] = run_sim_sanafe(cores, TIMESTEPS)
            nemo_runtimes[i] = run_sim_nemo(cores, TIMESTEPS, debug=False)

        print(sanafe_runtimes)
        with open(CSV_RESULTS_FILENAME, "w") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(("cores", "SANA-FE", "NeMo"))
            for i, cores in enumerate(core_counts):
                writer.writerow((cores, sanafe_runtimes[i], nemo_runtimes[i]))

        print("Saved results to file")

    if plot_experiments:
        plot_results()
