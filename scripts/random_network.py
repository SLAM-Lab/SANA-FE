"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.
"""
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
import csv
import yaml

# SANA-FE libraries
import sys
sys.path.insert(0, '/home/usr1/jboyle/neuro/sana-fe')
import utils

# Use a dumb seed to get consistent results
random.seed(1)

# Global experiment parameters
NETWORK_FILENAME = "runs/random.net"
ARCH_FILENAME = "loihi.arch"
LOIHI_CORES = 128
LOIHI_CORES_PER_TILE = 4
LOIHI_TILES = int(LOIHI_CORES / LOIHI_CORES_PER_TILE)
TIMESTEPS = 1
#TIMESTEPS = 4

def create_random_network(cores, neurons_per_core, messages_per_neuron,
                          spikes_per_message):
    network = utils.Network()
    compartments = utils.init_compartments(LOIHI_TILES, LOIHI_CORES_PER_TILE,
                                           1024)

    neurons = cores * neurons_per_core
    mappings = []
    for i in range(0, cores):
        m = (i/4, i%4)
        mappings.extend((m,) * neurons_per_core)

    print("Creating neuron population")
    population = utils.create_layer(network, neurons,
                                    compartments, 0, 0, 1, 0.0, 0.0, 0.0,
                                    mappings=mappings)

    print("Generating randomized network connections")
    weight = 1.0
    print(f"Cores: {cores}, messages per neuron: {messages_per_neuron}")
    print(f"Spikes per message: {spikes_per_message}")
    for n in range(0, neurons):
        src = population.neurons[n]
        # All neurons with outgoing connections should fire every timestep
        src.add_bias(1.0)
        if (n % 1024) == 0:
            print(f"Generating synaptic connections for neuron {n}")

        dest_core = random.sample(range(0, cores), messages_per_neuron)
        for c in dest_core:
            dest_neurons = random.sample(range(0, 1024), spikes_per_message)
            for d in dest_neurons:
                dest_id = (c * 1024) + d
                assert(dest_id < neurons)
                dest = population.neurons[dest_id]
                src.add_connection(dest, weight)

    network.save(NETWORK_FILENAME)


# Run the simulation on SANA-FE, generating the network and immediately using it
#  Return the total runtime measured by Python, including setup and processing
#  time.
def run_sim(timesteps, cores, messages_per_core, spikes_per_message):
    create_random_network(cores, 1024, messages_per_core, spikes_per_message)
    run_command = ("./sim", ARCH_FILENAME, NETWORK_FILENAME,
                   "{0}".format(timesteps))
    print("sana-fe command: {0}".format(" ".join(run_command)))
    subprocess.call(run_command)

    with open("stats.yaml", "r") as summary_file:
        summary = yaml.safe_load(summary_file)

    return summary


def plot_results():
    df = pd.read_csv("runs/sanafe_perf.csv")
    plt.rcParams.update({'font.size': 12, 'lines.markersize': 5})

    plt.figure(figsize=(4.5, 4.5))
    for cores in df["cores"].unique():
        plt.plot(df.loc[(df["spikes_per_message"] == 1) & (df["cores"] == cores), "messages"],
                 df.loc[(df["spikes_per_message"] == 1) & (df["cores"] == cores), "runtime"],
                 "o-")

    legend_str = [f"{c} cores" for c in df["cores"].unique()]
    plt.legend(legend_str, reverse=True)
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0,0))
    plt.xlabel("Spike Messages")
    plt.ylabel("Run-time (s)")
    plt.tight_layout()
    plt.savefig("runs/sanafe_perf.png")
    plt.close()

    plt.figure(figsize=(4.5, 4.5))
    for spikes_per_message in df["spikes_per_message"].unique():
        print(spikes_per_message)
        print(df.loc[(df["spikes_per_message"] == spikes_per_message) & (df["cores"] == cores)])
        plt.plot(df.loc[(df["spikes_per_message"] == spikes_per_message) & (df["cores"] == cores), "messages"],
                 df.loc[(df["spikes_per_message"] == spikes_per_message) & (df["cores"] == cores), "runtime"],
                 "o-")

    legend_str = [f"{s} spikes per message" for s in df["spikes_per_message"].unique()]
    plt.legend(legend_str)
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0,0))
    plt.xlabel("Spike Messages")
    plt.ylabel("Run-time (s)")
    plt.tight_layout()
    plt.savefig("runs/sanafe_perf2.png")
    plt.close()

    return


if __name__ == "__main__":
    run_experiments = False
    plot = True
    if run_experiments:
        cores = (16, 32, 64, 128)
        messages_per_neuron = (1, 2, 4, 8, 16, 32, 64, 128)
        spikes_per_message = (1, 4, 8)

        # Part one is we try different numbers of cores and different numbers
        #  of messages per core
        data_points = []
        with open("runs/sanafe_perf.csv", "w") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(("cores", "messages_per_neuron", "spikes_per_message",
                             "messages", "spikes", "runtime"))
            for c in cores:
                for m in messages_per_neuron:
                    if m <= c:
                        for s in spikes_per_message:
                            results = run_sim(TIMESTEPS, c, m, s)
                            row = (c, m, s, results["total_packets"],
                                   results["total_spikes"],
                                   results["wall_time"])
                            writer.writerow(row)
        print("Saved results to file")

    if plot:
        plot_results()
