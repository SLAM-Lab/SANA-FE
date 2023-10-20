"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Run Latin Square solver benchmark (CSP solver)
"""
import matplotlib
matplotlib.use('Agg')

import csv
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
import sim

ARCH_FILENAME = "arch/loihi_latin.yaml"
LOIHI_CORES = 128
LOIHI_CORES_PER_TILE = 4
LOIHI_TILES = int(LOIHI_CORES / LOIHI_CORES_PER_TILE)
LOIHI_COMPARTMENTS = 1024
TIMESTEPS = 10240

def calculate_graph_index(N, row, col, digit):
    return ((row*N + col)*N) + digit

def latin_square(N, tiles=LOIHI_TILES, cores_per_tile=LOIHI_CORES_PER_TILE,
                 neurons_per_core=LOIHI_COMPARTMENTS):
    network = sim.Network(save_mappings=True)
    compartments = sim.init_compartments(tiles, cores_per_tile,
                                         neurons_per_core)
    print(f"Creating WTA networks for {N} digits")
    #G = nx.DiGraph()
    #G.add_nodes_from(range(0, N**3))

    # For every position in the square, create a WTA layer representing all
    #  possible digit choices
    square = []
    for i in range(0, N):
        row = []
        for j in range(0, N):
            wta = sim.create_layer(network, N, compartments,
                                   log_spikes=False,
                                   log_potential=False,
                                   force_update=False,
                                   threshold=64.0,
                                   reset=0.0,
                                   leak=1,
                                   reverse_threshold=-2**7 + 1.0,
                                   reverse_reset_mode="saturate")
            for neuron in wta.neurons:
                neuron.add_bias(1 * 2**7)
            row.append(wta)
        square.append(row)

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
                        i = calculate_graph_index(N, row, col, digit)
                        j = calculate_graph_index(N, row, col, d)
                        # G.add_edge(i, j, weights=-1)
                        connections += 1

                for r in range(0, N):
                    if r != row:
                        # Add inhibiting connection for this digit at all other
                        #  rows
                        dest = square[r][col]
                        post_neuron = dest.neurons[digit]
                        pre_neuron.add_connection(post_neuron, -1)
                        j = calculate_graph_index(N, r, col, digit)
                        # G.add_edge(i, j, weights=-1)
                        connections += 1

                for c in range(0, N):
                    if c != col:
                        # Add inhibiting connection for this digit at other cols
                        dest = square[row][c]
                        post_neuron = dest.neurons[digit]
                        pre_neuron.add_connection(post_neuron, -1)
                        j = calculate_graph_index(N, row, c, digit)
                        # G.add_edge(i, j, weights=-1)
                        connections += 1

    print(f"Latin square network has {connections} connections")
    network_filename = os.path.join("runs", "dse", f"latin_square_N{N}.net")
    network_path = os.path.join(PROJECT_DIR, network_filename)
    network.save(network_path)


def plot_results(N, network_path):
    if N < 4:
        pos = nx.nx_agraph.graphviz_layout(G)
        nx.draw_networkx(G, pos)
        plt.savefig(os.path.join(PROJECT_DIR, "runs/latin/latin_net.png"))

    # Now execute the network using SANA-FE and extract the spike timings
    arch_path = os.path.join(PROJECT_DIR, ARCH_FILENAME)
    sim.run(arch_path, network_path, TIMESTEPS,
            spike_trace=True, potential_trace=True)

    # Use spiking data to create the grid solution produced by the Loihi run
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


def run_experiment(network_filename):
    arch_path = os.path.join(PROJECT_DIR, ARCH_FILENAME)
    network_path = os.path.join(PROJECT_DIR, network_filename)
    results = sim.run(arch_path, network_path, TIMESTEPS,
                      spike_trace=False, potential_trace=False)

    return results


if __name__ == "__main__":
    run_experiments = False
    plot_experiment = True

    if run_experiments:
        if (os.path.isfile(os.path.join(PROJECT_DIR, "runs", "latin",
                           "loihi_latin.csv"))):
            open(os.path.join(PROJECT_DIR, "runs", "latin", "sim_latin.csv"),
                 "w")
            with open(os.path.join(PROJECT_DIR, "runs", "latin",
                                    "loihi_latin.csv")) as latin_squares_file:
                reader = csv.DictReader(latin_squares_file)

                for line in reader:
                    # Each line of loihi_latin.csv is another experiment,
                    #  containing the network to run and the results measured
                    #  on Loihi
                    results = run_experiment(line["network"])
                    time = results["time"] / TIMESTEPS
                    energy = results["energy"] / TIMESTEPS
                    row = (line["N"], line["network"], energy, time)
                    with open(os.path.join(PROJECT_DIR, "runs/latin/sim_latin.csv"),
                              "a") as csv_file:
                        writer = csv.writer(csv_file)
                        writer.writerow(row)

    if plot_experiment:
        sim_df = pd.read_csv(os.path.join(PROJECT_DIR,
                                          "runs", "latin", "sim_latin.csv"))
        loihi_df = pd.read_csv(os.path.join(PROJECT_DIR,
                                            "runs", "latin", "loihi_latin.csv"))
        df = pd.merge(sim_df, loihi_df)
        sim_energy = df["sim_energy"].values * 1.0e6
        loihi_energy = df["loihi_energy"].values * 1.0e6
        sim_latency = df["sim_latency"].values * 1.0e6
        loihi_latency = df["loihi_latency"].values * 1.0e6

        # Plot the simulated vs measured energy
        plt.rcParams.update({"font.size": 7, "lines.markersize": 5})
        plt.figure(figsize=(1.7, 1.7))
        plt.minorticks_on()
        plt.gca().set_box_aspect(1)

        plt.plot(sim_energy, loihi_energy, "x", mew=1.5)
        plt.plot(np.linspace(min(sim_energy), max(sim_energy)),
                 np.linspace(min(sim_energy), max(sim_energy)), "k--")
        plt.xlabel("Simulated Energy ($\mu$J)")
        plt.ylabel("Measured Energy ($\mu$J)")
        plt.xticks(np.arange(0, 1.1, 0.4))
        plt.yticks(np.arange(0, 1.1, 0.4))
        plt.tight_layout(pad=0.3)

        plt.savefig(os.path.join(PROJECT_DIR, "runs", "latin",
                                 "latin_energy.pdf"))
        plt.savefig(os.path.join(PROJECT_DIR, "runs", "latin",
                                 "latin_energy.png"))

        # Plot the simulated vs measured latency
        plt.figure(figsize=(1.7, 1.7))
        plt.minorticks_on()
        plt.gca().set_box_aspect(1)

        plt.plot(sim_latency, loihi_latency, "x", mew=1.5)
        plt.plot(np.linspace(min(sim_latency), max(sim_latency)),
                 np.linspace(min(sim_latency), max(sim_latency)), "k--")
        plt.xlabel("Simulated Latency ($\mu$s)")
        plt.ylabel("Measured Latency ($\mu$s)")
        plt.xticks(np.arange(0, 41, 20))
        plt.yticks(np.arange(0, 41, 20))
        plt.tight_layout(pad=0.3)

        plt.savefig(os.path.join(PROJECT_DIR, "runs", "latin",
                                 "latin_latency.pdf"))
        plt.savefig(os.path.join(PROJECT_DIR, "runs", "latin",
                                 "latin_latency.png"))

        absolute_latency_error = np.abs(loihi_latency - sim_latency) / loihi_latency
        absolute_energy_error = np.abs(loihi_energy - sim_energy) / loihi_energy

        print(f"latency absolute mean error: {np.mean(absolute_latency_error) * 100.0}")
        print(f"energy absolute mean {np.mean(absolute_energy_error) * 100.0}")

        total_latency_error = (np.sum(loihi_latency) - np.sum(sim_latency)) / np.sum(loihi_latency)
        total_energy_error = (np.sum(loihi_energy) - np.sum(sim_energy)) / np.sum(loihi_energy)

        print(f"total latency error: {total_latency_error * 100.0}%")
        print(f"total energy error: {total_energy_error * 100.0}%")
