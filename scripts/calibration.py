"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Calibrating SANA-FE against real-world hardware

Use several partitions of a small benchmark to calibrate
the simulator.
"""
import matplotlib
matplotlib.use('Agg')

import csv
from matplotlib import pyplot as plt
import pickle
import sys
import os
import numpy as np
import pandas as pd
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))

sys.path.insert(0, PROJECT_DIR)
import sim

MAX_TILES = 32
MAX_CORES = 4
MAX_COMPARTMENTS = 1024
NETWORK_FILENAME = "runs/calibration/connected_layers.net"
ARCH_FILENAME = "arch/loihi.yaml"

import random
def fully_connected(layer_neuron_count, spiking=True, force_update=False,
                    connection_probability=1.0):
    # Two layers, fully connected
    network = sim.Network(save_mappings=True)
    loihi_compartments = sim.init_compartments(32, 4, 1024)

    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neuron_count

    reset = 0
    log_spikes = False
    log_potential = False
    force_update = False

    # Create layers
    layer_1 = sim.create_layer(network, layer_neuron_count, loihi_compartments,
                               log_spikes=False, log_potential=False,
                               force_update=False, threshold=threshold,
                               reset=reset)
    layer_2 = sim.create_layer(network, layer_neuron_count, loihi_compartments,
                               log_spikes=False, log_potential=False,
                               force_update=False, threshold=threshold,
                               reset=reset)

    # Create connections
    weight = 1.0
    for src in layer_1.neurons:
        for dest in layer_2.neurons:
            if random.random() < connection_probability:
                # Same weight for all connections
                src.add_connection(dest, weight)

    return network


def connected_layers(weights, spiking=True, mapping="luke"):
    network = sim.Network(save_mappings=True)
    loihi_compartments = sim.init_compartments(32, 4, 1024)

    layer_neuron_count = len(weights)
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neuron_count

    neurons_per_core = [0, 0, 0, 0]
    if mapping == "luke" or mapping == "l2_split" or mapping == "fixed":
        neurons_per_core[0] = layer_neuron_count
    elif (mapping == "split_2" or mapping == "l1_split" or mapping == "split_4"
          or mapping == "split_4_diff_tiles"):
        neurons_per_core[0] = ((layer_neuron_count + 1) // 2)
        neurons_per_core[1] = layer_neuron_count // 2
    else:
        print("Error: mapping not supported")
        exit(1)

    layer_mapping = [(0, 0) for _ in range(0, neurons_per_core[0])]
    for _ in range(0, neurons_per_core[1]):
        layer_mapping.append((0, 1))
    for _ in range(0, neurons_per_core[2]):
        layer_mapping.append((0, 2))
    for _ in range(0, neurons_per_core[3]):
        layer_mapping.append((0, 3))

    layer_1 = sim.create_layer(network, layer_neuron_count, loihi_compartments,
                               log_spikes=False, log_potential=False,
                               force_update=False, threshold=threshold,
                               reset=0.0, leak=1.0, mappings=layer_mapping)

    neurons_per_core = [0, 0, 0, 0, 0]
    if mapping == "luke":
        neurons_per_core[0] = min(layer_neuron_count, 1024 - layer_neuron_count)
        neurons_per_core[1] = layer_neuron_count - neurons_per_core[0]
    elif mapping == "fixed":
        neurons_per_core[1] = layer_neuron_count
    elif mapping == "split_2":
        neurons_per_core[0] = ((layer_neuron_count + 1) // 2)
        neurons_per_core[1] = layer_neuron_count // 2
    elif mapping == "l2_split":
        neurons_per_core[1] = ((layer_neuron_count + 1) // 2)
        neurons_per_core[2] = layer_neuron_count // 2
    elif mapping == "split_4":
        neurons_per_core[2] = ((layer_neuron_count + 1) // 2)
        neurons_per_core[3] = layer_neuron_count // 2
    elif mapping == "l1_split":
        neurons_per_core[2] = layer_neuron_count
    elif mapping == "split_4_diff_tiles":
        neurons_per_core[2] = ((layer_neuron_count + 1) // 2)
        neurons_per_core[3] = 0
        neurons_per_core[4] = layer_neuron_count // 2

    layer_mapping = [(0, 0) for _ in range(0, neurons_per_core[0])]
    for _ in range(0, neurons_per_core[1]):
        layer_mapping.append((0, 1))
    for _ in range(0, neurons_per_core[2]):
        layer_mapping.append((0, 2))
    for _ in range(0, neurons_per_core[3]):
        layer_mapping.append((0, 3))
    for _ in range(0, neurons_per_core[4]):
        layer_mapping.append((1, 0))

    layer_2 = sim.create_layer(network, layer_neuron_count,
                               loihi_compartments, log_spikes=False,
                               log_potential=False, force_update=False,
                               threshold=threshold, reset=0.0, leak=1.0,
                               mappings=layer_mapping)

    for src in layer_1.neurons:
        # Add bias to force neuron to fire
        src.add_bias(1.0)
        for dest in layer_2.neurons:
            # Take the ID of the neuron in the 2nd layer
            weight = float(weights[src.id][dest.id]) / 256
            if abs(weight) >= (1.0 / 256):
                # Zero weights are pruned i.e. removed
                src.add_connection(dest, weight)

    return network


def run_spiking_experiment(mapping, cores_blocking, tiles_blocking,
                           max_size=30):
    with open("runs/sandia_data/weights_loihi.pkl", "rb") as weights_file:
        weights = pickle.load(weights_file)

    # Setup the correct blocking and save to a temporary arch file
    with open(os.path.join(PROJECT_DIR, "arch", "loihi.yaml"), "r") as arch_file:
        loihi_arch = yaml.safe_load(arch_file)

    tiles = loihi_arch["architecture"]["tile"]
    for t in tiles:
        t["attributes"]["blocking"] = tiles_blocking
        cores = t["core"]
        for c in cores:
            c["attributes"]["blocking"] = cores_blocking

    generated_arch_filename = os.path.join(PROJECT_DIR, "runs", "calibration",
                                 "calibrated_loihi.arch")
    with open(generated_arch_filename, "w") as arch_file:
        yaml.safe_dump(loihi_arch, arch_file)

    timesteps = 1
    for i in range(1, max_size):
        # Sweep across range of network sizes
        layer_neurons = i*i
        snn = connected_layers(weights[i-1].transpose(), spiking=True,
                                    mapping=mapping)
        snn.save(NETWORK_FILENAME)

        print("Testing network with {0} neurons".format(2*layer_neurons))
        results = sim.run(generated_arch_filename, NETWORK_FILENAME,
                            timesteps)

        with open(os.path.join(PROJECT_DIR, "runs",
                               "calibration", "sim_spiking.csv"),
                  "a") as spiking_csv:
            spiking_writer = csv.DictWriter(spiking_csv,
                                            ("neuron_counts", "energy", "time",
                                             "mapping", "cores_blocking",
                                             "tiles_blocking"))
            spiking_writer.writerow({"neuron_counts": layer_neurons*2,
                                      "time": results["time"],
                                      "energy": results["energy"],
                                      "mapping": mapping,
                                      "cores_blocking": cores_blocking,
                                      "tiles_blocking": tiles_blocking})

    return


mappings = ("fixed", "luke", "split_2")
if __name__ == "__main__":
    run_experiments = False
    plot_experiments = True

    times = {0: [], 256: [], 512: [], 768: [], 1024: []}
    energy = {0: [], 256: [], 512: [], 768: [], 1024: []}

    # This experiment looks at two fully connected layers, spiking
    if run_experiments:
        with open(f"runs/calibration/sim_spiking.csv", "w") as spiking_csv:
            spiking_writer = csv.DictWriter(spiking_csv,
                                       ("neuron_counts", "energy", "time",
                                        "mapping", "cores_blocking",
                                        "tiles_blocking"))
            spiking_writer.writeheader()
        for mapping in mappings:
            run_spiking_experiment(mapping, cores_blocking=False,
                                   tiles_blocking=False, max_size=30)
            run_spiking_experiment(mapping, cores_blocking=True,
                                   tiles_blocking=False, max_size=30)
            run_spiking_experiment(mapping, cores_blocking=True,
                                   tiles_blocking=True, max_size=30)
        with open("runs/sandia_data/weights_loihi.pkl", "rb") as weights_file:
            weights = pickle.load(weights_file)

        neuron_counts = []
        spiking_times = []
        spiking_energy = []

    # **************************************************************************
    # Read Loihi measurement data from csv, this is only available to me locally
    #  since this is restricted data!
    if plot_experiments:
        loihi_times_spikes = {"fixed": [], "luke": [], "split_2": []}
        loihi_energy_spikes = []

        spiking_energy = []
        with open("runs/calibration/sim_spiking.csv", "r") as spiking_csv:
            df = pd.read_csv(spiking_csv)

        with open("runs/sandia_data/loihi_spiking.csv", "r") as spiking_csv:
            spiking_reader = csv.DictReader(spiking_csv)
            for row in spiking_reader:
                for mapping in mappings:
                    loihi_times_spikes[mapping].append(float(row[mapping]))
                loihi_energy_spikes.append(float(row["dynamic energy"]))

        plt.rcParams.update({'font.size': 7, 'lines.markersize': 4})
        # First plot results for the simple fixed mapping, where one layer is on
        #  one core and the second layer is on another
        spiking_frame = df.loc[(df["mapping"] == "luke") &
                               (df["cores_blocking"] == False) &
                               (df["tiles_blocking"] == False)]
        spiking_times = np.array(spiking_frame["time"])
        spiking_energy = np.array(spiking_frame["energy"])
        neuron_counts = np.array(spiking_frame["neuron_counts"])

        plt.figure(figsize=(2.2, 2.2))
        plt.plot(neuron_counts[:-7],
                 np.array(loihi_times_spikes["luke"][:-7]) * 1.0e3, "-")
        plt.plot(neuron_counts[:-7], spiking_times[:-7] * 1.0e3, "x")
        ax = plt.gca()
        ax.set_box_aspect(1)
        print(ax.get_position())
        #plt.gca().set_aspect("equal", adjustable="box")
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Time-step Latency (ms)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.xticks(np.arange(0, neuron_counts[-7]+1, 500))
        plt.legend(("Measured", "Simulated"), fontsize=7)
        plt.tight_layout(pad=0.3)
        ax_pos = ax.get_position()
        plt.savefig("runs/calibration/calibration_time_core.pdf")
        plt.savefig("runs/calibration/calibration_time_core.png")

        plt.figure(figsize=(2.2, 2.2))
        plt.plot(neuron_counts[:-7], np.array(loihi_energy_spikes[:-7]) * 1.0e6, "-")
        plt.plot(neuron_counts[:-7], np.array(spiking_energy[:-7]) * 1.0e6, "x")
        ax = plt.gca()
        ax.set_box_aspect(1)
        #plt.gca().set_aspect("equal", adjustable="box")
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Energy ($\mu$J)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.xticks(np.arange(0, neuron_counts[-7]+1, 500))
        plt.legend(("Measured", "Simulated"), fontsize=7)
        plt.tight_layout(pad=0.3)
        ax.set_position(ax_pos)
        plt.savefig("runs/calibration/calibration_energy.pdf")
        plt.savefig("runs/calibration/calibration_energy.png")

        plt.rcParams.update({'font.size': 7, 'lines.markersize': 4})
        ## Plot the effect of cores blocking
        spiking_frame = df.loc[(df["mapping"] == "luke") &
                               (df["tiles_blocking"] == False)]
        cores_nonblocking = spiking_frame.loc[spiking_frame["cores_blocking"] == False]
        cores_blocking = spiking_frame.loc[spiking_frame["cores_blocking"] == True]

        plt.figure(figsize=(2.1, 2.1))
        plt.plot(neuron_counts, np.array(loihi_times_spikes["luke"]) * 1.0e3, "-")
        plt.plot(neuron_counts, np.array(cores_nonblocking["time"] * 1.0e3), "x")
        plt.plot(neuron_counts, np.array(cores_blocking["time"]) * 1.0e3, "ko",
                 fillstyle="none")

        #plt.figure(figsize=(2.5, 2.5))
        #plt.plot(neuron_counts, np.array(cores_blocking["time"]) * 1.0e3, "-o")
        #plt.plot(neuron_counts, np.array(loihi_times_spikes["luke"]) * 1.0e3, "-x")
        #plt.legend(("Simulated", "Measured on Loihi"))

        plt.gca().set_box_aspect(1)
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Time-step Latency (ms)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.legend(("Measured", "No blocking", "Cores blocking"), fontsize=7)
        plt.tight_layout(pad=0.3)
        plt.savefig("runs/calibration/calibration_time_partition_2.pdf")
        plt.savefig("runs/calibration/calibration_time_partition_2.png")

        # Plot the effect of network tiles blocking
        spiking_frame = df.loc[(df["mapping"] == "split_2") &
                               (df["cores_blocking"] == True)]
        tiles_nonblocking = spiking_frame.loc[spiking_frame["tiles_blocking"] == False]
        tiles_blocking = spiking_frame.loc[spiking_frame["tiles_blocking"] == True]

        plt.figure(figsize=(2.1, 2.1))
        plt.plot(neuron_counts, np.array(loihi_times_spikes["split_2"]) * 1.0e3, "-")
        plt.plot(neuron_counts, np.array(tiles_nonblocking["time"]) * 1.0e3, "x")
        plt.plot(neuron_counts, np.array(tiles_blocking["time"]) * 1.0e3, "ko",
                 fillstyle="none"),
        plt.gca().set_box_aspect(1)
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Time-step Latency (ms)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.legend(("Measured", "Cores blocking", "Tiles blocking"),
                    fontsize=7)
        plt.tight_layout(pad=0.3)
        plt.savefig("runs/calibration/calibration_time_partition_3.pdf")
        plt.savefig("runs/calibration/calibration_time_partition_3.png")
