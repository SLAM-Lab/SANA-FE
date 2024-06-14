"""
Copyright (c) 2024 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Biologically-inspired power benchmark for neuromorphic hardware

Implements several partitions of a small power benchmark taken from
[Parker, 2021].
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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))

sys.path.insert(0, PROJECT_DIR)
import utils as sf
import sanafecpp as kernel

MAX_TILES = 32
MAX_CORES = 4
MAX_COMPARTMENTS = 1024
ARCH_FILENAME = "arch/loihi.yaml"

def connected_layers(arch, weights, spiking=True, mapping="l2_split",
                     copy_network=False):

    net = kernel.Network()
    layer_neuron_count = len(weights)
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neuron_count

    duplicate = (2 if copy_network else 1)
    for i in range(0, duplicate):
        neurons_per_core = [0, 0, 0, 0]
        if (mapping == "luke" or mapping == "l2_split" or mapping == "fixed" or
            mapping == "split_2_diff_tiles"):
            neurons_per_core[0] = layer_neuron_count
        elif (mapping == "split_2" or mapping == "l1_split" or mapping == "split_4"
            or mapping == "split_4_diff_tiles"):
            neurons_per_core[0] = ((layer_neuron_count + 1) // 2)
            neurons_per_core[1] = layer_neuron_count // 2
        else:
            print("Error: mapping not supported")
            exit(1)

        layer1_mapping = [(0, i) for _ in range(0, neurons_per_core[0])]
        for _ in range(0, neurons_per_core[1]):
            layer1_mapping.append((0, 1))
        for _ in range(0, neurons_per_core[2]):
            layer1_mapping.append((0, 2))
        for _ in range(0, neurons_per_core[3]):
            layer1_mapping.append((0, 3))

        layer1_attributes = {
            "log_spikes": 1,
            "log_potential": 1,
            "force_update": 1,
            "threshold": threshold,
            "reset": 0.0,
            "leak_decay": 0.0,
            "soma_hw_name": "loihi_lif",
            "synapse_hw_name": "loihi_dense_synapse",
        }
        layer_1 = net.create_neuron_group(layer_neuron_count,
                                          layer1_attributes)

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
        elif mapping == "split_2_diff_tiles":
            neurons_per_core[4] = layer_neuron_count

        layer2_mapping = [(0, 0) for _ in range(0, neurons_per_core[0])]
        for _ in range(0, neurons_per_core[1]):
            layer2_mapping.append((0, 1))
        for _ in range(0, neurons_per_core[2]):
            layer2_mapping.append((0, 2))
        for _ in range(0, neurons_per_core[3]):
            layer2_mapping.append((0, 3))
        for _ in range(0, neurons_per_core[4]):
            layer2_mapping.append((1, i))

        layer2_attributes = {
            "log_spikes": 1,
            "log_potential": 1,
            "force_update": 0,
            "threshold": threshold,
            "soma_hw_name": "loihi_lif",
            "synapse_hw_name": "loihi_dense_synapse",
            # Add bias to force neuron to fire
            "bias": 1.0,
        }
        layer_2 = net.create_neuron_group(layer_neuron_count,
                                          layer2_attributes)

        for src in layer_1.neurons:
            for dest in layer_2.neurons:
                # Take the ID of the neuron in the 2nd layer
                weight = float(weights[src.get_id()][dest.get_id()]) / 256
                if abs(weight) >= (1.0 / 256):
                    # Zero weights are pruned i.e. removed
                    src.connect_to_neuron(dest, { "weight": weight } )

        # TODO: we have to map the neurons after defining them and their
        #  connections. Add flags in the kernel to make sure this happens.
        #  Currently when we save the net file this won't be adhered to.. do
        #  we need to fix the save_net code or allow us to add connections
        #  after mapping to cores
        sf.map_neuron_group(layer_1, arch, layer1_mapping)
        sf.map_neuron_group(layer_2, arch, layer2_mapping)

    return net


def run_spiking_experiment(mapping, max_size=30):
    with open(os.path.join(PROJECT_DIR, "runs", "sandia_data",
                           "weights_loihi.pkl"), "rb") as weights_file:
        weights = pickle.load(weights_file)

    timesteps = 1
    for i in range(1, max_size):
        # Sweep across range of network sizes
        layer_neurons = i*i
        copy_network = (True if mapping == "split_2_diff_tiles" else False)
        arch = kernel.load_arch_description(ARCH_FILENAME)

        snn = connected_layers(arch, weights[i-1].transpose(), spiking=True,
                               mapping=mapping, copy_network=copy_network)
        network_filename = os.path.join(PROJECT_DIR, "runs", "power", "snn",
                                        f"connected_layers_N{layer_neurons}_map_{mapping}.net")
        snn.save_net_description(network_filename)

        print("Testing network with {0} neurons".format(2*layer_neurons))
        sim = kernel.Simulation(arch, snn)
        results = sim.run(timesteps)

        with open(os.path.join(PROJECT_DIR, "runs",
                               "power", "sim_spiking.csv"),
                  "a") as spiking_csv:
            spiking_writer = csv.DictWriter(spiking_csv,
                                            ("neuron_counts", "energy", "time",
                                             "mapping"))
            neuron_counts = layer_neurons * 2
            if copy_network:
                neuron_counts *= 2
            spiking_writer.writerow({"neuron_counts": neuron_counts,
                                      "time": results["sim_time"],
                                      "energy": results["energy"],
                                      "mapping": mapping})
            print(f"Run results:{results}")

    return


mappings = ("fixed", "l2_split", "split_2", "luke", "split_4")
#mappings = ("split_2_diff_tiles",)
if __name__ == "__main__":
    run_experiments = True
    plot_experiments = True

    times = {0: [], 256: [], 512: [], 768: [], 1024: []}
    energy = {0: [], 256: [], 512: [], 768: [], 1024: []}

    # This experiment looks at two fully connected layers, spiking
    if run_experiments:
        with open(os.path.join("runs", "power",
                               "sim_spiking.csv"), "w") as spiking_csv:
            spiking_writer = csv.DictWriter(spiking_csv,
                                       ("neuron_counts", "energy", "time",
                                        "mapping"))
            spiking_writer.writeheader()
        for mapping in mappings:
            run_spiking_experiment(mapping, max_size=30)
        with open(os.path.join("runs", "sandia_data",
                               "weights_loihi.pkl"), "rb") as weights_file:
            weights = pickle.load(weights_file)

        neuron_counts = []
        spiking_times = []
        spiking_energy = []

    # **************************************************************************
    # Read Loihi measurement data from csv, this is only available to me locally
    #  since this is restricted data!
    if plot_experiments:
        loihi_times_spikes = {"fixed": [], "l2_split": [], "split_2": [],
                              "luke": [], "split_4": []}
        loihi_energy_spikes = []

        spiking_energy = []
        with open(os.path.join("runs", "power",
                               "sim_spiking.csv"), "r") as spiking_csv:
            df = pd.read_csv(spiking_csv)

        with open(os.path.join("runs", "sandia_data",
                               "loihi_spiking.csv"), "r") as spiking_csv:
            spiking_reader = csv.DictReader(spiking_csv)
            for row in spiking_reader:
                for mapping in mappings:
                    loihi_times_spikes[mapping].append(float(row[mapping]))
                loihi_energy_spikes.append(float(row["dynamic energy"]))

        plt.rcParams.update({'font.size': 6, 'lines.markersize': 3})
        # First plot results for the simple fixed mapping, where one layer is on
        #  one core and the second layer is on another
        spiking_frame = df.loc[(df["mapping"] == "luke")]
        spiking_times = np.array(spiking_frame["time"])
        spiking_energy = np.array(spiking_frame["energy"])
        neuron_counts = np.array(spiking_frame["neuron_counts"])

        plt.figure(figsize=(1.5, 1.5))
        plt.plot(neuron_counts[6:-7],
                 np.array(loihi_times_spikes["luke"][6:-7]) * 1.0e3, "-")
        plt.plot(neuron_counts[6:-7], spiking_times[6:-7] * 1.0e3, "ko",
                 fillstyle="none")
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
        plt.legend(("Measured", "Simulated"), fontsize=6)
        plt.tight_layout(pad=0.3)
        ax_pos = ax.get_position()
        plt.savefig(os.path.join("runs", "power", "power_time_core.pdf"))
        plt.savefig(os.path.join("runs", "power", "power_time_core.png"))

        energy_error = 100 * np.abs(np.array(loihi_energy_spikes[6:]) - spiking_energy[6:]) / np.array(loihi_energy_spikes[6:])
        print(energy_error)
        energy_error = np.mean(np.abs(np.array(loihi_energy_spikes[6:]) - spiking_energy[6:]) / np.array(loihi_energy_spikes[6:]))
        latency_error = np.mean(np.abs(np.array(loihi_times_spikes["luke"][6:]) - spiking_times[6:]) / np.array(loihi_times_spikes["luke"][6:]))
        print(f"Energy error %: {energy_error*100}")
        print(f"Latency error %: {latency_error*100}")

        plt.figure(figsize=(1.6, 1.6))
        plt.plot(neuron_counts[6:], np.array(loihi_energy_spikes[6:]) * 1.0e6, "-")
        plt.plot(neuron_counts[6:], np.array(spiking_energy[6:]) * 1.0e6, "ko",
                 fillstyle="none", mew=0.8)
        ax = plt.gca()
        ax.set_box_aspect(1)
        #plt.gca().set_aspect("equal", adjustable="box")
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Energy ($\mu$J)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.xticks(np.arange(0, neuron_counts[-1]+1, 500))
        plt.legend(("Measured", "Simulated"), fontsize=6)
        plt.tight_layout(pad=0.1)
        #ax.set_position(ax_pos)
        plt.savefig(os.path.join("runs", "power", "power_energy.pdf"))
        plt.savefig(os.path.join("runs", "power", "power_energy.png"))

        plt.rcParams.update({'font.size': 6, 'lines.markersize': 3,})
        ## Plot the effect of cores blocking
        spiking_frame = df.loc[(df["mapping"] == "l2_split")]

        plt.figure(figsize=(1.6, 1.6))
        plt.plot(neuron_counts[6:], np.array(loihi_times_spikes["l2_split"][6:]) * 1.0e3, "-")
        plt.plot(neuron_counts[6:], np.array(spiking_frame["time"][6:]) * 1.0e3, "ko",
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
        plt.legend(("Measured", "Simulated"), fontsize=6)
        plt.tight_layout(pad=0.3)
        plt.savefig(os.path.join("runs", "power", "power_time_partition_2.pdf"))
        plt.savefig(os.path.join("runs", "power", "power_time_partition_2.png"))

        # Plot the effect of network tiles blocking
        spiking_frame = df.loc[(df["mapping"] == "split_4")]

        plt.figure(figsize=(1.6, 1.6))
        plt.plot(neuron_counts, np.array(loihi_times_spikes["split_4"]) * 1.0e3, "-")
        plt.plot(neuron_counts, np.array(spiking_frame["time"]) * 1.0e3, "ko",
                 fillstyle="none")
        plt.gca().set_box_aspect(1)
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Time-step Latency (ms)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.legend(("Measured", "Simulated"),
                    fontsize=6)
        plt.tight_layout(pad=0.3)
        plt.savefig(os.path.join("runs", "power", "power_time_partition_3.pdf"))
        plt.savefig(os.path.join("runs", "power", "power_time_partition_3.png"))

        spiking_frame = df.loc[(df["mapping"] == "luke")]
        plt.figure(figsize=(1.5, 1.5))
        plt.plot(neuron_counts, np.array(loihi_times_spikes["luke"]) * 1.0e3, "-")
        plt.plot(neuron_counts, np.array(spiking_frame["time"]) * 1.0e3, "ko",
                 fillstyle="none")
        plt.gca().set_box_aspect(1)
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Time-step Latency (ms)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.legend(("Measured", "Simulated"), fontsize=6)
        plt.tight_layout(pad=0.3)
        plt.savefig(os.path.join("runs", "power", "power_time_partition_luke.pdf"))
        plt.savefig(os.path.join("runs", "power", "power_time_partition_luke.png"))

        # Combine onto one plot
        plt.figure(figsize=(1.6, 1.6))
        plt.plot(neuron_counts[6:], np.array(loihi_times_spikes["luke"][6:]) * 1.0e3, "-")
        plt.plot(neuron_counts[6:], np.array(loihi_times_spikes["l2_split"][6:]) * 1.0e3, "--")
        plt.plot(neuron_counts[6:], np.array(loihi_times_spikes["split_4"][6:]) * 1.0e3, "-.")

        spiking_frame = df.loc[(df["mapping"] == "luke")]
        plt.plot(neuron_counts[6:], np.array(spiking_frame["time"][6:]) * 1.0e3, "o",
                 fillstyle="none", mew=0.8, color="#1f77b4")
        spiking_frame = df.loc[(df["mapping"] == "l2_split")]
        plt.plot(neuron_counts[10:], np.array(spiking_frame["time"][10:]) * 1.0e3, "s",
                 fillstyle="none", mew=0.8, color="#ff7f0e")
        spiking_frame = df.loc[(df["mapping"] == "split_4")]
        plt.plot(neuron_counts[10:], np.array(spiking_frame["time"][10:]) * 1.0e3, "^",
                 fillstyle="none", mew=0.8, color="#2ca02c")
        plt.gca().set_box_aspect(1)
        plt.yscale("linear")
        plt.xscale("linear")
        plt.ylabel("Time-step Latency (ms)")
        plt.xlabel("Neurons")
        plt.minorticks_on()
        plt.legend(("Default", "3 Cores", "4 Cores"),
                    fontsize=6)
        plt.tight_layout(pad=0.1)
        plt.savefig(os.path.join("runs", "power", "power_time.pdf"))
        plt.savefig(os.path.join("runs", "power", "power_time.png"))
