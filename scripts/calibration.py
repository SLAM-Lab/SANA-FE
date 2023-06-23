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
    network = sim.Network()
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
    layer_1 = sim.create_layer(network, layer_neuron_count,
                                     loihi_compartments,
                                     log_spikes, log_potential, force_update,
                                     threshold, reset)
    layer_2 = sim.create_layer(network, layer_neuron_count, loihi_compartments,
                                     log_spikes, log_potential, force_update,
                                     threshold, reset)

    # Create connections
    weight = 1.0
    for src in layer_1.neurons:
        for dest in layer_2.neurons:
            if random.random() < connection_probability:
                src.add_connection(dest, weight)  # Same weight for all connections

    return network


def connected_layers(weights, spiking=True, mapping="luke"):
    network = sim.Network()
    loihi_compartments = sim.init_compartments(32, 4, 1024)

    layer_neuron_count = len(weights)
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neuron_count

    reset = 0
    force_update = True
    log_spikes = False
    log_potential = False
    leak = 1.0

    neurons_per_core = [0, 0, 0, 0]
    if mapping == "luke" or mapping == "l2_split" or mapping == "fixed":
        neurons_per_core[0] = layer_neurons
    elif (mapping == "split_2" or mapping == "l1_split" or mapping == "split_4"
          or mapping == "split_4_diff_tiles"):
        neurons_per_core[0] = ((layer_neurons + 1) // 2)
        neurons_per_core[1] = layer_neurons // 2
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

    layer_1 = sim.create_layer(network, layer_neuron_count,
                                     loihi_compartments, log_spikes,
                                     log_potential, force_update, threshold,
                                     reset, leak, mappings=layer_mapping)

    force_update = False
    neurons_per_core = [0, 0, 0, 0, 0]
    if mapping == "luke":
        neurons_per_core[0] = min(layer_neurons, 1024 - layer_neurons)
        neurons_per_core[1] = layer_neurons - neurons_per_core[0]
    elif mapping == "fixed":
        neurons_per_core[1] = layer_neurons
    elif mapping == "split_2":
        neurons_per_core[0] = ((layer_neurons + 1) // 2)
        neurons_per_core[1] = layer_neurons // 2
    elif mapping == "l2_split":
        neurons_per_core[1] = ((layer_neurons + 1) // 2)
        neurons_per_core[2] = layer_neurons // 2
    elif mapping == "split_4":
        neurons_per_core[2] = ((layer_neurons + 1) // 2)
        neurons_per_core[3] = layer_neurons // 2
    elif mapping == "l1_split":
        neurons_per_core[2] = layer_neurons
    elif mapping == "split_4_diff_tiles":
        neurons_per_core[2] = ((layer_neurons + 1) // 2)
        neurons_per_core[3] = 0
        neurons_per_core[4] = layer_neurons // 2

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
                                     loihi_compartments, log_spikes,
                                     log_potential, force_update, threshold,
                                     reset, leak, mappings=layer_mapping)

    for src in layer_1.neurons:
        for dest in layer_2.neurons:
            # Take the ID of the neuron in the 2nd layer
            weight = float(weights[src.id][dest.id]) / 256
            if abs(weight) >= (1.0 / 256):
                # Zero weights are pruned i.e. removed
                src.add_connection(dest, weight)

    return network


if __name__ == "__main__":
    #core_count = [1, 2, 4, 8, 16, 32, 64, 128]
    times = {0: [], 256: [], 512: [], 768: [], 1024: []}
    energy = {0: [], 256: [], 512: [], 768: [], 1024: []}
    #mapping = "split_4_diff_tiles"

    mapping = "fixed"
    #mapping = "luke"
    #mapping = "split_2"
    """
    for cores in core_count:
        for compartments in range(0, MAX_COMPARTMENTS+1, 256):
            n = compartments * cores
            network = empty(n, compartments)
            results = run_sim(network, cores)

            times[compartments].append(results["time"])
            energy[compartments].append(results["energy"])

    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(5.5, 5.5))
    for compartments in range(0, MAX_COMPARTMENTS+1, 256):
        plt.plot(core_count, times[compartments], "-o")
    plt.ylabel("Time (s)")
    plt.xlabel("Cores Used")
    plt.legend(("0", "256", "512", "768", "1024"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("empty_times.png")

    plt.figure(figsize=(5.5, 5.5))
    for compartments in range(0, MAX_COMPARTMENTS+1, 256):
        plt.plot(core_count, energy[compartments], "-o")
    plt.ylabel("Energy (J)")
    plt.xlabel("Cores Used")
    plt.legend(("0", "256", "512", "768", "1024"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("empty_energy.png")
    """
    # This experiment looks at two fully connected layers, spiking

    with open("runs/sandia_data/weights_loihi.pkl", "rb") as weights_file:
        weights = pickle.load(weights_file)

    neuron_counts = []
    spiking_times = []
    spiking_energy = []

    timesteps = 1
    for i in range(1, 30):
    #for i in range(1, 4):
        layer_neurons = i*i

        #network = fully_connected(layer_neurons, spiking=True, probability=connection_probabilities[i-1])
        snn = connected_layers(weights[i-1].transpose(), spiking=True,
                                    mapping=mapping)
        print("Testing network with {0} neurons".format(2*layer_neurons))
        snn.save(NETWORK_FILENAME)
        results = sim.run(ARCH_FILENAME, NETWORK_FILENAME, timesteps)

        neuron_counts.append(layer_neurons*2)
        spiking_times.append(results["time"])
        spiking_energy.append(results["energy"])

    # Write all the simulation data to csv
    with open("runs/calibration/sim_spiking.csv", "w") as spiking_csv:
        spiking_writer = csv.DictWriter(spiking_csv,
                                        ("neuron_counts", "energy", "time"))
        spiking_writer.writeheader()
        for count, time, energy_val in zip(neuron_counts, spiking_times,
                                                  spiking_energy):
            spiking_writer.writerow({"neuron_counts": count,
                                     "energy": energy_val,
                                     "time": time})

    neuron_counts = []
    nonspiking_times = []
    nonspiking_energy = []

    """
    # The second experiment looks at two fully connected layers, not spiking
    for i in range(1, 31):
        layer_neurons = i*i

        commands = fully_connected(layer_neurons, spiking=False,
                                   force_update=True)
        print("Testing network with {0} neurons".format(layer_neurons*2))
        results = run_sim(commands, timesteps)

        neuron_counts.append(layer_neurons*2)
        nonspiking_times.append(results["time"] / timesteps)
        nonspiking_energy.append(results["energy"] / timesteps)

    with open("runs/sim_nonspiking.csv", "w") as nonspiking_csv:
        nonspiking_writer = csv.DictWriter(nonspiking_csv,
                                           ("neurons", "energy", "time"))
        nonspiking_writer.writeheader()
        for neuron_count, time, energy_val in zip(neuron_counts, nonspiking_times,
                                                  nonspiking_energy):
            nonspiking_writer.writerow({"neurons": neuron_count,
                                        "energy": energy_val,
                                        "time": time})
    """
    # **************************************************************************
    # Read Loihi measurement data from csv, this is only available to me locally
    #  since this is restricted data!
    neuron_counts = []
    loihi_times_spikes = []
    loihi_energy_spikes = []

    spiking_energy = []
    spiking_times = []
    with open("runs/calibration/sim_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            spiking_times.append(float(row["time"]))
            spiking_energy.append(float(row["energy"]))
            neuron_counts.append(int(row["neuron_counts"]))

    nonspiking_times = []
    nonspiking_energy = []
    with open("runs/calibration/sim_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            nonspiking_times.append(float(row["time"]))
            nonspiking_energy.append(float(row["energy"]))

    with open("runs/sandia_data/loihi_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            loihi_times_spikes.append(float(row[mapping]))
            loihi_energy_spikes.append(float(row["dynamic energy"]))

    loihi_times_no_spikes = []
    loihi_energy_no_spikes = []
    with open("runs/sandia_data/loihi_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            loihi_times_no_spikes.append(float(row["time"]))
            loihi_energy_no_spikes.append(float(row["energy"]))

    plt.rcParams.update({'font.size': 8, 'lines.markersize': 3})

    plt.figure(figsize=(2.5, 2.5))
    plt.plot(neuron_counts, spiking_times, "-o")
    plt.plot(neuron_counts, loihi_times_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time-step Latency (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured on Loihi"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.tight_layout()
    plt.savefig("runs/calibration/connected_spiking_time.pdf")
    plt.savefig("runs/calibration/connected_spiking_time.png")

    plt.figure(figsize=(2.5, 2.5))
    plt.plot(neuron_counts, spiking_energy, "-o")
    plt.plot(neuron_counts, loihi_energy_spikes, "--x", color="orange")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured on Loihi"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.tight_layout()
    plt.savefig("runs/calibration/connected_spiking_energy.pdf")
    plt.savefig("runs/calibration/connected_spiking_energy.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, nonspiking_energy[0:-1], "-o")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/calibration/connected_spiking_energy_sim_only.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, nonspiking_times[0:-1], "-o")
    plt.plot(neuron_counts, loihi_times_no_spikes[0:-1], "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured on Loihi"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/calibration/connected_not_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, nonspiking_energy[0:-1], "-o")
    plt.plot(neuron_counts, loihi_energy_no_spikes[0:-1], "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured on Loihi"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/calibration/connected_not_spiking_energy.png")

    """
    # Some additional plots to highlight trends
    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neuron_counts, loihi_times_spikes, "--x", color="orange")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Measured on Loihi",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("runs/connected_spiking_time_loihi_only.png")
    """
    #plt.show()
