"""First attempt at simulation run script.

Run a few basic experiments, show how we might interface a python
script with the simulator kernel.
"""

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt
import pickle

MAX_COMPARTMENTS = 1024
MAX_CORES = 128
NETWORK_FILENAME = "connected_layer.csv"
TECH_FILENAME = "loihi.tech"
ARCH_FILENAME = "loihi.list"

def run_sim(network):
    fields = ["Neuron ID", "Compartment ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
    with open(NETWORK_FILENAME, "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(fields)
        writer.writerows(network)

    timesteps = 2
    command = ("./sim", TECH_FILENAME, ARCH_FILENAME,
               NETWORK_FILENAME, "{0}".format(timesteps))
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    with open("results.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)

    return results


def connected_layer(weights, spiking=True):
    network = []

    layer_neurons = len(weights)
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neurons

    reset = 0
    force_update = True
    for n in range(0, layer_neurons):
        neuron = [n, n, threshold, reset, 0, 0, int(force_update)]
        for dest in range(0, layer_neurons):
            # Take the ID of the neuron in the 2nd layer
            weight = float(weights[n][dest]) / 255
            if weight != 0:
                # Zero weights are pruned i.e. removed
                neuron.extend((dest+layer_neurons, weight))
        network.append(neuron)

    for n in range(layer_neurons, 2*layer_neurons):
        neuron = [n, n, threshold, reset, 0, 0, int(force_update)]
        network.append(neuron)

    return network

import random
def fully_connected(layer_neurons, spiking=True, probability=1.0):
    # Two layers, fully connected
    network = []
    if spiking:  # always spike
        threshold = -1.0
    else:  # never spike
        threshold = 2*layer_neurons

    weight = 1.0
    reset = 0
    force_update = True
    for n in range(0, layer_neurons):
        neuron = [n, n, threshold, reset, 0, 0, int(force_update)]
        for dest in range(layer_neurons, 2*layer_neurons):
            if random.random() < probability:
                neuron.extend((dest, weight))  # Same weight for all connections
        network.append(neuron)

    for n in range(layer_neurons, 2*layer_neurons):
        neuron = [n, n, threshold, reset, 0, 0, int(force_update)]
        network.append(neuron)

    return network


def empty(neurons, max_compartments=MAX_COMPARTMENTS):
    network = []
    threshold = -1.0  # never spike
    reset = 0.0

    compartment = 0
    # Map a number of neurons onto cores, but we may not necessarily use all
    #  compartments of each core
    # TODO: how to do 0 compartments yet still activate the core?
    # Maybe the simulator can take number of cores to simulate as an arg
    #  then we assume all simulated cores are powered
    for n in range(0, neurons):
        neuron = [n, n, threshold, reset, 0, 0, 0]
        network.append(neuron)
        compartment += 1

    return network


if __name__ == "__main__":
    #core_count = [1, 2, 4, 8, 16, 32, 64, 128]
    times = {0: [], 256: [], 512: [], 768: [], 1024: []}
    energy = {0: [], 256: [], 512: [], 768: [], 1024: []}
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

    with open("sandia_data/weights.pkl", "rb") as weights_file:
        weights = pickle.load(weights_file)

    neurons = []
    spiking_times = []
    spiking_energy = []

    for i in range(1, 31):
        layer_neurons = i*i

        #network = fully_connected(layer_neurons, spiking=True, probability=connection_probabilities[i-1])
        network = connected_layer(weights[i-1], spiking=True)
        print("Testing network with {0} neurons".format(len(network)))
        results = run_sim(network)

        neurons.append(len(network))
        spiking_times.append(results["time"])
        spiking_energy.append(results["energy"])

    # Write all the simulation data to csv
    with open("sim_spiking.csv", "w") as spiking_csv:
        spiking_writer = csv.DictWriter(spiking_csv,
                                        ("neurons", "energy", "time"))
        spiking_writer.writeheader()
        for neuron_count, time, energy_val in zip(neurons, spiking_times,
                                                  spiking_energy):
            spiking_writer.writerow({"neurons": neuron_count,
                                     "energy": energy_val,
                                     "time": time})

    neurons = []
    nonspiking_times = []
    nonspiking_energy = []

    # The second experiment looks at two fully connected layers, not spiking
    for i in range(1, 31):
        layer_neurons = i*i

        network = fully_connected(layer_neurons, spiking=False)
        print("Testing network with {0} neurons".format(len(network)))
        results = run_sim(network)

        neurons.append(len(network))
        nonspiking_times.append(results["time"])
        nonspiking_energy.append(results["energy"])

    with open("sim_nonspiking.csv", "w") as nonspiking_csv:
        nonspiking_writer = csv.DictWriter(nonspiking_csv,
                                           ("neurons", "energy", "time"))
        nonspiking_writer.writeheader()
        for neuron_count, time, energy_val in zip(neurons, nonspiking_times,
                                                  nonspiking_energy):
            nonspiking_writer.writerow({"neurons": neuron_count,
                                        "energy": energy_val,
                                        "time": time})
    # Read Loihi measurement data from csv, this is only available to me locally
    #  since this is restricted data!
    neurons = []
    loihi_times_spikes = []
    loihi_energy_spikes = []

    spiking_energy = []
    spiking_times = []
    with open("sim_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            spiking_times.append(float(row["time"]))
            spiking_energy.append(float(row["energy"]))
            neurons.append(int(row["neurons"]))

    nonspiking_times = []
    nonspiking_energy = []
    with open("sim_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            nonspiking_times.append(float(row["time"]))
            nonspiking_energy.append(float(row["energy"]))

    with open("loihi_spiking.csv", "r") as spiking_csv:
        spiking_reader = csv.DictReader(spiking_csv)
        for row in spiking_reader:
            loihi_times_spikes.append(float(row["time"]))
            loihi_energy_spikes.append(float(row["energy"]))

    loihi_times_no_spikes = []
    loihi_energy_no_spikes = []
    with open("loihi_nonspiking.csv", "r") as nonspiking_csv:
        nonspiking_reader = csv.DictReader(nonspiking_csv)
        for row in nonspiking_reader:
            loihi_times_no_spikes.append(float(row["time"]))
            loihi_energy_no_spikes.append(float(row["energy"]))

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, spiking_times, "-o")
    plt.plot(neurons, loihi_times_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, spiking_energy, "-o")
    plt.plot(neurons, loihi_energy_spikes, "--x", color="orange")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_energy.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, nonspiking_energy, "-o")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_energy_sim_only.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, nonspiking_times, "-o")
    plt.plot(neurons, loihi_times_no_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_not_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, nonspiking_energy, "-o")
    plt.plot(neurons, loihi_energy_no_spikes, "--x")
    plt.yscale("linear")
    plt.xscale("linear")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_not_spiking_energy.png")

    # Some additional plots to highlight trends
    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, loihi_times_spikes, "--x", color="orange")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Measured",))
    #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_time_loihi_only.png")

    plt.show()
