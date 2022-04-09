"""First attempt at simulation run script.

Run a few basic experiments, show how we might interface a python
script with the simulator kernel.
"""

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt

MAX_COMPARTMENTS = 1024
MAX_CORES = 128
NETWORK_FILENAME = "connected_layer.csv"


def run_sim(network, core_count):
    fields = ["Neuron ID", "Core ID", "Threshold", "Reset", "Is Input",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
    with open(NETWORK_FILENAME, "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(fields)
        writer.writerows(network)

    timesteps = 1
    command = ("./sim", "{0}".format(timesteps), "{0}".format(core_count),
               NETWORK_FILENAME)
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    with open("results.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)

    return results


def fully_connected(layer_neurons, spiking=True):
    # Two layers, fully connected
    network = []
    if spiking:
        threshold = -1.0
        is_input = 1
    else:  # never spike
        threshold = 2*layer_neurons
        is_input = 0

    weight = 1.0
    reset = 0
    for n in range(0, layer_neurons):
        core_id = n / MAX_COMPARTMENTS

        neuron = [n, core_id, threshold, reset, is_input, 0, 0]
        for dest in range(layer_neurons, 2*layer_neurons):
            neuron.extend((dest, weight))  # Same weight for all connections
        network.append(neuron)

    for n in range(layer_neurons, 2*layer_neurons):
        core_id = n / MAX_COMPARTMENTS
        neuron = [n, core_id, threshold, reset, is_input, 0, 0]
        network.append(neuron)

    core_count = core_id + 1

    return network, core_count


def empty(neurons, max_compartments=MAX_COMPARTMENTS):
    network = []
    threshold = -1.0  # never spike
    reset = 0.0
    is_input = 1

    core_id = 0
    compartment = 0
    # Map a number of neurons onto cores, but we may not necessarily use all
    #  compartments of each core
    # TODO: how to do 0 compartments yet still activate the core?
    # Maybe the simulator can take number of cores to simulate as an arg
    #  then we assume all simulated cores are powered
    for n in range(0, neurons):
        if compartment == max_compartments:
            core_id += 1
            compartment = 0

        neuron = [n, core_id, threshold, reset, is_input, 0, 0]
        network.append(neuron)
        compartment += 1

    core_count = core_id + 1

    return network, core_count

if __name__ == "__main__":
    core_count = [1, 2, 4, 8, 16, 32, 64, 128]
    times = {0: [], 256: [], 512: [], 768: [], 1024: []}
    energy = {0: [], 256: [], 512: [], 768: [], 1024: []
    """
    for cores in core_count:
        for compartments in range(0, MAX_COMPARTMENTS+1, 256):
            n = compartments * cores
            network, _ = empty(n, compartments)
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
    neurons = []
    times = []
    energy = []

    for i in range(1, 31):
        layer_neurons = i*i

        network, core_count = fully_connected(layer_neurons)
        print("Testing network with {0} neurons".format(len(network)))
        results = run_sim(network, core_count)

        neurons.append(len(network))
        times.append(results["time"])
        energy.append(results["energy"])

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, times, "-o")
    plt.plot(neurons, loihi_times_spikes, "--x")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, energy, "-o")
    plt.plot(neurons, loihi_energy_spikes, "--x", color="orange")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_energy.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, energy, "-o")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated",))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_energy_sim_only.png")

    neurons = []
    times = []
    energy = []

    # The second experiment looks at two fully connected layers, not spiking
    for i in range(1, 31):
        layer_neurons = i*i

        network, core_count = fully_connected(layer_neurons, spiking=False)
        print("Testing network with {0} neurons".format(len(network)))
        results = run_sim(network, core_count)

        neurons.append(len(network))
        times.append(results["time"])
        energy.append(results["energy"])

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, times, "-o")
    plt.plot(neurons, loihi_times_no_spikes, "--x")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_not_spiking_time.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, energy, "-o")
    plt.plot(neurons, loihi_energy_no_spikes, "--x")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated", "Measured"))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_not_spiking_energy.png")

    # Some additional plots to highlight trends
    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, energy, "-o")
    plt.ylabel("Energy (J)")
    plt.xlabel("Neurons")
    plt.legend(("Simulated",))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_not_spiking_energy_sim_only.png")

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(neurons, loihi_times_spikes, "--x", color="orange")
    plt.ylabel("Time (s)")
    plt.xlabel("Neurons")
    plt.legend(("Measured",))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("connected_spiking_time_loihi_only.png")

    plt.show()
