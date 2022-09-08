"""First attempt at simulation run script.

Run a few basic experiments, show how we might interface a python
script with the simulator kernel.
"""
#import matplotlib
#matplotlib.use('Agg')

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt
import pandas as pd

ARCH_FILENAME = "loihi.arch"
NETWORK_FILENAME = "examples/dvs_gesture.net"

def run_sim():
    # Create the network to run
    fields = ["Neuron ID", "Core ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
    timesteps = 128
    command = ("./sim", ARCH_FILENAME,
               NETWORK_FILENAME, "{0}".format(timesteps),)
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    # Parse the detailed perf statistics
    print("Reading performance data")
    stats = pd.read_csv("perf.csv")
    analysis = parse_stats(stats)

    return analysis


def parse_stats(stats):
    print("Parsing statistics")
    total = stats.sum()
    print(total)
    update_keys, synapse_keys, spike_gen_energy, network_energy = [], [], [], []
    for k in total.keys():
        if "+[" in k: update_keys.append(k)
        elif "s[" in k: synapse_keys.append(k)
        elif "o[" in k: spike_gen_energy.append(k)
        elif "t[" in k: network_energy.append(k)

    analysis = {}
    analysis["update_energy"] = total[update_keys].sum()
    analysis["synapse_energy"] = total[synapse_keys].sum()
    analysis["spike_gen_energy"] = total[spike_gen_energy].sum()
    analysis["network_energy"] = total[network_energy].sum()

    return analysis


if __name__ == "__main__":
    neurons = []
    spiking_times = []
    spiking_update_energy = []
    spiking_spike_gen_energy = []
    spiking_synapse_energy = []
    spiking_network_energy = []

    # Use a pre-generated network for a realistic use case i.e. dvs-gesture
    analysis = run_sim()
    spiking_update_energy.append(analysis["update_energy"])
    spiking_spike_gen_energy.append(analysis["spike_gen_energy"])
    spiking_synapse_energy.append(analysis["synapse_energy"])
    spiking_network_energy.append(analysis["network_energy"])

    # Write all the simulation data to csv
    #with open("sim_network.csv", "w") as spiking_csv:
    #    spiking_writer = csv.DictWriter(spiking_csv,
    #                                    ("neurons", "update_energy", "spike_gen_energy",
    #                                     "synapse_energy", "network_energy"))
    #    spiking_writer.writeheader()
    #    for neuron_count, update_energy, spike_gen_energy, synapse_energy, network_energy in zip(neurons, spiking_update_energy,
    #                                              spiking_spike_gen_energy, spiking_synapse_energy,
    #                                              spiking_network_energy):
    #        spiking_writer.writerow({"neurons": neuron_count,
    #                                 "update_energy": update_energy,
    #                                 "spike_gen_energy": spike_gen_energy,
    #                                 "synapse_energy": synapse_energy,
    #                                 "network_energy": network_energy})

    plt.figure(figsize=(5.5, 5.5))
    plt.bar(1, spiking_update_energy)
    plt.bar(2, spiking_spike_gen_energy, color="orange")
    plt.bar(3, spiking_synapse_energy, color="red")
    plt.bar(4, spiking_network_energy, color="green")
    plt.xticks([1,2,3,4], ["Update", "Spike Gen", "Synapse", "Network"])
    plt.ylabel("Energy (J)")
    plt.xlabel("Operation Type")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("energy_breakdown.png")
    plt.show()
