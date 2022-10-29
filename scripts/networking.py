"""Run DVS Gesture and extract some network based statistics

This is currently being used to explore dvs gesture rather than networking..
I should rename
"""
#import matplotlib
#matplotlib.use('Agg')

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

ARCH_FILENAME = "loihi.arch"
NETWORK_FILENAME = "examples/dvs_gesture_mini.net"

def run_sim(timesteps):
    # Create the network to run
    fields = ["Neuron ID", "Core ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
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
    pd.set_option('display.max_rows', None)
    print(total)
    update_keys, synapse_keys, spike_gen_energy, network_energy = [], [], [], []
    for k in total.keys():
        if "energy" in k:
            if "+[" in k: update_keys.append(k)
            elif "s[" in k: synapse_keys.append(k)
            elif "o[" in k: spike_gen_energy.append(k)
            elif "t[" in k: network_energy.append(k)

    analysis = {}
    analysis["update_energy"] = total[update_keys].sum()
    analysis["synapse_energy"] = total[synapse_keys].sum()
    analysis["spike_gen_energy"] = total[spike_gen_energy].sum()
    analysis["network_energy"] = total[network_energy].sum()
    analysis["times"] = stats.loc[:, "time"]

    return analysis


if __name__ == "__main__":
    neurons = []
    spiking_times = []
    spiking_update_energy = []
    spiking_spike_gen_energy = []
    spiking_synapse_energy = []
    spiking_network_energy = []
    times = []
    timesteps = 127

    # Use a pre-generated network for a realistic use case i.e. dvs-gesture
    analysis = run_sim(timesteps)
    spiking_update_energy.append(analysis["update_energy"])
    spiking_spike_gen_energy.append(analysis["spike_gen_energy"])
    spiking_synapse_energy.append(analysis["synapse_energy"])
    spiking_network_energy.append(analysis["network_energy"])
    times = analysis["times"]

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

    loihi_data = pd.read_csv("runs/loihi_gesture.csv")
    loihi_times = loihi_data.loc[:, "spiking"] / 1.0e6
    plt.figure(figsize=(15, 4))
    plt.plot(np.arange(1, timesteps+1), times, marker='x')
    plt.plot(np.arange(1, timesteps+1), loihi_times[0:timesteps], marker='x')
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.legend(("Simulated", "Measured on Loihi"))
    plt.ylabel("Latency (s)")
    plt.xlabel("Timestep")
    plt.savefig("dvs_gesture_sim.png")

    voltage_data = pd.read_csv("probe_potential.csv")
    voltages = voltage_data.loc[:, "2.0"]
    plt.figure(figsize=(5.5, 5.5))
    plt.plot(np.arange(1, timesteps+1), voltages)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.savefig("probe_potential.png")

    with open("stats.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)
    #network_percentage = (results["network_time"] /
    #                                    results["time"]) * 100.0
    #print("Percentage of time used for only network: {0}".format(
    #      network_percentage))

    #plt.show()
