"""
Copyright (c) 2023 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Run DVS Gesture and extract some network based statistics
"""
#import matplotlib
#matplotlib.use('Agg')

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os

#ARCH_FILENAME = "loihi_big.arch"
#NETWORK_FILENAME = "examples/dvs_gesture_big.net"

ARCH_FILENAME = "loihi.arch"
NETWORK_FILENAME = "examples/dvs_gesture_32x32_i16.net" # Input 16, 32x32 net
LOIHI_TIME_DATA_FILENAME = "runs/loihi_gesture_32x32_i16_time.csv"
LOIHI_ENERGY_DATA_FILENAME = "runs/loihi_gesture_32x32_i16_energy.csv"

def run_sim(timesteps):
    # Create the network to run
    fields = ["Neuron ID", "Core ID", "Threshold", "Reset",
              "Log Spikes", "Log Voltage", "Synapse Info..."]
    #timesteps = 100000
    command = ("./sim", ARCH_FILENAME,
               NETWORK_FILENAME, "{0}".format(timesteps),)
    print("Command: {0}".format(" ".join(command)))
    subprocess.call(command)

    return


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
    #analysis["synapse_energy"] = total[synapse_keys].sum()
    analysis["spike_gen_energy"] = total[spike_gen_energy].sum()
    analysis["network_energy"] = total[network_energy].sum()
    # Sum across all keys for each timestep
    #analysis["total_energies"] = (stats.loc[:, update_keys].sum(axis=1) +
    #    stats.loc[:, synapse_keys].sum(axis=1) +
    #    stats.loc[:, spike_gen_energy].sum(axis=1) +
    #    stats.loc[:, network_energy].sum(axis=1))
    #print(analysis["total_energies"])

    analysis["times"] = stats.loc[:, "time"]
    analysis["hops"] = stats.loc[:, "hops"]
    analysis["fired"] = stats.loc[:, "fired"]
    analysis["packets"] = stats.loc[:, "packets"]
    analysis["total_energies"] = stats.loc[:, "total_energy"]

    return analysis


def parse_loihi_spiketrains(total_timesteps):
    # Parse the CSV generated from the DVS gesture runs on Loihi
    #  The format is - first line is the neuron ID
    #  Second line is the timestep
    files = ("inputs.csv", "0Conv2D_15x15x16.csv", "1Conv2D_13x13x32.csv",
             "2Conv2D_11x11x64.csv", "3Conv2D_9x9x11.csv", "5Dense_11.csv")

    neurons = []
    timesteps = []

    for i in range(0, len(files)):
        f = files[i]
        path = "runs/spiketrains/" + f

        with open(path, "r") as spiketrain:
            reader = csv.reader(spiketrain)

            neurons += next(reader)
            timesteps += next(reader)

    spiketrain = {}
    for t in range(0, total_timesteps+1):
        spiketrain[t] = []

    for n, t in zip(neurons, timesteps):
        # Why are we -2 out of sync?
        #  Need to subtract 1, because the SNN toolbox starts from timestep 1,
        #   whereas my simulator goes from timestep 0
        t = int(t) - 2
        n = int(n)
        spiketrain[t].append(int(n))

    return spiketrain

if __name__ == "__main__":
    neurons = []
    spiking_times = []
    spiking_update_energy = []
    spiking_spike_gen_energy = []
    spiking_synapse_energy = []
    spiking_network_energy = []
    times = []
    energies = []
    timesteps = 126

    loihi_spiketrains = parse_loihi_spiketrains(timesteps)

    # Use a pre-generated network for a realistic use case i.e. dvs-gesture
    run_sim(timesteps)
    # Parse the detailed perf statistics
    print("Reading performance data")
    stats = pd.read_csv("perf.csv")
    analysis = parse_stats(stats)


    #spiking_update_energy.append(analysis["update_energy"])
    #spiking_spike_gen_energy.append(analysis["spike_gen_energy"])
    #spiking_synapse_energy.append(analysis["synapse_energy"])
    #spiking_network_energy.append(analysis["network_energy"])
    times = analysis["times"]
    energies = analysis["total_energies"]

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
    plt.rcParams.update({'font.size': 8, 'lines.markersize': 3})
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

    # Plot the latency
    loihi_data = pd.read_csv(LOIHI_TIME_DATA_FILENAME)
    loihi_times = loihi_data.loc[:, "spiking"] / 1.0e6
    plt.figure(figsize=(7, 8))
    plt.subplot(311)
    plt.plot(np.arange(1, timesteps+1), times, marker='x')
    plt.plot(np.arange(1, timesteps+1), loihi_times[0:timesteps], marker='x')
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.legend(("Simulated", "Measured on Loihi"))
    plt.ylabel("Latency (s)")
    plt.xlabel("Timestep")
    # Also plot underneath the different activity on the chip for both my
    #  simulator and Loihi
    plt.subplot(312)
    plt.plot(np.arange(1, timesteps+1), analysis["packets"], marker='x')
    plt.plot(np.arange(1, timesteps+1), analysis["hops"], marker='x')
    plt.legend(("Packets Sent", "Total Hops"))
    plt.xlabel("Timestep")


    plt.subplot(313)
    plt.plot(np.arange(1, timesteps+1), analysis["fired"], marker='x')
    # Figure out how many neurons fired in the Loihi data
    fired_count = [len(loihi_spiketrains[i]) for
                   i in range(0, len(loihi_spiketrains))]
    plt.plot(np.arange(1, timesteps+1), fired_count[0:timesteps], marker='x')
    plt.ylabel("Neurons Fired")
    plt.xlabel("Timestep")
    plt.legend(("Simulated", "Measured on Loihi"))

    print("diff = {}".format(
        analysis["fired"] - np.array(fired_count[0:timesteps])))
    plt.savefig("dvs_gesture_sim_time.png")

    # Plot the latency

    loihi_data = pd.read_csv(LOIHI_TIME_DATA_FILENAME)
    loihi_times = loihi_data.loc[:, "spiking"] / 1.0e6
    plt.figure(figsize=(5.0, 2.5))
    plt.plot(np.arange(1, timesteps+1), times, "-o")
    plt.plot(np.arange(1, timesteps+1), loihi_times[0:timesteps], "--x")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.legend(("Simulated", "Measured on Loihi"))
    plt.ylabel("Time-step Latency (s)")
    plt.xlabel("Time-step")
    plt.legend(("Simulated", "Measured on Loihi"))
    plt.tight_layout()
    plt.savefig("dvs_gesture_sim_time.pdf")

    # Plot the correlation between simulated and measured time-step latency
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(times, loihi_times[0:timesteps], "x")
    plt.plot(np.linspace(min(times), max(times)), np.linspace(min(times), max(times)), "k--")
    plt.xticks((1.0e-5, 1.5e-5, 2.0e-5, 2.5e-5, 3.0e-5))
    plt.yticks((1.0e-5, 1.5e-5, 2.0e-5, 2.5e-5, 3.0e-5))
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0,0))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.ylabel("Measured Latency (s)")
    plt.xlabel("Simulated Latency (s)")
    plt.tight_layout()
    plt.savefig("dvs_gesture_sim_correlation.pdf")
    plt.savefig("dvs_gesture_sim_correlation.png")

    # Calculate total error
    relative_error = abs(loihi_times[0:timesteps] - times) / loihi_times[0:timesteps]
    print(relative_error)
    mean_error = sum(relative_error) / len(relative_error)
    print("Time Absolute Mean error: {0} ({1} %)".format(mean_error, mean_error * 100))

    total_error =  (sum(loihi_times[0:timesteps]) - sum(times)) / sum(loihi_times[0:timesteps])
    print("Time Total error: {0} ({1} %)".format(total_error, total_error * 100))


    # TODO: I'm simulating dynamic energy, but the measurements are for static
    #  energy consumption, which pretty much just mirror the time simulation...
    loihi_data = pd.read_csv(LOIHI_ENERGY_DATA_FILENAME)
    loihi_energies = loihi_data.loc[:, "spiking"] / 1.0e6
    print(loihi_energies)
    plt.figure(figsize=(5.0, 2.5))
    plt.plot(np.arange(1, timesteps+1), energies, '-o')
    plt.plot(np.arange(1, timesteps+1), loihi_energies[0:timesteps], '--x')
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    plt.legend(("Simulated", "Measured on Loihi"))
    plt.ylabel("Energy (J)")
    plt.xlabel("Timestep")
    plt.savefig("dvs_gesture_sim_energy.png")


    relative_error = abs(loihi_energies[0:timesteps] - energies) / loihi_energies[0:timesteps]
    print(relative_error)
    mean_error = sum(relative_error) / len(relative_error)
    print("Energy Absolute Mean error: {0} ({1} %)".format(mean_error, mean_error * 100))

    total_error =  (sum(loihi_energies[0:timesteps]) - sum(energies)) / sum(loihi_energies[0:timesteps])
    print("Energy Total error: {0} ({1} %)".format(total_error, total_error * 100))



    # Plot the potential probes from simulation
    layers = ("inputs", "0Conv2D_15x15x16", "1Conv2D_13x13x32",
             "2Conv2D_11x11x64", "3Conv2D_9x9x11", "5Dense_11")
    layer_sizes = (1024, 3600, 5408, 7744, 891, 11)
    thresholds = (255, 293, 486, 510, 1729, 473)
    voltage_data = pd.read_csv("probe_potential.csv")

    plt.rcParams.update({'font.size': 12, 'lines.markersize': 5})
    plot_neurons = {"inputs": [], "0Conv2D_15x15x16": [50, 67], "1Conv2D_13x13x32": [], "2Conv2D_11x11x64": [], "3Conv2D_9x9x11": [], "5Dense_11":[]}
    for layer_id, layer in enumerate(layers):
        layer_path = "runs/potentials/{0}".format(layer)
        if not os.path.exists(layer_path):
            os.makedirs(layer_path)

        for neuron_id in plot_neurons[layer]:
            voltages = voltage_data.loc[:, "{0}.{1}".format(
                layer_id, neuron_id)]
            plt.figure(figsize=(5.5, 5.5))
            plt.plot(np.arange(1, timesteps+1), voltages*64, "-x")
            plt.plot(np.arange(1, timesteps+1),
                     np.ones(timesteps) * thresholds[layer_id] * 64)
            #plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
            plt.ylabel("Membrane Potential")
            plt.xlabel("Time-step")
            #plt.ylim((-0.2, 1.0))
            plt.tight_layout()
            plt.savefig("{0}/probe_voltage_{1}.png".format(layer_path,
                                                           neuron_id))
            plt.close()

    with open("stats.yaml", "r") as results_file:
       results = yaml.safe_load(results_file)
    #network_percentage = (results["network_time"] /
    #                                    results["time"]) * 100.0
    #print("Percentage of time used for only network: {0}".format(
    #      network_percentage))

    #plt.show()
