"""
Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Experiment to simulate Loihi running DVS gesture at varying levels of detail.

This adapts code from tcad2025/dvs_gesture.py for a different set of latency plots
"""
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import yaml
import copy

import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
sys.path.insert(0, PROJECT_DIR)
import sanafe

ARCH_FILENAME = "loihi.yaml"
NETWORK_FILENAME = "dvs_gesture_32x32.net.tagged"
LOIHI_TIME_DATA_FILENAME = "loihi_gesture_32x32_time.csv"
LOIHI_ENERGY_DATA_FILENAME = "loihi_gesture_32x32_energy.csv"
SIM_TIME_DATA_FILENAME = "sim_gesture_32x32_time.csv"
SIM_ENERGY_DATA_FILENAME = "sim_gesture_32x32_energy.csv"

NETWORK_DIR = os.path.join(PROJECT_DIR, "runs", "dvs")
INPUTS_DIR = os.path.join(PROJECT_DIR, "runs", "dvs", "loihi_gesture_32x32")
RUN_DIR = os.path.join(PROJECT_DIR, "runs", "detail")

ARCH_PATH = os.path.join(PROJECT_DIR, "arch", ARCH_FILENAME)
NETWORK_PATH = os.path.join(NETWORK_DIR, NETWORK_FILENAME)
LOIHI_TIME_DATA_PATH = os.path.join(RUN_DIR, LOIHI_TIME_DATA_FILENAME)
LOIHI_ENERGY_DATA_PATH = os.path.join(RUN_DIR, LOIHI_ENERGY_DATA_FILENAME)
SIM_TIME_DATA_PATH = os.path.join(RUN_DIR, SIM_TIME_DATA_FILENAME)
SIM_ENERGY_DATA_PATH = os.path.join(RUN_DIR, SIM_ENERGY_DATA_FILENAME)


def parse_stats(stats):
    print("Parsing statistics")
    analysis = {}
    analysis["hops"] = stats.loc[:, "hops"]
    analysis["fired"] = stats.loc[:, "fired"]
    analysis["packets"] = stats.loc[:, "packets"]
    analysis["times"] = stats.loc[:, "sim_time"]
    analysis["total_energy"] = sum(stats.loc[:, "total_energy"])

    print("Finished parsing statistics")
    return analysis


def run_dvs(inputs, detailed_network=False, detailed_synapse=False):
    if detailed_synapse:
        arch_filename = "loihi_detailed.yaml"
    elif detailed_network:
        arch_filename = "loihi_net.yaml"
    else:
        arch_filename = "loihi_baseline.yaml"

    arch = sanafe.load_arch(os.path.join(RUN_DIR, arch_filename))
    net = sanafe.load_net(NETWORK_PATH, arch, use_netlist_format=True)
    chip = sanafe.SpikingChip(arch)
    chip.load(net)

    for frame in range(0, frames):
        print(f"Running for input: {frame}")
        dvs_inputs = inputs[frame, :]
        mapped_inputs = chip.mapped_neuron_groups["0"]
        for id, mapped_neuron in enumerate(mapped_inputs):
            mapped_neuron.set_model_attributes(
                model_attributes={"bias": dvs_inputs[id]})
        is_first_frame = (frame == 0)
        chip.sim(timesteps, perf_trace="perf.csv", message_trace="messages.csv",
                    write_trace_headers=is_first_frame,
                    timing_model="cycle" if detailed_network else "detailed",
                    processing_threads=16)
        chip.reset()

    # Parse the detailed perf statistics
    print("Reading performance data")
    stats = pd.read_csv(os.path.join(PROJECT_DIR, "perf.csv"))
    analysis = parse_stats(stats)
    return analysis["times"]

if __name__ == "__main__":
    run_experiments = False
    plot_experiments = True

    neurons = []
    spiking_times = []
    spiking_update_energy = []
    spiking_spike_gen_energy = []
    spiking_synapse_energy = []
    spiking_network_energy = []
    timesteps = 128
    frames = 100
    #frames = 1

    if run_experiments:
        # Clear data file
        open(SIM_TIME_DATA_PATH, "w")

        # Write temporary architecture files for different models
        with open("arch/loihi.yaml", "r") as loihi_file:
            loihi_yaml = yaml.safe_load(loihi_file)

        # *************
        loihi_baseline = copy.deepcopy(loihi_yaml)
        loihi_baseline["architecture"]["attributes"]["latency_sync"] = 0.0
        with open(os.path.join(RUN_DIR, "loihi_baseline.yaml"), "w") as baseline_yaml:
            yaml.dump(loihi_baseline, baseline_yaml)

        # ************
        loihi_net = copy.deepcopy(loihi_yaml)
        with open(os.path.join(RUN_DIR, "loihi_net.yaml"), "w") as cycle_yaml:
            yaml.dump(loihi_net, cycle_yaml)

        loihi_detailed = copy.deepcopy(loihi_yaml)
        tile_section = loihi_detailed["architecture"]["tile"][0]
        core_section = tile_section["core"][0]
        synapse_section = None
        for section in core_section["synapse"]:
            if section["name"] == "loihi_conv_synapse":
                synapse_section = section
                break
        synapse_section["attributes"]["model"] = "loihi"
        synapse_section["attributes"]["plugin"] = "/home/usr1/jboyle/neuro/sana-fe/plugins/libloihi_synapse.so" # On detailed-loihi-model INRC branch
        synapse_section["attributes"]["max_parallel_accesses"] = 4
        del synapse_section["attributes"]["latency_process_spike"]
        with open(os.path.join(RUN_DIR, "loihi_detailed.yaml"), "w") as detailed_yaml:
            yaml.dump(loihi_detailed, detailed_yaml)

        input_csv_filename = os.path.join(INPUTS_DIR, "inputs.csv")
        with open(input_csv_filename, "r") as input_csv:
            inputs = np.loadtxt(input_csv, delimiter=",", skiprows=1)

        latencies = {}
        latencies["full"] = run_dvs(inputs, detailed_network=True, detailed_synapse=True)
        latencies["cycle_accurate"] = run_dvs(inputs, detailed_network=True, detailed_synapse=False)
        latencies["baseline"] = run_dvs(inputs, detailed_network=False, detailed_synapse=False)

        # Convert dict of arrays to DataFrame and save to CSV
        df = pd.DataFrame(latencies)
        df.to_csv(SIM_TIME_DATA_PATH, index=False)
        print("Finished running experiments")

    if plot_experiments:
        # Plot the latency
        plt.rcParams.update({'font.size': 7, 'lines.markersize': 4})
        sim_data = pd.read_csv(SIM_TIME_DATA_PATH, delimiter=",")
        sim_times_baseline = sim_data["baseline"].to_numpy()
        sim_times_cycle = sim_data["cycle_accurate"].to_numpy()
        sim_times_detailed = sim_data["full"].to_numpy()

        print("Reading Loihi data")
        loihi_data = pd.read_csv(LOIHI_TIME_DATA_PATH)
        loihi_times = np.array(loihi_data.loc[:, :] / 1.0e6)
        total_loihi_times = np.sum(loihi_times[0:128,:], axis=0)
        print(f"Total Loihi: {total_loihi_times}")
        print(f"Max Total Loihi: {np.max(total_loihi_times)}")
        print(f"Min Total Loihi: {np.min(total_loihi_times)}")

        # There is a weird effect, that the first sample of all inputs > 1 is
        #  a 0 value. Just ignore the entries for both arrays (so we have
        #  timestep-1)

        # Plot the latency
        print("Plotting latency")
        loihi_data = pd.read_csv(LOIHI_TIME_DATA_PATH)
        loihi_times = np.array(loihi_data.loc[:, :] / 1.0e6)
        sim_times_baseline = np.delete(sim_times_baseline,
                list(range(timesteps-1, timesteps*frames, timesteps)))
        sim_times_cycle = np.delete(sim_times_cycle,
                list(range(timesteps-1, timesteps*frames, timesteps)))
        sim_times_detailed = np.delete(sim_times_detailed,
                list(range(timesteps-1, timesteps*frames, timesteps)))
        loihi_times = loihi_times[0:timesteps-1,:]

        total_sim_times_baseline = np.zeros(frames)
        total_sim_times_cycle = np.zeros(frames)
        total_sim_times_detailed = np.zeros(frames)
        total_hops = np.zeros(frames)
        loihi_total_times = np.zeros(frames)
        for i in range(0, frames):
            total_sim_times_baseline[i] = np.sum(sim_times_baseline[i*(timesteps-1)+1:(i+1)*(timesteps-1)])
            total_sim_times_cycle[i] = np.sum(sim_times_cycle[i*(timesteps-1)+1:(i+1)*(timesteps-1)])
            total_sim_times_detailed[i] = np.sum(sim_times_detailed[i*(timesteps-1)+1:(i+1)*(timesteps-1)])
            loihi_total_times[i] = np.sum(loihi_times[0:timesteps-2, i])

        plt.figure(figsize=(7.0, 1.6))
        plt.rcParams.update({'font.size': 7, 'lines.markersize': 2, "font.family": "sans-serif", "font.sans-serif": "Arial", "pdf.fonttype": 42})
        start_frame = 0
        plt.plot(np.arange(1, timesteps-1), sim_times_baseline[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, "-", alpha=0.8, color="#0b3a55")
        #plt.plot(np.arange(1, timesteps-1), sim_times_cycle[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, "-.")#0b3a55
        plt.plot(np.arange(1, timesteps-1), sim_times_detailed[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, "-", alpha=0.8, color="#56B4E9")
        plt.plot(np.arange(1, timesteps-1), loihi_times[0:timesteps-2, start_frame] * 1.0e6, "k--")

        plt.legend(("Less detail", "Most detail", "Measured on Loihi"), fontsize=6)
        plt.ylabel("Time-step Latency (μs)")
        plt.xlabel("Time-step")
        plt.yticks(np.arange(0, 61, 10))
        plt.minorticks_on()
        plt.tight_layout(pad=0.3)
        plt.savefig("runs/detail/dvs_gesture_time_series.pdf")
        plt.savefig("runs/detail/dvs_gesture_time_series.png", dpi=300)

        # Plot the correlation between simulated and measured time-step latency
        plt.figure(figsize=(1.5, 1.5))
        plt.minorticks_on()
        plt.gca().set_box_aspect(1)

        # TODO: plot the total test-case/inference latency, this will be less confusing next to the time-series! (Will be measured in ms)
        average_sim_times_baseline = total_sim_times_baseline / 127
        average_sim_times_cycle = total_sim_times_cycle / 127
        average_sim_times_detailed = total_sim_times_detailed / 127
        loihi_average_times = loihi_total_times / 127
        plt.plot(average_sim_times_baseline[0:frames]*1.0e6, loihi_average_times[0:frames]*1.0e6, "o", alpha=0.8, markeredgewidth=0, markerfacecolor="#0b3a55")[0]
        #plt.plot(average_sim_times_cycle[0:frames]*1.0e6, loihi_average_times[0:frames]*1.0e6, "o", alpha=0.7, markeredgewidth=0, markerfacecolor="#009E73")[0]
        plt.plot(average_sim_times_detailed[0:frames]*1.0e6, loihi_average_times[0:frames]*1.0e6, "o", alpha=0.8, markeredgewidth=0, markerfacecolor="#56B4E9")[0]
        plt.plot(np.linspace(10, 30), np.linspace(10, 30), "k--", alpha=0.8)

        plt.legend(("Less detail", "Most detail"),
           fontsize=6, markerscale=1.5, handlelength=0.6, handletextpad=0.05, loc="upper left")
        plt.ylabel("Measured Latency (μs)")
        plt.xlabel("Simulated Latency (μs)")
        plt.xlim((10, 31))
        plt.ylim((10, 31))
        plt.xticks(np.arange(10, 31, 10))
        plt.yticks(np.arange(10, 31, 10))
        plt.tight_layout(pad=0.1)
        plt.savefig("runs/detail/dvs_gesture_compare_plugin.pdf")
        plt.savefig("runs/detail/dvs_gesture_compare_plugin.png", dpi=300)

        # Calculate total error
        print("Calculating errors")
        relative_error_baseline = np.abs(loihi_total_times - total_sim_times_baseline) / loihi_total_times
        relative_error_cycle = np.abs(loihi_total_times - total_sim_times_cycle) / loihi_total_times
        relative_error_detailed = np.abs(loihi_total_times - total_sim_times_detailed) / loihi_total_times


        mean_error_baseline = np.sum(relative_error_baseline) / len(relative_error_baseline)
        mean_error_cycle = np.sum(relative_error_cycle) / len(relative_error_cycle)
        mean_error_detailed = np.sum(relative_error_detailed) / len(relative_error_detailed)
        print("Time Absolute Mean error (Baseline): {0} ({1} %)".format(mean_error_baseline, mean_error_baseline * 100))
        print("Time Absolute Mean error (Cycle-accurate): {0} ({1} %)".format(mean_error_cycle, mean_error_cycle * 100))
        print("Time Absolute Mean error (Cycle+Detailed synapse): {0} ({1} %)".format(mean_error_detailed, mean_error_detailed * 100))

        total_error_baseline =  (np.sum(loihi_total_times) - np.sum(total_sim_times_baseline)) / np.sum(loihi_total_times)
        print("Time Total error (Baseline): {0} ({1} %)".format(total_error_baseline, total_error_baseline * 100))

        # Calculate correlation
        frame_correlation_baseline = np.corrcoef(loihi_average_times[0:frames], average_sim_times_baseline[0:frames])[0,1]
        frame_correlation_cycle = np.corrcoef(loihi_average_times[0:frames], average_sim_times_cycle[0:frames])[0,1]
        frame_correlation_detailed = np.corrcoef(loihi_average_times[0:frames], average_sim_times_detailed[0:frames])[0,1]
        print(f"Pearson correlation for baseline frame: {frame_correlation_baseline}")
        print(f"Pearson correlation for cycle-accurate frame: {frame_correlation_cycle}")
        print(f"Pearson correlation for detailed synapse frame: {frame_correlation_detailed}")

        #plt.show()
        print("Time simulations finished")

        """
        # Do some per time-step analysis since this should be more revealing
        predicted_latency_baseline = np.array(0)
        loihi_latency = np.array(0)
        for i in range(0, frames):
            predicted_latency_baseline = np.append(predicted_latency_baseline, sim_times_baseline[i*(timesteps-1)+1:(i+1)*(timesteps-1)])
            loihi_latency = np.append(loihi_latency, loihi_times[0:timesteps-2, i])
        per_timestep_correlation =  np.corrcoef(predicted_latency_baseline, loihi_latency)[0,1]
        print(f"Pearson correlation for ts: {per_timestep_correlation}")

        mask = loihi_latency != 0
        error = np.abs(loihi_latency[mask] - predicted_latency_baseline[mask]) / loihi_latency[mask]
        error = np.sum(error) / len(error)
        print(f"Timestep MAPE (baseline): {error} ({error * 100} %)")


        plt.figure(figsize=(3.0, 3.0))
        scatter = plt.plot(predicted_latency_baseline[0:1000]*1.0e6, loihi_latency[0:1000]*1.0e6, "o", alpha=0.8)[0]
        plt.plot(np.linspace(min(loihi_latency)*1.0e6, max(loihi_latency)*1.0e6),
                    np.linspace(min(loihi_latency)*1.0e6, max(loihi_latency)*1.0e6), "k--")
        plt.ylabel("Measured Latency ($\mu$s)")
        plt.xlabel("Simulated Latency ($\mu$s)")
        plt.tight_layout(pad=0.3)
        plt.savefig("runs/detail/dvs_gesture_sim_correlation_ts.pdf")
        plt.savefig("runs/detail/dvs_gesture_sim_correlation_ts.png", dpi=300)
        """
        # print("Showing plots")
