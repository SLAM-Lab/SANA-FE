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


def run_dvs(inputs, detail_level="detailed", detailed_sync=False):
    if detailed_sync:
        arch_filename = "loihi_detailed.yaml"
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
                    timing_model=detail_level,
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
        tile_section = loihi_baseline["architecture"]["tile"][0]
        core_section = tile_section["core"][0]
        synapse_section = None
        for section in core_section["synapse"]:
            if section["name"] == "loihi_conv_synapse":
                synapse_section = section
                break
        loihi_baseline["architecture"]["attributes"]["latency_sync"] = 0.0
        synapse_section["attributes"]["model"] = "loihi"
        synapse_section["attributes"]["plugin"] = "/home/usr1/jboyle/neuro/sana-fe/plugins/libloihi_synapse.so" # On detailed-loihi-model INRC branch
        synapse_section["attributes"]["max_parallel_accesses"] = 4
        del synapse_section["attributes"]["latency_process_spike"]
        del synapse_section["attributes"]["energy_process_spike"]
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
        del synapse_section["attributes"]["energy_process_spike"]
        with open(os.path.join(RUN_DIR, "loihi_detailed.yaml"), "w") as detailed_yaml:
            yaml.dump(loihi_detailed, detailed_yaml)

        input_csv_filename = os.path.join(INPUTS_DIR, "inputs.csv")
        with open(input_csv_filename, "r") as input_csv:
            inputs = np.loadtxt(input_csv, delimiter=",", skiprows=1)

        latencies = {}
        latencies["cycle"] = run_dvs(inputs, detail_level="cycle", detailed_sync=True)
        latencies["detailed"] = run_dvs(inputs, detail_level="detailed", detailed_sync=False)
        latencies["simple"] = run_dvs(inputs, detail_level="simple", detailed_sync=False)

        # Convert dict of arrays to DataFrame and save to CSV
        df = pd.DataFrame(latencies)
        df.to_csv(SIM_TIME_DATA_PATH, index=False)
        print("Finished running experiments")

    if plot_experiments:
        # Plot the latency
        plt.rcParams.update({'font.size': 7, 'lines.markersize': 4})
        sim_data = pd.read_csv(SIM_TIME_DATA_PATH, delimiter=",")
        sim_times_cycle = sim_data["cycle"].to_numpy()
        sim_times_detailed = sim_data["detailed"].to_numpy()
        sim_times_simple = sim_data["simple"].to_numpy()

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
        sim_times_cycle = np.delete(sim_times_cycle,
                list(range(timesteps-1, timesteps*frames, timesteps)))
        sim_times_detailed = np.delete(sim_times_detailed,
                list(range(timesteps-1, timesteps*frames, timesteps)))
        sim_times_simple = np.delete(sim_times_simple,
                list(range(timesteps-1, timesteps*frames, timesteps)))
        loihi_times = loihi_times[0:timesteps-1,:]

        total_sim_times_cycle = np.zeros(frames)
        total_sim_times_detailed = np.zeros(frames)
        total_sim_times_simple = np.zeros(frames)
        total_hops = np.zeros(frames)
        loihi_total_times = np.zeros(frames)
        for i in range(0, frames):
            total_sim_times_cycle[i] = np.sum(sim_times_cycle[i*(timesteps-1)+1:(i+1)*(timesteps-1)])
            total_sim_times_detailed[i] = np.sum(sim_times_detailed[i*(timesteps-1)+1:(i+1)*(timesteps-1)])
            total_sim_times_simple[i] = np.sum(sim_times_simple[i*(timesteps-1)+1:(i+1)*(timesteps-1)])
            loihi_total_times[i] = np.sum(loihi_times[0:timesteps-2, i])

        plt.figure(figsize=(3.5, 1.0))
        plt.rcParams.update({'font.size': 7, 'lines.markersize': 2, "font.family": "sans-serif", "font.sans-serif": "Arial", "pdf.fonttype": 42})
        start_frame = 0
        # plt.plot(np.arange(1, timesteps-1), sim_times_cycle[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, "-.", linewidth=1.0, alpha=0.8, color="#002c59", zorder=5)
        # plt.plot(np.arange(1, timesteps-1), sim_times_detailed[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, "-", linewidth=1.0, alpha=0.8, color="#3983b8")
        # plt.plot(np.arange(1, timesteps-1), sim_times_simple[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, ":", linewidth=1.0, alpha=0.8, color="gray", zorder=10) # "#56b4e9"
        # plt.plot(np.arange(1, timesteps-1), loihi_times[0:timesteps-2, start_frame] * 1.0e6, "k--", linewidth=1.0)

        # plt.legend(("Cycle-accurate", "Semi-analytical", "Analytical"), fontsize=6)
        # plt.ylabel("Time-step Latency (μs)")
        # plt.xlabel("Time-step")
        # plt.yticks(np.arange(0, 61, 20))
        # #plt.minorticks_on()
        # plt.tight_layout(pad=0.1)
        # plt.savefig("runs/detail/dvs_gesture_time_series.pdf")
        # plt.savefig("runs/detail/dvs_gesture_time_series.png", dpi=300)

        plt.figure(figsize=(3.2, 1.2))  # Slightly increased height for legend
        plt.rcParams.update({'font.size': 7, 'lines.markersize': 2, "font.family": "sans-serif", "font.sans-serif": "Arial", "pdf.fonttype": 42})
        start_frame = 0
        plt.plot(np.arange(1, timesteps-1), sim_times_cycle[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, "-.", linewidth=1.0, alpha=0.8, color="#002c59", zorder=5)
        plt.plot(np.arange(1, timesteps-1), sim_times_detailed[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, "-", linewidth=1.0, alpha=0.8, color="#3983b8")
        plt.plot(np.arange(1, timesteps-1), sim_times_simple[start_frame*(timesteps-1)+1:(start_frame+1)*(timesteps-1)] * 1.0e6, ":", linewidth=1.0, alpha=0.8, color="gray", zorder=10)
        plt.plot(np.arange(1, timesteps-1), loihi_times[0:timesteps-2, start_frame] * 1.0e6, "k--", linewidth=1.0)

        plt.legend(("Cycle-accurate", "Semi-analytical", "Analytical", "Measured on Loihi"),
                   fontsize=6,
                   loc='upper center',
                   bbox_to_anchor=(0.5, 1.25),
                   ncol=4,
                   frameon=False,
                   handlelength=1.5,
                   handletextpad=0.1,
                   columnspacing=1.0)
        plt.ylabel("Time-step Latency (μs)")
        plt.xlabel("Time-step")
        plt.ylim(0, 70)  # Increased ylim to fit legend
        plt.yticks(np.arange(0, 61, 20))
        plt.tight_layout(pad=0.1)
        plt.subplots_adjust(top=0.82)  # Make room for legend at top
        plt.savefig("runs/detail/dvs_gesture_time_series.pdf")
        plt.savefig("runs/detail/dvs_gesture_time_series.png", dpi=300)

        # Plot the correlation between simulated and measured time-step latency
        plt.figure(figsize=(1.5, 1.5))
        plt.minorticks_on()
        plt.gca().set_box_aspect(1)

        # TODO: plot the total test-case/inference latency, this will be less confusing next to the time-series! (Will be measured in ms)
        average_sim_times_cycle = total_sim_times_cycle # / 127
        average_sim_times_detailed = total_sim_times_detailed # / 127
        average_sim_times_simple = total_sim_times_simple # / 127
        loihi_average_times = loihi_total_times # / 127
        # plt.plot(average_sim_times_detailed[0:frames]*1.0e6, loihi_average_times[0:frames]*1.0e6, "o", alpha=0.8, markeredgewidth=0, markerfacecolor="#0b3a55")[0]
        # plt.plot(average_sim_times_cycle[0:frames]*1.0e6, loihi_average_times[0:frames]*1.0e6, "o", alpha=0.8, markeredgewidth=0, markerfacecolor="#56B4E9")[0]

        plt.plot(average_sim_times_simple[0:frames]*1.0e3, loihi_average_times[0:frames]*1.0e3, "^", alpha=0.95, zorder=10, markeredgewidth=0.3, markeredgecolor='black', markerfacecolor='none')[0] # 0e4767
        plt.plot(average_sim_times_detailed[0:frames]*1.0e3, loihi_average_times[0:frames]*1.0e3, "o",  alpha=0.95, markeredgewidth=0, markerfacecolor="#56B4E9", zorder=0)[0]
        plt.plot(average_sim_times_cycle[0:frames]*1.0e3, loihi_average_times[0:frames]*1.0e3, "x", alpha=0.8, markersize=2.5, markeredgewidth=0.6, markeredgecolor='#0e4767')[0]
        plt.plot(np.linspace(1.5, 3.9), np.linspace(1.5, 3.9), "k--", alpha=0.7, linewidth=1.0)

        plt.legend(("Analytical", "Semi-analytical", "Cycle-accurate"),
           fontsize=5, markerscale=1.5, handlelength=0.6, handletextpad=0.05, loc="lower right")
        plt.ylabel("Measured Latency (ms)")
        plt.xlabel("Simulated Latency (ms)")
        plt.xlim((1.5, 3.9))
        plt.ylim((1.5, 3.9))
        plt.xticks(np.arange(1.5, 3.9, 0.5))
        plt.yticks(np.arange(1.5, 3.9, 0.5))
        plt.tight_layout(pad=0.1)
        plt.savefig("runs/detail/dvs_gesture_compare_plugin.pdf")
        plt.savefig("runs/detail/dvs_gesture_compare_plugin.png", dpi=300)

        # Calculate total error
        print("Calculating errors")
        relative_error_cycle = np.abs(loihi_total_times - total_sim_times_cycle) / loihi_total_times
        relative_error_detailed = np.abs(loihi_total_times - total_sim_times_detailed) / loihi_total_times
        relative_error_simple = np.abs(loihi_total_times - total_sim_times_simple) / loihi_total_times

        mean_error_simple = np.sum(relative_error_simple) / len(relative_error_simple)
        mean_error_cycle = np.sum(relative_error_cycle) / len(relative_error_cycle)
        mean_error_detailed = np.sum(relative_error_detailed) / len(relative_error_detailed)
        print("Time Absolute Mean error (Simple): {0} ({1} %)".format(mean_error_simple, mean_error_simple * 100))
        print("Time Absolute Mean error (Cycle-accurate): {0} ({1} %)".format(mean_error_cycle, mean_error_cycle * 100))
        print("Time Absolute Mean error (Detailed): {0} ({1} %)".format(mean_error_detailed, mean_error_detailed * 100))

        #total_error_baseline =  (np.sum(loihi_total_times) - np.sum(total_sim_times_baseline)) / np.sum(loihi_total_times)
        #print("Time Total error (Baseline): {0} ({1} %)".format(total_error_baseline, total_error_baseline * 100))

        # Calculate correlation
        frame_correlation_simple = np.corrcoef(loihi_average_times[0:frames], average_sim_times_simple[0:frames])[0,1]
        frame_correlation_cycle = np.corrcoef(loihi_average_times[0:frames], average_sim_times_cycle[0:frames])[0,1]
        frame_correlation_detailed = np.corrcoef(loihi_average_times[0:frames], average_sim_times_detailed[0:frames])[0,1]
        #print(f"Pearson correlation for baseline frame: {frame_correlation_baseline}")
        print(f"Pearson correlation for cycle-accurate frame: {frame_correlation_cycle}")
        print(f"Pearson correlation for detailed synapse frame: {frame_correlation_detailed}")

        import matplotlib
        matplotlib.rcParams['hatch.linewidth'] = 0.8
        # Create bar charts showing the simulator speed
        # Data for first bar chart (Scheduler model)
        scheduler_labels = ['Analytical', 'Semi-analytical', 'Cycle-accurate']
        scheduler_values = [4.12e8, 7.55e6, 1.59e5]
        scheduler_colors = ['#009E73', '#E69F00', '#56B4E9']  # Blue, Orange, Green
        schedule_alphas = [0.7, 0.8, 0.8]
        schedule_hatch = ['/////', None , 'xxxxxx']

        # Data for second bar chart (Architecture)
        arch_labels = ['Loihi', 'Loihi-Ind', 'Loihi-IMAC', 'Loihi-DD', 'Loihi-AD']
        arch_values = [7.55e6, 1.44e6, 7.12e5, 5.12e6, 2.39e6]
        arch_colors = ['#56B4E9', '#D55E00', '#009E73', '#CC79A7', '#E69F00']
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(1.5, 1.5))

        # Adjust spacing
        #plt.subplots_adjust(wspace=0.4, left=0.15, right=0.95, top=0.85, bottom=0.15)

        # First bar chart - Scheduler model (LOGARITHMIC)
        n1 = len(scheduler_labels)
        bar_width = 1.0  # No space between bars in cluster
        x1 = np.arange(n1) * bar_width  # Bars touching each other

        bars1 = ax1.bar(x1, np.array(scheduler_values)/1e6, color=['white', '#56B4E9', 'white'], alpha=0.8, width=bar_width, edgecolor='black', linewidth=0.5)
        # Overlay hatches with white edge (no fill)
        ax1.bar(x1, np.array(scheduler_values)/1e6, color='none', alpha=0.8, hatch=schedule_hatch, width=bar_width, edgecolor='#56B4E9', linewidth=0)

        # Add name labels on top of bars (vertical text, same color as bar)
        for bar, label, color, alpha, hatch in zip(bars1, scheduler_labels, scheduler_colors, schedule_alphas, schedule_hatch):
            if label == "Analytical":
                height = 30
            else:
                height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    label, ha='center', va='bottom', rotation=90,
                    color='black', fontsize=6, alpha=1.0)

        ax1.set_xlabel('Timing Model', fontsize=7)
        ax1.set_ylabel('SANA-FE Throughput (Mspike/s)', fontsize=7, va='center')
        ax1.set_xticks([])  # Remove x ticks since labels are on bars
        ax1.tick_params(axis='both', labelsize=6)

        # Logarithmic scale with no minor ticks
        ax1.set_yscale('log')
        #ax1.yaxis.set_minor_locator()

        # Second bar chart - Architecture (LINEAR)
        n2 = len(arch_labels)
        x2 = np.arange(n2) * bar_width  # Bars touching each other

        bars2 = ax2.bar(x2, np.array(arch_values)/1e6, color=arch_colors, width=bar_width, edgecolor='black', linewidth=0.5)

        # Add name labels on top of bars (vertical text, same color as bar)
        for bar, label, color in zip(bars2, arch_labels, arch_colors):
            height = bar.get_height() + 0.1
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    label,
                    ha='center', va='bottom', rotation=90,
                    color=color, fontsize=6)

        ax2.set_xlabel('Arch.', fontsize=7)
        ax2.set_xticks([])  # Remove x ticks since labels are on bars
        ax2.set_ylim([0, 10])
        ax2.tick_params(axis='both', labelsize=6)

        # Linear scale for architecture
        # No set_yscale needed (linear is default)

        # Remove top and right spines for cleaner look
        for ax in [ax1, ax2]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        plt.tight_layout(pad=0.1)
        plt.savefig("runs/detail/compare_models.pdf")
        plt.savefig("runs/detail/compare_models.png")


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
