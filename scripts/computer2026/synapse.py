"""
Copyright (c) 2026 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Explore synaptic concurrency within a Loihi-based architecture
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import time
import matplotlib
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
# Setup paths
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    sys.path.insert(0, PROJECT_DIR)
    import sanafe

parser = argparse.ArgumentParser(
        prog="Synaptic-concurrency",
        description="Explore impact of synaptic concurrency on Loihi.")
parser.add_argument("data_path", nargs="?",
                    default=os.path.abspath((os.path.join(PROJECT_DIR, "runs", "synapse"))))
parser.add_argument("run_path", nargs="?",
                    default=os.path.abspath((os.path.join(PROJECT_DIR, "runs", "synapse"))))
parser.add_argument("--run", action="store_true")
# TODO: Both of these are ignored for now, to keep interfaces consistent
parser.add_argument("--quick", action="store_true")
parser.add_argument("--plot", action="store_true")

args = parser.parse_args()

# Configuration
ARCH_FILENAME = "loihi.yaml"
NETWORK_FILENAME = "dvs_gesture_32x32.net.tagged"
DATA_PATH = args.data_path
RUN_DIR = args.run_path
RUN_EXPERIMENTS = args.run

ARCH_PATH = os.path.join(PROJECT_DIR, "arch", ARCH_FILENAME)
GENERATED_NETWORK_PATH = os.path.join(DATA_PATH, NETWORK_FILENAME)

# Simulation parameters
TIMESTEPS = 128
PARALLEL_ACCESS_VALUES = list(range(1, 17))  # [1, 2, 3, ..., 16]
PARALLELIZATION_MODIFIERS = (0, 25, 30, 50, 60, 75, 90, 100)  # %

def create_modified_arch_file(original_arch_path, max_parallel_accesses,
                              run_dir, parallelization_modifier):
    """
    Create a modified architecture YAML file with specified max_parallel_accesses
    """
    with open(original_arch_path, 'r') as f:
        arch_data = yaml.safe_load(f)

    # Navigate to the synapse section and update max_parallel_accesses
    for tile in arch_data['architecture']['tile']:
        for core in tile['core']:
            for synapse in core['synapse']:
                if synapse['name'] == 'loihi_conv_synapse':
                    synapse['attributes']['plugin'] = os.path.join(
                        PROJECT_DIR, "plugins", "libloihi_synapse.so")
                    synapse['attributes']['model'] = "loihi"
                    synapse['attributes']['max_parallel_accesses'] = max_parallel_accesses
                    synapse['attributes']['parallelizable_ratio'] = float(parallelization_modifier) / 100
                    synapse['attributes']['duplication_ratio'] = float(parallelization_modifier) / 100
                    del synapse["attributes"]["latency_process_spike"]
                    del synapse["attributes"]["energy_process_spike"]

    # Create new filename
    new_arch_filename = f"loihi_parallel_{max_parallel_accesses}_{parallelization_modifier}.yaml"
    new_arch_path = os.path.join(run_dir, new_arch_filename)

    # Save modified architecture file
    with open(new_arch_path, 'w') as f:
        yaml.dump(arch_data, f, default_flow_style=False)

    return new_arch_path


def run_single_experiment(max_parallel_accesses, parallelization_modifier):
    """
    Run a single experiment with specified max_parallel_accesses value
    """
    print(f"Starting experiment with max_parallel_accesses={max_parallel_accesses}")

    try:
        # Create modified architecture file
        modified_arch_path = create_modified_arch_file(
            ARCH_PATH, max_parallel_accesses, RUN_DIR, parallelization_modifier
        )

        # Load architecture and network
        arch = sanafe.load_arch(modified_arch_path)
        net = sanafe.load_net(GENERATED_NETWORK_PATH, arch, use_netlist_format=True)
        chip = sanafe.SpikingChip(arch)
        chip.load(net)

        # Setup performance tracking files
        perf_filename = f"perf_parallel_{max_parallel_accesses}_{parallelization_modifier}.csv"
        messages_filename = f"messages_parallel_{max_parallel_accesses}_{parallelization_modifier}.csv"
        perf_path = os.path.join(RUN_DIR, perf_filename)
        messages_path = os.path.join(RUN_DIR, messages_filename)

        start_time = time.perf_counter()
        # Run simulation
        chip.sim(
            TIMESTEPS,
            perf_trace=perf_path,
            message_trace=messages_path,
            write_trace_headers=True,
            timing_model="detailed",
            processing_threads=8,
            scheduler_threads=1
        )
        chip.reset()
        runtime = time.perf_counter() - start_time

        # Parse results
        stats = pd.read_csv(perf_path)
        total_time = stats.loc[:, "sim_time"].sum()
        total_energy = stats.loc[:, "total_energy"].sum()
        total_hops = stats.loc[:, "hops"].sum()
        total_packets = stats.loc[:, "packets"].sum()
        total_fired = stats.loc[:, "fired"].sum()
        total_spikes = stats.loc[:, "spikes"].sum()

        result = {
            'max_parallel_accesses': max_parallel_accesses,
            'parallelization_modifier': parallelization_modifier,
            'total_time': total_time,
            'total_energy': total_energy,
            'total_hops': total_hops,
            'total_packets': total_packets,
            'total_fired': total_fired,
            'total_spikes': total_spikes,
            'avg_time_per_timestep': total_time / TIMESTEPS,
            'avg_energy_per_timestep': total_energy / TIMESTEPS,
            'runtime_seconds': runtime
        }

        print(f"Completed experiment with max_parallel_accesses={max_parallel_accesses}")
        print(f"  Total time: {total_time:.2e}s")
        print(f"  Total energy: {total_energy:.2e}J")

        return result

    except Exception as e:
        print(f"Error in experiment with max_parallel_accesses={max_parallel_accesses}: {str(e)}")
        return {
            'max_parallel_accesses': max_parallel_accesses,
            'error': str(e)
        }


def main():
    """
    Main function to run parallel experiments
    """
    print("Starting DVS Gesture parallel architecture exploration")
    print(f"Running experiments with max_parallel_accesses values: {PARALLEL_ACCESS_VALUES}")
    print(f"Running experiments with {PARALLELIZATION_MODIFIERS}% parallelization")
    results_path = os.path.join(RUN_DIR, "parallel_exploration_results.csv")

    # Ensure run directory exists
    os.makedirs(RUN_DIR, exist_ok=True)

    sweep_start = time.perf_counter()
    if RUN_EXPERIMENTS:
        # Run experiments in parallel
        results = []
        for p in PARALLELIZATION_MODIFIERS:  # % of synaptic energy/latency
            for i, max_parallel_accesses in enumerate(PARALLEL_ACCESS_VALUES):
                print(f"\nProgress: {i+1}/{len(PARALLEL_ACCESS_VALUES)}")
                result = run_single_experiment(max_parallel_accesses, p)
                results.append(result)

        # Sort results by max_parallel_accesses for easier analysis
        results.sort(key=lambda x: x['max_parallel_accesses'])

        # Print summary
        print("\nSummary of results:")
        print("Max Parallel | Total Time (s) | Total Energy (J) | Avg Time/Step (s)")
        print("-" * 70)
        for result in results:
            if 'error' not in result:
                print(f"{result['max_parallel_accesses']:11d} | "
                    f"{result['total_time']:13.2e} | "
                    f"{result['total_energy']:15.2e} | "
                    f"{result['avg_time_per_timestep']:16.2e}")
            else:
                print(f"{result['max_parallel_accesses']:11d} | ERROR: {result['error']}")

        # Sort results by max_parallel_accesses for easier analysis
        #results.sort(key=lambda x: x['max_parallel_accesses'])

        # Save results to CSV
        results_df = pd.DataFrame(results)
        results_df.to_csv(results_path, index=False)

        print(f"\nAll experiments completed!")
        print(f"Results saved to: {results_path}")

        sweep_end = time.perf_counter()
        sweep_runtime = sweep_end - sweep_start
        print(f"Total sweep runtime: {sweep_runtime} s")


    # Plotting
    results_df = pd.read_csv(results_path)
    # Analyze and plot latency results
    #concurrency, total_latency, avg_latency = analyze_and_plot_latency(results)

    plt.rcParams.update({
        "font.size": 7,
        "font.family": "sans-serif",
        "font.sans-serif": "Arial",
        "pdf.fonttype": 42
    })

    # Split by calibration modifier
    df_0 = results_df[results_df["parallelization_modifier"] == 0].sort_values("parallelization_modifier")
    df_25 = results_df[results_df["parallelization_modifier"] == 25].sort_values("parallelization_modifier")
    df_30 = results_df[results_df["parallelization_modifier"] == 30].sort_values("parallelization_modifier")
    df_50 = results_df[results_df["parallelization_modifier"] == 50].sort_values("parallelization_modifier")
    df_60 = results_df[results_df["parallelization_modifier"] == 60].sort_values("parallelization_modifier")
    df_75 = results_df[results_df["parallelization_modifier"] == 75].sort_values("parallelization_modifier")
    df_90 = results_df[results_df["parallelization_modifier"] == 90].sort_values("parallelization_modifier")
    df_100 = results_df[results_df["parallelization_modifier"] == 100].sort_values("parallelization_modifier")
    x = df_100["max_parallel_accesses"]

    # --- Latency (Frames/s) plot ---
    fig = plt.figure(figsize=(1.5, 1.5))
    y_0 = 1.0 / df_0["total_time"].values
    y_25 = 1.0 / df_25["total_time"].values
    y_50 = 1.0 / df_50["total_time"].values
    y_75 = 1.0 / df_75["total_time"].values
    y_100 = 1.0 / df_100["total_time"].values

    n_blues = 3
    blues_colors = [matplotlib.colormaps['Blues'](v) for v in np.linspace(0.6, 0.9, n_blues)]
    colors = blues_colors
    colors.extend(['black',])

    plt.fill_between(x, y_25, y_100, alpha=0.25, color='#56B4E9', linewidth=0)
    #plt.plot(x, y_0, 'o-', linewidth=1, markersize=2, color=colors[0])
    plt.plot(x, y_25, 'o-', linewidth=1, markersize=2, color=colors[3])
    plt.plot(x, y_50, 'o-', linewidth=1, markersize=2, color=colors[2])
    plt.plot(x, y_75, 'o-', linewidth=1, markersize=2, color=colors[1])
    plt.plot(x, y_100, 'o-', linewidth=1, markersize=2, color=colors[0])

    x_val = 4
    y_val = y_100[x.values == x_val][0]
    plt.plot([0, x_val], [y_val, y_val], ':', color='gray', linewidth=1.5)
    plt.plot([x_val, x_val], [plt.ylim()[0], y_val], ':', color='gray', linewidth=1.5)

    ax = fig.axes[0]
    ax.ticklabel_format(style='plain', useOffset=False, axis='y')
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    plt.xlabel('Max. Parallel Reads')
    plt.ylabel('Simulated Frames/s')
    plt.grid(True, alpha=0.3)
    plt.xticks(range(0, 17, 4))
    plt.xlim((0, 17))
    plt.ylim((50, 650))

    # Annotate each line at the right end with colored text
    plt.text(17, y_25[-1]-35, '$p=0.25$', color=colors[3], fontsize=6, va='center', ha='right')
    plt.text(17, y_50[-1]+16, '$p=0.50$', color=colors[2], fontsize=6, va='center', ha='right')
    plt.text(17, y_75[-1]+25, '$p=0.75$', color=colors[1], fontsize=6, va='center', ha='right')
    plt.text(17, y_100[-1]+25, '$p=1.00$', color=colors[0], fontsize=6, va='center', ha='right')

    plt.tight_layout(pad=0.1)
    plt.savefig(os.path.join(RUN_DIR, "fig_5a_concurrency_loihi_latency.png"), dpi=300)
    plt.savefig(os.path.join(RUN_DIR, "fig_5a_concurrency_loihi_latency.pdf"))

    # --- Energy plot ---
    n_blues = 2
    blues_colors = [matplotlib.colormaps['Blues'](v) for v in np.linspace(0.6, 0.9, n_blues)]
    colors = blues_colors
    colors.extend(['black',])  # S=0 gets black


    fig = plt.figure(figsize=(1.5, 1.5))
    e_0 = df_0["total_energy"].values * 1.0e3
    e_30 = df_30["total_energy"].values * 1.0e3
    e_60 = df_60["total_energy"].values * 1.0e3
    e_90 = df_90["total_energy"].values * 1.0e3

    plt.fill_between(x, np.minimum(e_30, e_90), np.maximum(e_30, e_90),
                    alpha=0.25, color='#56B4E9', linewidth=0)
    # plt.plot(x, e_0, 'o-', linewidth=1, markersize=2, color="gray")
    plt.plot(x, e_90, 'o-', linewidth=1, markersize=2, color=colors[2])
    plt.plot(x, e_60, 'o-', linewidth=1, markersize=2, color=colors[1])
    plt.plot(x, e_30, 'o-', linewidth=1, markersize=2, color=colors[0])
    e_val = e_30[x.values == x_val][0]
    plt.plot([0, x_val], [e_val, e_val], ':', color='gray', linewidth=1.5)
    plt.plot([x_val, x_val], [plt.ylim()[0], e_val], ':', color='gray', linewidth=1.5)

    ax = fig.axes[0]
    ax.ticklabel_format(style='plain', useOffset=False, axis='y')
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    plt.xlabel('Max. Parallel Reads')
    plt.ylabel('H/W Energy Usage (mJ)')
    plt.grid(True, alpha=0.3)
    plt.xticks(range(0, 17, 4))
    plt.xlim((0, 18))
    plt.ylim((0.2, 0.75))

    plt.text(18, e_90[-1]+0.022, '$s=E_{\\mathrm{static}}/E_{\\mathrm{total}}=0.9$', color=colors[2], fontsize=6, va='center', ha='right')
    plt.text(18, e_60[-1]+0.022, '$s=0.6$', color=colors[1], fontsize=6, va='center', ha='right')
    plt.text(18, e_30[-1]+0.022, '$s=0.3$', color=colors[0], fontsize=6, va='center', ha='right')

    plt.tight_layout(pad=0.1)
    plt.savefig(os.path.join(RUN_DIR, "fig_5b_concurrency_loihi_energy.png"), dpi=300)
    plt.savefig(os.path.join(RUN_DIR, "fig_5b_concurrency_loihi_energy.pdf"))

    # plt.show()
    print(results_df)


if __name__ == "__main__":
    main()
