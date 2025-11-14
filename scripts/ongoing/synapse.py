"""
DVS Gesture Parallel Runner - Streamlined version for parallel architecture exploration
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import time


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
# Setup paths
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    sys.path.insert(0, PROJECT_DIR)
    import sanafe

# Configuration
ARCH_FILENAME = "loihi.yaml"
NETWORK_FILENAME = "dvs_gesture_32x32.net.tagged"
RUN_DIR = os.path.join(PROJECT_DIR, "runs", "synapse")

ARCH_PATH = os.path.join(PROJECT_DIR, "arch", ARCH_FILENAME)
GENERATED_NETWORK_PATH = os.path.join(RUN_DIR, NETWORK_FILENAME)

# Simulation parameters
TIMESTEPS = 128
PARALLEL_ACCESS_VALUES = list(range(1, 17))  # [1, 2, 3, ..., 16]

def create_modified_arch_file(original_arch_path, max_parallel_accesses, run_dir):
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
                    synapse['attributes']['plugin'] = "/home/usr1/jboyle/neuro/sana-fe/plugins/libloihi_synapse.so"
                    synapse['attributes']['model'] = "loihi"
                    synapse['attributes']['max_parallel_accesses'] = max_parallel_accesses
                    del synapse["attributes"]["latency_process_spike"]

    # Create new filename
    new_arch_filename = f"loihi_parallel_{max_parallel_accesses}.yaml"
    new_arch_path = os.path.join(run_dir, new_arch_filename)

    # Save modified architecture file
    with open(new_arch_path, 'w') as f:
        yaml.dump(arch_data, f, default_flow_style=False)

    return new_arch_path


def run_single_experiment(max_parallel_accesses):
    """
    Run a single experiment with specified max_parallel_accesses value
    """
    print(f"Starting experiment with max_parallel_accesses={max_parallel_accesses}")

    try:
        # Create modified architecture file
        modified_arch_path = create_modified_arch_file(
            ARCH_PATH, max_parallel_accesses, RUN_DIR
        )

        # Load architecture and network
        arch = sanafe.load_arch(modified_arch_path)
        net = sanafe.load_net(GENERATED_NETWORK_PATH, arch, use_netlist_format=True)
        chip = sanafe.SpikingChip(arch)
        chip.load(net)

        # Setup performance tracking files
        perf_filename = f"perf_parallel_{max_parallel_accesses}.csv"
        messages_filename = f"messages_parallel_{max_parallel_accesses}.csv"
        perf_path = os.path.join(RUN_DIR, perf_filename)
        messages_path = os.path.join(RUN_DIR, messages_filename)

        start_time = time.perf_counter()
        # Run simulation
        chip.sim(
            TIMESTEPS,
            perf_trace=perf_path,
            message_trace=messages_path,
            write_trace_headers=True,
            timing_model="cycle",
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

        result = {
            'max_parallel_accesses': max_parallel_accesses,
            'total_time': total_time,
            'total_energy': total_energy,
            'total_hops': total_hops,
            'total_packets': total_packets,
            'total_fired': total_fired,
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
    run_experiments = False
    print("Starting DVS Gesture parallel architecture exploration")
    print(f"Running experiments with max_parallel_accesses values: {PARALLEL_ACCESS_VALUES}")
    results_path = os.path.join(RUN_DIR, "parallel_exploration_results.csv")

    # Ensure run directory exists
    os.makedirs(RUN_DIR, exist_ok=True)

    sweep_start = time.perf_counter()
    if run_experiments:
        # Run experiments in parallel
        results = []
        for i, max_parallel_accesses in enumerate(PARALLEL_ACCESS_VALUES):
            print(f"\nProgress: {i+1}/{len(PARALLEL_ACCESS_VALUES)}")
            result = run_single_experiment(max_parallel_accesses)
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

    print(results_df)
    baseline = results_df[results_df["max_parallel_accesses"] == 1]["total_time"].to_numpy()
    results_df["speedup"] = baseline / results_df["total_time"]
    fig = plt.figure(figsize=(1.5, 1.5))
    plt.plot(results_df["max_parallel_accesses"], results_df["speedup"], 'o-',
             linewidth=1, markersize=2, color='#56B4E9')
    ax = fig.axes[0]
    ax.ticklabel_format(style='plain', useOffset=False, axis='y')
    plt.xlabel('Max. Parallel Reads')
    plt.ylabel('Hardware Speed-up')
    plt.grid(True, alpha=0.3)
    plt.xticks(range(0, 17, 4))
    plt.xlim((0, 17))
    plt.tight_layout(pad=0.1)
    plt.savefig(os.path.join(RUN_DIR, "latency_concurrency.png"), dpi=300)
    plt.savefig(os.path.join(RUN_DIR, "latency_concurrency.pdf"))
    plt.show()


if __name__ == "__main__":
    main()
