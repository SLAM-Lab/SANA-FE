"""
DVS Gesture Parallel Runner - Streamlined version for parallel architecture exploration
"""

import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path

# Setup paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
sys.path.insert(0, PROJECT_DIR)
import sanafe

# Configuration
ARCH_FILENAME = "loihi.yaml"
NETWORK_FILENAME = "dvs_gesture_32x32.net.tagged"
DVS_RUN_DIR = os.path.join(PROJECT_DIR, "runs", "synapse")

ARCH_PATH = os.path.join(PROJECT_DIR, "arch", ARCH_FILENAME)
GENERATED_NETWORK_PATH = os.path.join(DVS_RUN_DIR, NETWORK_FILENAME)

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
                    synapse['attributes']['max_parallel_accesses'] = max_parallel_accesses

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
            ARCH_PATH, max_parallel_accesses, DVS_RUN_DIR
        )

        # Load architecture and network
        arch = sanafe.load_arch(modified_arch_path)
        net = sanafe.load_net(GENERATED_NETWORK_PATH, arch, use_netlist_format=True)
        chip = sanafe.SpikingChip(arch)
        chip.load(net)

        # Setup performance tracking files
        perf_filename = f"perf_parallel_{max_parallel_accesses}.csv"
        messages_filename = f"messages_parallel_{max_parallel_accesses}.csv"
        perf_path = os.path.join(DVS_RUN_DIR, perf_filename)
        messages_path = os.path.join(DVS_RUN_DIR, messages_filename)

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
            'avg_energy_per_timestep': total_energy / TIMESTEPS
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
    print(f"Using {len(PARALLEL_ACCESS_VALUES)} parallel threads")

    # Ensure run directory exists
    os.makedirs(DVS_RUN_DIR, exist_ok=True)

    # Run experiments in parallel
    results = []
    for i, max_parallel_accesses in enumerate(PARALLEL_ACCESS_VALUES):
        print(f"\nProgress: {i+1}/{len(PARALLEL_ACCESS_VALUES)}")
        result = run_single_experiment(max_parallel_accesses)
        results.append(result)

    # Sort results by max_parallel_accesses for easier analysis
    results.sort(key=lambda x: x['max_parallel_accesses'])

    # Save results to CSV
    results_df = pd.DataFrame(results)
    results_path = os.path.join(DVS_RUN_DIR, "parallel_exploration_results.csv")
    results_df.to_csv(results_path, index=False)

    print(f"\nAll experiments completed!")
    print(f"Results saved to: {results_path}")

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


    import matplotlib.pyplot as plt
    # Sort results by max_parallel_accesses for easier analysis
    results.sort(key=lambda x: x['max_parallel_accesses'])

    # Save results to CSV
    results_df = pd.DataFrame(results)
    results_path = os.path.join(DVS_RUN_DIR, "parallel_exploration_results.csv")
    results_df.to_csv(results_path, index=False)

    print(f"\nAll experiments completed!")
    print(f"Results saved to: {results_path}")

    # Analyze and plot latency results
    #concurrency, total_latency, avg_latency = analyze_and_plot_latency(results)

    # plt.figure(figsize=(6, 4))
    # plt.plot(concurrency, avg_latency * 1e6, 'o-', linewidth=2, markersize=8, color='blue')
    # plt.xlabel('Max Parallel Accesses', fontsize=12)
    # plt.ylabel('Average Timestep Latency (μs)', fontsize=12)
    # plt.title('Hardware Latency vs Concurrency Level', fontsize=14)
    # plt.grid(True, alpha=0.3)
    # plt.xticks(range(1, max(concurrency)+1, 2))
    # plt.tight_layout()
    # plt.savefig(os.path.join(DVS_RUN_DIR, "latency_concurrency_clean.png"), dpi=300, bbox_inches='tight')
    # plt.savefig(os.path.join(DVS_RUN_DIR, "latency_concurrency_clean.pdf"), bbox_inches='tight')
    # plt.show()



if __name__ == "__main__":
    main()
