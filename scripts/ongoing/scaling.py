"""
Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
sys.path.insert(0, PROJECT_DIR)
import sanafe

# Okabe-Ito color palette (colorblind-friendly)
okabe_ito_colors = [
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermillion
    '#CC79A7',  # Reddish Purple
    '#000000'   # Black
]

def run_performance_experiments():
    """Run simulator performance experiments with varying thread counts."""

    # Configuration
    ARCH_PATH = "arch/loihi_large.yaml"
    #NETWORK_PATH = "dvs_gesture_32x32_sep23.yaml"
    NETWORK_PATH = "snn/fly.net"

    # Experiment parameters
    processing_threads_range = range(1, 49, 1)
    #processing_threads_range = range(32, 49, 1)
    scheduler_threads_range = [0, 1, 2] + list(range(4, 33, 4))   # Different scheduler thread counts
    #scheduler_threads_range = [2,] # Different scheduler thread counts
    num_repeats = 5
#    num_repeats = 1
    #timesteps = 10000
    timesteps = 1000

    # Results storage
    results = []

    print("Loading architecture and network...")
    arch = sanafe.load_arch(ARCH_PATH)
    net = sanafe.load_net(NETWORK_PATH, arch, use_netlist_format=True)

    # Load input data (first frame only)
    print("Starting performance experiments...")

    total_experiments = len(processing_threads_range) * len(scheduler_threads_range) * num_repeats
    experiment_count = 0

    for repeat in range(num_repeats):
        print(f"  Repeat {repeat + 1}/{num_repeats} ({experiment_count}/{total_experiments})")
        for proc_threads in processing_threads_range:
            for sched_threads in scheduler_threads_range:
                print(f"Testing: {proc_threads} processing threads, {sched_threads} scheduler threads")

                experiment_count += 1
                # Create fresh chip instance for each run
                chip = sanafe.SpikingChip(arch)
                chip.load(net)

                # Time the simulation
                start_time = time.perf_counter()
                chip.sim(timesteps,
                        processing_threads=proc_threads,
                        scheduler_threads=sched_threads)
                end_time = time.perf_counter()

                runtime = end_time - start_time
                print(f"runtime: {runtime}")

                # Store results
                results.append({
                    "processing_threads": proc_threads,
                    "scheduler_threads": sched_threads,
                    "repeat": repeat,
                    "runtime_seconds": runtime
                })

                # Reset chip for next run
                chip.reset()

    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Save to CSV
    output_filename = "runs/scaling/simulator_performance_results.csv"
    df.to_csv(output_filename, index=False)
    print(f"Results saved to {output_filename}")

    return df


def plot_results(df, show_plot=True, save_plot=True):
    """Create performance plots from the results DataFrame."""

    df = df[df["scheduler_threads"] <= 20]

    # Calculate mean and std for each configuration
    summary = df.groupby(["processing_threads", "scheduler_threads"])["runtime_seconds"].agg(["mean", "std"]).reset_index()

    plt.figure(figsize=(3, 3))

    #line_styles = (("-", "--", "-.", ":", "-", "-"))
    #markers = (("o", "x", "^", "s", "o", "o"))

    # Plot lines for each scheduler thread count
    scheduler_counts = sorted(summary['scheduler_threads'].unique())
    for i, sched_threads in enumerate(scheduler_counts):
        data = summary[summary['scheduler_threads'] == sched_threads]
        plt.plot(data['processing_threads'], data['mean'],
                    #yerr=data['std'],
                    #marker=markers[i],
                    color=okabe_ito_colors[i % len(okabe_ito_colors)],
                    linewidth=2,
                    #linestyle=line_styles[i],
                    markersize=6)

    plt.xlabel('Processing Threads (U)', fontsize=10)
    plt.ylabel('Runtime (seconds)', fontsize=10)
    #plt.title('Simulator Performance vs Thread Configuration', fontsize=10)
    plt.legend(("Schedule on main thread",  "1 scheduler thread", "2 scheduler threads",
                "4 scheduler threads", "8 scheduler threads",
                "12 scheduler threads", "16 scheduler threads",
                "20 scheduler threads"), fontsize=8)
    plt.grid(True, alpha=0.3)
    plt.xticks(range(0, 49, 8))

    # Set y-axis to start from 0 for better comparison
    plt.ylim(bottom=0)
    plt.tight_layout()

    if save_plot:
        plt.savefig('runs/scaling/simulator_performance.png', dpi=300)
        print("Plot saved as simulator_performance.png")


   # Calculate mean runtime for each configuration
    summary = df.groupby(['processing_threads', 'scheduler_threads'])['runtime_seconds'].mean().reset_index()

    # Get baseline runtime (1 processing thread, 1 scheduler thread)
    baseline = summary[(summary['processing_threads'] == 1) &
                      (summary['scheduler_threads'] == 0)]['runtime_seconds'].iloc[0]

    # Calculate speedup
    summary['speedup'] = baseline / summary['runtime_seconds']

    plt.figure(figsize=(3, 3))

    # Plot lines for each scheduler thread count
    scheduler_counts = sorted(summary['scheduler_threads'].unique())

    for i, sched_threads in enumerate(scheduler_counts):
        data = summary[summary['scheduler_threads'] == sched_threads]
        plt.plot(data['processing_threads'], data['speedup'],
                #marker=markers[i],
                linewidth=2,
                markersize=6,
                label=f'{sched_threads} scheduler threads',
                #linestyle=line_styles[i],
                color=okabe_ito_colors[i % len(okabe_ito_colors)])

    # Add ideal speedup reference line
    #plt.plot(range(1, 17), range(1, 17), '--', color='gray', alpha=0.5, label='Ideal speedup')

    plt.xlabel("Processing Threads")
    plt.ylabel('Speedup')
    plt.grid(True, alpha=0.3)
    plt.xticks(range(0, 49, 4))
    plt.xlim(0, 48)
    plt.ylim(bottom=0)

    plt.tight_layout()
    if save_plot:
        plt.savefig('runs/scaling/simulator_speedup.png', dpi=300)
        print("Plot saved as simulator_speedup.png")

    if show_plot:
        plt.show()


if __name__ == "__main__":
    # Configuration flags
    run_experiments = False
    plot_results_flag = True

    if run_experiments:
        print("Running performance experiments...")
        results_df = run_performance_experiments()

    # Load existing results
    results_df = pd.read_csv("runs/scaling/simulator_performance_results.csv")
    print("Loaded existing results from simulator_performance_results.csv")

    # Generate plots
    if plot_results_flag:
        print("Generating plots...")
        plot_results(results_df, show_plot=True, save_plot=True)

    print("Analysis complete!")
