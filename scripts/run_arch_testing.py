#!/usr/bin/env python3
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import sanafe
from sanafe.data.traces import TraceData
from sanafe.viz.raster import raster_plot
from sanafe.viz.potential import potential_plot, potential_heatmap, potential_subplots
from sanafe.viz.performance import (
    energy_breakdown_plot, throughput_plot,
    latency_histogram, latency_comparison,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = REPO / "tmp_arch_testing"
OUT.mkdir(exist_ok=True)

#1. Load architecture 
print("Loading architecture...")
arch = sanafe.load_arch(str(REPO / "sanafe" / "data" / "arch_testing.yaml"))

#2. Build SNN programmatically 
snn = sanafe.Network()

# bias 
inp = snn.create_neuron_group("input", 4,
    model_attributes={"threshold": 1.0, "reset": 0.0},
    log_spikes=True, log_potential=True)

# 4 neurons
hid = snn.create_neuron_group("hidden", 4,
    model_attributes={"threshold": 0.8, "reset": 0.0},
    log_spikes=True, log_potential=True)

# 2 neurons
out = snn.create_neuron_group("output", 2,
    model_attributes={"threshold": 1.5, "reset": 0.0},
    log_spikes=True, log_potential=True)

for i_n in inp:
    for h_n in hid:
        i_n.connect_to_neuron(h_n, {"weight": 0.5})

for h_n in hid:
    for o_n in out:
        h_n.connect_to_neuron(o_n, {"weight": 0.7})

inp[0].set_attributes(model_attributes={"bias": 1.2})
inp[1].set_attributes(model_attributes={"bias": 0.9})
inp[2].set_attributes(model_attributes={"bias": 1.1})
inp[3].set_attributes(model_attributes={"bias": 0.6})

# mapped neurons to cores on different tiles 
for i, n in enumerate(inp):
    n.map_to_core(arch.tiles[0].cores[i])
for i, n in enumerate(hid):
    n.map_to_core(arch.tiles[1].cores[i])
for i, n in enumerate(out):
    n.map_to_core(arch.tiles[2].cores[i])

# 3. Simulate 
print("Creating SpikingChip...")
chip = sanafe.SpikingChip(arch)
print("Loading SNN...")
chip.load(snn)
print("Running simulation (100 timesteps)...")

results = chip.sim(
    100,
    spike_trace=True,
    potential_trace=True,
    perf_trace=True,
    message_trace=True,
)

print(f"  Total energy: {results['energy']['total']:.2e} J")
print(f"  Neurons fired: {results['neurons_fired']}")
print(f"  Spikes: {results['spikes']}")

#  4. Build TraceData directly from in-memory results 
tr = TraceData.from_sim_results(results)
print(f"  TraceData: {tr}")

#  5. Export CSVs using proper serialization 
spike_csv = OUT / "spikes.csv"
tr.spikes_to_dataframe().to_csv(spike_csv, index=False)
print(f"Saved {spike_csv}")

pot_csv = OUT / "potentials.csv"
tr.potentials_to_dataframe().to_csv(pot_csv, index=True)
print(f"Saved {pot_csv}")

#  6. Generate all visualizations 
print("\nGenerating plots...")

fig, ax = raster_plot(tr, title="arch_testing — Spike Raster")
fig.savefig(OUT / "raster.png", dpi=150)
print("  Saved raster.png")

fig, ax = potential_plot(tr, title="arch_testing — Membrane Potentials")
fig.savefig(OUT / "potentials.png", dpi=150)
print("  Saved potentials.png")

fig, ax = potential_heatmap(tr, title="arch_testing — Potential Heatmap")
fig.savefig(OUT / "potential_heatmap.png", dpi=150)
print("  Saved potential_heatmap.png")

fig, axes = potential_subplots(tr, ncols=2)
fig.savefig(OUT / "potential_subplots.png", dpi=150)
print("  Saved potential_subplots.png")

fig, ax = energy_breakdown_plot(tr, title="arch_testing — Energy Breakdown")
fig.savefig(OUT / "energy_breakdown.png", dpi=150)
print("  Saved energy_breakdown.png")

fig, ax = energy_breakdown_plot(tr, mode="stacked_bar", normalize=True,
                                title="arch_testing — Energy % Breakdown")
fig.savefig(OUT / "energy_normalized.png", dpi=150)
print("  Saved energy_normalized.png")

fig, ax = throughput_plot(tr, metrics=["fired", "spikes", "hops"],
                         title="arch_testing — Throughput")
fig.savefig(OUT / "throughput.png", dpi=150)
print("  Saved throughput.png")

fig, ax = latency_histogram(tr, metric="generation_delay",
                            title="arch_testing — Generation Delay")
fig.savefig(OUT / "latency_gen.png", dpi=150)
print("  Saved latency_gen.png")

fig, ax = latency_comparison(tr,
                             metrics=["generation_delay", "receive_delay",
                                      "network_delay", "blocked_delay"],
                             title="arch_testing — Latency Comparison")
fig.savefig(OUT / "latency_comparison.png", dpi=150)
print("  Saved latency_comparison.png")

plt.close("all")
print(f"\nOutputs saved to {OUT}/")
