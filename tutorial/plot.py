"""
Copyright (c) 2024 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

plot.py: Plot figures for NICE 2024 Tutorial
"""
import csv
import os
import sys
# matplotlib is required for plotting script, in addition to pyyaml and numpy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
NETWORK_PATH = os.path.join(PROJECT_DIR, "tutorial", "solutions",
                            "snn_completed.net")
ARCH_PATH = os.path.join(PROJECT_DIR, "tutorial", "solutions",
                         "arch_completed.yaml")
TIMESTEPS = 10

sys.path.insert(0, PROJECT_DIR)
import sim

results = sim.run(ARCH_PATH, NETWORK_PATH, TIMESTEPS,
                  out_dir=os.path.join(PROJECT_DIR, "tutorial"),
                  spike_trace=True, potential_trace=True)

plt.rcParams.update({'font.size': 14, 'axes.linewidth': 1})

## Create the spike raster plot
timesteps = []
spikes = []
with open(os.path.join(PROJECT_DIR, "tutorial", "spikes.csv")) as spikes_file:
    spike_reader = csv.DictReader(spikes_file)
    for line in spike_reader:
        (nid, timestep) = line["neuron"], line["timestep"]
        group, neuron = nid.split(".")
        group, neuron, timestep = int(group), int(neuron), int(timestep)
        timesteps.append(timestep)
        spikes.append(neuron)

fig, ax = plt.subplots(2, figsize=(5, 4))
colors = matplotlib.colors.ListedColormap(("#ff7f0e", "#1f77b4"))
ax[0].scatter(timesteps, spikes, marker="|", s=700, linewidths=3, c=spikes, cmap=colors)
plt.xlabel("Time-step")

plt.minorticks_on()
ax[0].yaxis.set_tick_params(which='minor', bottom=False)
ax[0].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))

ax[0].set_ylim((-0.5, 1.5))
ax[0].set_yticks((0, 1))
ax[0].set_yticklabels(("Neuron 0.0", "Neuron 0.1"))

voltages = np.zeros((11, 2))
with open(os.path.join(PROJECT_DIR, "tutorial", "potential.csv")) as v_file:
    v_reader = csv.DictReader(v_file)
    for line in v_reader:
        (timestep, v1, v2) = (int(line["timestep"]), float(line["neuron 0.0"]),
                              float(line["neuron 0.1"]))
        voltages[timestep, 0] = v1
        voltages[timestep, 1] = v2

#plt.figure(figsize=(6, 2.5))
#plt.title("Voltages Recorded")
ax[1].plot(voltages[:, 1], '--^')
ax[1].plot(voltages[:, 0], '-o')
ax[1].set_ylabel("Neuron Potential (V)")
plt.legend(("Neuron 0.1", "Neuron 0.0"))

ax[1].yaxis.set_tick_params(which='minor', bottom=False)
ax[1].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))

ax[0].set_xlim((0, 10.1))
ax[1].set_xlim((0, 10.1))
ax[0].set_xticks((0, 2, 4, 6, 8, 10))
ax[1].set_xticks((0, 2, 4, 6, 8, 10))

plt.tight_layout()
plt.savefig("tutorial/neuron_traces.pdf")
plt.savefig("tutorial/neuron_traces.png")

fig, ax = plt.subplots(2, figsize=(5, 4))
fired = np.zeros((11))
energy = np.zeros((11))
with open(os.path.join(PROJECT_DIR, "tutorial", "perf.csv")) as perf_file:
    reader = csv.DictReader(perf_file)
    for line in reader:
        (timestep, f, e) = (int(line["timestep"]), int(line["fired"]),
                              float(line["total_energy"]))
        fired[timestep] = f
        energy[timestep] = e

#plt.figure(figsize=(6, 2.5))
#plt.title("Voltages Recorded")
ax[0].plot(fired, '-o')
ax[0].set_ylabel("Neurons Fired")
ax[1].plot(energy * 1.0e12, '-o')
ax[1].set_ylabel("Total Energy (pJ)")

plt.minorticks_on()
#ax[1].yaxis.set_tick_params(which='minor', bottom=False)
ax[1].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))

ax[0].set_xlim((0, 10.1))
ax[1].set_xlim((0, 10.1))
ax[0].set_xticks((0, 2, 4, 6, 8, 10))
ax[1].set_xticks((0, 2, 4, 6, 8, 10))
ax[1].set_xlabel("Time-step")

plt.tight_layout()
plt.savefig("tutorial/hw_traces.pdf")
plt.savefig("tutorial/hw_traces.png")


plt.show()