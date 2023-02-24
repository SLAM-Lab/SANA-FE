import matplotlib
matplotlib.use('Agg')

import csv
import subprocess
import yaml
from matplotlib import pyplot as plt

neurons = []
energies = []
with open("../runs/neuron_characterization.csv", "r") as nonspiking_csv:
    reader = csv.DictReader(nonspiking_csv)
    for row in reader:
        neurons.append(float(row["Neurons"]))
        energies.append(float(row["Energy (nJ)"]))

plt.rcParams.update({'axes.linewidth': 1.5, "scatter.marker": 'x'})
plt.rcParams.update({'font.size': 18, 'lines.markersize': 6})

plt.figure(figsize=(4.5, 4.5))
plt.scatter(neurons, energies, marker='x', s=30, lw=2)
plt.yscale("linear")
plt.xscale("linear")
plt.ylabel("Energy (nJ)")
plt.xlabel("Neuron Updates")
plt.ticklabel_format(style="sci", axis="x", scilimits=(0,0))
plt.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
plt.tight_layout()
plt.savefig("../runs/neuron_characterization.pdf")
plt.savefig("../runs/neuron_characterization.png")
