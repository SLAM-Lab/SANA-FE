"""
Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Compress spike train data from SAnA-FE for a chosen layer to snntoolbox's format
"""
import csv

layer = '1'
spikes = []

with open("spikes.csv", "r") as csvfile:
    reader = csv.DictReader(csvfile)
    neurons = []
    timesteps = []
    for row in reader:
        l, neuron = row["neuron"].split(".")
        if l == layer:
            spikes.append((neuron, int(row["timestep"])))

# Sort the list based on the neuron id rather than timestep
spikes.sort(key=lambda x:x[1])

spike_neurons = [spike[0] for spike in spikes]
spike_times = [spike[1] for spike in spikes]

# Now print in the new format
with open("spiketrain.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(spike_neurons)
    writer.writerow(spike_times)

print("Finished converting spike train format.")
