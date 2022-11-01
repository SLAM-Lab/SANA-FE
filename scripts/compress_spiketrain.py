"""Compress spike train data for a chosen layer, same as snntoolbox format"""
import csv

#layer = 'i'
layer = 0
spikes = []

with open("probe_spikes.csv", "r") as csvfile:
    reader = csv.DictReader(csvfile)
    neuron_ids = reader.fieldnames
    layer_neurons = [n for n in neuron_ids if "{0}.".format(layer) in n]
    print("{0} inputs / neurons found in network".format(len(layer_neurons)))
    timestep = 1

    for row in reader:
        for n in layer_neurons:
            spiked = int(row[n])
            if spiked != 0:
                neuron_id = n.split('.')[1]
                spikes.append((timestep, int(neuron_id)))

        timestep += 1

# Sort the list based on the neuron id rather than timestep
spikes.sort(key=lambda x:x[1])

spike_times = [spike[0] for spike in spikes]
spike_neurons = [spike[1] for spike in spikes]

# Now print in the new format
with open("spiketrain.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(spike_neurons)
    writer.writerow(spike_times)

print("Finished converting spike train format.")
