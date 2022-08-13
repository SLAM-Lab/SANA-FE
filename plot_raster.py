import matplotlib.pyplot as plt
import numpy as np
import csv


plt.figure()
with open("probe_spikes.csv") as spike_csv:
    spike_data = csv.DictReader(spike_csv)
    timesteps = 0
    neuron_ids = spike_data.fieldnames
    plt.ylim((0, len(neuron_ids)))

    for neuron_spikes in spike_data:
        values = list(neuron_spikes.values())
        values = values[0:-1]  # Trim empty field
        s = [int(v) for v in values] 
        spike_array = np.asarray(s)
        spikes = np.where(spike_array >= 1)[0] 
        
        plt.scatter([timesteps] * len(spikes), spikes.tolist(), c='b', s=3,
        marker='.', linewidths=0.1)
        timesteps += 1

plt.xlabel("Time-step")
plt.ylabel("Neuron")

plt.show()

