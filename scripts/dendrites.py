import snntorch as snn

from snntorch import spikeplot as splt
from snntorch import spikegen
import torch
import torch.nn as nn

import numpy as np
import matplotlib.pyplot as plt

network = snn.torch.load("runs/dendrites/dend.pt")
spike_recording = []

num_steps = 100
batch_size = 1
data_in = torch.rand(num_steps, batch_size, 1, 28, 28)

#for step in range(num_steps):
#    spike, state = network(data_in[step])
#    spike_recording.append(spike)

for param in network.parameters():
    print(param)