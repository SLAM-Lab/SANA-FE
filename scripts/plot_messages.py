import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd

MESSAGE_TRACE_FILENAME = "dvs_messages.trace"

def tile_idx(x, y):
    return (x * 4) + y


def track_hops(src_x, src_y, dest_x, dest_y):
    # Loihi is 8x4 grid (8 on the x axis)
    spikes_processed_per_router = np.zeros(8*4, dtype=int)

    # Add spike to tile that generates the spike
    spikes_processed_per_router[tile_idx(src_x, src_y)] += 1

    # Account for all the hops in between the src and dest tile
    while src_x < dest_x:
        src_x += 1
        spikes_processed_per_router[tile_idx(src_x, src_y)] += 1

    while src_x > dest_x:
        src_x -= 1
        spikes_processed_per_router[tile_idx(src_x, src_y)] += 1

    while src_y < dest_y:
        src_y += 1
        spikes_processed_per_router[tile_idx(src_x, src_y)] += 1

    while src_y > dest_y:
        src_y -= 1
        spikes_processed_per_router[tile_idx(src_x, src_y)] += 1

    return spikes_processed_per_router


hops = np.zeros(32, dtype=int)
with open(MESSAGE_TRACE_FILENAME) as trace:
    reader = csv.DictReader(trace)
    for line in reader:
        src_hw = line["src_hw"]
        dest_hw = line["dest_hw"]

        src_tile, src_core = src_hw.split('.')
        dest_tile, dest_core = dest_hw.split('.')

        src_x = int(src_tile) // 4
        src_y = int(src_tile) % 4

        dest_x = int(dest_tile) // 4
        dest_y = int(dest_tile) // 4
        hops += track_hops(src_x, src_y,
                           dest_x, dest_y)

plt.figure()

# Plot a heat map
x, y = np.meshgrid(np.linspace(0, 7, 8), np.linspace(0, 3, 4))
# Plot the circles

plt.scatter(x, y, c=hops, cmap="hot")
plt.colorbar()

df = pd.read_csv(MESSAGE_TRACE_FILENAME)
#df = pd.read_csv("latin_messages.trace")
#df = pd.read_csv("messages.trace")
plt.figure()
plt.hist(df["hops"], bins=50)

plt.figure()
plt.hist(df["generation_delay"], bins=50)

plt.figure()
plt.hist(df["processing_latency"], bins=50)

plt.show()
