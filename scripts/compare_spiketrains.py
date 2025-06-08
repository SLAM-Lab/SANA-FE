"""
Copyright (c) 2024 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.
"""
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("first")
parser.add_argument("second")

args = parser.parse_args()
first_spiketrain = args.first
second_spiketrain = args.second


with open(first_spiketrain, "r") as first:
    data = csv.reader(first)
    first_neurons = next(data)
    first_timesteps = next(data)

with open(second_spiketrain, "r") as second:
    data = csv.reader(second)
    second_neurons = next(data)
    second_timesteps = next(data)

same = True
for (n1, t1, n2, t2) in zip(first_neurons, first_timesteps, second_neurons, second_timesteps):
    if n1 != n2:
        print("Neuron {0} != {1}".format(n1, n2))
        same = False
        #break
    if t1 != t2:
        print("Timestep {0} != {1} (nid:{2})".format(t1, t2, n2))
        same = False
        #break

if same:
    print("Spike trains are the same")
