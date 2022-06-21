Performance Simulation for Neuromorphic Architectures

A framework to help with neuromorphic codesign providing performance (energy,
time, power) modeling.

# To Build
This project uses a Makefile based build.  To build use:
`make all`

Make sure that GCC is supported and working on your system.

# To Run an Example
`./sim loihi.tech 20 examples/random_network.csv`
This simulates 20 time-steps of a small randomly connected network neurons.

General usage:
`./sim <tech file> <N timesteps> <neuron list.csv>`

# Input Format
The simulator uses csv (comma separated values) to define the
spiking network (for now).  The simulator takes the number of timesteps to
simulate and a csv file containing the network parameters.  Each line of the
csv defines a neuron.  The line starts with neuron parameters then followed by a
variable number of synapse parameters (synapse 1, synapse 2, synapse 3, ...)

The csv files can optionally have a header row.  To see an example see
`examples/random_network.csv` 

One slightly awkward thing is the columns must exactly match what the simulator
is expecting; to see the column formats for both neurons / synapses see:

# Project Code
This project has been written in C.

main.c
sim.c
network.c
tech.c

# Contact
James Boyle: james.boyle@utexas.edu
