Simulating Performance of Neuromorphic Architectures

A framework to help with neuromorphic codesign, modeling performance (energy,
time, power).

# To Build

This project uses a Makefile based build.  To build use:
`make all`

For debug / profiling builds you can run
`make debug`

# To Run an Example

`./sim loihi.arch examples/random_network.net 20`

This simulates 20 time-steps of a tiny connected spiking network.

General usage:

`./sim  <architecture description> <neural network> <N timesteps>`

# Input Format

The simulator executes a series of commands.
The network, architecture and inputs are specified using a common command
interface.
Each command is a single line of plain-text, with space separated
values.

A command does one specific thing e.g. define a neuron in an SNN, or a tile
in a chip, or steps the simulation by one timestep.

`parse_arch.py` converts between a YAML description of a neuromorphic
architecture to a set of architecture commands for the simulator.

`scripts/connected_layers_experiment.py` is an example of a complete flow -
creating an architecture, network and simulating a network for number of
timesteps.

## YAML Architecture Description

Different architectures are defined using a hierarchical description using YAML.
This tool models neuromorphic designs with several assumptions, in order to
simplify the tool.

1) The chip is time-step based. A time-step is a small discrete amount of time.
    This is as opposed to a purely event driven simulation e.g. ROSS.
2) The neural cores follow some common design patterns

The top level is always the "architecture". This defines anything at the chip
level, including the NoC interconnect.
A chip contains one or more tiles, which with the interconnect form the NoC.
Each tile contains one or more cores, where a core performs computation.
Each neuromorphic core contains a subset of certain operations in a hardware
pipeline.
It is assumed that tiles and cores are all parallel processing elements.

Neuromorphic cores are assumed to share some common design principles.
Here we define the pipeline used for all cores in this simulator.
A pipeline here is a series of niche processors, which can each handle a part of
the computation of one or more neurons.
Those processors could be digital logic, analog circuits or a novel device.
Where multiple processors are defined for the same stage, we assume
these computations can be done in parallel.

Input axons receive spike packets from the network.
The synaptic processor looks up connectivity for incoming spikes and updates
the relevant synaptic currents.
The dendritic processor combines currents based on a tree structure and a set
of operations.
The soma processor updates membrane potentials based on dendritic currents
If the firing criteria is met, it generates a spike for that neuron.
The output axons send spikes from the soma out to the network.

# Commands

See `command.c`. Each command starts with a single character, followed by
a list of parameters.

## (Spiking) Network
* g: Add neuron group (neurons with the same properties e.g. a layer)
* n: Define neuron
* e: Add group of n external inputs to the spiking network
* \<: Define external input

## Hardware / Architecture
* @: Add network on chip interconnect
* t: Define tile
* c: Define core
* &: Map neuron to H/W units
* i: Add axon input
* s: Add synapse processor
* +: Add soma processor
* o: Add axon output

## SImulator control
* $: input vector, either spikes or spike rates
* \*: Step simulator
* l: Load commands from a file

# Output Format

Right now all outputs are hard fixed, either csv or yaml files. The plan is to
combine all data into perf.csv

`probe_spikes.csv` The spikes for each time-step on probed neurons

`probe_potential.csv` The potentials for each time-step on probed neurons

`perf.csv` Detailed statistics for each timestep and each hardware unit

`stats.yaml` High level statistics for the simulation e.g. runtime

# Project Code

This project has been written in C and Python. See header files for more
description.

`main.c`
`sim.c`
`network.c`
`arch.c`
`command.c`
`parse_arch.py`
`run.py`
`spikeperf.py`

This C code has been written based on the C99 standard.

# Contact
James Boyle: james.boyle@utexas.edu
