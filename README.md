Performance Simulation for Neuromorphic Architectures

A framework to help with neuromorphic codesign providing performance (energy,
time, power) modeling.

# To Build

This project uses a Makefile based build.  To build use:
`make all`

For debug / profiling builds you can run
`make debug`

# To Run an Example

`./sim loihi.arch examples/random_network.net 20`

This simulates 20 time-steps of a small randomly connected spiking network.

General usage:

`./sim  <architecture description> <spiking network> <N timesteps>`

# Input Format

Right now the simulator executes a series of commands, like a custom program.
Both configuration and input is a series of commands to the simulator. A command
is a single line of input, in future I'll add an interactive option.

Each command does one specific thing e.g. define a neuron in an SNN, or a tile
in a chip, or steps the simulation by one timestep. The complete list is given
below. The order only matters that things are defined before they get used. E.g.
you define the network and architecture before you map the two.

Probably what you want to do is create a script to create the input commands and
then kick off a simulation. An entire experiment can be captured in a single
file or program.

Note that `parse_arch.py` converts between a YAML description of a neuromorphic
architecture to a set of architecture commands for the simulator.

`scripts/connected_layers_experiment.py` is an example of a complete flow -
creating an architecture, network and simulating a network for number of
timesteps.

# Commands

See `command.c` in case I didn't update this document.

The plan is to make all of these callable from a binding. Then the user never
has to bother with this low level stuff. In other words, commands will become an
equivalent call to some binding from C to another language e.g. Python

Right now though my hacky way of doing it is by using a single character
followed by a variable number of fields.

## (Spiking) network related
* g: Define neuron group
* n: Define neuron
* e: Add group of n inputs to the spiking network
* \<: Define external input

## Hardware related
* @: Add network on chip interconnect
* t: Define tile
* c: Define core
* &: Map neuron group to H/W
* i: Add axon input
* s: Add synapse processor
* +: Add soma processor
* o: Add axon output

## Other
* $: input vector, either spikes or spike rates
* \*: Step simulator
* l: Load commands from a file


# Output Format

Right now all outputs are hard fixed, either csv or yaml files. At some point
the csv files will be combined into one.

`probe_spikes.csv` The spikes for each time-step on probed neurons

`probe_potential.csv` The potentials for each time-step on probed neurons

`perf.csv` Detailed statistics for each timestep, perf breakdown

`stats.yaml` High level statistics for the simulation e.g. total runtime

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

This C code has been written to compile under the C99 standard. The style is
using a mix ideas from the Google guide and Linux style guide. Python code
is for Python 3.x and should follow PEP8. At some point I'll add linting for
both.

## Software Testing / Continuous Integration

The code has been run using Valgrind to verify there are no memory leaks.
After making any changes (especially with memory allocation) make sure no leaks
are introduced using:

`valgrind ./sim loihi.arch examples/random_network.net 1`

Or something similar. Also clang static code analysis has been used, and should
ideally stay error free. In a clean build directory run: `scan-build make`.

# Contact
James Boyle: james.boyle@utexas.edu
