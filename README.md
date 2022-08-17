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
This simulates 20 time-steps of a small randomly connected network of neurons.

General usage:
`./sim  <architecture description> <spiking network> <N timesteps>`

TODO: in future it will just take whatever files you want to give it
TODO: add `--steps <N timesteps>` and `--output` options

# Input Format
For now the simulator runs a series of commands - all input is formed as
commands. A command is a single line of input, for now given as a bunch of
files, in future I'll add an interactive option. Each command does one specific
thing e.g. define a neuron in an SNN, or a tile in a chip, or steps the
simulation by one timestep. The complete list is given below. The order only
matters that things are defined before they get used. E.g. you define the
network and architecture before you map the two.

Probably what you want to do is create a script to create the input commands and
then kick off a simulation.

Note that `parse_arch.py` converts between a YAML description of a neuromorphic
architecture to a set of architecture commands for the simulator.

TODO: Create requirements for a venv / conda environment
TODO: Create bindings for other languages e.g. Python

# Output Format

Right now all outputs are hard fixed, either csv or yaml files.
`probe_spikes.csv`: The spikes for each time-step on probed neurons
`probe_potential.csv`: The potentials for each time-step on probed neurons
`stats.yaml`: High level statistics for the simulation e.g. total runtime
`perf.csv`: Detailed statistics for each timestep, perf breakdown

# Project Code
This project has been written in C and Python.

`main.c`
`sim.c`: Simulator kernel, simulates all models and estimates performance
`network.c`: (Spiking) network related functionality
`arch.c`: Architecture i.e., H/W design related functionality
`command.c`: Implements all commands as set of API function calls

`parse_arch.py`: Parse architecture description YAML file as set of commands

# Commands
See `command.c` in case I didn't update this

TODO: define the syntax for each one
The plan is to make this callable from a binding so the user never has
to bother with this low level stuff. I.e. commands will become a call to
some binding in another scripting language

Network related
* g: Define neuron group
* n: Define neuron
Hardware related
* @: Add network on chip interconnect
* t: Define tile
* c: Define core
* &: Map neuron group to H/W
* i: Add axon input
* s: Add synapse processor
* +: Add soma processor
* o: Add axon output
TODO - Other
* external input node
* Input spikes or spike rates
* Step simulator
* Load commands from a file

# Contact
James Boyle: james.boyle@utexas.edu
