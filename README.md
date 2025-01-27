Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Simulating Advanced Neuromorphic Architectures for Fast Exploration (SANA-FE)

A framework to model energy and performance of neuromorphic hardware.

# To Build

This project uses CMake as its build system and dependency manager. To setup
compilation, run the following command in the project directory:
`cmake .`

Then compile SANA-FE by using the command:
`make`

Optionally, enable multithreaded builds using the optional flag `-j`

## Dependencies

Building this project requires `cmake`, `make`, and a compiler that supports the
C++17 standard (e.g., GCC >= 8, Clang >= 5). This project makes use of the
opensource project RapidYAML for all YAML file parsing. To build with Python
interfaces (currently mandatory), you must also have Python 3.x installed with
PyBind11. Assuming Python3 is available, you can install PyBind11 using:

`pip install pybind11`

# To Run an Example

`./sim arch/example.yaml snn/example.yaml 10`

This simulates 10 time-steps of a tiny connected spiking neural network (SNN).

General usage:

`./sim <architecture description> <SNN description> <N timesteps>`

In addition to the standlone simulator, SANA-FE can also be scripted using a
Python API. For an example of how this can be done, see the Jupyter
notebook-based tutorial: `tutorial/tutorial.ipynb`.

Additional examples and experiments may be found in the `scripts/` directory.

# Simulator Inputs

SANA-FE takes command line arguments, an architecture description file (YAML)
and an SNN description file (YAML). The description files both use custom
file formats. Examples for architectures may be found in `arch/`. Examples for
SNNs may be found in `snn/`.

Optional command line flags can be used to enable simulation traces. Note that
after enabling traces globally, you will still have to create probes at the
neuron level to get trace output.

Flags:
* `-v`: Enable potential (v) traces to `potential.csv`
* `-s`: Enable spike traces to `spikes.csv`
* `-p`: Record the simulated performance of each timestep to `perf.csv`
* `-m`: Enable message traces to `messages.csv`
* `-n`: Use the (legacy) netlist format for SNNs, instead of YAML.

## SNN Description

The SNN description format is based on the YAML file format.

Different mapped SNNs can be defined flexibly and generally using sections for
neuron `groups`, `edges`, and hardware `mappings`. Each section allows for
custom `attributes` to be defined, which are converted to model parameters
within the simulator. While the keywords for sections are fixed, attributes
allow for custom user-defined parameters to be associated with neurons and
connections.

The SNN must be defined under the main `network` section. All other top-level
sections are ignored. Then, we have `groups`, `edges`, and `mappings`
sub-sections.

Groups of neurons are one or more neurons that may share some common attributes.
This is similar to how other frameworks may define populations, or layers of
similar neurons. How neurons are grouped is up to the user, but they can be
useful for sharing common attributes or connections. Under the `groups`
subsection, you must create a list of named neuron groups. Within each group is
an `attributes` section and a `neurons` section.

In each `neurons` subsection, list all sets of neurons belonging to the group as
an ID:attributes pair. For conciseness it is possible define ranges of neurons
at once using the range (..) notation. Following each neuron, give an ordered or
unordered list of attributes e.g.,

- 0..2: [attribute1: value1]
- 3: {attribute1: value1}

In the `edges` section, define neuron to neuron connections or group to group
hyper-edges, including any edge attributes. The edge format uses a notation
similar to the graph DOT format e.g.,

- layer1.0 -> layer2.1: [weight: 1]
- layer1 -> layer2: [weight: 1]

Finally, in the `mappings` section we map neurons to hardware cores. Under the
section heading is a list of mappings, with an example of one mapping as
follows:

- layer1.0..1: [core: 0.0]

Similar to before, neurons may give as a range for brevity. This maps two
neurons to tile 0, core 0 (the first core in the chip).

As long as valid YAML syntax is used, SANA-FE does not distinguish between
different styles (block style, flow style, or a mix of the two). For one example
of a simple SNN, see `snn/example.yaml`.

## Architecture Description

The architecture description format is also based on the YAML file format.

Different architectures are defined using a hierarchical description.
This tool models neuromorphic designs with several assumptions, in order to
simplify the tool.

1) The chip is time-step based. A time-step is a small discrete amount of time.
    This is as opposed to a purely event driven simulation e.g. ROSS.
2) The neural cores adhere to some common design patterns

At the top level, the description begins with the "architecture" keyword. Any
other top-level sections will be ignored. This defines anything at the chip
level, including the NoC interconnect.

A chip contains one or more network tiles, representing some shared network
resources e.g., a router. Each `tile` contains one or more cores, where a core
performs computation. Each neuromorphic `core` contains a fixed spike processing
hardware pipeline. It is assumed that tiles and cores are all parallel processing
elements.

Each core is assumed to have a neuromorphic pipeline which processes the updates
for one or more neurons. The pipeline is a fixed sequence of niche hardware
units. Those hardware units could contain digital logic, analog circuits or
even novel devices.

The pipeline contains the following units:

* The input axon unit receive spike packets from the network and generate
  synaptic addresses for memory lookups.

* The synaptic unit looks up connectivity for incoming spikes and updates any
  relevant synaptic currents.

* The dendritic unit combines currents based some internal structure and a set
  of operations.

* The soma unit updates membrane potentials based on the dendritic current and
  neuron model. If the firing criteria is met, it generates a spike for that
  neuron.

* The output axon unit send spikes from the soma out to the network to go to
  other cores' pipelines.

For an example, see `arch/loihi.yaml`. There are a nested series of
keywords, where keywords define required hardware blocks. Each block must be
contain a `name` keyword, which may optionally specify the number of instances.
Blocks are duplicated the number specified in the range, for example:

    # Define 8 cores, 0 through 7
    -name: neuromorphic_core[0..7]

Blocks much also have both an `attributes` section and the next hardware blocks
in the hierarchy. The attributes section will generate one or more parameters
that are passed to be parsed by the simulator and the relevant hardware models
implemented either internally (`models.cpp`) or externally (plugins).

# Simulator Outputs

If corresponding traces are enabled, output is saved to trace files with
hard-coded names using either csv or yaml extensions.

`spikes.csv`: The spikes for each time-step on probed neurons

`potential.csv`: The potentials for each time-step on probed neurons

`perf.csv`: Detailed statistics for each timestep and each hardware unit

`messages.csv`: Information on spike messages for each time-step

`run_summary.yaml`: High-level statistics for the simulation e.g. runtime

# Simulator Kernel

SANA-FE uses a user-provided spiking architecture, a mapped SNN, and run-time
configuration to simulate a spiking chip as it executes a spiking application.
SANA-FE uses the `Architecture` to compile a `SpikingChip`, which it then loads
the mapped SNN. SANA-FE then rapidly simulates the design at a time-step
granularity.

During each time-step SANA-FE models custom spike-processing pipelines executing
within each core, modeling the processing of neurons and spike messages. Using
our spiking hardware template, we enable custom hardware blocks to be
incorporated for axonal, synaptic, dendritic and somatic hardware. Each hardware
unit is implemented using a `model` - you can take the built-in hardware unit
models provided in `models.cpp`, or implement models externally as hardware unit
plugins using the fixed base class interfaces. The SANA-FE kernel coordinates
all on-chip activity, makes calls to the models and tracks the total energy and
latency across the chip.

SANA-FE includes efficient but detailed semi-analytical timing models. This
takes aggregated information about all spike messages generated in a time-step
and calls a custom scheduler in `schedule.cpp`. The on chip schedule ultimately
gives you a reasonably accurate prediction of the chip timings, accounting for
effects such as blocking in the NoC and custom latency simulations within
hardware units.

# Plugins

As part of SANA-FE, the user can implement different hardware models using
custom plugins. Models for synapses, dendrites and somas are all supported.
The plugin mechanism defines common interfaces for you to create .so files.
There are some example plugins provided in the `/plugins` folder, including an
implementation of a Hodgkin-Huxley neuron.

To include your plugin, several steps need to be taken:
1. The path to the plugin needs to be specified in the architecture yaml file,
inside a corresponding `synapse`, `dendrite` or `soma` hardware unit.
You must provide the plugin path using the attribute `plugin: <name>`. You must
also specify the name of the model inside the plugin using the attribute
`model: <name>`.
2. The hardware unit is mapped to in the net file with the keyword
`soma_hw_name`.
3. The `/plugins`folder contains the plugin compiled into a file with the
format `[plugin name].so`.

Each plugin must implement a class derived from the base classes `SynapseUnit`,
`DendriteUnit` or `SomaUnit` classes. On every creation of a neuron or
connection, a new instance of the specific plugin class will be created.

## Creating a New Plugin

SANA-FE's plugin model makes it easy to create new models and
integrate them directly into your SNN execution. There are, however,
two main interfacing requirements to allow plugins to smoothly execute.

1. The plugin class must extend the `BaseSoma` located in `plugins.hpp`.
This includes implementing the base constructor, `update()`, `reset()` and
`set_attribute()` functions.
2. One class factory functions needs to be created. This has to be in the
format `create_<plugin_name>`. This function is used to get an instance of your
new class. The objects are automatically destructed.

It is recommended new users read through the hodgkin_huxley.cpp file to see what
an example soma class looks like. The `update_soma` function will be passed
the input current spike as a double and returns `NeuronStatus`, an enum
representing whether the neuron is idle, updated, or fired. The `parameters`
function takes in a struct of `attributes` with a specified int length.
These parameters are arbitrary keyword=value pairs that can pass in any
information to the class. This is where parameters specified in the network
file will be passed to.

The .so should be kept, the .o file is unnecessary and can be deleted.

## Legacy SNN Description Format (netlist)

Version 1 of SANA-FE (written in C) defined a simpler, less capable SNN
description format (compared to the current YAML-based format). For
back-compatability, the netlist-style format is still supported. To use this
format, use the command-line flag (-n).

In the netlist format each line defines a new entry, which may either be a
neuron group (g), neuron (n), edge (e), or mapping (&). Each line starts with
the type of entry followed by one required field and then any number of named
attributes. Fields are separated by one or more spaces.

Attributes are defined using the syntax: `<attribute>=<value>`. Note, there
is no space before or after the equals. The attribute `soma_hw_name` is required
to be set for every neuron or neuron group.

A neuron group helps reduce the number of repeated, shared parameters for a
population of neurons e.g., for a layer of neurons in a deep SNN.

`g <number of neurons> <common attributes>`

Neurons are addressed (using the group number followed by neuron number), and
then all attributes are specified. Note the group must be defined first.

`n group_id.neuron_id <unique attributes>`

An edge connects one source neuron (presynaptic) to one destination neuron
(postsynaptic). The edge may also have attributes such as synaptic weight.

`e src_group_id.src_neuron_id->dest_group_id.dest_neuron_id <edge attributes>`

Finally, mappings place predefined neurons on a hardware core. Here we specify
the neuron and the core.

`& group_id.neuron_id@tile_id.core_id`

An example of how to use the netlist format is given in `snn/example.net`

# Project Code

This project has been written in C++ and Python. All code is in the `/src`
directory. See header files for more detail on supported classes and functions.

`arch.cpp`
`chip.cpp`
`description.cpp`
`main.cpp`
`network.cpp`
`param.hpp`
`plugins.cpp`
`print.cpp`
`pymodule.cpp`
`schedule.cpp`

C++ code has been written using the C++17 standard.

# Contact
James Boyle: james.boyle@utexas.edu
