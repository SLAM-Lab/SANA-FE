<p align="center">
<img src="https://raw.githubusercontent.com/SLAM-Lab/SANA-FE/main/sana_fe_logo.svg" alt="SANA-FE" width="400" style="max-width: 100%; height: auto;">
</p>

Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

Simulating Advanced Neuromorphic Architectures for Fast Exploration (SANA-FE)

A framework for modeling the energy usage and performance of different
neuromorphic hardware.

# Citation

We hope that you find this project useful. If you use SANA-FE in your work,
please cite our paper:
J. A. Boyle, M. Plagge, S. G. Cardwell, F. S. Chance and A. Gerstlauer,
"SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration,"
in IEEE Transactions on Computer-Aided Design of Integrated Circuits and
Systems, doi: 10.1109/TCAD.2025.3537971.

    @article{boyle2025sanafe,
      title={SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration},
      author={Boyle, James A and Plagge, Mark and Cardwell, Suma George and Chance, Frances S and Gerstlauer, Andreas},
      journal={IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems},
      year={2025},
    }

# To Build

This project uses CMake as its build system and dependency manager. To setup
compilation, run the following command in the project directory:
`cmake .`

Then compile SANA-FE by using the command:
`make`

Optionally, enable multithreaded builds using the optional flag `-j <nthreads>`

## Dependencies

Building this project requires `cmake`, `make`, and a compiler that supports the
C++17 standard (e.g., GCC >= 8, Clang >= 5). This project uses RapidYAML for all
YAML file parsing, and Booksim 2 for optional cycle-accurate NoC modeling. To
build the Python interfaces, you must also have Python >= 3.8 installed with
PyBind11. You can install PyBind11 using:

`pip install pybind11`

Booksim 2 requires the `bison` and `flex` packages for config parsing.
For example, in Ubuntu these can be installed using `apt`.

`apt install bison flex`

# To Run an Example

`./sim arch/example.yaml snn/example.yaml 100`

This simulates 100 time-steps of a tiny connected spiking neural network (SNN).

General usage:

`./sim [optional flags] <architecture description> <SNN description> <N timesteps>`

In addition to the standlone simulator, SANA-FE can also be scripted using a
Python API. For an example of how this can be done, see the Jupyter
notebook-based tutorials in the `tutorial/` directory.

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
* `-m`: Enable message traces to `messages.csv`
* `-n`: Use the (legacy) netlist format for SNNs, instead of YAML.
* `-o`: Output directory
* `-p`: Record the simulated performance of each timestep to `perf.csv`
* `-s`: Enable spike traces to `spikes.csv`
* `-t [simple/detailed/cycle]`: Specify the timing model  (default=`detailed`)
* `-v`: Enable potential (voltage) traces to `potential.csv`
* `-N`: Number of neuron/message processing threads (default=1)
* `-S`: Number of scheduling threads (default=0, use main thread)

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

In each `neurons` subsection, list all sets of neurons belonging to the group.
For conciseness we support specifying multiple neurons using the range (..)
notation. Following each neuron, give an ordered or unordered list of
attributes e.g.,

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
SANA-FE supports a base hardware model base class, with which it implements
all of its synaptic, dendritic and somatic hardware models. Using SANA-FE's
`PipelineUnit` base class, you can implement your own models as hardware
plugins.

## Using Plugins

There is one example already provided in the `/plugins` folder implementing a
Hodgkin-Huxley neuron model (`hodgkin_huxley.cpp`). There are a few steps
required to use plugins in SANA-FE:

1. Specify the plugin path in the architecture yaml file, in the corresponding
`synapse`, `dendrite` or `soma` hardware section. Specify the plugin path using
the attribute `plugin: <pathname>`.
2. Specify the model name using the attribute `model: <name>`.
3. Map neurons to the hardware unit as usual with the attribute: `soma_hw_name`.

For example, for the Hodgkin-Huxley example provided with SANA-FE, you could use
it as follows:

    # Rest of arch description
    ...
    soma:
    - name: plugin_example_soma
      attributes:
        plugin: plugins/hodgkin_huxley.cpp
        model: HodgkinHuxley
    ...

## Creating a New Plugin

SANA-FE can run any models provided as user plugins. The plugin must be compiled
as a shared library containing one or more hardware models.
Models can execute arbitrary code, but interfaces must be derived either from
the  general `PipelineUnit` class, or one of the specialized `SynapseUnit`,
`DendriteUnit` or `SomaUnit` base classes.

SANA-FE's plugin mechanism makes it easy to integrate plugins with your
architectural simulations. However, a few steps are needed to get plugins
running:

1. You must make sure your plugin has been built as a shared library (`.so`),
either by updating the plugin CMake file or providing your own build scripts.
2. Your new plugin must implement a hardware model class with the hardware
functionality you want. The model class you implement must be derived from
`PipelineUnit` in `chip.hpp`, which defines the required interfaces. These are
enforced by pure virtual methods, including attribute parsing methods update
methods. For examples of different hardware models, see either `models.cpp` or
the `plugins` folder.
3. Finally, provide a class factory function that returns a new instance of
your model class. This has to be in the format `create_<modelname>`. For
example, for a `HodgkinHuxley` model, we would specify the following code in
the plugin C++ file:

```
extern "C" sanafe::PipelineUnit *create_HodgkinHuxley()
{
    return (sanafe::PipelineUnit *) new HodgkinHuxley();
}
```

It is recommended new users look through the rest of the `hodgkin_huxley.cpp`
file to see what an example plugin looks like.

## Legacy SNN Description (Netlist) Format

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

C++ code has been written using the C++17 standard.

# References

J. A. Boyle, M. Plagge, S. G. Cardwell, F. S. Chance and A. Gerstlauer,
"SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration,"
in IEEE Transactions on Computer-Aided Design of Integrated Circuits and
Systems, 2025,
[doi:10.1109/TCAD.2025.3537971](https://doi.org/10.1109/TCAD.2025.3537971).

J. A. Boyle, M. Plagge, S. G. Cardwell, F. S. Chance and A. Gerstlauer,
"Tutorial: Large-Scale Spiking Neuromorphic Architecture Exploration using
SANA-FE," in 2024 International Conference on Hardware/Software Codesign and
System Synthesis (CODES+ISSS), Raleigh, NC, USA,
[doi:10.1109/CODES-ISSS60120.2024.00007](https://doi.org/10.1109/CODES-ISSS60120.2024.00007).

J. A. Boyle, M. Plagge, S. G. Cardwell, F. S. Chance and A. Gerstlauer,
"Performance and Energy Simulation of Spiking Neuromorphic Architectures for
Fast Exploration," in 2023 International Conference on Neuromorphic Systems
(ICONS), Santa Fe, NM, USA,
[doi:10.1145/3589737.3605970](https://doi.org/10.1145/3589737.3605970).

# Contact
James Boyle: james.boyle@utexas.edu
