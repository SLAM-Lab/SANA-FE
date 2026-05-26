<p align="center">
<img src="https://raw.githubusercontent.com/SLAM-Lab/SANA-FE/main/sana_fe_logo.svg" alt="SANA-FE" width="400" style="max-width: 100%; height: auto;">
</p>

[![PyPI version](https://img.shields.io/pypi/v/sanafe.svg)](https://pypi.org/project/sanafe/)
[![Documentation Status](https://readthedocs.org/projects/sana-fe/badge/?version=latest)](https://sana-fe.readthedocs.io/en/latest/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1109%2FTCAD.2025.3537971-blue)](https://doi.org/10.1109/TCAD.2025.3537971)

Simulating Advanced Neuromorphic Architectures for Fast Exploration (SANA-FE)
is a framework for modeling the energy usage and performance of different
neuromorphic hardware. Given a description of a neuromorphic chip and a mapped
spiking neural network (SNN), SANA-FE simulates execution at time-step
granularity and reports energy, latency, and detailed per-unit performance
statistics.

# Quickstart

SANA-FE has two interfaces: A Python API (recommended for most users) and a
standalone command-line simulator.

## Python

```python
import sanafe;

arch, net = sanafe.load_example()
chip = sanafe.SpikingChip(arch)
chip.load(sanafe.load_net('snn/example_snn.yaml', arch))
results = chip.sim(100)
print(results)
```

## Standalone Simulator

`./sim arch/example_chip.yaml snn/example_snn.yaml 100`

This simulates 100 time-steps of a small connected spiking neural network (SNN).
For an example simulation of a real-world architecture (Intel's Loihi) running
DVS gesture classification with spike traces, run:

`./sim -s arch/loihi.yaml snn/dvs.yaml 1000`

General usage:

`./sim [optional flags] <architecture description> <SNN description> <N timesteps>`

### Command-line Flags
* `-m`: Enable message traces to `messages.csv`
* `-n`: Use the (legacy) netlist format for SNNs, instead of YAML.
* `-o`: Output directory
* `-p`: Record the simulated performance of each timestep to `perf.csv`
* `-s`: Enable spike traces to `spikes.csv`
* `-t [simple/detailed/cycle]`: Specify the timing model  (default=`detailed`)
* `-v`: Enable potential (voltage) traces to `potential.csv`
* `-x`: Enable extra user neuron traces to `neurons.csv`
* `-N`: Number of neuron/message processing threads (default=1)
* `-S`: Number of scheduling threads (default=0, use main thread)

# Installation

The latest Python release is available on PyPI and can be installed with:
`pip install sanafe`

## Custom Builds

The project can also be manually installed from source. This project uses CMake
as its build system and dependency manager. To setup compilation, first create
a temporary build directory:

`mkdir build && cd build`

Run the following command in this build directory:

`cmake ..`

Then compile SANA-FE and copy it to the project directory by running the command:

`make -j ${NPROC} && make install && cd ..`

The option `-j` indicates the number of parallel build threads.

## Dependencies

Building this project requires `cmake`, `make`, and a compiler that supports the
C++17 standard (e.g., GCC >= 8, Clang >= 5). This project uses RapidYAML for all
YAML file parsing, and Booksim 2 for optional cycle-accurate NoC modeling. To
build the Python interfaces, you must also have Python >= 3.10 installed with
PyBind11. You can install PyBind11 using:

`pip install pybind11`

# Simulator Overview

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

# Simulator Outputs

If corresponding traces are enabled, output is saved to trace files with
hard-coded names using either csv or yaml extensions.

`spikes.csv`: The spikes for each time-step on probed neurons

`potential.csv`: The potentials for each time-step on probed neurons

`neurons.csv`: Any extra neuron state traces for each time-step (if supported).

`perf.csv`: Detailed statistics for each timestep and each hardware unit

`messages.csv`: Information on spike messages for each time-step

`run_summary.yaml`: High-level statistics for the simulation e.g. runtime

This project has been written in C++ and Python C++ code has been written using
the C++17 standard.

For more details, see:

- [Python API reference](https://sana-fe.readthedocs.io/en/latest/) — full API documentation on Read the Docs
- [SNN file format](snn/README.md) — YAML format for describing mapped spiking networks
- [Architecture file format](arch/README.md) — YAML format for describing neuromorphic hardware
- [Plugin authoring guide](plugins/README.md) — implementing custom hardware models in C++
- [Tutorials](tutorial/) — Jupyter notebooks walking through end-to-end use

# Citation

We hope that you find this project useful. If you use SANA-FE in your work,
please cite our paper:

James A. Boyle, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration," in
IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems (TCAD), vol. 44, no. 8, pp. 3165–3178, 2025,
[doi:10.1109/TCAD.2025.3537971](https://doi.org/10.1109/TCAD.2025.3537971).

```bibtex
@article{boyle2025sanafe,
  title={SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration},
  author={James A. Boyle and Mark Plagge and Suma George Cardwell and Frances S. Chance and Andreas Gerstlauer},
  journal={IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems (TCAD)},
  volume={44},
  number={8},
  pages={3165--3178},
  year={2025},
  doi={10.1109/TCAD.2025.3537971}
}
```

# Projects using SANA-FE

We are currently inviting any SANA-FE users to send their projects! We would
like to include references here.

# References

James A. Boyle, Jason Ho, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"Exploring Dendrites in Large-Scale Neuromorphic Architectures,"
in International Conference on Neuromorphic Systems (ICONS), Seattle, WA, USA, 2025,
[doi:10.1109/ICONS69015.2025.00018](https://doi.org/10.1109/ICONS69015.2025.00018).

James A. Boyle, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration,"
in IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems (TCAD), vol. 44, no. 8, pp. 3165–3178, 2025,
[doi:10.1109/TCAD.2025.3537971](https://doi.org/10.1109/TCAD.2025.3537971).

James A. Boyle, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"Tutorial: Large-Scale Spiking Neuromorphic Architecture Exploration using SANA-FE,"
in International Conference on Hardware/Software Codesign and System Synthesis (CODES+ISSS), Raleigh, NC, USA, 2024,
[doi:10.1109/CODES-ISSS60120.2024.00007](https://doi.org/10.1109/CODES-ISSS60120.2024.00007).

James A. Boyle, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"Performance and Energy Simulation of Spiking Neuromorphic Architectures for Fast Exploration,"
in International Conference on Neuromorphic Systems (ICONS), Santa Fe, NM, USA, 2023,
[doi:10.1145/3589737.3605970](https://doi.org/10.1145/3589737.3605970).

# License and Acknowledgements

Copyright (c) 2026 - The University of Texas at Austin (GPL 3)

This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.

# Contact
James Boyle: james.boyle@utexas.edu
