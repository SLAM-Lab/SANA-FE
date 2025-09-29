Welcome to SANA-FE's Documentation
=======================================

Architecting new neuromorphic chips involves several design decisions that can
affect overall system power consumption and performance. Architecture models can
be used to estimate the impact of different design decisions and thus inform
decision-making and exploration. SANA-FE (Simulating Advanced Neuromorphic
Architectures for Fast Exploration) is an open-source tool developed in a
collaboration between UT Austin and Sandia National Laboratories (SNL) to
rapidly and accurately model and simulate the energy and performance of
different neuromorphic hardware platforms.

.. image:: sana_fe.svg
   :alt: SANA-FE diagram
   :align: center


SANA-FE requires a description of the architecture you want to simulate, a
specification of the spiking neural network (SNN) application you wish to max
and execute. Optionally, you can also provide shared library hardware plugins
and simulation configuration. SANA-FE will simulate the given architecture and
give performane and energy estimates. Optionally, SANA-FE will generate detailed
hardware traces.

The SANA-FE kernel is primarily written in C++, but uses PyBind11 to generate
Python interfaces. SANA-FE can alternatively be used as a standalone C++
simulator. For more information on using SANA-FE without Python, visit the
project repository homepage at: https://github.com/SLAM-Lab/SANA-FE
The SANA-FE Python interface supports creating mapped SNNs, launching
simulations, and getting prediction data out for analysis.

Quick-start
===========

To install SANA-FE from PyPI, run the following

.. code-block:: bash

   python
   pip install sanafe

If you would rather install the latest SANA-FE from source, run

.. code-block:: bash

   git clone https://github.com/SLAM-Lab/SANA-FE sana-fe
   cd sana-fe
   pip install .

You can test your SANA-FE installation by running the following command

.. code-block:: bash

   python -c "import sanafe; \
   arch=sanafe.load_arch('arch/example.yaml'); c=sanafe.SpikingChip(arch); \
   c.load(sanafe.load_net('snn/example.yaml',arch)); c.sim(100)"

If this ran successfully, you should see the following print on the console

.. code-block:: bash

   Executed steps: [100/100]

For more examples of how to use SANA-FE, we have provided a set of Jupyter
notebooks at:
https://github.com/SLAM-Lab/SANA-FE/blob/main/tutorial

These tutorials can be run locally or using Google Colab.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api
   layers

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citation
========

We hope that you find this project useful. If you use SANA-FE in your work,
please cite our paper:

.. code-block:: bibtex

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


References
==========

James A. Boyle, Jason Ho, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"Exploring Dendrites in Large-Scale Neuromorphic Architectures,"
in International Conference on Neuromorphic Systems (ICONS), Seattle, WA, USA, 2025,
`doi: <https://doi.org/>`_.

James A. Boyle, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration,"in IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems (TCAD), vol. 44, no. 8, pp. 3165–3178, 2025,
`doi:10.1109/TCAD.2025.3537971 <https://doi.org/10.1109/TCAD.2025.3537971>`_.

James A. Boyle, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"Tutorial: Large-Scale Spiking Neuromorphic Architecture Exploration using SANA-FE,"
in International Conference on Hardware/Software Codesign and System Synthesis (CODES+ISSS), Raleigh, NC, USA, 2024,
`doi:10.1109/CODES-ISSS60120.2024.00007 <https://doi.org/10.1109/CODES-ISSS60120.2024.00007>`_.

James A. Boyle, Mark Plagge, Suma George Cardwell, Frances S. Chance, and Andreas Gerstlauer,
"Performance and Energy Simulation of Spiking Neuromorphic Architectures for Fast Exploration,"
in International Conference on Neuromorphic Systems (ICONS), Santa Fe, NM, USA, 2023,
`doi:10.1145/3589737.3605970 <https://doi.org/10.1145/3589737.3605970>`_.


Acknowledgement
===============

Copyright (c) 2025 - The University of Texas at Austin
This work was produced under contract #2317831 to National Technology and
Engineering Solutions of Sandia, LLC which is under contract
No. DE-NA0003525 with the U.S. Department of Energy.