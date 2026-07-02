# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
Example SNNs and Architectures.
"""
from importlib.resources import files
from sanafecpp import load_arch, load_net

def load_example():
    """Load a bundled example architecture and SNN."""
    base = files("sanafe.examples")
    arch = load_arch(base / "example_chip.yaml")
    net  = load_net(base / "example_snn.yaml", arch)
    return arch, net

def load_loihi():
    """Load the Loihi architecture file"""
    base = files("sanafe.examples")
    loihi_arch = load_arch(base / "loihi.yaml")
    return loihi_arch

def load_truenorth():
    """Load the TrueNorth architecture file"""
    base = files("sanafe.examples")
    truenorth_arch = load_arch(base / "truenorth.yaml")
    return truenorth_arch
