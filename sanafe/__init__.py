"""Top-level SANA-FE module"""
# Import pybind11 (C++) kernel under top-level
from sanafecpp import *

# Import Python submodules for convenient access
from ._examples import load_example, load_loihi, load_truenorth
