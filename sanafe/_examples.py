from importlib.resources import files
from . import load_arch, load_net

def load_example():
    """Load a bundled example architecture and SNN."""
    base = files("sanafe.examples")
    arch = load_arch(base / "example_chip.yaml")
    net  = load_net(base / "example_snn.yaml", arch)
    return arch, net

def load_loihi():
    """Loihi the Loihi architecture file"""
    base = files("sanafe.examples")
    loihi_arch = load_arch(base / "loihi.yaml")
    return loihi_arch

def load_truenorth():
    """Loihi the Loihi architecture file"""
    base = files("sanafe.examples")
    truenorth_arch = load_arch(base / "truenorth.yaml")
    return truenorth_arch
