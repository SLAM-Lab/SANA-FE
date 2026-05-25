from importlib.resources import files
from . import load_arch, load_net

def load_example():
    """Load a bundled example architecture and SNN."""
    base = files("sanafe.examples")
    arch = load_arch(base / "example_chip.yaml")
    net  = load_net(base / "example_snn.yaml", arch)
    return arch, net
