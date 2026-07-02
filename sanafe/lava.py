# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
# Implemented by Lance Lui as part of the capstone senior design project
"""
Lava to SANA-FE Backend

This module converts LAVA processes or serialization
object to the SNN representation runnable on SANA-FE

keep this file structure-
lava
├── src
│   ├── lava
│   │   ├── utils
│   │   │   ├── sanafe.py
│   │   │   │...
│   │   │...
│   │...
│...
SANAFE
├── sim.py
│...
"""
# As part of this backend, we have to access private Lava class details
# pylint: disable=protected-access

import lava
import sanafe


# Structural keys from Lava proc_params that describe the process shape
# rather than neuron model behavior. These are consumed during network
# construction and should NOT be forwarded as model_attributes.
_STRUCTURAL_KEYS = {"shape"}


def _extract_params(process):
    """Extract shape and model attributes from a Lava process.

    Splits ``proc_params._parameters`` into structural information (the
    neuron count derived from *shape*) and everything else, which is
    passed through to SANA-FE as ``model_attributes``.

    Args:
        process: A Lava ``AbstractProcess`` instance.

    Returns:
        tuple: (neuron_count, model_attributes_dict)
    """
    params = process.proc_params._parameters
    dim = params["shape"]
    neuron_count = dim[0] * (dim[1] if len(dim) > 1 else 1)

    model_attributes = {
        k: v for k, v in params.items() if k not in _STRUCTURAL_KEYS
    }
    return neuron_count, model_attributes


def _greedy_map_to_arch(network, arch):
    """Map every neuron in *network* to cores in *arch* using a greedy policy.

    Neurons are assigned to cores in the order they appear across groups.
    Each core is filled up to its ``max_neurons_supported`` limit before
    moving on to the next core.

    Args:
        network (sanafe.Network): The network whose neurons will be mapped.
        arch (sanafe.Architecture): Target architecture (e.g. from
            ``sanafe.load_loihi()``).

    Raises:
        RuntimeError: If the architecture runs out of cores before all
            neurons have been placed.
    """
    # Build a flat list of all cores across all tiles
    all_cores = []
    for tile in arch.tiles:
        for core in tile.cores:
            all_cores.append(core)

    core_idx = 0
    neurons_on_current_core = 0
    max_neurons = all_cores[core_idx].max_neurons_supported \
        if hasattr(all_cores[core_idx], "max_neurons_supported") else 1024

    for group_name in network.groups:
        group = network.groups[group_name]
        for neuron in group:
            if neurons_on_current_core >= max_neurons:
                core_idx += 1
                if core_idx >= len(all_cores):
                    raise RuntimeError(
                        f"Architecture has only {len(all_cores)} cores, "
                        f"which cannot hold all neurons in the network."
                    )
                neurons_on_current_core = 0
                max_neurons = (
                    all_cores[core_idx].max_neurons_supported
                    if hasattr(all_cores[core_idx], "max_neurons_supported")
                    else 1024
                )
            neuron.map_to_core(all_cores[core_idx])
            neurons_on_current_core += 1


def serial_to_sanafe(
    filename,
    arch=None,
    save_path="converted_network.yaml",
):
    """Load a serialized Lava process file and convert it to a SANA-FE network.

    If the serialized object is a single ``AbstractProcess``, it is
    forwarded to :func:`process_to_sanafe`.  Otherwise each element in
    the serialized list is turned into a neuron group and consecutive
    groups are connected pairwise.

    Args:
        filename (str): Path to the serialized Lava process file.
        arch (sanafe.Architecture, optional): Target architecture for
            hardware mapping.  Defaults to ``sanafe.load_loihi()``.
        save_path (str, optional): Where to save the resulting SANA-FE
            network.  Defaults to ``"converted_network.yaml"``.

    Returns:
        sanafe.Network: The constructed (and mapped) network.
    """
    p = lava.utils.serialization.load(filename)

    # Single process shortcut
    if isinstance(p[0], lava.magma.core.process.process.AbstractProcess):
        return process_to_sanafe(p[0], arch=arch, save_path=save_path)

    if arch is None:
        arch = sanafe.load_loihi()

    network = sanafe.Network()
    groups = []

    for i, proc in enumerate(p[0]):
        neuron_count, model_attrs = _extract_params(proc)
        group = network.create_neuron_group(
            f"layer_{i}",
            neuron_count,
            model_attributes=model_attrs,
        )
        groups.append(group)

        # Connect to the previous group (sequential topology)
        if len(groups) > 1:
            groups[-2].connect_neurons_dense(groups[-1], {})

    _greedy_map_to_arch(network, arch)
    network.save(save_path)
    return network


def process_to_sanafe(
    process,
    arch=None,
    save_path="converted_network.yaml",
):
    """Convert a live Lava ``AbstractProcess`` to a SANA-FE network.

    The process shape is interpreted as (*rows*, *cols*); each column
    becomes a separate neuron group and consecutive groups are connected.

    Args:
        process (AbstractProcess): A Lava process instance.
        arch (sanafe.Architecture, optional): Target architecture for
            hardware mapping.  Defaults to ``sanafe.load_loihi()``.
        save_path (str, optional): Where to save the resulting SANA-FE
            network.  Defaults to ``"converted_network.yaml"``.

    Returns:
        sanafe.Network: The constructed (and mapped) network.
    """
    if arch is None:
        arch = sanafe.load_loihi()

    network = sanafe.Network()
    params = process.proc_params._parameters
    dim = params["shape"]
    rows = dim[0]
    cols = dim[1] if len(dim) > 1 else 1

    model_attrs = {
        k: v for k, v in params.items() if k not in _STRUCTURAL_KEYS
    }

    groups = []
    for col_idx in range(cols):
        group = network.create_neuron_group(
            f"group_{col_idx}",
            rows,
            model_attributes=model_attrs,
        )
        groups.append(group)

        if len(groups) > 1:
            groups[-2].connect_neurons_dense(groups[-1], {})

    _greedy_map_to_arch(network, arch)
    network.save(save_path)
    return network
