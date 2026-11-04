# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
SANA-FE backend for Sango networks on the NeuroScale chip.

Extends the generic SimSANAFE backend with what the NeuroScale benchmarks need:

  * its own model registry (the RTL-exact `neuroscale_plif` soma and the
    `neuroscale_synapse`), written for the installed Sango's registry layout
  * explicit placement: {Sango node name: ((x, y, core), neuron index)}, mapped in
    neuron-index order (the RTL updates neurons in index order, which fixes the
    order of LFSR draws within a core)
  * per-core LFSR seeding from the first neuron mapped to each core
  * chip-variant settings: floor at zero (lif), LFSR width, refractory period

Usage:

    from sanafe.sango.neuroscale import SimNeuroScale

    sim = SimNeuroScale(net)               # net is a built Sango network
    sim.compile(arch="arch/neuroscale3.yaml", placements=placements,
                floor_at_zero=True)
    sim.run(100)
    spikes = sim.get_spikes()              # per node, spike times (0-based ticks)
    print(sim.latency)                     # total simulated time (s)
"""

import copy
import os
import tempfile
from dataclasses import dataclass, field

import pandas as pd
import sanafe

from .sango import SimSANAFE


def _field(dsl, default):
    # 'dsl' is read by SimSANAFE, 'mapfrom' by the shared Sango Backend base class
    return {"dsl": dsl, "mapfrom": dsl, "default": default}


_NEURON = {
    "graph_type": "neuron",
    "hw_type": "soma",
    "hw_name": "neuroscale_plif",
    "model_type": "neuroscale_plif",
    "param": {},
    "state": {
        "potential": _field("voltage", 0.0),
        "threshold": _field("threshold", 1.0),
        "bias": _field("bias", 0.0),
        # Sango leak 1.0 (full leak) is the chip's decay bit; it has nothing in between
        "decay": _field("leak", 1.0),
        "spike_probability": _field("prob", 1.0),
        "floor_at_zero": _field(None, False),
        "refractory_delay": _field(None, 0),
    },
}

NEUROSCALE_REGISTRY = {
    "LIF": copy.deepcopy(_NEURON),
    "pLIF": copy.deepcopy(_NEURON),
    "RefractoryLIF": copy.deepcopy(_NEURON),
    "PSP": {
        "graph_type": "synapse",
        "hw_type": "synapse",
        "model_type": "neuroscale_synapse",
        "param": {},
        "state": {"weight": _field("weight", 1.0)},
    },
}


def default_seed(arch, core):
    """The chip variants' LFSR seed for a core: its logical index + 1 (0 would lock the LFSR)."""
    x, y, core_id = core
    tiles = {(t.x, t.y): t for t in arch.tiles}
    grid_x = max(t.x for t in arch.tiles) + 1
    return (y * grid_x + x) * len(tiles[(x, y)].cores) + core_id + 1


def seed_core(first_neuron, lfsr_bits, seed):
    """One LFSR per core: any neuron mapped to it carries the width and seed. The width
    is set first, since the seed is masked to it."""
    first_neuron.set_attributes(model_attributes={"lfsr_bits": lfsr_bits})
    first_neuron.set_attributes(model_attributes={"lfsr_seed": seed})


class SimNeuroScale(SimSANAFE):
    """SANA-FE backend for Sango networks mapped to the NeuroScale chip."""

    def __init__(self, net, debug=False, verbose=False):
        super().__init__(net, debug=debug, verbose=verbose)
        self.placements = None
        self.lfsr_bits = 8
        self.seed_for_core = None
        self.core_times = None

    def import_registry(self):
        # The Sango release this was written for looks for model/*.py, the branch's
        #  generic registry for registry/*.py; this backend carries its own
        return copy.deepcopy(NEUROSCALE_REGISTRY)

    def compile(self, placements=None, floor_at_zero=False, lfsr_bits=None,
                refractory_period=0, seed_for_core=None, **kwargs):
        """
        placements: {node name: ((x, y, core), neuron index)}; None maps greedily
        floor_at_zero: the lif chip floors the potential at zero, the rw chips do not
        lfsr_bits: LFSR width of the chip's probabilistic neurons (8, 12 or 16)
        refractory_period: ticks skipped after a spike by RefractoryLIF neurons
        seed_for_core: f((x, y, core)) -> LFSR seed; default (y * grid_x + x) *
            cores_per_tile + core + 1, the chip variants' rule
        """
        self.placements = placements
        self.seed_for_core = seed_for_core
        if lfsr_bits is not None:
            self.lfsr_bits = lfsr_bits
        for name in ("LIF", "pLIF", "RefractoryLIF"):
            state = self.model_registry[name]["state"]
            state["floor_at_zero"]["default"] = bool(floor_at_zero)
        self.model_registry["RefractoryLIF"]["state"]["refractory_delay"]["default"] = \
            int(refractory_period)
        super().compile(**kwargs)

    # -- mapping --------------------------------------------------------------
    def _map_to_cores(self):
        if self.placements is None:
            super()._map_to_cores()
            return
        tiles = {(t.x, t.y): t for t in self.arch.tiles}
        names = {n: name for name, n in self.node_index.items()}
        by_core = {}
        for n in range(self.num_nodes):
            core, index = self.placements[names[n]]
            by_core.setdefault(tuple(core), []).append((index, n))
        for core, members in by_core.items():
            x, y, core_id = core
            target = tiles[(x, y)].cores[core_id]
            members.sort()
            for _, n in members:
                self.sanafe_neurons[n].map_to_core(target)
            seed = (self.seed_for_core(core) if self.seed_for_core
                    else default_seed(self.arch, core))
            seed_core(self.sanafe_neurons[members[0][1]], self.lfsr_bits, seed)

    # -- run ------------------------------------------------------------------
    def run_backend(self, timesteps=1, spike_trace=None, **kwargs):
        message_trace = os.path.join(tempfile.mkdtemp(), "messages.csv")
        super().run_backend(timesteps, spike_trace=spike_trace, timing_model="simple",
                            message_trace=message_trace, **kwargs)
        messages = pd.read_csv(message_trace)
        # Time at which each core finished each tick
        self.core_times = (messages.groupby(["src_hw", "timestep"]).sent_timestamp
                           .first().unstack("timestep"))
        self.latency = float(self.core_times[timesteps].max())
        return self.run_data


@dataclass
class ScaffoldRun:
    spikes: dict                 # {node name: [tick, ...]}, 0-based ticks
    latency: float               # total simulated time (s)
    core_times: pd.DataFrame     # time (s) each core finished each tick: index "tile.core"
    arch: object = None
    cores: dict = field(default_factory=dict)   # {node name: (x, y, core)}


def run_scaffold(graph, ticks, arch, placements, floor_at_zero=False, lfsr_bits=None,
                 seed_for_core=None, spike_trace=None):
    """Simulate a Fugu scaffold graph (as lowered to the chip) in SANA-FE.

    Use this instead of SimNeuroScale when the chip network is not the Sango network
    itself, e.g. after delays are expanded into relay neurons. Nodes carry threshold,
    bias, voltage, decay and optionally p (spike probability); edges carry weight.
    placements: {node name: ((x, y, core), neuron index)}.
    """
    arch = sanafe.load_arch(arch) if isinstance(arch, str) else arch
    names = list(graph.nodes)
    net = sanafe.Network()
    group = net.create_neuron_group("n", len(names), soma_hw_name="neuroscale_plif")
    neurons = {name: group[i] for i, name in enumerate(names)}
    offset = {name: i for i, name in enumerate(names)}
    for name, a in graph.nodes(data=True):
        attributes = {
            "potential": float(a.get("voltage", 0.0)),
            "threshold": float(a["threshold"]),
            "bias": float(a.get("bias", 0.0)),
            "decay": float(a.get("decay", 0.0)),
            "spike_probability": float(a.get("p", 1.0)),
            "floor_at_zero": bool(floor_at_zero),
        }
        compartment = a.get("compartment") or {}
        if compartment.get("name") == "NeuroScaleRefractory":
            attributes["refractory_delay"] = int(compartment["period"])
        neurons[name].set_attributes(model_attributes=attributes, log_spikes=True)
    for source, target, e in graph.edges(data=True):
        neurons[source].connect_to_neuron(neurons[target], {"weight": float(e["weight"])})

    tiles = {(t.x, t.y): t for t in arch.tiles}
    by_core = {}
    for name in names:
        core, index = placements[name]
        by_core.setdefault(tuple(core), []).append((index, name))
    cores = {}
    for core, members in by_core.items():
        x, y, core_id = core
        target = tiles[(x, y)].cores[core_id]
        members.sort()
        for _, name in members:
            neurons[name].map_to_core(target)
            cores[name] = core
        if lfsr_bits is not None:
            seed = seed_for_core(core) if seed_for_core else default_seed(arch, core)
            seed_core(neurons[members[0][1]], lfsr_bits, seed)

    chip = sanafe.SpikingChip(arch)
    chip.load(net)
    work = tempfile.mkdtemp()
    spike_path = spike_trace or os.path.join(work, "spikes.csv")
    message_path = os.path.join(work, "messages.csv")
    with open(spike_path, "w", encoding="utf-8") as trace:
        chip.sim(ticks, spike_trace=trace, timing_model="simple", message_trace=message_path)
    spikes = {name: [] for name in names}
    with open(spike_path, encoding="utf-8") as trace:
        trace.readline()
        for line in trace:
            label, timestep = line.strip().rsplit(",", 1)
            spikes[names[int(label.rsplit(".", 1)[1])]].append(int(float(timestep)) - 1)
    messages = pd.read_csv(message_path)
    core_times = messages.groupby(["src_hw", "timestep"]).sent_timestamp.first().unstack("timestep")
    return ScaffoldRun(spikes, float(core_times[ticks].max()), core_times, arch, cores)
