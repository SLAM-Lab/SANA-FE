# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
SANA-FE backend for Sango networks on the NeuroScale chip.

Extends the generic SimSANAFE backend with what the NeuroScale benchmarks need:

  * its own model registry: the RTL-exact soma for the chip's neuron kind
    (`lif`, `rw` random walker or `lca` trace neuron) and the `neuroscale_synapse`
  * explicit placement: {Sango node name: ((x, y, core), neuron index)}, mapped in
    neuron-index order (the RTL updates neurons in index order, which fixes the
    order of LFSR draws within a core)
  * per-core LFSR seeding from the first neuron mapped to each core
  * chip-variant settings: neuron kind, LFSR width, trace shift, refractory period

Usage:

    from sanafe.sango.neuroscale import SimNeuroScale

    sim = SimNeuroScale(net)               # net is a built Sango network
    sim.compile(arch="arch/neuroscale.yaml", placements=placements,
                kind="lif")
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


# One soma per RTL neuron. The chip variant decides which one: every neuron of a
# chip is of its kind (neuron_models/lif, rw_*, lca_*).
KINDS = {
    "lif": dict(hw_name="neuroscale_lif", model="neuroscale_lif"),
    "rw": dict(hw_name="neuroscale_rw", model="neuroscale_rw"),
    "lca": dict(hw_name="neuroscale_lca", model="neuroscale_lca"),
}


def arch_for_kind(arch_dir, kind):
    """The architecture file for a neuron kind: neuroscale.yaml (lif, fig4's chip),
    neuroscale_rw.yaml, neuroscale_lca.yaml. Each has the one soma unit of its kind."""
    return os.path.join(arch_dir, "neuroscale.yaml" if kind == "lif" else f"neuroscale_{kind}.yaml")


def kind_for_variant(variant):
    """The soma kind of a chip variant (benchmarks.variants.Variant)."""
    if "trace" in variant.features:
        return "lca"
    return "rw" if "prob" in variant.features else "lif"


def _neuron_entry(kind):
    state = {
        "potential": _field("voltage", 0.0),
        "threshold": _field("threshold", 1.0),
        "bias": _field("bias", 0.0),
        "refractory_delay": _field(None, 0),
    }
    if kind == "rw":
        # Sango leak 1.0 (full leak) is the chip's decay bit; it has nothing in between
        state["decay"] = _field("leak", 1.0)
        state["spike_probability"] = _field("prob", 1.0)
    if kind == "lca":
        state["trace_shift"] = _field(None, 5)
        state["trace_frac_bits"] = _field(None, 8)
    return {"graph_type": "neuron", "hw_type": "soma", "hw_name": KINDS[kind]["hw_name"],
            "model_type": KINDS[kind]["model"], "param": {}, "state": state}


def neuroscale_registry(kind):
    """Sango registry for a chip of the given neuron kind ("lif", "rw" or "lca")."""
    neuron = _neuron_entry(kind)
    return {
        "LIF": copy.deepcopy(neuron),
        "pLIF": copy.deepcopy(neuron),
        "RefractoryLIF": copy.deepcopy(neuron),
        "PSP": {
            "graph_type": "synapse",
            "hw_type": "synapse",
            "model_type": "neuroscale_synapse",
            "param": {},
            "state": {"weight": _field("weight", 1.0)},
        },
    }


def send_times(messages, labels):
    """{(key, 0-based tick): seconds at which the neuron's spike message was sent}, for the
    neurons in ``labels`` = {"group.index": key}. A neuron without a fan-out sends no message."""
    first = messages.groupby(["src_neuron", "timestep"]).sent_timestamp.min()
    return {(labels[label], int(tick) - 1): float(when)
            for (label, tick), when in first.items() if label in labels}


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
        self.kind = "rw"
        self.lfsr_bits = 8
        self.seed_for_core = None
        self.core_times = None

    def import_registry(self):
        # The generic registry (model/registry.py) is for arch/sango.yaml; this
        #  backend's registry depends on the chip's neuron kind, set in compile()
        return neuroscale_registry("rw")

    def compile(self, placements=None, kind="rw", lfsr_bits=None, refractory_period=0,
                seed_for_core=None, trace_shift=None, trace_frac_bits=None, **kwargs):
        """
        placements: {node name: ((x, y, core), neuron index)}; None maps greedily
        kind: the chip's neuron: "lif" (floor at 0), "rw" (probabilistic, LFSR) or
            "lca" (trace neuron); see kind_for_variant
        lfsr_bits: LFSR width of the rw chips' neurons (8, 12 or 16)
        refractory_period: ticks skipped after a spike by RefractoryLIF neurons
        seed_for_core: f((x, y, core)) -> LFSR seed (rw); default (y * grid_x + x) *
            cores_per_tile + core + 1, the chip variants' rule
        trace_shift, trace_frac_bits: the lca chip's trace (TRACE_SHIFT, TRACE_BITS - 16)
        """
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {sorted(KINDS)}, not {kind!r}")
        self.kind = kind
        self.model_registry = neuroscale_registry(kind)
        self.placements = placements
        self.seed_for_core = seed_for_core
        if lfsr_bits is not None:
            self.lfsr_bits = lfsr_bits
        for name in ("LIF", "pLIF", "RefractoryLIF"):
            state = self.model_registry[name]["state"]
            if kind == "lca":
                if trace_shift is not None:
                    state["trace_shift"]["default"] = int(trace_shift)
                if trace_frac_bits is not None:
                    state["trace_frac_bits"]["default"] = int(trace_frac_bits)
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
            if self.kind == "rw":
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
        self.messages = messages
        return self.run_data

    def spike_times(self):
        """{(node index, 0-based tick): seconds the spike left its neuron}, see send_times."""
        labels = {f"{self.node_data[n]['group_name']}.{self.local_index[n]}": n
                  for n in range(self.num_nodes)}
        return send_times(self.messages, labels)


@dataclass
class ScaffoldRun:
    spikes: dict                 # {node name: [tick, ...]}, 0-based ticks
    latency: float               # total simulated time (s)
    core_times: pd.DataFrame     # time (s) each core finished each tick: index "tile.core"
    arch: object = None
    cores: dict = field(default_factory=dict)   # {node name: (x, y, core)}
    send_times: dict = field(default_factory=dict)   # {(node name, tick): seconds the spike left}


def run_scaffold(graph, ticks, arch, placements, kind="rw", lfsr_bits=None,
                 seed_for_core=None, trace_shift=None, trace_frac_bits=None, spike_trace=None):
    """Simulate a Fugu scaffold graph (as lowered to the chip) in SANA-FE.

    Use this instead of SimNeuroScale when the chip network is not the Sango network
    itself, e.g. after delays are expanded into relay neurons. Nodes carry threshold,
    bias, voltage, decay and optionally p (spike probability); edges carry weight.
    placements: {node name: ((x, y, core), neuron index)}. kind, lfsr_bits, trace_shift and
    trace_frac_bits are as in SimNeuroScale.compile.
    """
    arch = sanafe.load_arch(arch) if isinstance(arch, str) else arch
    names = list(graph.nodes)
    net = sanafe.Network()
    group = net.create_neuron_group("n", len(names), soma_hw_name=KINDS[kind]["hw_name"])
    neurons = {name: group[i] for i, name in enumerate(names)}
    offset = {name: i for i, name in enumerate(names)}
    for name, a in graph.nodes(data=True):
        attributes = {
            "potential": float(a.get("voltage", 0.0)),
            "threshold": float(a["threshold"]),
            "bias": float(a.get("bias", 0.0)),
        }
        if kind == "rw":
            attributes["decay"] = float(a.get("decay", 0.0))
            attributes["spike_probability"] = float(a.get("p", 1.0))
        if kind == "lca":
            attributes["trace_shift"] = int(trace_shift)
            attributes["trace_frac_bits"] = int(trace_frac_bits)
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
        if kind == "rw" and lfsr_bits is not None:
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
    times = send_times(messages, {f"n.{i}": name for i, name in enumerate(names)})
    return ScaffoldRun(spikes, float(core_times[ticks].max()), core_times, arch, cores, times)
