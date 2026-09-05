# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
Minimal SANA-FE backend for Sango.

Usage:

    from sanafe.sango import SimSANAFE

    sim = SimSANAFE(net)
    sim.compile(arch="arch/sango.yaml")
    sim.run(10)
    spikes = sim.get_spikes()
    print(sim.energy, sim.latency)
"""

import os
import tempfile
from collections import defaultdict

try:
    import sanafe
except ImportError:  # keep the import optional, as the other backends do
    sanafe = None

from sango.backend import Backend

# Which group-level keyword a registry entry's 'hw_type' selects
HW_KEYWORDS = {"soma": "soma_hw_name",
               "dendrite": "default_dendrite_hw_name",
               "synapse": "default_synapse_hw_name"}

# pylint: disable=invalid-name,no-self-argument,attribute-defined-outside-init
class SimSANAFE(Backend):
    """SANA-FE hardware-performance backend."""

    def __init__(self, net, debug=False, verbose=False):
        if sanafe is None:
            raise ImportError("sanafe package is required for SimSANAFE "
                              "(pip install sanafe)")
        # Populates self.model_registry from registry/*.py
        super().__init__(net, debug=debug, verbose=verbose)

        self.dt = 1.0
        self.arch = None
        self.arch_name = None
        self.chip = None
        self.sanafe_net = None
        self.sanafe_groups = {}   # {group_name: sanafe.NeuronGroup}
        self.sanafe_neurons = {}  # {global node index: sanafe.Neuron}
        self.run_data = None
        self.energy = None
        self.latency = None

        self.soma_hw_name = None  # None means pick first soma unit
        self.input_hw_name = "sango_inputs"
        self.dendrite_hw_name = None
        self.max_neurons_per_core = 1024

    def to_backend(self, arch="arch/sango.yaml", record="all", dt=None,
                   soma_hw_name=None, input_hw_name=None,
                   dendrite_hw_name=None, max_neurons_per_core=None, **kwargs):
        """Translate the processed Sango graph into a SANA-FE network."""
        if dt is not None:
            self.dt = dt
        if soma_hw_name is not None:
            self.soma_hw_name = soma_hw_name
        if input_hw_name is not None:
            self.input_hw_name = input_hw_name
        if dendrite_hw_name is not None:
            self.dendrite_hw_name = dendrite_hw_name
        if max_neurons_per_core is not None:
            self.max_neurons_per_core = max_neurons_per_core

        self.arch_name = arch
        self.arch = sanafe.load_arch(arch)
        self.record_spec = record

        self._build_network()
        self._map_to_cores()

    @staticmethod
    def _represent(entry, value):
        """Apply a registry "rep" hint, the way Brian's backend applies "unit"."""
        if entry.get("rep") == "tick":  # i.e. timesteps
            return int(round(value))
        return value

    def rekey_model(self, data):
        """As the shared base class, but honouring "rep"."""
        entry = self.model_registry[data["model"]]
        for key, value in entry.get("state", {}).items():
            raw = (data.pop(value["dsl"]) if value["dsl"] is not None
                   else value["default"])
            data[key] = self._represent(value, raw)
        return data

    def rekey_param(self, data):
        """As the shared base class, but honouring "rep"."""
        entry = self.model_registry[data["model"]]
        keys, params = [], []
        for key, value in entry.get("param", {}).items():
            raw = (data.pop(value["dsl"]) if value["dsl"] is not None
                   else value["default"])
            keys.append(key)
            params.append(self._represent(value, raw))
        return tuple(keys), tuple(params)

    def _is_input(self, model):
        return self.model_registry[model]["graph_type"] == "input"

    def _record_set(self):
        """Which node indices should be traced."""
        if self.record_spec == "all":
            return set(range(self.num_nodes))
        if not self.record_spec:
            return set()

        # A list of Sango path names, as the Fugu backend accepts
        return {self.node_index[name] for name in self.record_spec}

    def _spike_raster(self, times):
        """Translate Sango spike times to dense raster array for SANA-FE."""
        if len(times) == 0:
            return []
        raster = [False] * (int(max(times) / self.dt) + 1)
        for t in times:
            raster[int(round(t / self.dt))] = True
        return raster

    def _build_network(self):
        self.sanafe_net = sanafe.Network()
        self.sanafe_groups = {}
        self.sanafe_neurons = {}
        to_record = self._record_set()
        node_names = {n: name for name, n in self.node_index.items()}

        # Sango has already grouped neurons, so reuse this for SANA-FE's neuron
        #  groups
        for group_name, count in self.group_count.items():
            model = self.group_models[group_name]
            entry = self.model_registry[model]
            # Sango's shared() parameters map onto SANA-FE group attributes
            hw = {"model_attributes": dict(self.group_params.get(group_name, {}))}
            # Input somas are per-core indexed units
            if not self._is_input(model) and self.soma_hw_name:
                hw["soma_hw_name"] = self.soma_hw_name
            if self.dendrite_hw_name:
                hw["default_dendrite_hw_name"] = self.dendrite_hw_name

            hw_name = entry.get("hw_name")
            if isinstance(hw_name, str):
                hw[HW_KEYWORDS[entry["hw_type"]]] = hw_name
            elif hw_name is not None:
                for hw_type, name in hw_name.items():
                    hw[HW_KEYWORDS[hw_type]] = name
            self.sanafe_groups[group_name] = \
                self.sanafe_net.create_neuron_group(group_name, count, **hw)
            if self.debug:
                print(f"[sanafe] group {group_name}: {count} x {model}")

        # Per-neuron attributes, forwarded as-is
        for n, data in enumerate(self.node_data):
            neuron = self.sanafe_groups[data["group_name"]][self.local_index[n]]
            self.sanafe_neurons[n] = neuron

            entry = self.model_registry[data["model"]]
            attributes = {key: data[key] for key in entry["state"]}
            if self._is_input(data["model"]):
                attributes["spikes"] = \
                    self._spike_raster(self.input_data[node_names[n]])

            neuron.set_attributes(model_attributes=attributes)
            if n in to_record:
                neuron.set_attributes(log_spikes=True, log_potential=True)

        # Edges and edge attributes, likewise forwarded as-is
        for s in range(self.num_nodes):
            for t, data in self.edge_data[s].items():
                entry = self.model_registry[data['model']]
                attributes = {key: data[key] for key in entry['state']}
                self.sanafe_neurons[s].connect_to_neuron(
                    self.sanafe_neurons[t], attributes)

    def _map_to_cores(self):
        """Fill cores in order, one Sango group at a time using greedy approach."""
        cores = self.arch.cores()
        core_id = 0
        neurons_on_core = 0
        inputs_on_core = 0

        for n in range(self.num_nodes):
            if neurons_on_core >= self.max_neurons_per_core:
                core_id += 1
                neurons_on_core = 0
                inputs_on_core = 0
            if core_id >= len(cores):
                raise RuntimeError(
                    f"Network needs more than {len(cores)} cores on "
                    f"'{self.arch_name}'")

            neuron = self.sanafe_neurons[n]
            if self._is_input(self.node_data[n]["model"]):
                # Input somas are replicated per-core units and need an index
                neuron.set_attributes(
                    soma_hw_name=f"{self.input_hw_name}[{inputs_on_core}]")
                inputs_on_core += 1

            neuron.map_to_core(cores[core_id])
            neurons_on_core += 1

    def run_backend(self, timesteps=1, spike_trace=None, **kwargs):
        self.timesteps = timesteps
        self.chip = sanafe.SpikingChip(self.arch)
        self.chip.load(self.sanafe_net)

        self._trace_path = spike_trace or os.path.join(
            tempfile.mkdtemp(), "spikes.csv")
        with open(self._trace_path, "w", encoding="utf-8") as trace:
            self.run_data = self.chip.sim(timesteps, spike_trace=trace,
                                          **kwargs)

        self.energy = self.run_data.get("energy")
        self.latency = self.run_data.get("sim_time")
        return self.run_data

    def read_spikes(self):
        """Parse the SANA-FE spike trace back into Sango node order."""
        spikes = defaultdict(list)
        with open(self._trace_path, encoding="utf-8") as trace:
            trace.readline()  # header, "neuron,timestep"
            for line in trace:
                line = line.strip()
                if not line:
                    continue
                name, timestep = line.rsplit(",", 1)
                group_name, offset = name.rsplit(".", 1)
                index = (self.group_offset[self.group_index[group_name]]
                         + int(offset))
                # Account for offset (1) between SANA-FE and Sango timesteps
                spikes[index].append((float(timestep) - 1.0) * self.dt)

        self.spike_list = [sorted(spikes[i]) for i in range(self.num_nodes)]
        return self.spike_list