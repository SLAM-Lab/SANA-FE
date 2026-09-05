# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
Minimal SANA-FE backend for Sango.

Usage:

    from sanafe.sango import SimSANAFE

    sim = SimSANAFE(net)
    sim.compile(arch="arch/loihi.yaml")
    sim.run(10)
    spikes = sim.get_spikes()
    print(sim.energy, sim.latency)
"""

# TODO: there's a bug with the delay being provided as a float by sango when it should be integer

import os
import tempfile
from collections import defaultdict

try:
    import sanafe
except ImportError:  # keep the import optional, as the other backends do
    sanafe = None

from sango.backend import Backend

# Graph attributes that are Sango specific
RESERVED_ATTRIBUTES = ("model", "group_name", "times")

# pylint: disable=invalid-name,no-self-argument,attribute-defined-outside-init
class SimSANAFE(Backend):
    """SANA-FE hardware-performance backend."""

    def __init__(self, net, debug=False, verbose=False):
        if sanafe is None:
            raise ImportError("sanafe package is required for SimSANAFE "
                              "(pip install sanafe)")
        super().__init__(net, debug=debug, verbose=verbose)

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

    def to_backend(self, arch="arch/sango.yaml", record='all',
                   soma_hw_name=None, input_hw_name=None,
                   dendrite_hw_name=None, max_neurons_per_core=None, **kwargs):
        """Translate the processed Sango graph into a SANA-FE network."""
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

    def _is_input(self, model):
        return self.model_registry[model]["graph_type"] == "input"

    def _record_set(self):
        """Which node indices should be traced."""
        if self.record_spec == 'all':
            return set(range(self.num_nodes))
        if not self.record_spec:
            return set()
        return None

    def _spike_raster(self, times):
        """Translate Sango spike times to dense raster array for SANA-FE."""
        if not len(times):
            return []
        steps = int(max(times)) + 1
        raster = [False] * steps
        for t in times:
            raster[int(t)] = True
        return raster

    def _build_network(self):
        self.sanafe_net = sanafe.Network()
        self.sanafe_groups = {}
        self.sanafe_neurons = {}
        to_record = self._record_set()
        node_names = {n: name for name, n in self.node_index.items()}

        # Sango has already formed neuron groups that SANA-FE needs
        for group_name, count in self.group_count.items():
            model = self.group_models[group_name]
            hw = {}

            if not self._is_input(model) and self.soma_hw_name:
                hw['soma_hw_name'] = self.soma_hw_name
            if self.dendrite_hw_name:
                hw['default_dendrite_hw_name'] = self.dendrite_hw_name
            self.sanafe_groups[group_name] = \
                self.sanafe_net.create_neuron_group(group_name, count, **hw)
            if self.debug:
                print(f"[sanafe] group {group_name}: {count} x {model}")

        # Per-neuron attributes, forwarded as-is
        for n, data in enumerate(self.node_data):
            neuron = self.sanafe_groups[data['group_name']][self.local_index[n]]
            self.sanafe_neurons[n] = neuron

            attributes = self._attributes(data)
            if self._is_input(data['model']):
                attributes['spikes'] = \
                    self._spike_raster(self.input_data[node_names[n]])

            neuron.set_attributes(model_attributes=attributes)
            if n in to_record:
                neuron.set_attributes(log_spikes=True, log_potential=True)

        # Connect edges and forward edge attributes
        for s in range(self.num_nodes):
            for t, data in self.edge_data[s].items():
                self.sanafe_neurons[s].connect_to_neuron(
                    self.sanafe_neurons[t], self._attributes(data))


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
            model = self.node_data[n]["model"]
            if self.model_registry[model]["graph_type"] == "input":
                # Input somas are replicated per-core units and need an index
                neuron.set_attributes(
                    soma_hw_name=f"{self.input_hw_name}[{inputs_on_core}]")
                inputs_on_core += 1

            neuron.map_to_core(cores[core_id])
            neurons_on_core += 1

    def _attributes(self, data):
        """Forward Sango attributes as-is, in the types the models declare."""
        types = self.model_registry[data['model']].get('attribute_types', {})
        return {key: types[key](value) if key in types else value
                for key, value in data.items()
                if key not in RESERVED_ATTRIBUTES}

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
            header = trace.readline()  # "neuron,timestep", ignore header
            for line in trace:
                line = line.strip()
                if not line:
                    continue
                name, timestep = line.rsplit(",", 1)
                group_name, offset = name.rsplit(".", 1)
                index = (self.group_offset[self.group_index[group_name]]
                         + int(offset))
                # Account for offset (1) between SANA-FE and Sango timesteps
                spikes[index].append(float(timestep) - 1.0)

        self.spike_list = [sorted(spikes[i]) for i in range(self.num_nodes)]
        return self.spike_list
