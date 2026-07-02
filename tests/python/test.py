#!/usr/bin/env python3
"""Regression unit tests for the SANA-FE Python API (`import sanafe`).

These tests exercise the public, documented API surface: building SNNs,
connecting neurons (direct / dense / sparse / conv2d), constructing and
loading architectures, compiling and simulating a SpikingChip, attribute
round-tripping, neuron indexing/slicing semantics, repr output, and pickling.

Run with:
    python -m unittest test_sanafe -v
"""

import math
import pickle
import unittest

import sanafe


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def load_reference_arch():
    """Return a ready-made ``(arch, snn)`` from the built-in example."""
    arch, snn = sanafe.load_example()
    return arch, snn


def first_core(arch):
    """Return the first mappable core of an architecture using the documented
    ``arch.tiles[0].cores[0]`` access path."""
    return arch.tiles[0].cores[0]


# --------------------------------------------------------------------------
# Module-level surface
# --------------------------------------------------------------------------

class TestModuleSurface(unittest.TestCase):
    def test_core_symbols_present(self):
        for sym in ("Network", "Architecture", "SpikingChip", "NeuronGroup",
                    "Neuron", "Connection", "Tile", "Core", "MappedNeuron"):
            self.assertTrue(hasattr(sanafe, sym), f"missing symbol: {sym}")

    def test_loader_functions_present(self):
        for fn in ("load_arch", "load_net"):
            self.assertTrue(callable(getattr(sanafe, fn, None)),
                            f"{fn} should be callable")

    def test_framework_and_model_attributes_exposed(self):
        self.assertTrue(hasattr(sanafe, "framework_attributes"))
        self.assertTrue(hasattr(sanafe, "model_attributes"))

    def test_builtin_models_include_lif_and_input(self):
        models = sanafe.model_attributes
        names = set(models.keys())
        # Documented built-in models that the test architecture relies on.
        self.assertIn("leaky_integrate_fire", names)
        self.assertIn("input", names)

    def test_buffer_position_enum(self):
        bp = getattr(sanafe, "BufferPosition", None)
        self.assertIsNotNone(bp)
        for member in ("buffer_before_dendrite_unit", "buffer_before_soma_unit",
                       "buffer_before_axon_out_unit", "buffer_positions"):
            self.assertTrue(hasattr(bp, member), f"enum missing {member}")

    def test_hardware_mapping_error_exposed(self):
        self.assertTrue(hasattr(sanafe, "HardwareMappingError"))
        self.assertTrue(issubclass(sanafe.HardwareMappingError, Exception))


# --------------------------------------------------------------------------
# Network / NeuronGroup construction
# --------------------------------------------------------------------------

class TestNetworkConstruction(unittest.TestCase):
    def setUp(self):
        self.net = sanafe.Network()

    def test_empty_network_repr(self):
        self.assertIsInstance(repr(self.net), str)

    def test_create_neuron_group(self):
        group = self.net.create_neuron_group("layer0", 8)
        self.assertEqual(group.get_name(), "layer0")
        self.assertEqual(len(group), 8)

    def test_groups_property_contains_created_group(self):
        self.net.create_neuron_group("layer0", 4)
        self.assertIn("layer0", self.net.groups)

    def test_getitem_returns_group(self):
        self.net.create_neuron_group("layerX", 3)
        self.assertEqual(self.net["layerX"].get_name(), "layerX")

    def test_getitem_missing_raises_indexerror(self):
        with self.assertRaises(IndexError):
            _ = self.net["does_not_exist"]

    def test_create_group_with_lif_attributes(self):
        # Documented LIF attribute names.
        group = self.net.create_neuron_group(
            "lif", 2,
            model_attributes={"threshold": 1.0, "reset": 0.0,
                              "leak_decay": 0.9, "bias": 0.1})
        self.assertEqual(len(group), 2)

    def test_create_group_logging_flags(self):
        group = self.net.create_neuron_group(
            "logged", 2, {"threshold": 1.0},
            log_spikes=True, log_potential=True)
        self.assertEqual(len(group), 2)


# --------------------------------------------------------------------------
# Neuron access: indexing, slicing, iteration
# --------------------------------------------------------------------------

class TestNeuronAccess(unittest.TestCase):
    def setUp(self):
        self.net = sanafe.Network()
        self.group = self.net.create_neuron_group("g", 5)

    def test_len(self):
        self.assertEqual(len(self.group), 5)

    def test_integer_index(self):
        self.assertEqual(self.group[0].get_id(), 0)

    def test_negative_index(self):
        self.assertEqual(self.group[-1].get_id(), 4)

    def test_out_of_range_raises(self):
        with self.assertRaises(IndexError):
            _ = self.group[5]

    def test_slice_returns_list(self):
        sl = self.group[1:4]
        self.assertIsInstance(sl, list)
        self.assertEqual([n.get_id() for n in sl], [1, 2, 3])

    def test_strided_slice(self):
        sl = self.group[0:5:2]
        self.assertEqual([n.get_id() for n in sl], [0, 2, 4])

    def test_iteration(self):
        self.assertEqual([n.get_id() for n in self.group], [0, 1, 2, 3, 4])

    def test_neurons_view(self):
        view = self.group.neurons
        self.assertEqual(len(view), 5)
        self.assertEqual(view[2].get_id(), 2)
        self.assertEqual([n.get_id() for n in view], [0, 1, 2, 3, 4])
        self.assertIsInstance(repr(view), str)

    def test_neuron_repr(self):
        self.assertIsInstance(repr(self.group[0]), str)

    def test_bad_index_type_raises_typeerror(self):
        with self.assertRaises(TypeError):
            _ = self.group["not_an_index"]

    def test_set_attributes_on_neuron(self):
        # set_attributes accepts model_attributes (documented).
        self.group[0].set_attributes(model_attributes={"bias": 0.2})


# --------------------------------------------------------------------------
# Connectivity
# --------------------------------------------------------------------------

class TestConnectivity(unittest.TestCase):
    def setUp(self):
        self.net = sanafe.Network()
        self.src = self.net.create_neuron_group("src", 4)
        self.dst = self.net.create_neuron_group("dst", 4)

    def test_connect_to_neuron_returns_index(self):
        idx = self.src[0].connect_to_neuron(self.dst[0])
        self.assertIsInstance(idx, int)
        self.assertGreaterEqual(idx, 0)

    def test_connect_to_neuron_with_weight(self):
        self.src[0].connect_to_neuron(self.dst[1], {"weight": 2.5})
        self.assertGreaterEqual(len(self.src[0].edges_out), 1)

    def test_edges_out_populated(self):
        self.src[0].connect_to_neuron(self.dst[0])
        self.src[0].connect_to_neuron(self.dst[1])
        self.assertEqual(len(self.src[0].edges_out), 2)

    def test_connection_properties(self):
        self.src[2].connect_to_neuron(self.dst[3], {"weight": 1.0})
        con = self.src[2].edges_out[0]
        self.assertIsInstance(repr(con), str)
        _ = con.pre_neuron
        _ = con.post_neuron
        _ = con.synapse_hw_name

    def test_connect_neurons_dense(self):
        weights = [1.0] * (len(self.src) * len(self.dst))
        self.src.connect_neurons_dense(self.dst, {"weight": weights})
        total = sum(len(self.src[i].edges_out) for i in range(len(self.src)))
        self.assertEqual(total, len(self.src) * len(self.dst))

    def test_connect_neurons_sparse(self):
        pairs = [(0, 1), (1, 2), (3, 0)]
        self.src.connect_neurons_sparse(
            self.dst, {"weight": [1.0, 2.0, 3.0]}, pairs)
        total = sum(len(self.src[i].edges_out) for i in range(len(self.src)))
        self.assertEqual(total, len(pairs))

    def test_connect_neurons_sparse_rejects_bad_pair(self):
        with self.assertRaises((ValueError, TypeError, RuntimeError)):
            self.src.connect_neurons_sparse(
                self.dst, {"weight": [1.0]}, [(0, 1, 2)])

    def test_connect_neurons_conv2d(self):
        # 2x2x1 input, 2x2 kernel, 1 output channel -> 1x1x1 output.
        src = self.net.create_neuron_group("conv_in", 4)
        dst = self.net.create_neuron_group("conv_out", 1)
        src.connect_neurons_conv2d(
            dst, {"weight": [1.0, 1.0, 1.0, 1.0]},
            input_width=2, input_height=2, input_channels=1,
            kernel_width=2, kernel_height=2, kernel_count=1,
            stride_width=1, stride_height=1)
        total = sum(len(src[i].edges_out) for i in range(len(src)))
        self.assertGreater(total, 0)


# --------------------------------------------------------------------------
# Attribute round-tripping through the binding layer
# --------------------------------------------------------------------------

class TestAttributeRoundTrip(unittest.TestCase):
    def setUp(self):
        self.net = sanafe.Network()
        self.src = self.net.create_neuron_group("src", 2)
        self.dst = self.net.create_neuron_group("dst", 2)

    def test_scalar_types_survive(self):
        self.src[0].connect_to_neuron(
            self.dst[0], {"weight": 1.5, "delay": 7, "label": "syn"})
        attrs = self.src[0].edges_out[0].synapse_attributes
        self.assertIn("weight", attrs)
        self.assertAlmostEqual(float(attrs["weight"]), 1.5, places=5)
        self.assertEqual(int(attrs["delay"]), 7)
        self.assertEqual(attrs["label"], "syn")

    def test_list_attribute(self):
        self.src[0].connect_to_neuron(self.dst[1], {"taps": [0.1, 0.2, 0.3]})
        attrs = self.src[0].edges_out[0].synapse_attributes
        self.assertIn("taps", attrs)
        self.assertEqual(len(attrs["taps"]), 3)

    def test_nested_dict_attribute(self):
        self.src[1].connect_to_neuron(self.dst[0], {"cfg": {"a": 1, "b": 2.0}})
        attrs = self.src[1].edges_out[0].synapse_attributes
        self.assertIn("cfg", attrs)


# --------------------------------------------------------------------------
# Architecture construction (programmatic) + standalone components
# --------------------------------------------------------------------------

class TestArchitectureComponents(unittest.TestCase):
    def test_reference_arch_has_tile_and_core(self):
        arch, _ = load_reference_arch()
        self.assertGreaterEqual(len(arch.tiles), 1)
        self.assertGreaterEqual(len(arch.tiles[0].cores), 1)
        self.assertIsInstance(repr(arch), str)


# --------------------------------------------------------------------------
# NeuronAddress pickling
# --------------------------------------------------------------------------

class TestNeuronAddressPickle(unittest.TestCase):
    def _make_address(self):
        net = sanafe.Network()
        a = net.create_neuron_group("a", 1)
        b = net.create_neuron_group("b", 1)
        a[0].connect_to_neuron(b[0], {"weight": 1.0})
        return a[0].edges_out[0].pre_neuron

    def test_address_repr(self):
        self.assertIsInstance(repr(self._make_address()), str)

    def test_address_pickle_round_trip(self):
        addr = self._make_address()
        restored = pickle.loads(pickle.dumps(addr))
        self.assertEqual(restored.group_name, addr.group_name)
        self.assertEqual(restored.neuron_offset, addr.neuron_offset)


# --------------------------------------------------------------------------
# End-to-end: build, map, load, simulate
# --------------------------------------------------------------------------

class TestSimulation(unittest.TestCase):
    """Build a tiny deterministic SNN, map it to the reference architecture's
    first core, and simulate. Uses the documented `input` model with an
    explicit spike train so firing is fully predictable, plus a downstream LIF
    neuron with a low threshold so it fires from a single weighted input.
    """

    def _build_chip(self, timesteps=20):
        arch, _ = load_reference_arch()
        core = first_core(arch)

        net = sanafe.Network()
        # Driver: input model emitting a known spike train.
        spike_train = [bool(t % 2 == 0) for t in range(timesteps)]
        drivers = net.create_neuron_group(
            "drivers", 1,
            model_attributes={"spikes": spike_train},
            soma_hw_name="demo_input[0]")
        # Target: LIF that fires when V > threshold.
        targets = net.create_neuron_group(
            "targets", 1,
            model_attributes={"threshold": 0.5, "reset": 0.0,
                              "leak_decay": 1.0},
            log_spikes=True, log_potential=True)

        drivers[0].connect_to_neuron(targets[0], {"weight": 1.0})
        drivers[0].map_to_core(core)
        targets[0].map_to_core(core)

        chip = sanafe.SpikingChip(arch)
        chip.load(net)
        return chip

    def test_construct_and_load(self):
        chip = self._build_chip()
        self.assertIsInstance(chip, sanafe.SpikingChip)

    def test_load_example_smoke(self):
        # Mirrors the documented quick-start smoke test exactly.
        # load_example() is a guaranteed built-in, so call it directly.
        arch, snn = sanafe.load_example()
        chip = sanafe.SpikingChip(arch)
        chip.load(snn)
        result = chip.sim(10)
        self.assertEqual(result["timesteps_executed"], 10)

    def test_sim_returns_documented_keys(self):
        chip = self._build_chip()
        result = chip.sim(timesteps=10)
        for key in ("timesteps_executed", "energy", "sim_time", "spikes",
                    "packets_sent", "neurons_updated", "neurons_fired"):
            self.assertIn(key, result, f"sim result missing {key}")
        self.assertEqual(result["timesteps_executed"], 10)

    def test_sim_energy_breakdown_keys(self):
        chip = self._build_chip()
        energy = chip.sim(timesteps=5)["energy"]
        for key in ("total", "synapse", "dendrite", "soma", "network"):
            self.assertIn(key, energy)
            self.assertTrue(math.isfinite(float(energy[key])))

    def test_sim_trace_keys_present(self):
        chip = self._build_chip()
        # Note: the public sim() signature exposes spike/potential/perf/message
        # traces (no separate neuron_trace argument).
        result = chip.sim(timesteps=3, spike_trace=True, potential_trace=True,
                          perf_trace=True, message_trace=True)
        for key in ("spike_trace", "potential_trace", "perf_trace",
                    "message_trace"):
            self.assertIn(key, result)

    def test_input_drives_lif_spiking(self):
        # With a low-threshold non-leaky LIF and unit-weight input spikes,
        # the target should fire at least once over the run.
        chip = self._build_chip(timesteps=20)
        result = chip.sim(timesteps=20)
        self.assertGreater(result["spikes"], 0,
                           "expected at least one spike from driven LIF")

    def test_get_power(self):
        chip = self._build_chip()
        chip.sim(timesteps=2)
        self.assertTrue(math.isfinite(float(chip.get_power())))

    def test_reset(self):
        chip = self._build_chip()
        chip.sim(timesteps=2)
        chip.reset()  # should not raise

    def test_mapped_neuron_groups_after_load(self):
        chip = self._build_chip()
        mapped = chip.mapped_neuron_groups
        self.assertIn("targets", mapped)
        self.assertEqual(len(mapped["targets"]), 1)


# --------------------------------------------------------------------------
# Determinism / regression guard
# --------------------------------------------------------------------------

class TestDeterminism(unittest.TestCase):
    def _run(self):
        arch, _ = load_reference_arch()
        core = first_core(arch)
        net = sanafe.Network()
        drivers = net.create_neuron_group(
            "d", 1, {"spikes": [True, False] * 5}, soma_hw_name="demo_input[0]")
        targets = net.create_neuron_group(
            "t", 1, {"threshold": 0.5, "reset": 0.0, "leak_decay": 1.0})
        drivers[0].connect_to_neuron(targets[0], {"weight": 1.0})
        drivers[0].map_to_core(core)
        targets[0].map_to_core(core)
        chip = sanafe.SpikingChip(arch)
        chip.load(net)
        return chip.sim(timesteps=10)

    def test_two_identical_runs_match(self):
        r1 = self._run()
        r2 = self._run()
        self.assertEqual(r1["spikes"], r2["spikes"])
        self.assertEqual(r1["neurons_fired"], r2["neurons_fired"])
        self.assertAlmostEqual(
            float(r1["energy"]["total"]),
            float(r2["energy"]["total"]), places=9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
