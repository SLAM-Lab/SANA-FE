"""Sync processing latency, overlap and synchronization compatibility."""

import tempfile
import unittest
from pathlib import Path

import sanafe
import yaml

ROOT = Path(__file__).resolve().parents[2]


class SyncLatency(unittest.TestCase):
    def architecture(
        self,
        processing=None,
        neuron=10e-9,
        fixed=0.0,
        table=False,
        protocol="neuroscale",
    ):
        data = yaml.safe_load((ROOT / "arch/neuroscale.yaml").read_text())
        attrs = data["architecture"]["attributes"]
        attrs["energy_sync"] = {}
        attrs["energy_timestep_fixed"] = 0.0
        attrs["sync_protocol"] = protocol
        attrs.pop("latency_sync", None)
        attrs["latency_timestep_fixed"] = {0: fixed, 1: fixed} if table else fixed
        if processing is not None:
            attrs["latency_sync"] = processing
        tile = data["architecture"]["tile"][0]

        def clear(value):
            if isinstance(value, dict):
                for key, child in value.items():
                    if key.startswith(("energy_", "latency_")):
                        value[key] = 0.0
                    else:
                        clear(child)
            elif isinstance(value, list):
                for child in value:
                    clear(child)

        clear(tile)
        tile["core"][0]["soma"][0]["attributes"]["latency_access_neuron"] = neuron
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "arch.yaml"
            path.write_text(yaml.safe_dump(data))
            return sanafe.load_arch(str(path))

    def run_chip(
        self, ticks, links, processing=None, neuron=10e-9, timing="simple", **kwargs
    ):
        arch = self.architecture(processing, neuron, **kwargs)
        chip = sanafe.SpikingChip(arch)
        net = sanafe.Network()
        neurons = net.create_neuron_group("quiet", 8, {"threshold": 1000.0})
        for i in range(8):
            neurons[i].map_to_core(arch.tiles[0].cores[i // 2])
        for source, target in links:
            neurons[source].connect_to_neuron(neurons[target], {"weight": 1.0})
        chip.load(net)
        result = chip.sim(ticks, timing_model=timing)
        self.assertEqual(result["energy"]["total"], 0.0)
        return result["sim_time"]

    def costs(self):
        return dict(
            per_sending_core=100e-9,
            per_receiving_core=100e-9,
            per_message_sent=20e-9,
            per_message_received=20e-9,
        )

    def test_omitted_and_zero_costs_keep_existing_latency(self):
        links = [(0, 2), (2, 0)]
        for processing in (None, {}):
            self.assertAlmostEqual(
                self.run_chip(12, links, processing), 12 * 20e-9, delta=1e-18
            )

    def test_cycles_expose_receive_work_and_duplicate_links_are_counted_once(self):
        links = [(0, 2), (1, 3), (2, 0), (3, 1)]
        self.assertAlmostEqual(
            self.run_chip(12, links, self.costs()), 12 * 240e-9, delta=1e-18
        )

    def test_sync_delay_charges_only_work_not_hidden_by_neurons(self):
        for neuron_latency in (10e-9, 100e-9, 200e-9):
            with self.subTest(neuron_latency=neuron_latency):
                self.assertAlmostEqual(
                    self.run_chip(
                        12, [(0, 2), (2, 0)], self.costs(), neuron=neuron_latency
                    ),
                    12 * max(2 * neuron_latency, 240e-9),
                    delta=1e-18,
                )

    def test_acyclic_progress_within_advance_window(self):
        times = [
            self.run_chip(ticks, [(0, 2), (2, 4)], self.costs(), neuron=100e-9)
            for ticks in (12, 24)
        ]
        self.assertAlmostEqual(times[1] - times[0], 12 * 200e-9, delta=1e-18)

    def test_source_only_core_needs_no_active_synapse_model(self):
        costs = dict(per_sending_core=90e-9, per_message_sent=5e-9)
        self.assertAlmostEqual(self.run_chip(1, [(0, 2)], costs), 95e-9, delta=1e-18)

    def test_detailed_receiver_occupancy_without_spikes(self):
        self.assertAlmostEqual(
            self.run_chip(12, [(0, 2), (2, 0)], self.costs(), timing="detailed"),
            12 * 240e-9,
            delta=1e-18,
        )

    def test_scalar_and_table_fixed_sync_remain_independent(self):
        for table in (False, True):
            self.assertAlmostEqual(
                self.run_chip(
                    12, [(0, 2), (2, 0)], self.costs(), fixed=1e-6, table=table
                ),
                12 * (1e-6 + 240e-9),
                delta=1e-18,
            )

    def test_analytical_barrier_overlaps_receive_work_with_neurons(self):
        self.assertAlmostEqual(
            self.run_chip(12, [(0, 2), (2, 0)], self.costs(), protocol="barrier"),
            12 * 240e-9,
            delta=1e-18,
        )

    def test_invalid_costs_are_rejected(self):
        for field in self.costs():
            for value in (-1.0, float("inf"), float("nan")):
                with (
                    self.subTest(field=field, value=value),
                    self.assertRaises(ValueError),
                ):
                    sanafe.SpikingChip(self.architecture({field: value}))

    def test_external_network_timing_rejects_unsupported_processing_budget(self):
        with self.assertRaisesRegex(ValueError, "Sync processing latency"):
            self.run_chip(1, [(0, 2)], self.costs(), timing="cycle")


if __name__ == "__main__":
    unittest.main()
