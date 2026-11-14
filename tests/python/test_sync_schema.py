"""Consistent timestep/sync costs and compatibility with legacy architectures."""

import tempfile
import unittest
from pathlib import Path

import sanafe
import yaml

ROOT = Path(__file__).resolve().parents[2]


class SyncSchema(unittest.TestCase):
    def run_chip(self, attributes, positions=((0, 0),), links=(), neuron_latency=0.0):
        payload = yaml.safe_load((ROOT / "arch/neuroscale.yaml").read_text())
        attrs = payload["architecture"]["attributes"]
        for name in (
            "energy_timestep_fixed",
            "latency_timestep_fixed",
            "energy_sync",
            "latency_sync",
        ):
            attrs.pop(name, None)
        attrs.update(attributes)

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

        if positions and max(tile for tile, _ in positions) >= 16:
            attrs["width"] = 8
            attrs["height"] = 4
            payload["architecture"]["tile"][0]["name"] = "neuroscale_tile[0..31]"
        clear(payload["architecture"]["tile"])
        payload["architecture"]["tile"][0]["core"][0]["soma"][0]["attributes"][
            "latency_access_neuron"
        ] = neuron_latency
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "arch.yaml"
            path.write_text(yaml.safe_dump(payload))
            arch = sanafe.load_arch(str(path))
        net = sanafe.Network()
        nodes = (
            net.create_neuron_group("quiet", len(positions), {"threshold": 1000.0})
            if positions
            else []
        )
        for node, (tile, core) in zip(nodes, positions):
            node.map_to_core(arch.tiles[tile].cores[core])
        for source, destination in links:
            nodes[source].connect_to_neuron(nodes[destination], {"weight": 1.0})
        chip = sanafe.SpikingChip(arch)
        chip.load(net)
        return chip.sim(3, timing_model="simple")

    def test_tile_costs_accept_the_same_scalar_list_and_map_forms(self):
        for energy, latency in (
            (7e-9, 11e-9),
            ([3e-9, 7e-9], [5e-9, 11e-9]),
            ({0: 3e-9, 1: 7e-9}, {0: 5e-9, 1: 11e-9}),
        ):
            result = self.run_chip(
                dict(
                    energy_timestep_fixed=2e-9,
                    latency_timestep_fixed=3e-9,
                    energy_sync={"per_tile": energy},
                    latency_sync={"per_tile": latency},
                )
            )
            self.assertAlmostEqual(result["energy"]["total"], 27e-9, delta=1e-20)
            self.assertAlmostEqual(result["sim_time"], 42e-9, delta=1e-20)

    def test_lookup_counts_tiles_not_cores(self):
        costs = dict(
            energy_timestep_fixed=2e-9,
            latency_timestep_fixed=3e-9,
            energy_sync={"per_tile": {1: 7e-9, 2: 9e-9}},
            latency_sync={"per_tile": {1: 11e-9, 2: 13e-9}},
        )
        for positions, energy, latency in (
            (((0, 0), (0, 1)), 27e-9, 42e-9),
            (((0, 0), (1, 0)), 33e-9, 48e-9),
        ):
            result = self.run_chip(costs, positions=positions)
            self.assertAlmostEqual(result["energy"]["total"], energy, delta=1e-20)
            self.assertAlmostEqual(result["sim_time"], latency, delta=1e-20)

    def test_loihi_split_preserves_total_delays_including_zero_and_one_tile(self):
        attrs = yaml.safe_load((ROOT / "arch/loihi.yaml").read_text())["architecture"][
            "attributes"
        ]
        for tiles, delay in (
            (0, 0.6e-6),
            (1, 0.6e-6),
            (2, 1.0e-6),
            (3, 1.0e-6),
            (4, 1.4e-6),
            (28, 1.4e-6),
            (29, 1.8e-6),
            (32, 1.8e-6),
        ):
            result = self.run_chip(
                {
                    "latency_timestep_fixed": attrs["latency_timestep_fixed"],
                    "latency_sync": attrs["latency_sync"],
                },
                positions=tuple((i, 0) for i in range(tiles)),
            )
            self.assertAlmostEqual(result["sim_time"], 3 * delay, delta=1e-18)

    def test_tile_delay_is_added_after_neuron_processing(self):
        result = self.run_chip(
            dict(
                latency_timestep_fixed=0.6e-6,
                latency_sync={"per_tile": {0: 0.0, 2: 0.4e-6}},
            ),
            positions=((0, 0), (1, 0)),
            neuron_latency=1e-6,
        )
        self.assertAlmostEqual(result["sim_time"], 6e-6, delta=1e-18)

    def test_energy_and_latency_share_endpoint_field_semantics(self):
        result = self.run_chip(
            dict(
                energy_timestep_fixed=3e-9,
                latency_timestep_fixed=1e-9,
                energy_sync=dict(
                    per_sending_core=5e-9,
                    per_receiving_core=7e-9,
                    per_message_sent=11e-9,
                    per_message_received=13e-9,
                ),
                latency_sync=dict(
                    per_sending_core=2e-9,
                    per_receiving_core=3e-9,
                    per_message_sent=4e-9,
                    per_message_received=5e-9,
                ),
            ),
            positions=((0, 0), (0, 1)),
            links=((0, 1), (1, 0)),
        )
        self.assertAlmostEqual(result["energy"]["total"], 225e-9, delta=1e-20)
        self.assertAlmostEqual(result["sim_time"], 45e-9, delta=1e-20)

    def test_legacy_scalar_costs_still_load(self):
        result = self.run_chip(
            dict(energy_sync=2e-9, latency_sync=5e-9, sync_model="fixed")
        )
        self.assertAlmostEqual(result["energy"]["total"], 6e-9, delta=1e-20)
        self.assertAlmostEqual(result["sim_time"], 15e-9, delta=1e-20)

    def test_legacy_loihi_table_retains_lookup_behavior(self):
        for positions, delay in ((((0, 0),), 0.6e-6), (((0, 0), (1, 0)), 1.0e-6)):
            result = self.run_chip(
                dict(
                    sync_model="table", latency_sync={1: 0.6e-6, 2: 1.0e-6, 4: 1.4e-6}
                ),
                positions=positions,
            )
            self.assertAlmostEqual(result["sim_time"], 3 * delay, delta=1e-18)

    def test_legacy_fixed_and_message_energy_keep_their_meaning(self):
        result = self.run_chip(
            dict(energy_sync=dict(fixed=3e-9, per_sync_message=7e-9)),
            positions=((0, 0), (0, 1)),
            links=((0, 1),),
        )
        self.assertAlmostEqual(result["energy"]["total"], 30e-9, delta=1e-20)

    def test_previous_processing_map_can_be_migrated_independently(self):
        result = self.run_chip(
            dict(
                latency_sync=1e-9,
                latency_sync_processing=dict(per_message_received=5e-9),
            ),
            positions=((0, 0), (0, 1)),
            links=((0, 1), (1, 0)),
        )
        self.assertAlmostEqual(result["sim_time"], 18e-9, delta=1e-20)

    def test_invalid_fixed_costs_are_rejected(self):
        for field in ("energy_timestep_fixed", "latency_timestep_fixed"):
            for value in (-1.0, float("inf"), float("nan"), [], {}, [1.0], {0: 1.0}):
                with (
                    self.subTest(field=field, value=value),
                    self.assertRaises(ValueError),
                ):
                    self.run_chip({field: value})

    def test_invalid_tile_costs_are_rejected(self):
        for field in ("energy_sync", "latency_sync"):
            for value in (-1.0, float("inf"), float("nan"), [], {}):
                with (
                    self.subTest(field=field, value=value),
                    self.assertRaises(ValueError),
                ):
                    self.run_chip({field: {"per_tile": value}})

    def test_new_fixed_fields_override_legacy_baselines_without_double_counting(self):
        result = self.run_chip(
            dict(
                energy_timestep_fixed=2e-9,
                energy_sync=dict(fixed=99e-9),
                latency_timestep_fixed=5e-9,
                latency_sync=99e-9,
            )
        )
        self.assertAlmostEqual(result["energy"]["total"], 6e-9, delta=1e-20)
        self.assertAlmostEqual(result["sim_time"], 15e-9, delta=1e-20)

    def test_conflicting_latency_maps_are_rejected(self):
        with self.assertRaises(ValueError):
            self.run_chip(dict(latency_sync={}, latency_sync_processing={}))


if __name__ == "__main__":
    unittest.main()
