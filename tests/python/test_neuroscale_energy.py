"""Regression tests for NeuroScale synapse and chip synchronization energy."""

import copy
import tempfile
import unittest
from pathlib import Path

import yaml

import sanafe

SANAFE_ROOT = Path(__file__).resolve().parents[2]


class EnergyAccounting(unittest.TestCase):
    def architecture(
        self,
        sync=None,
        read=0.0,
        scan=0.0,
        handshake=0.0,
        handshake_base=0.0,
        extra_synapse=False,
    ):
        payload = yaml.safe_load((SANAFE_ROOT / "arch/neuroscale.yaml").read_text())
        attrs = payload["architecture"]["attributes"]
        attrs.pop("energy_sync", None)
        attrs.pop("energy_timestep_fixed", None)
        if handshake or handshake_base:
            attrs["energy_sync"] = {
                "per_sending_core": handshake_base,
                "per_receiving_core": handshake_base,
                "per_message_received": handshake,
            }
        if sync is not None:
            attrs["energy_timestep_fixed"] = sync
        tile = payload["architecture"]["tile"][0]
        core = tile["core"][0]

        # Isolate sync and synapse energy from all other hardware costs.
        def clear(value):
            if isinstance(value, dict):
                for key, child in value.items():
                    if key.startswith("energy_"):
                        value[key] = 0.0
                    else:
                        clear(child)
            elif isinstance(value, list):
                for child in value:
                    clear(child)

        clear(tile)
        core["synapse"][0]["attributes"].update(
            energy_read_synapse=read,
            energy_scan_axon=scan,
        )
        if extra_synapse:
            second = copy.deepcopy(core["synapse"][0])
            second["name"] = "extra_synapses"
            core["synapse"].append(second)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "arch.yaml"
            path.write_text(yaml.safe_dump(payload))
            arch = sanafe.load_arch(str(path))
        return arch

    def chip(self, sync=None, read=0.0, scan=0.0, spiking=False, cores=1):
        arch = self.architecture(sync=sync, read=read, scan=scan)
        net = sanafe.Network()
        quiet = (
            net.create_neuron_group("quiet", cores, {"threshold": 1000.0})
            if cores
            else None
        )
        for i in range(cores):
            quiet[i].map_to_core(arch.tiles[i].cores[0])
        if spiking:
            src = net.create_neuron_group("source", 1, {"threshold": 1.0, "bias": 2.0})
            src[0].map_to_core(arch.tiles[0].cores[0])
            src[0].connect_to_neuron(quiet[0], {"weight": 1.0, "axon": 3})
        chip = sanafe.SpikingChip(arch)
        chip.load(net)
        return chip

    def test_sync_is_once_per_chip_tick_and_repeated_calls(self):
        for timing in ("simple", "detailed"):
            for cores in (0, 1, 2):
                chip = self.chip(sync=2e-9, cores=cores)
                for ticks in (3, 4):
                    energy = chip.sim(ticks, timing_model=timing)["energy"]
                    self.assertAlmostEqual(energy["total"], ticks * 2e-9, delta=1e-20)
                    self.assertAlmostEqual(energy["sync"], energy["total"], delta=1e-20)

    def test_omitted_sync_defaults_to_zero(self):
        self.assertEqual(
            self.chip().sim(3, timing_model="simple")["energy"]["total"], 0.0
        )

    def test_read_and_persistent_scan_costs_are_linear(self):
        # Tick 1 emits a spike; ticks 2..5 each read one synapse and scan
        # axon IDs 0..3. Returned costs must never reuse the unit's total.
        energy = self.chip(read=2e-9, scan=3e-9, spiking=True).sim(
            5, timing_model="simple"
        )["energy"]
        self.assertAlmostEqual(energy["synapse"], 4 * (2e-9 + 4 * 3e-9), delta=1e-20)
        self.assertAlmostEqual(energy["total"], energy["synapse"], delta=1e-20)

    def handshake_network(self, arch, links, split_units=False):
        net = sanafe.Network()
        neurons = net.create_neuron_group("quiet", 8, {"threshold": 1000.0})
        for i in range(8):
            neurons[i].map_to_core(arch.tiles[0].cores[i // 2])
        for index, (source, target) in enumerate(links):
            neurons[source].connect_to_neuron(
                neurons[target],
                {"weight": 1.0},
                synapse_hw_name="extra_synapses" if split_units and index % 2 else None,
            )
        return net

    def test_handshake_counts_unique_directed_core_links(self):
        for links, count in [
            ([], 0),
            ([(0, 2), (1, 2), (0, 3)], 1),
            ([(0, 2), (2, 0)], 2),
            ([(0, 2), (4, 2)], 2),
            ([(0, 1), (1, 0)], 1),
        ]:
            with self.subTest(links=links):
                arch = self.architecture(handshake=2e-9)
                chip = sanafe.SpikingChip(arch)
                chip.load(self.handshake_network(arch, links))
                for ticks in (3, 4):
                    result = chip.sim(ticks, timing_model="simple")
                    self.assertEqual(result["neurons_fired"], 0)
                    self.assertAlmostEqual(
                        result["energy"]["total"], ticks * count * 2e-9, delta=1e-20
                    )
                    self.assertAlmostEqual(
                        result["energy"]["sync"],
                        result["energy"]["total"],
                        delta=1e-20,
                    )

    def test_handshake_is_shared_between_synapse_units(self):
        arch = self.architecture(handshake=2e-9, extra_synapse=True)
        chip = sanafe.SpikingChip(arch)
        chip.load(self.handshake_network(arch, [(0, 2), (1, 3)], split_units=True))
        self.assertAlmostEqual(
            chip.sim(3, timing_model="simple")["energy"]["total"], 6e-9, delta=1e-20
        )

    def test_reloading_network_clears_handshake_links(self):
        arch = self.architecture(handshake=2e-9)
        chip = sanafe.SpikingChip(arch)
        chip.load(self.handshake_network(arch, [(0, 2)]))
        self.assertAlmostEqual(
            chip.sim(3, timing_model="simple")["energy"]["total"], 6e-9, delta=1e-20
        )
        chip.load(self.handshake_network(arch, []))
        self.assertEqual(chip.sim(3, timing_model="simple")["energy"]["total"], 0.0)

    def test_handshake_base_counts_nonempty_tables(self):
        cases = [
            ([], 0, 0),
            ([(0, 2), (1, 3)], 2, 1),
            ([(0, 2), (0, 4)], 3, 2),
            ([(0, 4), (2, 4)], 3, 2),
            ([(0, 2), (2, 4)], 4, 2),
            ([(0, 2), (2, 0)], 4, 2),
            ([(0, 1), (1, 0)], 2, 1),
            ([(0, 2), (0, 4), (0, 6)], 4, 3),
        ]
        for links, tables, count in cases:
            for timing in ("simple", "detailed"):
                with self.subTest(links=links, timing=timing):
                    arch = self.architecture(handshake=2e-9, handshake_base=3e-9)
                    chip = sanafe.SpikingChip(arch)
                    chip.load(self.handshake_network(arch, links))
                    for ticks in (3, 4):
                        energy = chip.sim(ticks, timing_model=timing)["energy"]
                        self.assertAlmostEqual(
                            energy["total"],
                            ticks * (count * 2e-9 + tables * 3e-9),
                            delta=1e-20,
                        )

    def test_handshake_base_is_shared_between_synapse_units(self):
        for links, tables, count in [
            ([(0, 2), (1, 3)], 2, 1),
            ([(0, 4), (2, 5)], 3, 2),
        ]:
            with self.subTest(links=links):
                arch = self.architecture(
                    handshake=2e-9, handshake_base=3e-9, extra_synapse=True
                )
                chip = sanafe.SpikingChip(arch)
                chip.load(self.handshake_network(arch, links, split_units=True))
                self.assertAlmostEqual(
                    chip.sim(3, timing_model="simple")["energy"]["total"],
                    3 * (count * 2e-9 + tables * 3e-9),
                    delta=1e-20,
                )

    def test_reloading_clears_fixed_table_ownership(self):
        arch = self.architecture(handshake=2e-9, handshake_base=3e-9)
        chip = sanafe.SpikingChip(arch)
        for links, tables, count in [
            ([(0, 2)], 2, 1),
            ([], 0, 0),
            ([(0, 2), (0, 4)], 3, 2),
        ]:
            chip.load(self.handshake_network(arch, links))
            self.assertAlmostEqual(
                chip.sim(3, timing_model="simple")["energy"]["total"],
                3 * (count * 2e-9 + tables * 3e-9),
                delta=1e-20,
            )

    def test_reset_keeps_handshake_energy(self):
        arch = self.architecture(handshake=2e-9, handshake_base=3e-9)
        chip = sanafe.SpikingChip(arch)
        chip.load(self.handshake_network(arch, [(0, 2)]))
        expected = 3 * (2e-9 + 2 * 3e-9)
        self.assertAlmostEqual(
            chip.sim(3, timing_model="simple")["energy"]["total"], expected, delta=1e-20
        )
        chip.reset()
        self.assertAlmostEqual(
            chip.sim(3, timing_model="simple")["energy"]["total"], expected, delta=1e-20
        )

    def test_handshake_rejects_invalid_energy(self):
        for attribute in ("handshake", "handshake_base", "sync"):
            for cost in (-1.0, float("nan"), float("inf")):
                with (
                    self.subTest(attribute=attribute, cost=cost),
                    self.assertRaises(ValueError),
                ):
                    kwargs = {attribute: cost}
                    if attribute == "sync":
                        kwargs["handshake"] = 1e-9
                    sanafe.SpikingChip(self.architecture(**kwargs))


if __name__ == "__main__":
    unittest.main()
