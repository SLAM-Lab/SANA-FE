"""Regression tests for NeuroScale synapse and chip synchronization energy."""

import tempfile
import unittest
from pathlib import Path

import yaml

import sanafe

SANAFE_ROOT = Path(__file__).resolve().parents[2]


class EnergyAccounting(unittest.TestCase):
    def chip(self, sync=None, read=0.0, scan=0.0, spiking=False, cores=1):
        payload = yaml.safe_load((SANAFE_ROOT / "arch/neuroscale.yaml").read_text())
        attrs = payload["architecture"]["attributes"]
        attrs.pop("energy_sync", None)
        if sync is not None:
            attrs["energy_sync"] = sync
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
            energy_read_synapse=read, energy_scan_axon=scan
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "arch.yaml"
            path.write_text(yaml.safe_dump(payload))
            arch = sanafe.load_arch(str(path))
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
                    self.assertAlmostEqual(
                        energy["network"], energy["total"], delta=1e-20
                    )

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


if __name__ == "__main__":
    unittest.main()
