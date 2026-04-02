"""
Convert SANA-FE trace outputs into pandas DataFrames.

Each function accepts whatever shape is most convenient:
    - a path or Path to a CSV produced by chip.sim()
    - the dict returned by chip.sim()
    - the raw in-memory value held under the matching key
    - a DataFrame (returned as-is)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd


def _is_path(source: Any) -> bool:
    return isinstance(source, (str, Path))


def _maybe_unwrap(source: Any, key: str) -> Any:
    if isinstance(source, dict) and key in source and not _looks_like_perf_dict(source):
        return source[key]
    return source


def _looks_like_perf_dict(source: Any) -> bool:
    if not isinstance(source, dict):
        return False
    return "total_energy" in source and "timestep" in source and isinstance(
        source.get("timestep"), (list, tuple, np.ndarray))


def spikes_to_dataframe(source: Any) -> pd.DataFrame:
    if isinstance(source, pd.DataFrame):
        return source

    if _is_path(source):
        df = pd.read_csv(source)
        split = df["neuron"].str.rsplit(".", n=1, expand=True)
        df["group"] = split[0]
        df["neuron_offset"] = split[1].astype(int)
        df["neuron_id"] = df["neuron"]
        return df[["timestep", "group", "neuron_offset", "neuron_id"]]

    spikes = _maybe_unwrap(source, "spike_trace")
    if spikes is None:
        raise ValueError("No spike trace data in source")

    records = []
    for t, row in enumerate(spikes):
        for s in row:
            group = getattr(s, "group_name", None) or s["group_name"]
            offset = getattr(s, "neuron_offset", None)
            if offset is None:
                offset = s["neuron_offset"]
            records.append({
                "timestep": t,
                "group": group,
                "neuron_offset": offset,
                "neuron_id": f"{group}.{offset}",
            })
    return pd.DataFrame(records, columns=["timestep", "group",
                                          "neuron_offset", "neuron_id"])


def potentials_to_dataframe(source: Any,
                            neuron_ids: Sequence[str] | None = None) -> pd.DataFrame:
    if isinstance(source, pd.DataFrame):
        df = source.copy()
        if "timestep" in df.columns:
            df = df.set_index("timestep")
        return df

    if _is_path(source):
        df = pd.read_csv(source)
        if "timestep" in df.columns:
            df = df.set_index("timestep")
        df.columns = [c.replace("neuron ", "", 1) for c in df.columns]
        return df

    potentials = _maybe_unwrap(source, "potential_trace")
    if potentials is None or len(potentials) == 0:
        raise ValueError("No potential trace data in source")

    arr = np.asarray(potentials, dtype=float)
    n_neurons = arr.shape[1] if arr.ndim == 2 else 0
    if neuron_ids is None:
        neuron_ids = [f"Neuron {i}" for i in range(n_neurons)]
    elif len(neuron_ids) != n_neurons:
        raise ValueError(
            f"neuron_ids has {len(neuron_ids)} entries but trace has {n_neurons} columns")

    df = pd.DataFrame(arr, columns=list(neuron_ids))
    df.index.name = "timestep"
    return df


def performance_to_dataframe(source: Any) -> pd.DataFrame:
    if isinstance(source, pd.DataFrame):
        return source

    if _is_path(source):
        return pd.read_csv(source)

    if _looks_like_perf_dict(source):
        return pd.DataFrame(source)

    perf = _maybe_unwrap(source, "perf_trace")
    if perf is None or (isinstance(perf, dict) and not perf):
        raise ValueError("No performance trace data in source")

    return pd.DataFrame(perf)


_MESSAGE_RENAME = {
    "receive_delay": "processing_delay",
    "blocked_delay": "blocking_delay",
    "sent_timestamp": "send_timestamp",
}


def messages_to_dataframe(source: Any) -> pd.DataFrame:
    if isinstance(source, pd.DataFrame):
        return source

    if _is_path(source):
        return pd.read_csv(source)

    messages = _maybe_unwrap(source, "message_trace")
    if messages is None:
        raise ValueError("No message trace data in source")

    rows: list[dict] = []
    for ts_msgs in messages:
        for m in ts_msgs:
            row = dict(m)
            for old, new in _MESSAGE_RENAME.items():
                if old in row:
                    row[new] = row.pop(old)
            if "src_neuron_group_id" in row and "src_neuron_offset" in row:
                row.setdefault("src_neuron",
                               f"{row['src_neuron_group_id']}.{row['src_neuron_offset']}")
            if "src_tile_id" in row and "src_core_offset" in row:
                row.setdefault("src_hw",
                               f"{row['src_tile_id']}.{row['src_core_offset']}")
            if "dest_tile_id" in row and "dest_core_offset" in row:
                row.setdefault("dest_hw",
                               f"{row['dest_tile_id']}.{row['dest_core_offset']}")
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import sanafe

    arch = sanafe.load_arch("tutorial/arch.yaml")
    snn = sanafe.Network()
    g = snn.create_neuron_group(
        "in", 2, {"bias": 0.5, "threshold": 1.0, "reset": 0.0},
        log_spikes=True, log_potential=True)
    for n in g:
        n.map_to_core(arch.tiles[0].cores[0])
    chip = sanafe.SpikingChip(arch)
    chip.load(snn)
    r = chip.sim(5, spike_trace=True, potential_trace=True,
                 perf_trace=True, message_trace=True)

    spikes = spikes_to_dataframe(r)
    pots = potentials_to_dataframe(r)
    perf = performance_to_dataframe(r)
    msgs = messages_to_dataframe(r)

    assert list(spikes.columns) == ["timestep", "group", "neuron_offset", "neuron_id"]
    assert pots.index.name == "timestep"
    assert "total_energy" in perf.columns
    assert "processing_delay" in msgs.columns
    assert "receive_delay" not in msgs.columns
    print("OK", len(spikes), len(pots), len(perf), len(msgs))
