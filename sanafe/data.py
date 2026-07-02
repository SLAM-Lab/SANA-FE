# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
Convert SANA-FE trace outputs into useful formats.

Each function accepts whatever type is most convenient:
    - a path or Path to a CSV file produced by chip.sim()
    - the dict returned by chip.sim()
    - the raw in-memory value held under the matching key
    - a DataFrame (returned as-is)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

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


def spikes_to_raster(
    source: Any,
    groups: Sequence[str] | None = None,
    time_range: tuple[int, int] | None = None,
    n_timesteps: int | None = None,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Convert a spike trace into a dense 2D raster matrix.

    Args:
        source: The spike trace data. May be a :class:`pandas.DataFrame`, a
            path to a CSV file, the dict returned by ``chip.sim()``, or the raw
            ``spike_trace`` values i.e., a per-timestep list, where each element
            is a list of spiking neurons for that timestep.
        groups: Optional subset of neuron groups to include. Row order in the
            returned matrix follows the order given here; when omitted, all
            groups are included in sorted order.
        time_range (tuple(int, int), optional): Optional ``(start, stop)``
            half-open interval of timesteps to include. When omitted, spans from
            the minimum to maximum timestep present in the data.
        n_timesteps (int, optional): Optional explicit number of columns in the
            matrix. Useful when the trace ends before the simulation does (e.g.,
            no spikes in the final timesteps) and you want the matrix to span
            the full run. Ignored when ``time_range`` is given.

    Returns:
        tuple: A 3-tuple ``(matrix, neuron_ids, timesteps)`` where

        * ``matrix`` is a 2D boolean :class:`numpy.ndarray` of shape
          ``(n_neurons, n_timesteps)`` with ``True`` where a spike occurred.
        * ``neuron_ids`` is a list of ``"group.offset"`` strings labelling
          each row of the matrix.
        * ``timesteps`` is a 1D :class:`numpy.ndarray` giving the timestep
          index for each column.

    Raises:
        ValueError: If no spike trace data can be found in ``source``, or if
            ``groups`` contains unknown group names.
    """
    df = spikes_to_dataframe(source)

    all_groups = sorted(df["group"].unique())
    if groups is None:
        groups = all_groups
    else:
        unknown = set(groups) - set(all_groups)
        if unknown:
            raise ValueError(
                f"Unknown groups: {unknown}. Available: {all_groups}")
        df = df[df["group"].isin(groups)]

    # Determine time axis
    if time_range is not None:
        t_start, t_stop = time_range
        df = df[(df["timestep"] >= t_start) & (df["timestep"] < t_stop)]
    else:
        t_start = int(df["timestep"].min()) if len(df) else 0
        if n_timesteps is not None:
            t_stop = t_start + n_timesteps
        else:
            t_stop = int(df["timestep"].max()) + 1 if len(df) else t_start + 1

    timesteps = np.arange(t_start, t_stop)

    # Determine row order: stable per-group ordering by neuron_offset
    neuron_ids: list[str] = []
    row_of: dict[str, int] = {}
    for g in groups:
        offsets = sorted(df.loc[df["group"] == g, "neuron_offset"].unique())
        for off in offsets:
            nid = f"{g}.{int(off)}"
            row_of[nid] = len(neuron_ids)
            neuron_ids.append(nid)

    matrix = np.zeros((len(neuron_ids), len(timesteps)), dtype=bool)
    if len(df) and neuron_ids and timesteps:
        rows = df["neuron_id"].map(row_of).to_numpy()
        cols = df["timestep"].to_numpy() - t_start
        valid = (rows >= 0) & (cols >= 0) & (cols < len(timesteps))
        matrix[rows[valid], cols[valid]] = True

    return matrix, neuron_ids, timesteps


def spikes_to_dataframe(source: Any) -> pd.DataFrame:
    """Convert a spike trace into a pandas DataFrame.

    Args:
        source (pandas.DataFrame or Path or str or dict or list[list[NeuronAddress]]):
           The spike trace data. May be a :class:`pandas.DataFrame`, a path to a
           CSV file, the dict returned by ``chip.sim()``, or the raw
           ``spike_trace`` value. In memory, spike traces are a per-timestep
           list, with each entry a list of :class:`NeuronAddress` of neurons
           that fired that timestep.

    Returns:
        pandas.DataFrame: A DataFrame with columns ``timestep``, ``group``,
        ``neuron_offset``, and ``neuron_id``.

    Raises:
        ValueError: If no spike trace data can be found in ``source``.
    """
    if isinstance(source, pd.DataFrame):
        return source

    if _is_path(source):
        df = pd.read_csv(source, dtype={"neuron": str})
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
    """Convert a membrane-potential trace into a pandas DataFrame indexed by timestep.

    Args:
        source (pandas.DataFrame or str or Path or dict or list[list[float]]):
            The potential trace data. May be a :class:`pandas.DataFrame`, a
            path to a CSV file, the dict returned by ``chip.sim()``, or a 2D
            array of ``potential_trace`` values (i.e., a per-timestep list,
            containing per-neuron lists of potentials).
        neuron_ids (List[str]): Optional column labels for the neurons. When
            omitted, columns are named ``Neuron 0``, ``Neuron 1``, and so on.
            Only used when ``source`` is in-memory trace data.

    Returns:
        pandas.DataFrame: A DataFrame whose index is named ``timestep`` and
        whose columns correspond to individual neurons.

    Raises:
        ValueError: If no potential trace data is found, or if ``neuron_ids``
            length does not match the number of trace columns.
    """
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


def neuron_traces_to_dataframe(source: Any,
                              neuron_ids: Sequence[str] | None = None) -> pd.DataFrame:
    """Convert additional neuron traces into a pandas DataFrame indexed by timestep.

    Args:
        source (pandas.DataFrame or str or Path or dict or dict[str, list[list]]):
            The neuron trace data. May be a :class:`pandas.DataFrame`, a
            path to a CSV file, the dict returned by ``chip.sim()``, or a 2D
            array of neuron_traces (i.e., a per-timestep list, containing
            per-neuron lists of values).
        neuron_ids (List[str]): Optional column labels for the neurons. When
            omitted, columns are named ``Neuron 0``, ``Neuron 1``, and so on.
            Only used when ``source`` is in-memory trace data.

    Returns:
        pandas.DataFrame: A DataFrame whose index is named ``timestep`` and
        whose columns correspond to individual neurons.

    Raises:
        ValueError: If no potential trace data is found, or if ``neuron_ids``
            length does not match the number of trace columns.
    """
    if isinstance(source, pd.DataFrame):
        df = source.copy()
        index_cols = [c for c in ("timestep", "neuron_id") if c in df.columns]
        if index_cols:
            df = df.set_index(index_cols)
        return df

    if _is_path(source):
        # Drop pandas' phantom column from the trailing-comma rows.
        raw = pd.read_csv(source)


        if "timestep" not in raw.columns:
            raise ValueError(
                f"Neuron trace CSV {source!r} is missing a 'timestep' column")

        long_rows: list[dict] = []
        for col in raw.columns:
            if col == "timestep" or col.startswith("Unnamed"):
                continue
            # Parse "neuron <group>.<offset>/<variable>"
            stripped = col[len("neuron "):] if col.startswith("neuron ") else col
            if "/" not in stripped:
                raise ValueError(
                    f"Unrecognised neuron trace column {col!r}; "
                    "expected 'neuron <group>.<offset>/<variable>'")
            neuron_id, var = stripped.rsplit("/", 1)
            for t, v in zip(raw["timestep"], raw[col]):
                long_rows.append(
                    {"timestep": t, "neuron_id": neuron_id, "var": var, "value": v})

        df = pd.DataFrame(long_rows)
        # Pivot to (timestep, neuron_id) index, one column per variable.
        return df.pivot_table(
            index=["timestep", "neuron_id"], columns="var", values="value")

    traces = _maybe_unwrap(source, "neuron_trace")
    if traces is None or (isinstance(traces, dict) and not traces):
        raise ValueError("No neuron trace data in source")
    if not isinstance(traces, dict):
        raise ValueError(
            f"Expected neuron_trace to be a dict[str, list[list]], got {type(traces).__name__}")

    # Stack each variable into a (T, N) array and check shapes agree.
    arrays: dict[str, np.ndarray] = {}
    n_timesteps: int | None = None
    n_neurons: int | None = None
    for var, values in traces.items():
        arr = np.asarray(values, dtype=float)
        if arr.ndim != 2:
            raise ValueError(
                f"neuron_trace[{var!r}] is not 2D (got shape {arr.shape})")
        if n_timesteps is None:
            n_timesteps, n_neurons = arr.shape
        elif arr.shape != (n_timesteps, n_neurons):
            raise ValueError(
                f"neuron_trace[{var!r}] has shape {arr.shape}, "
                f"expected ({n_timesteps}, {n_neurons})")
        arrays[var] = arr

    if n_neurons == 0:
        raise ValueError("Neuron traces contain no neurons")

    if neuron_ids is None:
        neuron_ids = [f"Neuron {i}" for i in range(n_neurons)]
    elif len(neuron_ids) != n_neurons:
        raise ValueError(
            f"neuron_ids has {len(neuron_ids)} entries but trace has {n_neurons} columns")

    index = pd.MultiIndex.from_product(
        [range(n_timesteps), list(neuron_ids)],
        names=["timestep", "neuron_id"],
    )
    return pd.DataFrame(
        {var: arr.reshape(-1) for var, arr in arrays.items()},
        index=index,
    )


def performance_to_dataframe(source: Any) -> pd.DataFrame:
    """Convert a performance trace into a DataFrame.

    Args:
        source (pandas.DataFrame or Path or str or dict:
            The performance trace data. May be a :class:`pandas.DataFrame`,
            a path to a CSV file, the dict from chip.sim() or a dict containing
            performance data (dict[str, list], i.e. lists of per-timestep values
            for each metric).

    Returns:
        pandas.DataFrame: A DataFrame of per-timestep performance metrics.

    Raises:
        ValueError: If no performance trace data can be found in ``source``.
    """
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


def messages_to_dataframe(source: Any) -> pd.DataFrame:
    """Convert a message trace into a flat pandas DataFrame, one row per message.

    In addition, convenience identifier columns are derived when their
    components are present:

    * ``src_neuron`` from ``src_neuron_group_id`` and ``src_neuron_offset``
    * ``src_hw`` from ``src_tile_id`` and ``src_core_offset``
    * ``dest_hw`` from ``dest_tile_id`` and ``dest_core_offset``

    Args:
        source (pandas.DataFrame or Path or str or dict or list[list[dict]]):
            The message trace data. May be a :class:`pandas.DataFrame`, a path
            to a CSV file, the dict returned by ``chip.sim()``, or the raw
            ``message_trace`` values.

    Returns:
        pandas.DataFrame: A DataFrame with one row per message.

    Raises:
        ValueError: If no message trace data can be found in ``source``.
    """
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
