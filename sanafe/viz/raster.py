"""
Spike raster plot. Accepts whatever chip.sim() returned, a CSV path,
or a pre-built DataFrame.
"""

from __future__ import annotations

from typing import Any, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from sanafe.data import spikes_to_dataframe
from sanafe.viz.styles import (
    SANAFEStyle,
    create_figure,
    get_default_style,
    get_group_colors,
    style_axis,
)


def plot_raster(
    source: Any,
    groups: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    colors: Optional[dict] = None,
    group_spacing: float = 0.5,
    show_legend: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: str = "Neuron",
    **scatter_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot spike raster."""
    style = style or get_default_style()
    df = spikes_to_dataframe(source)

    if time_range is not None:
        df = df[(df["timestep"] >= time_range[0])
                & (df["timestep"] < time_range[1])]

    all_groups = sorted(df["group"].unique())
    if groups is None:
        groups = all_groups
    else:
        unknown = set(groups) - set(all_groups)
        if unknown:
            raise ValueError(
                f"Unknown groups: {unknown}. Available: {all_groups}")
        df = df[df["group"].isin(groups)]

    if colors is None:
        colors = get_group_colors(list(groups), style)

    y_of: dict[Tuple[str, int], float] = {}
    y_ticks: list[float] = []
    y_tick_labels: list[str] = []
    y_cursor = 0.0
    for g in groups:
        offsets = sorted(df.loc[df["group"] == g, "neuron_offset"].unique())
        if not offsets:
            continue
        start_y = y_cursor
        for off in offsets:
            y_of[(g, int(off))] = y_cursor
            y_cursor += 1.0
        y_ticks.append((start_y + y_cursor - 1.0) / 2)
        y_tick_labels.append(g)
        y_cursor += group_spacing

    xs = df["timestep"].values
    ys = np.array([y_of[(row.group, int(row.neuron_offset))]
                   for row in df.itertuples()])
    cs = [colors.get(g, "#333333") for g in df["group"]]

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    scatter_defaults = {
        "marker": style.raster_marker,
        "s": style.raster_marker_size,
        "linewidths": style.raster_line_width,
    }
    scatter_defaults.update(scatter_kwargs)

    if len(xs):
        ax.scatter(xs, ys, c=cs, **scatter_defaults)

    if y_ticks:
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)
    else:
        ax.set_ylabel(ylabel)

    t_max = df["timestep"].max() if len(df) else 0
    ax.set_xlim(-0.5, float(t_max) + 0.5)
    ax.set_ylim(-0.5, max(y_cursor - group_spacing, 0.5) + 0.5)

    style_axis(ax, style, xlabel=xlabel, title=title)

    if show_legend and len(groups) > 1:
        handles = [
            Line2D([0], [0], marker=style.raster_marker, color="w",
                   markerfacecolor=colors.get(g, "#333333"),
                   markeredgecolor=colors.get(g, "#333333"),
                   markersize=10, label=g)
            for g in groups
        ]
        ax.legend(handles=handles, loc="upper right", framealpha=0.9)

    if style.tight_layout:
        fig.tight_layout()

    return fig, ax


def raster_plot_matrix(
    spike_matrix: np.ndarray,
    neuron_ids: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: str = "Neuron",
    color: str = "#1f77b4",
    **scatter_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot spike raster given matrix of spikes."""
    style = style or get_default_style()

    if time_range is not None:
        spike_matrix = spike_matrix[time_range[0]:time_range[1], :]

    n_timesteps, n_neurons = spike_matrix.shape
    spike_times, spike_neurons = np.where(spike_matrix)

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    scatter_defaults = {
        "marker": style.raster_marker,
        "s": style.raster_marker_size,
        "linewidths": style.raster_line_width,
        "c": color,
    }
    scatter_defaults.update(scatter_kwargs)

    if len(spike_times) > 0:
        ax.scatter(spike_times, spike_neurons, **scatter_defaults)

    if neuron_ids is not None:
        if n_neurons > 20:
            step = max(1, n_neurons // 10)
            tick_positions = list(range(0, n_neurons, step))
            ax.set_yticks(tick_positions)
            ax.set_yticklabels([neuron_ids[i] for i in tick_positions])
        else:
            ax.set_yticks(range(n_neurons))
            ax.set_yticklabels(neuron_ids)

    ax.set_xlim(-0.5, n_timesteps + 0.5)
    ax.set_ylim(-0.5, n_neurons - 0.5)

    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)

    if style.tight_layout:
        fig.tight_layout()

    return fig, ax
