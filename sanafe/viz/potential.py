"""
Membrane potential plots. Default rendering is a heatmap (consistent
with the raster view for spikes); plot_potential_lines is available
for traditional per-neuron line plots.
"""

from __future__ import annotations

from typing import Any, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from sanafe.data import potentials_to_dataframe
from sanafe.viz.styles import (
    SANAFEStyle,
    DEFAULT_COLORS,
    create_figure,
    get_default_style,
    style_axis,
)


def plot_potential(
    source: Any,
    neuron_ids: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_colorbar: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: str = "Neuron",
    **imshow_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot neuron potential time series as heatmaps."""
    style = style or get_default_style()
    df = potentials_to_dataframe(source, neuron_ids=neuron_ids)

    if time_range is not None:
        df = df.iloc[time_range[0]:time_range[1]]

    _, n_neurons = df.shape
    timesteps = np.asarray(df.index)
    t_start = int(timesteps[0]) if len(timesteps) else 0
    t_end = int(timesteps[-1]) + 1 if len(timesteps) else 1

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    imshow_defaults = {
        "aspect": "auto",
        "origin": "lower",
        "cmap": cmap,
        "extent": [t_start - 0.5, t_end - 0.5, -0.5, n_neurons - 0.5],
    }
    if vmin is not None:
        imshow_defaults["vmin"] = vmin
    if vmax is not None:
        imshow_defaults["vmax"] = vmax
    imshow_defaults.update(imshow_kwargs)

    im = ax.imshow(df.values.T, **imshow_defaults)

    if show_colorbar:
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Membrane Potential")

    labels = list(df.columns)
    if n_neurons <= 20:
        ax.set_yticks(range(n_neurons))
        ax.set_yticklabels(labels)
    else:
        step = max(1, n_neurons // 10)
        tick_positions = list(range(0, n_neurons, step))
        ax.set_yticks(tick_positions)
        ax.set_yticklabels([labels[i] for i in tick_positions])

    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)
    if style.tight_layout:
        fig.tight_layout()
    return fig, ax


def plot_potential_lines(
    source: Any,
    neuron_ids: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    colors: Optional[Sequence[str]] = None,
    show_threshold: Optional[float] = None,
    threshold_color: str = "#d62728",
    threshold_linestyle: str = "--",
    show_legend: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: str = "Membrane Potential",
    **plot_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot potential time-series as line graphs."""
    style = style or get_default_style()
    df = potentials_to_dataframe(source, neuron_ids=neuron_ids)

    if time_range is not None:
        df = df.iloc[time_range[0]:time_range[1]]

    n_neurons = df.shape[1]
    if colors is None:
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(n_neurons)]

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    plot_defaults = {"linewidth": style.potential_line_width}
    if style.potential_marker:
        plot_defaults["marker"] = style.potential_marker
        plot_defaults["markersize"] = style.potential_marker_size
    plot_defaults.update(plot_kwargs)

    timesteps = np.asarray(df.index)
    for i, col in enumerate(df.columns):
        ax.plot(timesteps, df[col].values, color=colors[i],
                label=str(col), **plot_defaults)

    if show_threshold is not None:
        ax.axhline(y=show_threshold, color=threshold_color,
                   linestyle=threshold_linestyle,
                   linewidth=style.line_width * 0.8, label="Threshold", zorder=0)

    if len(timesteps):
        ax.set_xlim(timesteps[0] - 0.5, timesteps[-1] + 0.5)

    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)
    if show_legend:
        ax.legend(loc="upper right", framealpha=0.9)
    if style.tight_layout:
        fig.tight_layout()
    return fig, ax
