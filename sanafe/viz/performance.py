"""
Hardware performance plots: energy breakdown, throughput, latency.
All accept sim() results, a CSV path, a raw dict, or a DataFrame.
"""

from __future__ import annotations

from typing import Any, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from sanafe.data import (
    performance_to_dataframe,
    messages_to_dataframe,
)
from sanafe.viz.styles import (
    SANAFEStyle,
    DEFAULT_COLORS,
    create_figure,
    get_default_style,
    style_axis,
)

_ENERGY_COLUMNS = [
    "synapse_energy",
    "dendrite_energy",
    "soma_energy",
    "network_energy",
]

_LATENCY_COLUMNS = [
    "generation_delay",
    "processing_delay",
    "network_delay",
    "blocking_delay",
]

_ENERGY_UNITS = [(1e-15, "fJ"), (1e-12, "pJ"), (1e-9, "nJ"),
                 (1e-6, "µJ"), (1e-3, "mJ"), (1.0, "J")]
_TIME_UNITS = [(1e-12, "ps"), (1e-9, "ns"), (1e-6, "µs"),
               (1e-3, "ms"), (1.0, "s")]


def _auto_scale(values: np.ndarray, units) -> Tuple[float, str]:
    peak = float(np.nanmax(np.abs(values))) if len(values) else 0.0
    if peak == 0.0:
        return 1.0, units[-1][1]
    for factor, label in units:
        if peak / factor < 1000:
            return factor, label
    return units[-1][0], units[-1][1]


def plot_energy(
    source: Any,
    time_range: Optional[Tuple[int, int]] = None,
    mode: str = "stacked_area",
    components: Optional[Sequence[str]] = None,
    normalize: bool = False,
    show_total: bool = False,
    colors: Optional[Sequence[str]] = None,
    component_labels: Optional[Sequence[str]] = None,
    show_legend: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: Optional[str] = None,
    **kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot energy time series."""
    style = style or get_default_style()
    df = performance_to_dataframe(source)

    if components is None:
        components = list(_ENERGY_COLUMNS)
    missing = [c for c in components if c not in df.columns]
    if missing:
        raise ValueError(
            f"Columns not found in perf trace: {missing}. "
            f"Available: {list(df.columns)}")

    if time_range is not None:
        if "timestep" in df.columns:
            df = df[(df["timestep"] >= time_range[0])
                    & (df["timestep"] < time_range[1])]
        else:
            df = df.iloc[time_range[0]:time_range[1]]

    timesteps = (df["timestep"].values if "timestep" in df.columns
                 else np.arange(len(df)))

    raw = [df[c].values.astype(float) for c in components]
    factor, unit = _auto_scale(np.stack(raw).sum(axis=0), _ENERGY_UNITS)
    scaled = [arr / factor for arr in raw]

    if normalize:
        totals = sum(scaled)
        with np.errstate(invalid="ignore", divide="ignore"):
            scaled = [np.where(totals > 0, a / totals * 100.0, 0.0)
                      for a in scaled]
        unit = "%"

    if colors is None:
        colors = style.energy_colors[:len(components)]
    if component_labels is None:
        if len(components) <= len(style.energy_component_names):
            component_labels = style.energy_component_names[:len(components)]
        else:
            component_labels = [c.replace("_energy", "").replace("_", " ").title()
                                for c in components]

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    if mode == "stacked_area":
        bottom = np.zeros_like(scaled[0])
        for vals, label, color in zip(scaled, component_labels, colors):
            ax.fill_between(timesteps, bottom, bottom + vals,
                            label=label, color=color,
                            alpha=style.perf_fill_alpha + 0.4, **kwargs)
            bottom = bottom + vals
        if show_total and not normalize and "total_energy" in df.columns:
            ax.plot(timesteps, df["total_energy"].values / factor,
                    color="black", linewidth=style.perf_line_width,
                    linestyle="--", label="Total", zorder=5)
    elif mode == "stacked_bar":
        bottoms = np.zeros(len(timesteps))
        for vals, label, color in zip(scaled, component_labels, colors):
            ax.bar(timesteps, vals, bottom=bottoms, width=0.8,
                   label=label, color=color,
                   edgecolor=style.hist_edgecolor,
                   linewidth=style.hist_edgewidth, **kwargs)
            bottoms += vals
    elif mode == "bar":
        n = len(components)
        bw = 0.8 / n
        offsets = np.linspace(-(n - 1) * bw / 2, (n - 1) * bw / 2, n)
        for off, vals, label, color in zip(offsets, scaled, component_labels, colors):
            ax.bar(timesteps + off, vals, width=bw, label=label, color=color,
                   edgecolor=style.hist_edgecolor,
                   linewidth=style.hist_edgewidth, **kwargs)
    else:
        raise ValueError(
            f"Unknown mode '{mode}'. Use 'stacked_area', 'stacked_bar', or 'bar'.")

    if ylabel is None:
        ylabel = "Energy (%)" if normalize else f"Energy ({unit})"
    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)

    if show_legend:
        ax.legend(loc="upper right", framealpha=0.9)
    if style.tight_layout:
        fig.tight_layout()
    return fig, ax


def plot_throughput(
    source: Any,
    metrics: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    colors: Optional[Sequence[str]] = None,
    labels: Optional[Sequence[str]] = None,
    secondary_y: Optional[Sequence[str]] = None,
    show_legend: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: Optional[str] = None,
    **plot_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot throughput time series."""
    style = style or get_default_style()
    df = performance_to_dataframe(source)

    if metrics is None:
        metrics = [m for m in ("fired", "spikes", "hops") if m in df.columns]
    missing = [m for m in metrics if m not in df.columns]
    if missing:
        raise ValueError(
            f"Columns not found in perf trace: {missing}. "
            f"Available: {list(df.columns)}")

    secondary = set(secondary_y or [])

    if time_range is not None:
        if "timestep" in df.columns:
            df = df[(df["timestep"] >= time_range[0])
                    & (df["timestep"] < time_range[1])]
        else:
            df = df.iloc[time_range[0]:time_range[1]]

    timesteps = (df["timestep"].values if "timestep" in df.columns
                 else np.arange(len(df)))

    if colors is None:
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(len(metrics))]
    if labels is None:
        labels = [m.replace("_", " ").title() for m in metrics]

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    ax_right = None
    if secondary:
        ax_right = ax.twinx()
        ax.right_ax = ax_right

    plot_defaults = {"linewidth": style.perf_line_width}
    if style.perf_marker:
        plot_defaults["marker"] = style.perf_marker
        plot_defaults["markersize"] = style.perf_marker_size
    plot_defaults.update(plot_kwargs)

    left_lines, right_lines = [], []
    for metric, color, label in zip(metrics, colors, labels):
        values = df[metric].values.astype(float)
        scaled_label = label
        if "energy" in metric:
            factor, unit = _auto_scale(values, _ENERGY_UNITS)
            values = values / factor
            scaled_label = f"{label} ({unit})"
        elif "sim_time" in metric or "delay" in metric:
            factor, unit = _auto_scale(values, _TIME_UNITS)
            values = values / factor
            scaled_label = f"{label} ({unit})"

        target = ax_right if metric in secondary else ax
        (line,) = target.plot(timesteps, values, color=color,
                              label=scaled_label, **plot_defaults)
        (right_lines if metric in secondary else left_lines).append(line)

    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel or "Count", title=title)

    if ax_right is not None:
        ax_right.tick_params(axis="y", labelsize=style.tick_size)
        if show_legend:
            all_lines = left_lines + right_lines
            ax.legend(all_lines, [l.get_label() for l in all_lines],
                      loc="upper right", framealpha=0.9)
    elif show_legend:
        ax.legend(loc="upper right", framealpha=0.9)

    if style.tight_layout:
        fig.tight_layout()
    return fig, ax


def plot_message_latency(
    source: Any,
    metrics: Optional[Sequence[str]] = None,
    filter_placeholder: bool = True,
    time_range: Optional[Tuple[int, int]] = None,
    bins: Optional[int] = None,
    colors: Optional[Sequence[str]] = None,
    labels: Optional[Sequence[str]] = None,
    log_scale: bool = False,
    show_legend: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: str = "Count",
    **hist_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot spike message (packet) metrics."""
    style = style or get_default_style()
    df = messages_to_dataframe(source)

    if metrics is None:
        metrics = list(_LATENCY_COLUMNS)
    if isinstance(metrics, str):
        metrics = [metrics]
    missing = [m for m in metrics if m not in df.columns]
    if missing:
        raise ValueError(
            f"Columns not found in message trace: {missing}. "
            f"Available: {list(df.columns)}")

    if filter_placeholder and "mid" in df.columns:
        df = df[df["mid"] >= 0]

    if time_range is not None and "timestep" in df.columns:
        df = df[(df["timestep"] >= time_range[0])
                & (df["timestep"] < time_range[1])]

    all_vals = np.concatenate([df[m].values.astype(float) for m in metrics])
    all_vals = all_vals[np.isfinite(all_vals)]

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    if len(all_vals) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=style.label_size, color="#999")
        return fig, ax

    factor, unit = _auto_scale(all_vals, _TIME_UNITS)

    if colors is None:
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(len(metrics))]
    if labels is None:
        labels = [m.replace("_", " ").title() for m in metrics]
    if bins is None:
        bins = style.hist_bins

    hist_defaults = {
        "bins": bins,
        "alpha": style.hist_alpha * 0.85,
        "edgecolor": style.hist_edgecolor,
        "linewidth": style.hist_edgewidth,
    }
    hist_defaults.update(hist_kwargs)

    for metric, color, label in zip(metrics, colors, labels):
        vals = df[metric].values.astype(float)
        vals = vals[np.isfinite(vals)] / factor
        if len(vals) > 0:
            ax.hist(vals, color=color, label=label, **hist_defaults)

    if log_scale:
        ax.set_yscale("log")
    if xlabel is None:
        xlabel = f"Delay ({unit})"

    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel,
               title=title or ("Delay Distribution" if len(metrics) > 1 else None))

    if show_legend and len(metrics) > 1:
        ax.legend(loc="upper right", framealpha=0.9)
    if style.tight_layout:
        fig.tight_layout()
    return fig, ax
