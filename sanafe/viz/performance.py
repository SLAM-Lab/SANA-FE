"""
This module provides functions for visualizing hardware performance metrics
and message-level latency data from SANA-FE simulations. Includes energy
breakdown plots, throughput timeseries, and latency distribution histograms.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

from sanafe.data.traces import TraceData
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


def _coerce_to_trace_data(data: Union[TraceData, Dict[str, Any]]) -> TraceData:
    """Accept raw results dict or TraceData, always return TraceData."""
    if isinstance(data, dict):
        return TraceData.from_sim_results(data)
    return data


def _get_perf_df(data: TraceData) -> pd.DataFrame:
    """Extract perf DataFrame or raise a clear error."""
    if not data.has_performance():
        raise ValueError(
            "No performance trace data available. "
            "Pass perf_trace to chip.sim() or load a perf CSV."
        )
    return data.performance_to_dataframe()


def _get_message_df(data: TraceData) -> pd.DataFrame:
    """Extract message DataFrame or raise a clear error."""
    if not data.has_messages():
        raise ValueError(
            "No message trace data available. "
            "Pass message_trace to chip.sim() or load a message CSV."
        )
    return data.messages_to_dataframe()


def _auto_energy_scale(
    values: np.ndarray,
) -> Tuple[np.ndarray, str]:
    """
    Pick a human-readable energy unit and return scaled values + unit label.

    Chooses the unit so that peak values fall in roughly 0.1 - 999 range.
    """
    peak = np.nanmax(values) if len(values) > 0 else 0.0
    if peak == 0.0:
        return values, "J"

    scales = [
        (1e-15, "fJ"),
        (1e-12, "pJ"),
        (1e-9,  "nJ"),
        (1e-6,  "µJ"),
        (1e-3,  "mJ"),
        (1.0,   "J"),
    ]
    for factor, label in scales:
        if peak / factor < 1000:
            return values / factor, label

    return values, "J"


def _auto_time_scale(
    values: np.ndarray,
) -> Tuple[np.ndarray, str]:
    """Pick a human-readable time unit for delay values."""
    peak = np.nanmax(np.abs(values)) if len(values) > 0 else 0.0
    if peak == 0.0:
        return values, "s"

    scales = [
        (1e-12, "ps"),
        (1e-9,  "ns"),
        (1e-6,  "µs"),
        (1e-3,  "ms"),
        (1.0,   "s"),
    ]
    for factor, label in scales:
        if peak / factor < 1000:
            return values / factor, label

    return values, "s"


def energy_breakdown_plot(
    data: Union[TraceData, Dict[str, Any]],
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
    """
    Plot energy consumption breakdown over time from the performance trace.

    Visualizes per-component energy (synapse, dendrite, soma, network) at each
    timestep. Supports stacked area, stacked bar, and grouped bar modes.

    Args:
        data: TraceData or raw results dict from chip.sim().
        time_range: (start, end) timestep range to display.
        mode: Plot mode - one of:
            - "stacked_area": filled area chart with components stacked
            - "stacked_bar": stacked bar chart per timestep
            - "bar": grouped (side-by-side) bar chart
        components: Which energy columns to include. Defaults to
            ["synapse_energy", "dendrite_energy", "soma_energy",
             "network_energy"]. Can also be ["total_energy"] for a single line.
        normalize: If True, show each timestep as a percentage breakdown
            (100% stacked). Only applies to stacked modes.
        show_total: If True, overlay the total_energy line on top of
            the stacked plot.
        colors: List of colors for each component. If None, uses
            style.energy_colors.
        component_labels: Display names for each component.
        show_legend: Whether to show a legend.
        ax: Existing Axes to plot on.
        style: Style configuration.
        figsize: Figure size (width, height).
        title: Plot title.
        xlabel: X-axis label.
        ylabel: Y-axis label. Auto-generated with unit if None.
        **kwargs: Extra arguments forwarded to the underlying matplotlib call
            (fill_between for area, bar for bar chart).

    Returns: Tuple of (Figure, Axes).
    """
    data = _coerce_to_trace_data(data)
    df = _get_perf_df(data)

    if style is None:
        style = get_default_style()

    if components is None:
        components = list(_ENERGY_COLUMNS)
    missing = [c for c in components if c not in df.columns]
    if missing:
        raise ValueError(
            f"Columns not found in perf trace: {missing}. "
            f"Available: {list(df.columns)}"
        )

    if time_range is not None:
        start, end = time_range
        if "timestep" in df.columns:
            df = df[(df["timestep"] >= start) & (df["timestep"] < end)]
        else:
            df = df.iloc[start:end]

    if "timestep" in df.columns:
        timesteps = df["timestep"].values
    else:
        timesteps = np.arange(len(df))

    # Pull raw energy arrays and compute auto scale
    raw_arrays = [df[c].values.astype(float) for c in components]

    # Find a common unit from the total of all selected components
    combined_peak = np.stack(raw_arrays).sum(axis=0)
    _, unit = _auto_energy_scale(combined_peak)
    factor = {"fJ": 1e-15, "pJ": 1e-12, "nJ": 1e-9,
              "µJ": 1e-6, "mJ": 1e-3, "J": 1.0}.get(unit, 1.0)
    scaled = [arr / factor for arr in raw_arrays]

    if normalize:
        totals = sum(scaled)
        with np.errstate(invalid="ignore", divide="ignore"):
            scaled = [
                np.where(totals > 0, arr / totals * 100.0, 0.0)
                for arr in scaled
            ]
        unit = "%"

    # Resolve colors and labels
    if colors is None:
        colors = style.energy_colors[: len(components)]
    if component_labels is None:
        component_labels = (
            style.energy_component_names[: len(components)]
            if len(components) <= len(style.energy_component_names)
            else [c.replace("_energy", "").replace("_", " ").title()
                  for c in components]
        )

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    if mode == "stacked_area":
        cumulative = np.zeros_like(scaled[0])
        for i, (vals, label, color) in enumerate(
            zip(scaled, component_labels, colors)
        ):
            ax.fill_between(
                timesteps, cumulative, cumulative + vals,
                label=label, color=color,
                alpha=style.perf_fill_alpha + 0.4,
                **kwargs,
            )
            cumulative = cumulative + vals

        if show_total and not normalize:
            total_scaled = df["total_energy"].values.astype(float) / factor
            ax.plot(
                timesteps, total_scaled,
                color="black", linewidth=style.perf_line_width,
                linestyle="--", label="Total", zorder=5,
            )

    elif mode == "stacked_bar":
        bar_width = 0.8
        bottoms = np.zeros(len(timesteps))
        for vals, label, color in zip(scaled, component_labels, colors):
            ax.bar(
                timesteps, vals, bottom=bottoms,
                width=bar_width, label=label, color=color,
                edgecolor=style.hist_edgecolor,
                linewidth=style.hist_edgewidth,
                **kwargs,
            )
            bottoms += vals

    elif mode == "bar":
        n_comp = len(components)
        bar_width = 0.8 / n_comp
        offsets = np.linspace(
            -(n_comp - 1) * bar_width / 2,
            (n_comp - 1) * bar_width / 2,
            n_comp,
        )
        for offset, vals, label, color in zip(
            offsets, scaled, component_labels, colors
        ):
            ax.bar(
                timesteps + offset, vals,
                width=bar_width, label=label, color=color,
                edgecolor=style.hist_edgecolor,
                linewidth=style.hist_edgewidth,
                **kwargs,
            )
    else:
        raise ValueError(
            f"Unknown mode '{mode}'. Use 'stacked_area', 'stacked_bar', or 'bar'."
        )

    if ylabel is None:
        ylabel = f"Energy ({unit})" if not normalize else "Energy (%)"
    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)

    if show_legend:
        ax.legend(loc="upper right", framealpha=0.9)

    if style.tight_layout:
        fig.tight_layout()

    return fig, ax


def throughput_plot(
    data: Union[TraceData, Dict[str, Any]],
    metrics: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    colors: Optional[Sequence[str]] = None,
    labels: Optional[Sequence[str]] = None,
    show_legend: bool = True,
    secondary_y: Optional[Sequence[str]] = None,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: Optional[str] = None,
    **plot_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot throughput / activity metrics over time from the performance trace.

    Renders one or more per-timestep performance counters as line plots.
    Useful for visualizing how firing rate, packet count, and spike volume
    evolve across the simulation.

    Args:
        data: TraceData or raw results dict from chip.sim().
        metrics: List of column names to plot. Defaults to
            ["fired", "packets", "spikes"].
        time_range: (start, end) timestep range to display.
        colors: List of colors for each metric.
        labels: Display names for each metric. If None, derived from column
            names (underscores replaced, title-cased).
        show_legend: Whether to show a legend.
        secondary_y: List of metric names that should use a secondary y-axis
            on the right. Useful when combining counts (fired) with energy
            (total_energy) on one plot.
        ax: Existing Axes to plot on.
        style: Style configuration.
        figsize: Figure size (width, height).
        title: Plot title.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        **plot_kwargs: Extra arguments forwarded to ax.plot().

    Returns:
        Tuple of (Figure, Axes). When secondary_y is used, the secondary
        Axes is accessible via ``ax.right_ax`` (stored as an attribute).
    """
    data = _coerce_to_trace_data(data)
    df = _get_perf_df(data)

    if style is None:
        style = get_default_style()

    if metrics is None:
        metrics = ["fired", "packets", "spikes"]
    missing = [m for m in metrics if m not in df.columns]
    if missing:
        raise ValueError(
            f"Columns not found in perf trace: {missing}. "
            f"Available: {list(df.columns)}"
        )

    if secondary_y is None:
        secondary_y = []
    secondary_y_set = set(secondary_y)

    if time_range is not None:
        start, end = time_range
        if "timestep" in df.columns:
            df = df[(df["timestep"] >= start) & (df["timestep"] < end)]
        else:
            df = df.iloc[start:end]

    if "timestep" in df.columns:
        timesteps = df["timestep"].values
    else:
        timesteps = np.arange(len(df))

    if colors is None:
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
                  for i in range(len(metrics))]
    if labels is None:
        labels = [m.replace("_", " ").title() for m in metrics]

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    ax_right = None
    if secondary_y_set:
        ax_right = ax.twinx()
        ax.right_ax = ax_right  # stash for user access

    plot_defaults = {
        "linewidth": style.perf_line_width,
    }
    if style.perf_marker:
        plot_defaults["marker"] = style.perf_marker
        plot_defaults["markersize"] = style.perf_marker_size
    plot_defaults.update(plot_kwargs)

    right_lines = []
    left_lines = []

    for metric, color, label in zip(metrics, colors, labels):
        values = df[metric].values.astype(float)
        target_ax = ax_right if metric in secondary_y_set else ax

        # Auto-scale energy columns
        scaled_label = label
        if "energy" in metric:
            values, unit = _auto_energy_scale(values)
            scaled_label = f"{label} ({unit})"
        elif "sim_time" in metric:
            values, unit = _auto_time_scale(values)
            scaled_label = f"{label} ({unit})"

        line, = target_ax.plot(
            timesteps, values,
            color=color, label=scaled_label,
            **plot_defaults,
        )
        if metric in secondary_y_set:
            right_lines.append(line)
        else:
            left_lines.append(line)

    if ylabel is None:
        ylabel = "Count"
    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)

    if ax_right is not None:
        ax_right.tick_params(axis="y", labelsize=style.tick_size)
        # Build unified legend from both axes
        if show_legend:
            all_lines = left_lines + right_lines
            all_labels = [l.get_label() for l in all_lines]
            ax.legend(all_lines, all_labels, loc="upper right", framealpha=0.9)
    elif show_legend:
        ax.legend(loc="upper right", framealpha=0.9)

    if style.tight_layout:
        fig.tight_layout()

    return fig, ax



def latency_histogram(
    data: Union[TraceData, Dict[str, Any]],
    metric: str = "processing_delay",
    filter_placeholder: bool = True,
    time_range: Optional[Tuple[int, int]] = None,
    bins: Optional[int] = None,
    log_scale: bool = False,
    color: Optional[str] = None,
    show_stats: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: str = "Count",
    **hist_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot a histogram of message-level latency or delay values.

    Draws a distribution of a chosen delay metric from the message trace.
    Useful for understanding timing characteristics and identifying
    congestion-related outliers.

    Args:
        data: TraceData or raw results dict from chip.sim().
        metric: Which delay column to histogram. One of:
            - "generation_delay" : time from soma fire to message send
            - "processing_delay" : time from receive to processed
            - "network_delay"    : time in the NoC (send to receive)
            - "blocking_delay"   : extra delay due to congestion
            - "hops"             : number of NoC hops (integer)
            Or any other numeric column in the message trace.
        filter_placeholder: If True (default), exclude placeholder messages
            (mid == -1) that represent timesteps where no real packet was sent.
        time_range: (start, end) timestep range to include.
        bins: Number of histogram bins. If None, uses style.hist_bins.
        log_scale: If True, use a log scale on the y-axis.
        color: Bar color. If None, uses the first default color.
        show_stats: If True, annotate the plot with mean, median, and std.
        ax: Existing Axes to plot on.
        style: Style configuration.
        figsize: Figure size (width, height).
        title: Plot title. If None, auto-generated from metric name.
        xlabel: X-axis label. If None, auto-generated from metric name.
        ylabel: Y-axis label.
        **hist_kwargs: Extra arguments forwarded to ax.hist().

    Returns: Tuple of (Figure, Axes).

    """
    data = _coerce_to_trace_data(data)
    df = _get_message_df(data)

    if style is None:
        style = get_default_style()

    if metric not in df.columns:
        raise ValueError(
            f"Column '{metric}' not found in message trace. "
            f"Available: {list(df.columns)}"
        )

    if filter_placeholder and "mid" in df.columns:
        df = df[df["mid"] >= 0]

    if time_range is not None:
        start, end = time_range
        df = df[(df["timestep"] >= start) & (df["timestep"] < end)]

    values = df[metric].values.astype(float)

    mask = np.isfinite(values)
    values = values[mask]

    if len(values) == 0:
        # Nothing to plot — still create the figure and return it
        if ax is None:
            fig, ax = create_figure(figsize=figsize, style=style)
        else:
            fig = ax.get_figure()
        ax.text(
            0.5, 0.5, "No data",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=style.label_size, color="#999999",
        )
        style_axis(ax, style, xlabel=xlabel or metric, ylabel=ylabel,
                   title=title)
        return fig, ax

    # Auto-scale time-based metrics
    is_time_metric = "delay" in metric or "latency" in metric or "timestamp" in metric
    unit = ""
    if is_time_metric:
        values, unit = _auto_time_scale(values)

    if bins is None:
        bins = style.hist_bins
    if color is None:
        color = DEFAULT_COLORS[0]

    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

    hist_defaults = {
        "bins": bins,
        "color": color,
        "alpha": style.hist_alpha,
        "edgecolor": style.hist_edgecolor,
        "linewidth": style.hist_edgewidth,
    }
    hist_defaults.update(hist_kwargs)
    ax.hist(values, **hist_defaults)

    if log_scale:
        ax.set_yscale("log")

    if show_stats and len(values) > 0:
        mean_val = np.mean(values)
        median_val = np.median(values)
        std_val = np.std(values)
        stat_text = (
            f"mean = {mean_val:.3g}\n"
            f"median = {median_val:.3g}\n"
            f"std = {std_val:.3g}\n"
            f"n = {len(values)}"
        )
        ax.text(
            0.97, 0.95, stat_text,
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=style.tick_size,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#cccccc", alpha=0.9),
        )

    if xlabel is None:
        pretty = metric.replace("_", " ").title()
        xlabel = f"{pretty} ({unit})" if unit else pretty
    if title is None:
        title = f"Distribution of {metric.replace('_', ' ').title()}"

    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)

    if style.tight_layout:
        fig.tight_layout()

    return fig, ax


def latency_comparison(
    data: Union[TraceData, Dict[str, Any]],
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
    title: str = "Delay Component Comparison",
    xlabel: Optional[str] = None,
    ylabel: str = "Count",
    **hist_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Overlay multiple delay distributions on the same axes.

    Plots overlapping histograms for several delay components so their
    distributions can be visually compared.

    Args:
        data: TraceData or raw results dict from chip.sim().
        metrics: Delay columns to compare. Defaults to
            ["generation_delay", "processing_delay", "network_delay",
             "blocking_delay"].
        filter_placeholder: Exclude placeholder messages (mid == -1).
        time_range: (start, end) timestep range.
        bins: Number of bins for all histograms.
        colors: Color per metric.
        labels: Display names per metric.
        log_scale: Log y-axis.
        show_legend: Show legend.
        ax: Existing Axes.
        style: Style configuration.
        figsize: Figure size.
        title: Plot title.
        xlabel: X-axis label. Auto-scaled with time unit if None.
        ylabel: Y-axis label.
        **hist_kwargs: Extra arguments forwarded to ax.hist().

    Returns: Tuple of (Figure, Axes).

    """
    data = _coerce_to_trace_data(data)
    df = _get_message_df(data)

    if style is None:
        style = get_default_style()

    if metrics is None:
        metrics = list(_LATENCY_COLUMNS)
    missing = [m for m in metrics if m not in df.columns]
    if missing:
        raise ValueError(
            f"Columns not found in message trace: {missing}. "
            f"Available: {list(df.columns)}"
        )

    if filter_placeholder and "mid" in df.columns:
        df = df[df["mid"] >= 0]

    if time_range is not None:
        start, end = time_range
        df = df[(df["timestep"] >= start) & (df["timestep"] < end)]

    # Collect all finite values across metrics to find a common time unit
    all_values = np.concatenate([
        df[m].values.astype(float) for m in metrics
    ])
    all_values = all_values[np.isfinite(all_values)]
    if len(all_values) == 0:
        if ax is None:
            fig, ax = create_figure(figsize=figsize, style=style)
        else:
            fig = ax.get_figure()
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=style.label_size, color="#999")
        return fig, ax

    _, unit = _auto_time_scale(all_values)
    factor = {"ps": 1e-12, "ns": 1e-9, "µs": 1e-6,
              "ms": 1e-3, "s": 1.0}.get(unit, 1.0)

    # Resolve colors and labels
    if colors is None:
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
                  for i in range(len(metrics))]
    if labels is None:
        labels = [m.replace("_", " ").title() for m in metrics]

    if bins is None:
        bins = style.hist_bins

    # Create figure
    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()

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

    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)

    if show_legend:
        ax.legend(loc="upper right", framealpha=0.9)

    if style.tight_layout:
        fig.tight_layout()

    return fig, ax
