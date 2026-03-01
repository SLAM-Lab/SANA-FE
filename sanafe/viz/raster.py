"""
raster.py - Spike raster plot visualization for SANA-FE.

This module provides functions for creating spike raster plots, which display
neural spiking activity over time with neurons on the y-axis and timesteps
on the x-axis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from sanafe.data.traces import TraceData
from sanafe.viz.styles import (
    SANAFEStyle,
    create_figure,
    get_default_style,
    get_group_colors,
    style_axis,
)


def raster_plot(
    data: Union[TraceData, Dict[str, Any]],
    groups: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    colors: Optional[Dict[str, str]] = None,
    group_spacing: float = 0.5,
    show_group_labels: bool = True,
    show_legend: bool = True,
    ax: Optional[plt.Axes] = None,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Time-step",
    ylabel: str = "Neuron",
    **scatter_kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a spike raster plot from SANA-FE trace data.
    
    Displays neural spiking activity with neurons on the y-axis and timesteps
    on the x-axis. Each spike is shown as a vertical tick mark. Neurons can
    be colored by group.
    
    Args:
        data: Either a TraceData object or raw results dict from chip.sim()
        groups: List of group names to include. If None, includes all groups.
        time_range: Tuple of (start, end) timesteps to display. If None, shows all.
        colors: Dict mapping group names to color strings. If None, auto-assigned.
        group_spacing: Vertical spacing between neuron groups (in neuron units)
        show_group_labels: Whether to show group names on y-axis
        show_legend: Whether to show a legend for group colors
        ax: Existing Axes to plot on. If None, creates new figure.
        style: Style configuration. If None, uses default style.
        figsize: Figure size (width, height). If None, uses style default.
        title: Plot title. If None, no title is shown.
        xlabel: X-axis label
        ylabel: Y-axis label
        **scatter_kwargs: Additional arguments passed to plt.scatter()
        
    Returns: Tuple of (Figure, Axes) objects.
        
    """
    if isinstance(data, dict):
        data = TraceData.from_sim_results(data)
    
    if not data.has_spikes():
        raise ValueError("No spike trace data available in the provided data")
    
    if style is None:
        style = get_default_style()
    
    all_groups = data.get_neuron_groups()
    if groups is not None:
        # Validate groups exist
        invalid = set(groups) - set(all_groups)
        if invalid:
            raise ValueError(f"Unknown groups: {invalid}. Available: {all_groups}")
        groups = list(groups)
    else:
        groups = all_groups
    
    if time_range is not None:
        data = data.filter_by_time(start=time_range[0], end=time_range[1])
    
    if colors is None:
        colors = get_group_colors(groups, style)
    
    # Build neuron ordering: map (group, offset) -> y-position
    neurons_by_group: Dict[str, set] = {g: set() for g in groups}
    for spikes_at_t in data.spike_trace:
        for spike in spikes_at_t:
            if spike.group_name in neurons_by_group:
                neurons_by_group[spike.group_name].add(spike.neuron_offset)
    
    # Create y-position mapping with group spacing
    neuron_to_y: Dict[Tuple[str, int], float] = {}
    y_ticks: List[float] = []
    y_tick_labels: List[str] = []
    current_y = 0.0
    group_y_ranges: Dict[str, Tuple[float, float]] = {}
    
    for group in groups:
        offsets = sorted(neurons_by_group[group])
        if not offsets:
            continue
            
        group_start_y = current_y
        for offset in offsets:
            neuron_to_y[(group, offset)] = current_y
            current_y += 1.0
        group_end_y = current_y - 1.0
        
        group_y_ranges[group] = (group_start_y, group_end_y)
        
        group_center = (group_start_y + group_end_y) / 2
        y_ticks.append(group_center)
        y_tick_labels.append(group)
        
        current_y += group_spacing
    
    # Prepare scatter data
    timesteps = []
    y_positions = []
    spike_colors = []
    
    for timestep, spikes_at_t in enumerate(data.spike_trace):
        for spike in spikes_at_t:
            key = (spike.group_name, spike.neuron_offset)
            if key in neuron_to_y:
                timesteps.append(timestep)
                y_positions.append(neuron_to_y[key])
                spike_colors.append(colors.get(spike.group_name, "#333333"))
    
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
    
    # Plot spikes
    if timesteps:
        ax.scatter(timesteps, y_positions, c=spike_colors, **scatter_defaults)
    
    if show_group_labels and y_ticks:
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tick_labels)
    else:
        ax.set_ylabel(ylabel)
    
    time_max = data.timesteps
    y_max = current_y - group_spacing  # Remove last spacing
    
    ax.set_xlim(-0.5, time_max + 0.5)
    ax.set_ylim(-0.5, y_max + 0.5)
    
    style_axis(ax, style, xlabel=xlabel, title=title)
    
    if show_legend and len(groups) > 1:
        from matplotlib.lines import Line2D
        handles = [
            Line2D([0], [0], marker=style.raster_marker, color="w",
                   markerfacecolor=colors.get(g, "#333333"),
                   markeredgecolor=colors.get(g, "#333333"),
                   markersize=10, label=g)
            for g in groups if g in group_y_ranges
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
    """
    Create a spike raster plot from a binary spike matrix.
    
    Alternative to raster_plot() when you have spike data as a numpy array
    rather than TraceData.
    
    Args:
        spike_matrix: 2D binary array of shape (timesteps, neurons) where
            1 indicates a spike occurred
        neuron_ids: Optional list of neuron ID strings for y-axis labels
        time_range: Tuple of (start, end) timesteps to display
        ax: Existing Axes to plot on
        style: Style configuration
        figsize: Figure size (width, height)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        color: Color for all spikes
        **scatter_kwargs: Additional arguments passed to plt.scatter()
        
    Returns: Tuple of (Figure, Axes) objects.
    """
    if style is None:
        style = get_default_style()
    
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
        # Show subset of labels if too many neurons
        if n_neurons > 20:
            step = n_neurons // 10
            tick_positions = list(range(0, n_neurons, step))
            tick_labels = [neuron_ids[i] for i in tick_positions]
            ax.set_yticks(tick_positions)
            ax.set_yticklabels(tick_labels)
        else:
            ax.set_yticks(range(n_neurons))
            ax.set_yticklabels(neuron_ids)
    
    ax.set_xlim(-0.5, n_timesteps + 0.5)
    ax.set_ylim(-0.5, n_neurons - 0.5)
    
    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)
    
    if style.tight_layout:
        fig.tight_layout()
    
    return fig, ax
