"""
This module provides functions for creating potential timeseries plots,
which display membrane voltage traces for individual neurons over time.
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
    style_axis,
    DEFAULT_COLORS,
)


def potential_plot(
    data: Union[TraceData, Dict[str, Any]],
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
    """
    Create a membrane potential timeseries plot from SANA-FE trace data.
    Displays membrane voltage traces for neurons over time. Multiple neurons
    can be shown on the same plot with different colors.
    
    Args:
        data: Either a TraceData object or raw results dict from chip.sim()
        neuron_ids: List of neuron ID strings for legend labels. If None, neurons are labeled numerically (Neuron 0, Neuron 1, ...).
        time_range: Tuple of (start, end) timesteps to display. If None, shows all.
        colors: List of colors for each neuron trace. If None, auto-assigned.
        show_threshold: If set, draws a horizontal threshold line at this value.
        threshold_color: Color for the threshold line.
        threshold_linestyle: Line style for the threshold line.
        show_legend: Whether to show a legend.
        ax: Existing Axes to plot on. If None, creates new figure.
        style: Style configuration. If None, uses default style.
        figsize: Figure size (width, height). If None, uses style default.
        title: Plot title. If None, no title is shown.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        **plot_kwargs: Additional arguments passed to plt.plot().
        
    Returns: Tuple of (Figure, Axes) objects.
    """
    # Handle raw results dict
    if isinstance(data, dict):
        data = TraceData.from_sim_results(data)
    
    if not data.has_potentials():
        raise ValueError("No potential trace data available in the provided data")
    
    if style is None:
        style = get_default_style()
    
    # Get potential data as array
    potentials = data.potentials_to_array()
    n_timesteps, n_neurons = potentials.shape
    
    # Apply time range filter
    time_offset = 0
    if time_range is not None:
        start, end = time_range
        potentials = potentials[start:end, :]
        n_timesteps = potentials.shape[0]
        time_offset = start
    
    # Generate default neuron labels
    if neuron_ids is None:
        neuron_ids = [f"Neuron {i}" for i in range(n_neurons)]
    elif len(neuron_ids) != n_neurons:
        raise ValueError(
            f"Number of neuron_ids ({len(neuron_ids)}) doesn't match "
            f"number of traced neurons ({n_neurons})"
        )
    
    if colors is None:
        colors = DEFAULT_COLORS[:n_neurons]
        if len(colors) < n_neurons:
            # Extend with cycling if not enough colors
            colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(n_neurons)]
    
    # Create figure if needed
    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()
    
    # Time axis
    timesteps = np.arange(time_offset, time_offset + n_timesteps)
    
    plot_defaults = {
        "linewidth": style.potential_line_width,
    }
    if style.potential_marker:
        plot_defaults["marker"] = style.potential_marker
        plot_defaults["markersize"] = style.potential_marker_size
    plot_defaults.update(plot_kwargs)
    
    # Plot each neuron's potential trace
    lines = []
    for i in range(n_neurons):
        line, = ax.plot(
            timesteps,
            potentials[:, i],
            color=colors[i],
            label=neuron_ids[i],
            **plot_defaults,
        )
        lines.append(line)
    
    # Add threshold line if specified
    if show_threshold is not None:
        ax.axhline(
            y=show_threshold,
            color=threshold_color,
            linestyle=threshold_linestyle,
            linewidth=style.line_width * 0.8,
            label="Threshold",
            zorder=0,  # Behind the traces
        )
    
    ax.set_xlim(timesteps[0] - 0.5, timesteps[-1] + 0.5)
    
    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)
    
    # Legend
    if show_legend:
        ax.legend(loc="upper right", framealpha=0.9)
    
    if style.tight_layout:
        fig.tight_layout()
    
    return fig, ax


def potential_heatmap(
    data: Union[TraceData, Dict[str, Any], np.ndarray],
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
    """
    Heatmap visualization of membrane potentials over time.
    
    Args:
        data: TraceData, results dict, or numpy array of potentials
        neuron_ids: List of neuron ID strings for y-axis labels
        time_range: Tuple of (start, end) timesteps to display
        cmap: Colormap name for the heatmap
        vmin: Minimum value for color scaling
        vmax: Maximum value for color scaling
        show_colorbar: Whether to show a colorbar
        ax: Existing Axes to plot on
        style: Style configuration
        figsize: Figure size (width, height)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        **imshow_kwargs: Additional arguments passed to plt.imshow()
        
    Returns: Tuple of (Figure, Axes) objects.
    """
    # Handle different input types
    if isinstance(data, dict):
        data = TraceData.from_sim_results(data)
    
    if isinstance(data, TraceData):
        if not data.has_potentials():
            raise ValueError("No potential trace data available")
        potentials = data.potentials_to_array()
    else:
        potentials = np.asarray(data)
    
    if style is None:
        style = get_default_style()
    
    n_timesteps, n_neurons = potentials.shape
    
    # Apply time range filter
    time_offset = 0
    if time_range is not None:
        start, end = time_range
        potentials = potentials[start:end, :]
        n_timesteps = potentials.shape[0]
        time_offset = start
    
    if ax is None:
        fig, ax = create_figure(figsize=figsize, style=style)
    else:
        fig = ax.get_figure()
    
    # Create heatmap (transpose so neurons are on y-axis)
    imshow_defaults = {
        "aspect": "auto",
        "origin": "lower",
        "cmap": cmap,
        "extent": [time_offset - 0.5, time_offset + n_timesteps - 0.5,
                   -0.5, n_neurons - 0.5],
    }
    if vmin is not None:
        imshow_defaults["vmin"] = vmin
    if vmax is not None:
        imshow_defaults["vmax"] = vmax
    imshow_defaults.update(imshow_kwargs)
    
    im = ax.imshow(potentials.T, **imshow_defaults)
    
    if show_colorbar:
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Membrane Potential")
    
    # Configure y-axis labels
    if neuron_ids is not None:
        if n_neurons <= 20:
            ax.set_yticks(range(n_neurons))
            ax.set_yticklabels(neuron_ids)
        else:
            # Show subset of labels
            step = n_neurons // 10
            tick_positions = list(range(0, n_neurons, step))
            tick_labels = [neuron_ids[i] for i in tick_positions]
            ax.set_yticks(tick_positions)
            ax.set_yticklabels(tick_labels)
    
    style_axis(ax, style, xlabel=xlabel, ylabel=ylabel, title=title)
    
    if style.tight_layout:
        fig.tight_layout()
    
    return fig, ax


def potential_subplots(
    data: Union[TraceData, Dict[str, Any]],
    neuron_ids: Optional[Sequence[str]] = None,
    time_range: Optional[Tuple[int, int]] = None,
    ncols: int = 1,
    colors: Optional[Sequence[str]] = None,
    show_threshold: Optional[float] = None,
    sharex: bool = True,
    sharey: bool = True,
    style: Optional[SANAFEStyle] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **plot_kwargs,
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Grid of subplots, one per neuron potential trace.
    
    Args:
        data: TraceData or results dict with potential traces
        neuron_ids: List of neuron ID strings for subplot titles
        time_range: Tuple of (start, end) timesteps to display
        ncols: Number of columns in the subplot grid
        colors: List of colors for each subplot
        show_threshold: If set, draws threshold line on each subplot
        sharex: Whether subplots share x-axis
        sharey: Whether subplots share y-axis
        style: Style configuration
        figsize: Figure size (width, height)
        **plot_kwargs: Additional arguments passed to plt.plot()
        
    Returns:
        Tuple of (Figure, array of Axes) objects.
    """
    # Handle raw results dict
    if isinstance(data, dict):
        data = TraceData.from_sim_results(data)
    
    if not data.has_potentials():
        raise ValueError("No potential trace data available")
    
    if style is None:
        style = get_default_style()
    
    # Get potential data
    potentials = data.potentials_to_array()
    n_timesteps, n_neurons = potentials.shape
    
    time_offset = 0
    if time_range is not None:
        start, end = time_range
        potentials = potentials[start:end, :]
        n_timesteps = potentials.shape[0]
        time_offset = start
    
    # Generate default neuron labels
    if neuron_ids is None:
        neuron_ids = [f"Neuron {i}" for i in range(n_neurons)]
    
    if colors is None:
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i in range(n_neurons)]
    
    # Calculate grid dimensions
    nrows = (n_neurons + ncols - 1) // ncols
    
    # Determine figure size
    if figsize is None:
        width = style.figure_size[0] * min(ncols, 2)
        height = style.figure_size[1] * 0.6 * nrows
        figsize = (width, height)
    
    # Create subplots
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=figsize,
        sharex=sharex,
        sharey=sharey,
        squeeze=False,
    )
    axes = axes.flatten()
    
    timesteps = np.arange(time_offset, time_offset + n_timesteps)
    
    # Plot each neuron
    plot_defaults = {
        "linewidth": style.potential_line_width,
    }
    plot_defaults.update(plot_kwargs)
    
    for i in range(n_neurons):
        ax = axes[i]
        ax.plot(timesteps, potentials[:, i], color=colors[i], **plot_defaults)
        ax.set_title(neuron_ids[i], fontsize=style.label_size)
        
        if show_threshold is not None:
            ax.axhline(y=show_threshold, color="#d62728", linestyle="--",
                      linewidth=style.line_width * 0.6, alpha=0.7)
        
        # Add labels to edge subplots only
        if i >= n_neurons - ncols:
            ax.set_xlabel("Time-step", fontsize=style.tick_size)
        if i % ncols == 0:
            ax.set_ylabel("Potential", fontsize=style.tick_size)
    
    # Hide unused subplots
    for i in range(n_neurons, len(axes)):
        axes[i].set_visible(False)
    
    fig.tight_layout()
    
    return fig, axes[:n_neurons]
