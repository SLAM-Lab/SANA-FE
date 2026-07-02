"""
This module provides consistent styling across all SANA-FE plots, including
color palettes, figure sizes, and matplotlib configuration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

DEFAULT_COLORS = [
    "#1f77b4",  # Blue
    "#ff7f0e",  # Orange
    "#2ca02c",  # Green
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#8c564b",  # Brown
    "#e377c2",  # Pink
    "#7f7f7f",  # Gray
    "#bcbd22",  # Olive
    "#17becf",  # Cyan
]

# Extended palette
EXTENDED_COLORS = DEFAULT_COLORS + [
    "#aec7e8",  # Light blue
    "#ffbb78",  # Light orange
    "#98df8a",  # Light green
    "#ff9896",  # Light red
    "#c5b0d5",  # Light purple
    "#c49c94",  # Light brown
    "#f7b6d2",  # Light pink
    "#c7c7c7",  # Light gray
    "#dbdb8d",  # Light olive
    "#9edae5",  # Light cyan
]

NEUROMORPHIC_CMAP_COLORS = [
    "#0d0887",  # Deep purple
    "#46039f",
    "#7201a8",
    "#9c179e",
    "#bd3786",
    "#d8576b",
    "#ed7953",
    "#fb9f3a",
    "#fdca26",
    "#f0f921",  # Bright yellow
]


@dataclass
class SANAFEStyle:
    """
    Style configuration for plots

        colors: List of colors for different groups/series
        figure_size: Default (width, height) in inches
        dpi: Resolution for saved figures
        font_family: Font family for text
        font_size: Base font size in points
        title_size: Title font size in points
        label_size: Axis label font size in points
        tick_size: Tick label font size in points
        line_width: Default line width
        marker_size: Default marker size for scatter plots
        spine_width: Width of axis spines
        grid: Whether to show grid by default
        grid_alpha: Transparency of grid lines
        tight_layout: Whether to use tight_layout by default
    """
    colors: List[str] = field(default_factory=DEFAULT_COLORS.copy)
    figure_size: Tuple[float, float] = (8.0, 5.0)
    dpi: int = 100
    font_family: str = "sans-serif"
    font_size: float = 11.0
    title_size: float = 13.0
    label_size: float = 11.0
    tick_size: float = 10.0
    line_width: float = 1.5
    marker_size: float = 30.0
    spine_width: float = 1.0
    grid: bool = False
    grid_alpha: float = 0.3
    tight_layout: bool = True

    # Raster plot specific
    raster_marker: str = "|"
    raster_marker_size: float = 100.0
    raster_line_width: float = 1.5

    # Potential plot specific
    potential_line_width: float = 1.5
    potential_marker: Optional[str] = None
    potential_marker_size: float = 4.0

    # Performance plot specific
    perf_line_width: float = 1.5
    perf_marker: Optional[str] = "o"
    perf_marker_size: float = 3.0
    perf_fill_alpha: float = 0.3

    # Energy breakdown specific
    energy_colors: List[str] = field(default_factory=lambda: [
        "#1f77b4",  # Synapse - Blue
        "#2ca02c",  # Dendrite - Green
        "#ff7f0e",  # Soma - Orange
        "#d62728",  # Network - Red
    ])
    energy_component_names: List[str] = field(default_factory=lambda: [
        "Synapse", "Dendrite", "Soma", "Network",
    ])

    # Histogram specific
    hist_bins: int = 30
    hist_alpha: float = 0.7
    hist_edgecolor: str = "white"
    hist_edgewidth: float = 0.5

    def to_rc_params(self) -> Dict[str, Any]:
        """Convert style to matplotlib rcParams dictionary."""
        return {
            "figure.figsize": self.figure_size,
            "figure.dpi": self.dpi,
            "font.family": self.font_family,
            "font.size": self.font_size,
            "axes.titlesize": self.title_size,
            "axes.labelsize": self.label_size,
            "xtick.labelsize": self.tick_size,
            "ytick.labelsize": self.tick_size,
            "lines.linewidth": self.line_width,
            "lines.markersize": self.marker_size,
            "axes.linewidth": self.spine_width,
            "axes.grid": self.grid,
            "grid.alpha": self.grid_alpha,
            "figure.autolayout": self.tight_layout,
        }


PUBLICATION_STYLE = SANAFEStyle(
    figure_size=(6.0, 4.0),
    dpi=300,
    font_family="serif",
    font_size=10.0,
    title_size=11.0,
    label_size=10.0,
    tick_size=9.0,
    line_width=1.0,
    spine_width=0.8,
    raster_marker_size=80.0,
    raster_line_width=1.0,
)

PRESENTATION_STYLE = SANAFEStyle(
    figure_size=(12.0, 7.0),
    dpi=100,
    font_size=14.0,
    title_size=18.0,
    label_size=14.0,
    tick_size=12.0,
    line_width=2.5,
    spine_width=1.5,
    raster_marker_size=150.0,
    raster_line_width=2.5,
)

NOTEBOOK_STYLE = SANAFEStyle(
    figure_size=(10.0, 6.0),
    dpi=100,
    font_size=12.0,
    title_size=14.0,
    label_size=12.0,
    tick_size=11.0,
    grid=True,
    grid_alpha=0.2,
)

# Global default style
_default_style: SANAFEStyle = SANAFEStyle()


def get_default_style() -> SANAFEStyle:
    """Get default stylesheet"""
    return _default_style


def set_default_style(style: Optional[SANAFEStyle] = None) -> None:
    """
    Set the default style for all SANA-FE plots.
    """
    global _default_style  # pylint: disable=global-statement
    if style is None:
        _default_style = SANAFEStyle()
    else:
        _default_style = style


def apply_style(style: Optional[SANAFEStyle] = None) -> None:
    """
    Apply style settings to matplotlib's rcParams.
    """
    if style is None:
        style = _default_style

    for key, value in style.to_rc_params().items():
        mpl.rcParams[key] = value


def get_group_colors(
    groups: Sequence[str],
    style: Optional[SANAFEStyle] = None,
) -> Dict[str, str]:
    """
    Assigns colors from the style's color palette to each group name.
    Colors cycle if there are more groups than colors.

    Returns: Dictionary mapping group names to color strings.

    """
    if style is None:
        style = _default_style

    colors = style.colors if len(style.colors) >= len(groups) else EXTENDED_COLORS

    return {
        group: colors[i % len(colors)]
        for i, group in enumerate(groups)
    }


def get_colormap(
    name: str = "neuromorphic",
    n_colors: int = 256,
) -> LinearSegmentedColormap:
    """
    Get a colormap for SANA-FE visualizations.

    Args:
        name: Colormap name. Options:
            - "neuromorphic": Purple-to-yellow gradient (default)
            - "activity": Blue-to-red for activity levels
            - "energy": Green-to-red for energy consumption
            - Any matplotlib colormap name
        n_colors: Number of discrete colors in the colormap

    Returns: Matplotlib colormap object.
    """
    if name == "neuromorphic":
        return LinearSegmentedColormap.from_list(
            "neuromorphic",
            NEUROMORPHIC_CMAP_COLORS,
            N=n_colors,
        )
    if name == "activity":
        return LinearSegmentedColormap.from_list(
            "activity",
            ["#2166ac", "#f7f7f7", "#b2182b"],
            N=n_colors,
        )
    if name == "energy":
        return LinearSegmentedColormap.from_list(
            "energy",
            ["#1a9850", "#ffffbf", "#d73027"],
            N=n_colors,
        )

    return plt.get_cmap(name)


def create_figure(
    figsize: Optional[Tuple[float, float]] = None,
    style: Optional[SANAFEStyle] = None,
    **kwargs,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a figure with styling applied.

    Returns: Tuple of (Figure, Axes) objects.
    """
    if style is None:
        style = _default_style

    if figsize is None:
        figsize = style.figure_size

    fig, ax = plt.subplots(figsize=figsize, **kwargs)

    # Apply spine styling
    for spine in ax.spines.values():
        spine.set_linewidth(style.spine_width)

    # Apply grid if enabled
    if style.grid:
        ax.grid(True, alpha=style.grid_alpha)

    return fig, ax


def style_axis(
    ax: plt.Axes,
    style: Optional[SANAFEStyle] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
) -> plt.Axes:
    """Apply SANA-FE style to given axis."""
    if style is None:
        style = _default_style

    # Apply spine styling
    for spine in ax.spines.values():
        spine.set_linewidth(style.spine_width)

    # Apply grid if enabled
    if style.grid:
        ax.grid(True, alpha=style.grid_alpha)

    # Set labels and title
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=style.label_size)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=style.label_size)
    if title is not None:
        ax.set_title(title, fontsize=style.title_size)

    # Set limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # Set tick label sizes
    ax.tick_params(axis="both", labelsize=style.tick_size)

    return ax
