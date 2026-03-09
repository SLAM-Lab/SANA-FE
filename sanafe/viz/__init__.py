"""
Provides plotting utilities for visualizing SANA-FE simulation outputs,
including spike raster plots, potential timeseries, and performance metrics.
"""

from sanafe.viz.raster import raster_plot
from sanafe.viz.potential import potential_plot
from sanafe.viz.performance import (
    energy_breakdown_plot,
    throughput_plot,
    latency_histogram,
    latency_comparison,
)
from sanafe.viz.styles import (
    SANAFEStyle,
    get_group_colors,
    apply_style,
    set_default_style,
    PUBLICATION_STYLE,
    PRESENTATION_STYLE,
    NOTEBOOK_STYLE,
)

__all__ = [
    # SNN visualization
    "raster_plot",
    "potential_plot",
    # Performance visualization
    "energy_breakdown_plot",
    "throughput_plot",
    "latency_histogram",
    "latency_comparison",
    # Styling
    "SANAFEStyle",
    "get_group_colors",
    "apply_style",
    "set_default_style",
    "PUBLICATION_STYLE",
    "PRESENTATION_STYLE",
    "NOTEBOOK_STYLE",
]
