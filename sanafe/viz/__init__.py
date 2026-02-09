"""
Provides plotting utilities for visualizing SANA-FE simulation outputs,
including spike raster plots, potential timeseries, and performance metrics.
"""

# from sanafe.viz.raster import raster_plot
from sanafe.viz.potential import potential_plot
from sanafe.viz.styles import (
    SANAFEStyle,
    get_group_colors,
    apply_style,
    set_default_style,
)

__all__ = [
    # "raster_plot",
    "potential_plot",
    "SANAFEStyle",
    "get_group_colors",
    "apply_style",
    "set_default_style",
]
