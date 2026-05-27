"""
Convert SANA-FE trace outputs (in-memory or CSV) to pandas DataFrames.
"""

from sanafe.data.traces import (
    spikes_to_dataframe,
    potentials_to_dataframe,
    performance_to_dataframe,
    messages_to_dataframe,
)

__all__ = [
    "spikes_to_dataframe",
    "potentials_to_dataframe",
    "performance_to_dataframe",
    "messages_to_dataframe",
]
