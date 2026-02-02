"""
Provides the TraceData class for unified access to all trace types
produced by SANA-FE simulations, with conversion methods for pandas DataFrames
and numpy arrays.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd


@dataclass
class SpikeEvent:
    """
    Represents a single spike event.

        timestep: The simulation timestep when the spike occurred
        group_name: Name of the neuron group containing the spiking neuron
        neuron_offset: Index of the neuron within its group
    """
    timestep: int
    group_name: str
    neuron_offset: int
    
    @property
    def neuron_id(self) -> str:
        """Full neuron identifier in 'group.offset' format."""
        return f"{self.group_name}.{self.neuron_offset}"


@dataclass
class TraceData:
    """    
    Provides unified access to spike traces, potential traces, performance
    metrics, and message traces from SANA-FE simulations. Supports loading
    from both in-memory simulation results and CSV files.
    
        spike_trace: Raw spike trace data (list of lists of NeuronAddress)
        potential_trace: Raw potential trace data (list of lists of floats)
        perf_trace: Raw performance trace data (dict of metric lists)
        message_trace: Raw message trace data (list of lists of message dicts)
        neuron_labels: Optional mapping of (group, offset) to custom labels
        
        From simulation results
        From CSV files
        Convert to pandas
    """
    
    spike_trace: Optional[List[List[Any]]] = None
    potential_trace: Optional[List[List[float]]] = None
    perf_trace: Optional[Dict[str, List[Any]]] = None
    message_trace: Optional[List[List[Dict[str, Any]]]] = None
    neuron_labels: Dict[Tuple[str, int], str] = field(default_factory=dict)
    
    # Cache for converted data
    _spike_df_cache: Optional[pd.DataFrame] = field(default=None, repr=False)
    _potential_df_cache: Optional[pd.DataFrame] = field(default=None, repr=False)
    _perf_df_cache: Optional[pd.DataFrame] = field(default=None, repr=False)
    _message_df_cache: Optional[pd.DataFrame] = field(default=None, repr=False)
    
    @classmethod
    def from_sim_results(cls, results: Dict[str, Any]) -> "TraceData":
        """
        Create TraceData from chip.sim() return dictionary.
        
            Args: Dictionary returned by SpikingChip.sim(), containing keys like 'spike_trace', 'potential_trace', 'perf_trace', 'message_trace'.
            
            Returns: TraceData instance with loaded traces.

            results = chip.sim(100, spike_trace=True, perf_trace=True)
            traces = TraceData.from_sim_results(results)
        """
        return cls(
            spike_trace=results.get("spike_trace"),
            potential_trace=results.get("potential_trace"),
            perf_trace=results.get("perf_trace"),
            message_trace=results.get("message_trace"),
        )
    
    @classmethod
    def from_files(
        cls,
        spike_csv: Optional[Union[str, Path]] = None,
        potential_csv: Optional[Union[str, Path]] = None,
        perf_csv: Optional[Union[str, Path]] = None,
        message_csv: Optional[Union[str, Path]] = None,
    ) -> "TraceData":
        """
        Create TraceData from CSV trace files.

        spike_csv: Path to spike trace CSV file
        potential_csv: Path to potential trace CSV file
        perf_csv: Path to performance trace CSV file
        message_csv: Path to message trace CSV file
            
        Returns: TraceData instance with loaded traces.
        """
        instance = cls()
        
        if spike_csv is not None:
            instance._load_spike_csv(Path(spike_csv))
        if potential_csv is not None:
            instance._load_potential_csv(Path(potential_csv))
        if perf_csv is not None:
            instance._load_perf_csv(Path(perf_csv))
        if message_csv is not None:
            instance._load_message_csv(Path(message_csv))
            
        return instance
    
    def _load_spike_csv(self, path: Path) -> None:
        """Load spike trace from CSV file."""
        # SANA-FE spike CSV format: neuron,timestep
        self.spike_trace = []
        max_timestep = 0
        
        with open(path, "r") as f:
            reader = csv.DictReader(f)
            spikes_by_timestep: Dict[int, List[Any]] = {}
            
            for row in reader:
                timestep = int(row["timestep"])
                neuron_str = row["neuron"]
                group_name, offset_str = neuron_str.rsplit(".", 1)
                
                max_timestep = max(max_timestep, timestep)
                
                if timestep not in spikes_by_timestep:
                    spikes_by_timestep[timestep] = []
                    
                # Object mimicking NeuronAddress
                spike_info = type("NeuronAddress", (), {
                    "group_name": group_name,
                    "neuron_offset": int(offset_str)
                })()
                spikes_by_timestep[timestep].append(spike_info)
        
        # Convert to list format (indexed by timestep)
        self.spike_trace = [
            spikes_by_timestep.get(t, []) 
            for t in range(max_timestep + 1)
        ]
    
    def _load_potential_csv(self, path: Path) -> None:
        # Potential CSV format: columns are neuron IDs, rows are timesteps
        df = pd.read_csv(path)
        self.potential_trace = df.values.tolist()
    
    def _load_perf_csv(self, path: Path) -> None:
        df = pd.read_csv(path)
        self.perf_trace = df.to_dict(orient="list")
    
    def _load_message_csv(self, path: Path) -> None:
        self.message_trace = []
        
        with open(path, "r") as f:
            reader = csv.DictReader(f)
            messages_by_timestep: Dict[int, List[Dict]] = {}
            max_timestep = 0
            
            for row in reader:
                timestep = int(row["timestep"])
                max_timestep = max(max_timestep, timestep)
                
                if timestep not in messages_by_timestep:
                    messages_by_timestep[timestep] = []
                
                # Convert numeric fields
                msg = {}
                for key, value in row.items():
                    try:
                        if "." in value or "e" in value.lower():
                            msg[key] = float(value)
                        else:
                            msg[key] = int(value)
                    except (ValueError, AttributeError):
                        msg[key] = value
                        
                messages_by_timestep[timestep].append(msg)
        
        self.message_trace = [
            messages_by_timestep.get(t, [])
            for t in range(max_timestep + 1)
        ]
    
    def has_spikes(self) -> bool:
        return self.spike_trace is not None and len(self.spike_trace) > 0
    
    def has_potentials(self) -> bool:
        return self.potential_trace is not None and len(self.potential_trace) > 0
    
    def has_performance(self) -> bool:
        return self.perf_trace is not None and len(self.perf_trace) > 0
    
    def has_messages(self) -> bool:
        return self.message_trace is not None and len(self.message_trace) > 0
    
    @property
    def timesteps(self) -> int:
        """
        Get the number of timesteps in the trace data.
        
        Returns the length based on whichever trace is available.
        """
        if self.has_spikes():
            return len(self.spike_trace)
        if self.has_potentials():
            return len(self.potential_trace)
        if self.has_performance() and "timestep" in self.perf_trace:
            return len(self.perf_trace["timestep"])
        if self.has_messages():
            return len(self.message_trace)
        return 0
    
    def spikes_to_dataframe(self) -> pd.DataFrame:
        """
        Convert spike trace to pandas DataFrame.
        
        Returns:
            DataFrame with columns:
                - timestep: Simulation timestep
                - group: Neuron group name
                - neuron_offset: Index within group
                - neuron_id: Full identifier (group.offset)
        """
        if not self.has_spikes():
            raise ValueError("No spike trace data available")
        
        if self._spike_df_cache is not None:
            return self._spike_df_cache.copy()
        
        records = []
        for timestep, spikes_at_t in enumerate(self.spike_trace):
            for spike in spikes_at_t:
                records.append({
                    "timestep": timestep,
                    "group": spike.group_name,
                    "neuron_offset": spike.neuron_offset,
                    "neuron_id": f"{spike.group_name}.{spike.neuron_offset}",
                })
        
        self._spike_df_cache = pd.DataFrame(records)
        return self._spike_df_cache.copy()
    
    def spikes_to_events(self) -> List[SpikeEvent]:
        """
        Convert spike trace to list of SpikeEvent objects.
        
        Returns: List of SpikeEvent objects, sorted by timestep.
        """
        if not self.has_spikes():
            raise ValueError("No spike trace data available")
        
        events = []
        for timestep, spikes_at_t in enumerate(self.spike_trace):
            for spike in spikes_at_t:
                events.append(SpikeEvent(
                    timestep=timestep,
                    group_name=spike.group_name,
                    neuron_offset=spike.neuron_offset,
                ))
        return events
    
    def spikes_to_matrix(
        self,
        neuron_ids: Optional[Sequence[str]] = None,
    ) -> Tuple[np.ndarray, List[str]]:
        """
        Convert spike trace to binary spike matrix.
        
        neuron_ids: Optional list of neuron IDs to include. If None, all neurons that spiked at least once are included.
                
        Returns: Tuple of 2D numpy array of shape (timesteps, neurons) with 1s at spikes and list of neuron ID strings corresponding to columns
        """
        if not self.has_spikes():
            raise ValueError("No spike trace data available")
        
        # Collect all unique neuron IDs if not specified
        if neuron_ids is None:
            neuron_id_set = set()
            for spikes_at_t in self.spike_trace:
                for spike in spikes_at_t:
                    neuron_id_set.add(f"{spike.group_name}.{spike.neuron_offset}")
            neuron_ids = sorted(neuron_id_set)
        
        neuron_ids = list(neuron_ids)
        neuron_to_idx = {nid: idx for idx, nid in enumerate(neuron_ids)}
        
        n_timesteps = len(self.spike_trace)
        n_neurons = len(neuron_ids)
        
        matrix = np.zeros((n_timesteps, n_neurons), dtype=np.int8)
        
        for timestep, spikes_at_t in enumerate(self.spike_trace):
            for spike in spikes_at_t:
                neuron_id = f"{spike.group_name}.{spike.neuron_offset}"
                if neuron_id in neuron_to_idx:
                    matrix[timestep, neuron_to_idx[neuron_id]] = 1
        
        return matrix, neuron_ids
    
    def potentials_to_dataframe(
        self,
        neuron_ids: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        """
        Convert potential trace to pandas DataFrame.
        
        neuron_ids: Optional list of neuron ID strings to use as column names. If None, columns are named numerically (0, 1, 2, ...).
                
        Returns: DataFrame with timestep as index and neurons as columns.
        """
        if not self.has_potentials():
            raise ValueError("No potential trace data available")
        
        df = pd.DataFrame(self.potential_trace)
        df.index.name = "timestep"
        
        if neuron_ids is not None:
            if len(neuron_ids) != len(df.columns):
                raise ValueError(
                    f"Number of neuron_ids ({len(neuron_ids)}) doesn't match "
                    f"number of traced neurons ({len(df.columns)})"
                )
            df.columns = neuron_ids
        
        return df
    
    def potentials_to_array(self) -> np.ndarray:
        """
        Convert potential trace to numpy array.
        
        Returns: 2D numpy array of shape (timesteps, neurons).
        """
        if not self.has_potentials():
            raise ValueError("No potential trace data available")
        
        return np.array(self.potential_trace)
    
    def performance_to_dataframe(self) -> pd.DataFrame:
        """
        Convert performance trace to pandas DataFrame.
        
        Returns:
            DataFrame with columns for each performance metric:
                - timestep, fired, updated, packets, hops, spikes
                - sim_time, synapse_energy, dendrite_energy, soma_energy
                - network_energy, total_energy
        """
        if not self.has_performance():
            raise ValueError("No performance trace data available")
        
        if self._perf_df_cache is not None:
            return self._perf_df_cache.copy()
        
        self._perf_df_cache = pd.DataFrame(self.perf_trace)
        return self._perf_df_cache.copy()
    
    def messages_to_dataframe(self) -> pd.DataFrame:
        """
        Convert message trace to pandas DataFrame.
        
        Returns:
            DataFrame with columns for each message attribute:
                - timestep, mid, src_neuron, src_hw, dest_hw
                - hops, spikes, send_timestamp, received_timestamp
                - processed_timestamp, generation_delay, processing_delay
                - network_delay, blocking_delay, messages_along_route
        """
        if not self.has_messages():
            raise ValueError("No message trace data available")
        
        if self._message_df_cache is not None:
            return self._message_df_cache.copy()
        
        records = []
        for timestep, messages_at_t in enumerate(self.message_trace):
            for msg in messages_at_t:
                record = {"timestep": timestep}
                record.update(msg)
                records.append(record)
        
        self._message_df_cache = pd.DataFrame(records)
        return self._message_df_cache.copy()
    
    def get_neuron_groups(self) -> List[str]:
        """
        Get list of unique neuron group names from spike trace.
        
        Returns: Sorted list of group names.
        """
        if not self.has_spikes():
            return []
        
        groups = set()
        for spikes_at_t in self.spike_trace:
            for spike in spikes_at_t:
                groups.add(spike.group_name)
        
        return sorted(groups)
    
    def filter_by_groups(
        self,
        groups: Sequence[str],
    ) -> "TraceData":
        """
        Create a new TraceData with only spikes from specified groups.

        groups: List of group names to include.
            
        Returns: New TraceData instance with filtered spike data and other trace types are copied as-is.
        """
        groups_set = set(groups)
        
        filtered_spikes = None
        if self.has_spikes():
            filtered_spikes = []
            for spikes_at_t in self.spike_trace:
                filtered_at_t = [
                    s for s in spikes_at_t 
                    if s.group_name in groups_set
                ]
                filtered_spikes.append(filtered_at_t)
        
        return TraceData(
            spike_trace=filtered_spikes,
            potential_trace=self.potential_trace,
            perf_trace=self.perf_trace,
            message_trace=self.message_trace,
            neuron_labels=self.neuron_labels.copy(),
        )
    
    def filter_by_time(
        self,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> "TraceData":
        """
        Create a new TraceData filtered to a time range.
        
        start: Start timestep (inclusive). None means from beginning.
        end: End timestep (exclusive). None means to end.
            
        Returns: New TraceData instance with filtered time range.
        """
        start = start or 0
        end = end or self.timesteps
        
        filtered_spikes = None
        if self.has_spikes():
            filtered_spikes = self.spike_trace[start:end]
        
        filtered_potentials = None
        if self.has_potentials():
            filtered_potentials = self.potential_trace[start:end]
        
        filtered_messages = None
        if self.has_messages():
            filtered_messages = self.message_trace[start:end]
        
        # Performance trace needs special handling (dict of lists)
        filtered_perf = None
        if self.has_performance():
            filtered_perf = {}
            for key, values in self.perf_trace.items():
                filtered_perf[key] = values[start:end]
        
        return TraceData(
            spike_trace=filtered_spikes,
            potential_trace=filtered_potentials,
            perf_trace=filtered_perf,
            message_trace=filtered_messages,
            neuron_labels=self.neuron_labels.copy(),
        )
    
    def __repr__(self) -> str:
        parts = [f"TraceData(timesteps={self.timesteps}"]
        if self.has_spikes():
            total_spikes = sum(len(s) for s in self.spike_trace)
            parts.append(f"spikes={total_spikes}")
        if self.has_potentials():
            n_neurons = len(self.potential_trace[0]) if self.potential_trace else 0
            parts.append(f"potential_neurons={n_neurons}")
        if self.has_performance():
            parts.append("perf=True")
        if self.has_messages():
            total_msgs = sum(len(m) for m in self.message_trace)
            parts.append(f"messages={total_msgs}")
        return ", ".join(parts) + ")"
