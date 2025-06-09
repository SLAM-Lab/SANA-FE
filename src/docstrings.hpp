// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// docstrings.hpp - Documentation strings for Python bindings
#ifndef DOCSTRINGS_HEADER_INCLUDED_
#define DOCSTRINGS_HEADER_INCLUDED_

namespace docstrings
{

// Module level
constexpr const char *module_doc = R"pbdoc(
SANA-FE: Simulating Advanced Neuromorphic Architectures for Fast Exploration

A framework for modeling the energy usage and performance of different neuromorphic hardware.

Example:
    >>> import sanafe
    >>> arch = sanafe.load_arch("arch.yaml")
    >>> net = sanafe.load_net("snn.yaml")
    >>> chip = sanafe.SpikingChip(arch)
    >>> chip.load(net)
    >>> result = chip.sim(100)
)pbdoc";

// Network class
constexpr const char *network_doc = R"pbdoc(
Spiking Neural Network (SNN) container representing neurons and their connections.

A Network represents the graph of neurons and synapses, independent of any hardware
mapping. Use this to define your neural network (application) before mapping to
a neuromorphic hardware architecture. SANA-FE allows you to associate arbitrary
parameters with groups, neurons and edges, making it possible for you to pass
model attributes as user input.
)pbdoc";

constexpr const char *network_create_neuron_group_doc = R"pbdoc(
Create a new group of neurons with shared properties.

:param group_name: Unique identifier for this neuron group
:type group_name: str
:param neuron_count: Number of neurons to create in this group
:type neuron_count: int
:param model_attributes: Model parameters (e.g., threshold, leak). Default is None.
:type model_attributes: dict, optional
:param default_synapse_hw_name: Default synapse hardware type. Default is None.
:type default_synapse_hw_name: str, optional
:param default_dendrite_hw_name: Default dendrite hardware type. Default is None.
:type default_dendrite_hw_name: str, optional
:param force_dendrite_update: Force dendrite updates every timestep. Default is False.
:type force_dendrite_update: bool, optional
:param force_synapse_update: Force synapse updates every timestep. Default is False
:type force_synapase_update: bool, optional
:param force_soma_update: Force soma updates every timestep. Default is False.
:type force_soma_update: bool, optional
:param log_potential: Enable membrane potential logging. Default is False.
:type log_potential: bool, optional
:param log_spikes: Enable spike event logging. Default is False.
:type log_spikes: bool, optional
:param soma_hw_name: Soma hardware implementation name. Default is ""
:type soma_hw_name: str, optional

:return: The created neuron group
:rtype: NeuronGroup

Example:
    >>> group = net.create_neuron_group("layer1", 256,
    ...     model_attributes={"threshold": 1.0, "leak": 0.9})
)pbdoc";

constexpr const char *network_save_doc = R"pbdoc(
Save network to file in YAML or legacy netlist format.

:param path: Output file path
:type path: str
:param  use_netlist_format: Use legacy netlist format instead of YAML. Default is False.
:type use_netlist_format: bool, optional

)pbdoc";

// NeuronGroup class
constexpr const char *neuron_group_doc = R"pbdoc(
Collection of neurons with shared configuration and connectivity patterns.

NeuronGroups represent groups of similar or related neurons, akin to
'populations', or 'layers' on neurons in other SNN frameworks. NeuronGroups
provide efficient batch operations for creating connections and setting common
neuron parameters. Individual neurons can be accessed via indexing and slicing.
)pbdoc";

constexpr const char *neuron_group_connect_dense_doc = R"pbdoc(
Create fully-connected (dense) connections to destination group.

Connects every neuron in this group to every neuron in the destination
group. Connection attributes are provided as flattened arrays.

:param dest_group: Target neuron group
:type dest_group: NeuronGroup
:param attributes: Connection attributes with lists of values
    Length must equal source_neurons * dest_neurons
:type attributes: dict

Example:
    >>> weights = [random.random() for _ in range(src.neurons * dst.neurons)]
    >>> src.connect_neurons_dense(dst, {"weight": weights})
)pbdoc";

constexpr const char *neuron_group_connect_sparse_doc = R"pbdoc(
Create sparse connections using explicit source-destination pairs.

:param dest_group: Target neuron group
:type dest_group: NeuronGroup
:param attributes: Connection attributes with lists of values
:type attributes: dict
:param src_dest_id_pairs: List of (source_idx, dest_idx) tuples
:type src_dest_id_pairs: list

Example:
    >>> connections = [(0, 5), (1, 3), (2, 7)]
    >>> weights = [0.5, 0.8, 0.3]
    >>> src.connect_neurons_sparse(dst, {"weight": weights}, connections)
)pbdoc";

constexpr const char *neuron_group_connect_conv2d_doc = R"pbdoc(
Create 2D convolutional connections between neuron groups.

Implements CNN-style convolution with configurable kernel size, stride,
and channel dimensions. Input neurons are interpreted as flattened 2D/3D arrays.

:param dest_group: Output feature map neuron group
:type dest_group: NeuronGroup
:param attributes: Filter weights and other connection parameters
:type attributes: dict
:param input_width: Width of input feature map
:type input_width: int
:param input_height: Height of input feature map
:type input_height: int
:param input_channels: Number of input channels
:type input_channels: int
:param kernel_width: Convolution kernel width
:type kernel_width: int
:param kernel_height: Convolution kernel height
:type kernel_height: int
:param kernel_count: Number of output channels/filters. Default is 1.
:type kernel_count: int
:param stride_width: Horizontal stride. Default is 1.
:type stride_width: int, optional
:param stride_height: Vertical stride. Default is 1.
:type stride_height: int, optional

Example:
    >>> # 28x28x1 -> 26x26x32 convolution (3x3 kernels)
    >>> src.connect_neurons_conv2d(dst, {"weight": filter_weights},
    ...     28, 28, 1, 3, 3, 32, 1, 1)
)pbdoc";

// Neuron class
constexpr const char *neuron_doc = R"pbdoc(
Individual spiking neuron with configurable hardware mapping and parameters.

Neurons maintain their own state, connections, and hardware assignments.
Access via NeuronGroup indexing or iteration.
)pbdoc";

constexpr const char *neuron_set_attributes_doc = R"pbdoc(
Configure neuron-specific parameters and hardware assignments.

:param str soma_hw_name: Soma processing unit name. Default is None.
:param str default_synapse_hw_name: Default synapse type. Default is None.
:param str dendrite_hw_name: Dendrite processing unit name. Default is None.
:param bool log_spikes: Enable spike logging for this neuron. Default is False.
:param bool log_potential: Enable potential logging for this neuron. Default is False.
:param bool force_synapse_update: Force synapse updates every timestep. Default is False.
:param bool force_dendrite_update: Force dendrite updates every timestep. Default is False.
:param bool force_soma_update: Force soma updates every timestep. Default is False.
:param dict model_attributes: General model parameters. Default is None.
:param dict soma_attributes: Soma-specific parameters. Default is None.
:param dict dendrite_attributes: Dendrite-specific parameters. Default is None.
)pbdoc";

constexpr const char *neuron_map_to_core_doc = R"pbdoc(
Assign this neuron to a specific hardware core.

:param core_configuration: Target core for neuron placement
:type core_configuration: Core
)pbdoc";

constexpr const char *neuron_connect_to_neuron_doc = R"pbdoc(
Create direct connection to another neuron.

:param dest_neuron: Target neuron
:type dest_neuron: Neuron
:param attributes: Connection-specific parameters. Default is None.
:type attributes: dict

:return: Connection index for later reference
:rtype: int
)pbdoc";

// Architecture class
constexpr const char *architecture_doc = R"pbdoc(
Hardware architecture specification defining tiles, cores, and NoC topology.

Architecture describes the target neuromorphic hardware platform including
processing cores, network-on-chip parameters, and power/timing models.

Example:
    >>> arch = sanafe.Architecture("MyChip", noc_config)
    >>> tile = arch.create_tile("tile_0", energy_north_hop=1e-12)
    >>> core = arch.create_core("core_0", tile.id)
)pbdoc";

constexpr const char *architecture_create_tile_doc = R"pbdoc(
Add a new tile to the architecture with specified power metrics.

:param name: Unique tile identifier
:type name: str
:param energy_north_hop: Energy per northward hop (Joules). Default is 0.0.
:type energy_north_hop: float, optional
:param latency_north_hop: Latency per northward hop (seconds). Default is 0.0.
:type latency_north_hop: float, optional
:param energy_east_hop: Energy per eastward hop (Joules). Default is 0.0.
:type energy_east_hop: float, optional
:param latency_east_hop: Latency per eastward hop (seconds). Default is 0.0.
:type latency_east_hop: float, optional
:param energy_south_hop: Energy per southward hop (Joules). Default is 0.0.
:type energy_south_hop: float, optional
:param latency_south_hop: Latency per southward hop (seconds). Default is 0.0.
:type latency_south_hop: float, optional
:param energy_west_hop: Energy per westward hop (Joules). Default is 0.0.
:type energy_west_hop: float, optional
:param latency_west_hop: Latency per westward hop (seconds). Default is 0.0.
:type latency_west_hop: float, optional
:param log_energy: Enable tile energy logging. Default is False.
:type log_energy: bool, optional
:param log_latency: Enable tile latency logging. Default is False.
:type log_latency: bool, optional

:return: Reference to created tile configuration
:rtype: Tile
)pbdoc";

constexpr const char *architecture_create_core_doc = R"pbdoc(
Add a processing core to an existing tile.

:param name: Unique core identifier within tile
:type name: str
:param parent_tile_id: ID of containing tile
:type parent_tile_id: int
:param buffer_position: Pipeline buffer location. Default is before_soma.
:type buffer_position: BufferPosition, optional
:param buffer_inside_unit: Whether buffer is inside processing unit. Default is False.
:type buffer_inside_unit: bool, optional
:param max_neurons_supported: Maximum neurons per core. Default is None.
:type max_neurons_supported: int, optional
:param log_energy: Enable core energy logging. Default is False.
:type log_energy: bool, optional
:param log_latency: Enable core latency logging. Default is False.
:type log_latency: bool, optional

:return: Reference to created core configuration
:rtype: Core
)pbdoc";

// SpikingChip class
constexpr const char *spiking_chip_doc = R"pbdoc(
Simulated spiking chip executing a mapped SNN.

SpikingChip combines an Architecture with a mapped Network to provide
functional simulation with detailed energy and performance analysis. You create
a SpikingChip from an abstract Architecture object, and then load your mapped
Network application. Finally, you simulate the chip for one or more steps,
optionally with trace outputs.

Example:
    >>> chip = sanafe.SpikingChip(arch)
    >>> chip.load(network)
    >>> results = chip.sim(timesteps=100)
)pbdoc";

constexpr const char *spiking_chip_load_doc = R"pbdoc(
Map a neural network onto this chip's hardware architecture.

network (Network): Spiking network to map onto hardware
)pbdoc";

constexpr const char *spiking_chip_sim_doc = R"pbdoc(
Execute neuromorphic simulation for specified timesteps.

:param timesteps: Number of simulation timesteps. Default is 1.
:type timesteps: int, optional
:param timing_model: Timing model ("simple", "detailed", "cycle"). Default is "detailed".
:type timing_model: str, optional
:param processing_threads: Number of processing threads. Default is 1.
:type processing_threads: int, optional
:param scheduler_threads: Number of scheduler threads. Default is 0 (run in main thread).
:type scheduler_threads: int, optional
:param spike_trace: Spike trace output. Default is None.
:type spike_trace: str, file, bool, or None, optional
:param potential_trace: Potential trace output. Default is None.
:type potential_trace: str, file, bool, or None, optional
:param perf_trace: Performance metrics trace output. Default is None.
:type perf_trace: str, file, bool, or None, optional
:param message_trace: Message trace output. Default is None.
:type message_trace: str, file, bool, or None, optional
:param write_trace_headers: Write CSV headers to trace files. Default is True.
:type write_trace_headers: bool, optional

:return: Simulation results including energy, timing, and trace data:

    - timesteps_executed: Number of timesteps simulated
    - energy: Energy breakdown by component
    - sim_time: Simulated hardware time (seconds)
    - spikes: Total spike count
    - spike_trace: Spike data (if enabled)
    - potential_trace: Potential data (if enabled)
    - perf_trace: Performance metrics (if enabled)
    - message_trace: Network message data (if enabled)

:rtype: dict

Example:
    >>> results = chip.sim(timesteps=1000, spike_trace=True, perf_trace=True)
    >>> print(f"Total energy: {results['energy']['total']:.2e} J")
    >>> print(f"Spike count: {results['spikes']}")
)pbdoc";

constexpr const char *spiking_chip_reset_doc = R"pbdoc(
Reset all neuron states and hardware buffers to initial conditions.
)pbdoc";

constexpr const char *spiking_chip_get_power_doc = R"pbdoc(
Calculate average power consumption based on simulation results.

:return: Average power in Watts (total_energy / sim_time)
:rtype: float
)pbdoc";

// Connection class
constexpr const char *connection_doc = R"pbdoc(
Synaptic connection between two neurons with hardware and parameter mappings.
)pbdoc";

// NeuronAddress class
constexpr const char *neuron_address_doc = R"pbdoc(
Unique identifier for a neuron within the network.

Combines group name and offset to provide unambiguous neuron addressing
for connectivity and tracing.
)pbdoc";

// Utility functions
constexpr const char *load_arch_doc = R"pbdoc(
Load hardware architecture from YAML configuration file.

An architecture is an abstract hardware configuration. To simulate this
hardware, you must use the configuration to create one or more SpikingChip
objects.

:param path: Path to architecture YAML file
:type path: str

:return: Loaded architecture configuration
:rtype: Architecture

Example:
    >>> arch = sanafe.load_arch("loihi.yaml")
)pbdoc";

constexpr const char *load_net_doc = R"pbdoc(
Load neural network from file (YAML or legacy netlist format).

:param path: Path to network file
:type path: str
:param arch: Target architecture for validation
:type arch: Architecture
:param use_netlist_format: Parse as legacy netlist format. Default is False.
:type use_netlist_format: bool, optional

:return: Loaded spiking neural network
:rtype: Network

Example:
    >>> net = sanafe.load_net("snn.yaml", arch)
)pbdoc";

} // namespace docstrings

#endif