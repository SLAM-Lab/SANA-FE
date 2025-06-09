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

Args:
    group_name (str): Unique identifier for this neuron group
    neuron_count (int): Number of neurons to create in this group
    model_attributes (dict, optional): Model parameters (e.g., threshold, leak). Default is None.
    default_synapse_hw_name (str, optional): Default synapse hardware type. Default is None.
    default_dendrite_hw_name (str, optional): Default dendrite hardware type. Default is None.
    force_dendrite_update (bool, optional): Force dendrite updates every timestep. Default is False.
    force_synapse_update (bool, optional): Force synapse updates every timestep. Default is False.
    force_soma_update (bool, optional): Force soma updates every timestep. Default is False.
    log_potential (bool, optional): Enable membrane potential logging. Default is False.
    log_spikes (bool, optional): Enable spike event logging. Default is False.
    soma_hw_name (str, optional): Soma hardware implementation name. Default is None.

Returns:
    NeuronGroup: The created neuron group

Example:
    >>> group = net.create_neuron_group("layer1", 256,
    ...     model_attributes={"threshold": 1.0, "leak": 0.9})
)pbdoc";

constexpr const char *network_save_doc = R"pbdoc(
Save network to file in YAML or legacy netlist format.

Args:
    dest_group (NeuronGroup): Target neuron group
    attributes (dict): Connection attributes with lists of values
        Length must equal source_neurons * dest_neurons

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

Args:
    dest_group (NeuronGroup): Target neuron group
    attributes (dict): Connection attributes with lists of values
    src_dest_id_pairs (list): List of (source_idx, dest_idx) tuples

Example:
    >>> weights = [random.random() for _ in range(src.neurons * dst.neurons)]
    >>> src.connect_neurons_dense(dst, {"weight": weights})
)pbdoc";

constexpr const char *neuron_group_connect_sparse_doc = R"pbdoc(
Create sparse connections using explicit source-destination pairs.

Args:
    dest_group (NeuronGroup): Target neuron group
    attributes (dict): Connection attributes with lists of values
        Length must equal source_neurons * dest_neurons
Example:
    >>> connections = [(0, 5), (1, 3), (2, 7)]
    >>> weights = [0.5, 0.8, 0.3]
    >>> src.connect_neurons_sparse(dst, {"weight": weights}, connections)
)pbdoc";

constexpr const char *neuron_group_connect_conv2d_doc = R"pbdoc(
Create 2D convolutional connections between neuron groups.

Implements CNN-style convolution with configurable kernel size, stride,
and channel dimensions. Input neurons are interpreted as flattened 2D/3D arrays.

Args:
    dest_group (NeuronGroup): Output feature map neuron group
    attributes (dict): Filter weights and other connection parameters
    input_width (int): Width of input feature map
    input_height (int): Height of input feature map
    input_channels (int): Number of input channels
    kernel_width (int): Convolution kernel width
    kernel_height (int): Convolution kernel height
    kernel_count (int, optional): Number of output channels/filters. Default is 1.
    stride_width (int, optional): Horizontal stride. Default is 1.
    stride_height (int, optional): Vertical stride. Default is 1.

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

Args:
    soma_hw_name (str, optional): Soma processing unit name. Default is None.
    default_synapse_hw_name (str, optional): Default synapse type. Default is None.
    dendrite_hw_name (str, optional): Dendrite processing unit name. Default is None.
    log_spikes (bool, optional): Enable spike logging for this neuron. Default is False.
    log_potential (bool, optional): Enable potential logging for this neuron. Default is False.
    force_synapse_update (bool, optional): Force synapse updates every timestep. Default is False.
    force_dendrite_update (bool, optional): Force dendrite updates every timestep. Default is False.
    force_soma_update (bool, optional): Force soma updates every timestep. Default is False.
    model_attributes (dict, optional): General model parameters. Default is None.
    soma_attributes (dict, optional): Soma-specific parameters. Default is None.
    dendrite_attributes (dict, optional): Dendrite-specific parameters. Default is None.

)pbdoc";

constexpr const char *neuron_map_to_core_doc = R"pbdoc(
Assign this neuron to a specific hardware core.

:param core_configuration: Target core for neuron placement
:type core_configuration: Core
)pbdoc";

constexpr const char *neuron_connect_to_neuron_doc = R"pbdoc(
Create direct connection to another neuron.

Args:
    dest_neuron (Neuron): Target neuron
    attributes (dict, optional): Connection-specific parameters

Returns:
    int: Connection index for later reference

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

Args:
    name (str): Unique tile identifier
    energy_north_hop (float, optional): Energy per northward hop (Joules). Default is 0.0.
    latency_north_hop (float, optional): Latency per northward hop (seconds). Default is 0.0.
    energy_east_hop (float, optional): Energy per eastward hop (Joules). Default is 0.0.
    latency_east_hop (float, optional): Latency per eastward hop (seconds). Default is 0.0.
    energy_south_hop (float, optional): Energy per southward hop (Joules). Default is 0.0.
    latency_south_hop (float, optional): Latency per southward hop (seconds). Default is 0.0.
    energy_west_hop (float, optional): Energy per westward hop (Joules). Default is 0.0.
    latency_west_hop (float, optional): Latency per westward hop (seconds). Default is 0.0.
    log_energy (bool, optional): Enable tile energy logging. Default is False.
    log_latency (bool, optional): Enable tile latency logging. Default is False.

Returns:
    Tile: Reference to created tile configuration
)pbdoc";

constexpr const char *architecture_create_core_doc = R"pbdoc(
Add a processing core to an existing tile.


Args:
    name (str): Unique core identifier within tile
    parent_tile_id (int): ID of containing tile
    buffer_position (BufferPosition, optional): Pipeline buffer location. Default is before soma.
    buffer_inside_unit (bool, optional): Whether buffer is inside processing unit. Default is False.
    max_neurons_supported (int, optional): Maximum neurons per core. Default is None.
    log_energy (bool, optional): Enable core energy logging. Default is False.
    log_latency (bool, optional): Enable core latency logging. Default is False.

Returns:
    Core: Reference to created core configuration
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

Args:
    timesteps (int, optional): Number of simulation timesteps
    timing_model (str, optional): Timing model ("simple", "detailed", "cycle"). Default is "detailed".
    processing_threads (int, optional): Number of processing threads. Default is 1.
    scheduler_threads (int, optional): Number of scheduler threads. Default is 0 (run in main thread).
    spike_trace (object, optional): Spike trace output (file, True, or None). Default is None.
    potential_trace (object, optional): Potential trace output. Default is None.
    perf_trace (object, optional): Performance metrics trace output. Default is None.
    message_trace (object, optional): Message trace output. Default is None.
    write_trace_headers (bool, optional): Write CSV headers to trace files. Default is True.

Returns:
    dict: Simulation results including energy, timing, and trace data
        - timesteps_executed: Number of timesteps simulated
        - energy: Energy breakdown by component
        - sim_time: Simulated hardware time (seconds)
        - spikes: Total spike count
        - spike_trace: Spike data (if enabled)
        - potential_trace: Potential data (if enabled)
        - perf_trace: Performance metrics (if enabled)
        - message_trace: Network message data (if enabled)

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

Returns:
    float: Average power in Watts (total_energy / sim_time)
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

Args:
    path (str): Path to architecture YAML file

Returns:
    Architecture: Loaded architecture configuration

Example:
    >>> arch = sanafe.load_arch("loihi.yaml")
)pbdoc";

constexpr const char *load_net_doc = R"pbdoc(
Load neural network from file (YAML or legacy netlist format).

Args:
    path (str): Path to network file
    arch (Architecture): Target architecture for validation
    use_netlist_format (bool, optional): Parse as legacy netlist format

Returns:
    Network: Loaded spiking neural network

Example:
    >>> net = sanafe.load_net("snn.yaml", arch)
)pbdoc";

} // namespace docstrings

#endif