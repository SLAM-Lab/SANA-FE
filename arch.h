// arch.h: Create a neuromorphic design based on a set of commands
//  In this simulator an architecture is a represented as a set of different
//  hardware blocks. The design (chip) is a set of tiles, connected by NoC
//  interconnect. Within each tile is one or more cores. Each core contains
//  neuromorphic computation. The neuromorphic pipeline (which seems sufficient
//  for any design) is a series of elements:
/*
axon inputs->synapse processor->dendrite processor->soma processor->axon outputs
(spikes in) (spikes to current) (process input)    (membrane update)(spikes out)
*/
// Note importantly that a single processor might handle a bunch of neurons
//  i.e. it is *not* necessarily a 1-1 mapping. I will provide some different
//  implementations of each component. It is possible to extend though with
//  custom elements to do whatever you want
//  element.
#ifndef ARCH_HEADER_INCLUDED_
#define ARCH_HEADER_INCLUDED_

// TODO: for now hard define max numbers of hardware blocks. I did this so I
//  don't have quite as much dynamic allocation going on. If needed this can
//  fairly easily be removed and replaced with allocating arbitrary numbers of
//  elements off the heap e.g. in a linked list
#define ARCH_MAX_COMPARTMENTS
#define ARCH_MAX_TILES 128
#define ARCH_MAX_CORES 16
#define ARCH_MAX_LINKS 4

#define ARCH_INVALID_ID -1

enum buffer_positions
{
	BUFFER_AXON_IN = 0,	// Buffer incoming packets
	BUFFER_SYNAPSE,		// Buffer synaptic addresses
	BUFFER_DENDRITE,	// Buffer synaptic current
	BUFFER_SOMA,		// Buffer dendritic current
	BUFFER_AXON_OUT,	// Buffer axon addresses i.e. spikes out
	BUFFER_POSITIONS,
};

enum neuron_models
{
	NEURON_IF,
	NEURON_LIF,
};

struct axon_input
{
	struct tile *t;
	int packet_size, buffer_in, packets_buffer, spikes_buffer;
	long int packets_in;
	double energy, time;
};

struct synapse_processor
{
	int buffer_in, spikes_buffer;
	int weights_per_word, word_bits, weight_bits;
	long int total_spikes, memory_reads;
	double energy, time, busy_until;
	double energy_spike_op, energy_memory_access;
	double time_spike_op, time_memory_access;
};

struct dendrite_processor
{
	int buffer_in;
	double energy, time;
};

struct soma_processor
{
	int buffer_in, model;
	long int updates, spikes_sent;
	double energy, time;
	double energy_active_neuron_update, time_active_neuron_update;
	double energy_inactive_neuron_update, time_inactive_neuron_update;
	double energy_spiking, time_spiking;
};

struct axon_output
{
	struct tile *t;
	int buffer_in;
	long int packets_out;
	double energy, time;
	double energy_access, time_access;
};

struct core
{
	struct tile *t;
	struct core *next_timing;
	struct neuron **neurons;
	struct core **axon_map;
	double **message_processing_time;
	int *spikes_sent_per_core;
	struct axon_input axon_in;
	struct synapse_processor synapse;
	struct dendrite_processor dendrite;
	struct soma_processor soma;
	struct axon_output axon_out;
	double energy, time;
	int id, buffer_pos;
	int neuron_count, curr_neuron, neurons_left, status;
	int curr_axon;
};

struct tile
{
	struct core cores[ARCH_MAX_CORES];
	// TODO: maybe can associate energy and latency with each link, that
	//  will be the most general way to implement this!
	struct tile *links[ARCH_MAX_LINKS];
	double energy, time;
	double energy_east_west_hop, time_east_west_hop;
	double energy_north_south_hop, time_north_south_hop;
	double energy_spike_within_tile, time_spike_within_tile;
	double busy_until;
	int id, x, y, core_count;
	int max_dimensions, width; // For now just support 2 dimensions
};

struct architecture
{
	struct tile tiles[ARCH_MAX_TILES];
	double time_barrier;
	int tile_count, initialized;

	long int total_hops;
};

void arch_init(struct architecture *const arch);
void arch_free(struct architecture *const arch);
int arch_create_noc(struct architecture *const arch, const int width, const int height);
int arch_create_tile(struct architecture *const arch, const double energy_east_west_hop, const double energy_north_south_hop, const double time_east_west_hop, const double time_north_south_hop);
int arch_create_core(struct architecture *const arch, struct tile *const t);
void arch_create_axon_in(struct architecture *const arch, struct core *const c);
void arch_create_synapse(struct architecture *const arch, struct core *const c, const int weight_bits, const int word_bits, const double energy_spike_op, const double time_spike_op, const double energy_memory_access, const double time_memory_access);
void arch_create_soma(struct architecture *const arch, struct core *const c, int model, double energy_active_neuron_update, double time_active_neuron_update, double energy_inactive_neuron_update, double time_inactive_neuron_update, double energy_spiking, double time_spiking);
void arch_create_axon_out(struct architecture *const arch, struct core *const c, const double spike_energy, const double spike_time);

#endif
