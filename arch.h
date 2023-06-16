// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
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

#include <stdint.h>
#include "network.h"

// Hard define maximum defined h/w sizes
#define ARCH_MAX_COMPARTMENTS 1024
//#define ARCH_MAX_AXON_MAP 4096
#define ARCH_MAX_AXON_MAP 16384
#define ARCH_MAX_TILES 32
#define ARCH_MAX_CORES 4
#define ARCH_MAX_LINKS 4

#define ARCH_INVALID_ID -1
enum buffer_positions
{
	BUFFER_AXON_IN = 0,	// Buffer incoming packets
	BUFFER_SYNAPSE,		// Buffer synaptic addresses
	BUFFER_DENDRITE,	// Buffer synaptic current
	BUFFER_SOMA,		// Buffer dendritic current
	BUFFER_AXON_OUT,	// Buffer axon addresses i.e. spikes out
	BUFFER_NETWORK,		// Buffer messages to the network
	BUFFER_POSITIONS,
};

enum neuron_models
{
	NEURON_IF,
	NEURON_LIF,
	NEURON_TRUENORTH,
};

enum neuron_reset_modes
{
	NEURON_NO_RESET,
	NEURON_RESET_SOFT,
	NEURON_RESET_HARD,
	NEURON_RESET_SATURATE,
	NEURON_RESET_MODE_COUNT,
};

struct hardware_mapping
{
	struct core *core;
	struct synapse_processor *synapse_hw;
	struct dendrite_processor *dendrite_hw;
	struct soma_processor *soma_hw;
	struct axon_output *axon_out;
	struct axon_input *axon_in;
};

struct axon_map
{
	// List of all neuron connections to send spike to
	int connection_count, spikes_received;
	long int last_updated;
	double network_latency, receive_latency;
	struct connection **connections;
	struct neuron *pre_neuron;
};

struct axon_input
{
	struct tile *t;
	struct axon_map map[ARCH_MAX_AXON_MAP];
	long int packets_in;
	double energy, time;
	int packet_size, packets_buffer, spikes_buffer, map_count;
};

struct synapse_processor
{
	int spikes_buffer;
	int weights_per_word, word_bits, weight_bits;
	long int total_spikes, memory_reads;
	double energy, time;
	double energy_spike_op, energy_memory_access;
	double time_spike_op, time_memory_access;
};

struct dendrite_processor
{
	double energy, time;
};

struct soma_processor
{
	int model, leak_towards_zero, reset_mode, reverse_reset_mode;
	long int updates, spikes_sent;
	double energy, time;
	double energy_active_neuron_update, time_active_neuron_update;
	double energy_inactive_neuron_update, time_inactive_neuron_update;
	double energy_spiking, time_spiking;
};

struct axon_output
{
	// The axon output points to a number of axons, stored at the
	//  post-synaptic core. A neuron can point to a number of these
	struct axon_map *map_ptr[ARCH_MAX_AXON_MAP];
	struct tile *t;
	int map_count;

	long int packets_out;
	double energy, time;
	double energy_access, time_access;
};

struct core
{
	struct tile *t;
	struct core *next_timing;
	struct neuron **neurons;

	struct axon_input axon_in;
	struct synapse_processor synapse;
	struct dendrite_processor dendrite;
	struct soma_processor soma;
	struct axon_output axon_out;

	double energy, time, blocked_until;
	int id, buffer_pos, is_blocking;
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
	double blocked_until;
	int id, x, y, core_count, is_blocking;
	int max_dimensions, width; // For now just support 2 dimensions
};

struct architecture
{
	struct tile tiles[ARCH_MAX_TILES];
	double time_barrier;
	int tile_count, initialized;
};

struct architecture *arch_init(void);
void arch_free(struct architecture *const arch);
int arch_create_noc(struct architecture *const arch, const int width, const int height);
int arch_create_tile(struct architecture *const arch, const double energy_east_west_hop, const double energy_north_south_hop, const double time_east_west_hop, const double time_north_south_hop);
int arch_create_core(struct architecture *const arch, struct tile *const t);
void arch_create_axon_in(struct architecture *const arch, struct core *const c);
void arch_create_synapse(struct architecture *const arch, struct core *const c, const int weight_bits, const int word_bits, const double energy_spike_op, const double time_spike_op, const double energy_memory_access, const double time_memory_access);
void arch_create_soma(struct architecture *const arch, struct core *const c, int model, double energy_active_neuron_update, double time_active_neuron_update, double energy_inactive_neuron_update, double time_inactive_neuron_update, double energy_spiking, double time_spiking);
void arch_create_axon_out(struct architecture *const arch, struct core *const c, const double spike_energy, const double spike_time);
int arch_map_neuron(struct neuron *const n, const struct hardware_mapping);
void arch_create_axon_maps(struct architecture *const arch);
void arch_create_core_axon_map(struct core *const core);

#endif
