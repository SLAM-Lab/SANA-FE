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
axon input -> synapse --------> dendrite ------> soma -------> axon output
(spikes in) (spikes to current)(process input)(membrane update)(spikes out)
*/

#ifndef ARCH_HEADER_INCLUDED_
#define ARCH_HEADER_INCLUDED_

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "network.h"
#include "description.h"

// Hard define maximum defined h/w sizes
#define ARCH_MAX_COMPARTMENTS 16384
#define ARCH_MAX_CONNECTION_MAP (ARCH_MAX_COMPARTMENTS*4)
// TODO: better dynamically define or allocate these numbers, so that we can
//  support a range of architectures seamlessly. At the moment, a large amount
//  of memory is needed if we want to support lots of large cores
// TrueNorth
//#define ARCH_MAX_TILES 4096
//#define ARCH_MAX_CORES_PER_TILE 1
// Loihi
#define ARCH_MAX_TILES 256
//#define ARCH_MAX_TILES 32
#define ARCH_MAX_CORES_PER_TILE 4
#define ARCH_MAX_UNITS 3
#define ARCH_MAX_CORES (ARCH_MAX_TILES * ARCH_MAX_CORES_PER_TILE)

#define ARCH_MAX_LINKS 4
#define ARCH_MAX_DESCRIPTION_LINE 256
#define ARCH_MAX_ATTRIBUTES 256

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
	NEURON_LIF,
	NEURON_STOCHASTIC_LIF,
	NEURON_TRUENORTH,
};

enum synapse_models
{
	SYNAPSE_CUBA,
};


enum neuron_reset_modes
{
	NEURON_NO_RESET,
	NEURON_RESET_SOFT,
	NEURON_RESET_HARD,
	NEURON_RESET_SATURATE,
	NEURON_RESET_MODE_COUNT,
};

enum arch_description_blocks
{
	ARCH_DESCRIPTION_TOP = -1,
	ARCH_DESCRIPTION_ARCH,
	ARCH_DESCRIPTION_TILE,
	ARCH_DESCRIPTION_CORE,
	ARCH_DESCRIPTION_AXON_IN,
	ARCH_DESCRIPTION_SYNAPSE,
	ARCH_DESCRIPTION_DENDRITE,
	ARCH_DESCRIPTION_SOMA,
	ARCH_DESCRIPTION_AXON_OUT,
};

enum noise_type
{
	NOISE_NONE = -1,
	NOISE_FILE_STREAM,
	// TODO: implement different random noise generators
};

struct message
{
	struct neuron *src_neuron, *dest_neuron;
	struct message *next;
	double generation_latency, network_latency, receive_latency;
	double blocked_latency, sent_timestamp, processed_timestamp;
	long int timestep;
	int spikes, hops;
};

struct message_fifo
{
	int count;
	struct message *head, *tail;
	struct message_fifo *next;  // For priority queue of core fifos
};

struct connection_map
{
	// List of all neuron connections to send spike to
	struct connection **connections;
	struct message *message;
	struct neuron *pre_neuron;
	long int last_updated;
	int connection_count, spikes_received, active_synapses;
};

struct axon_input
{
	char name[MAX_FIELD_LEN];
	struct tile *t;
	struct connection_map map[ARCH_MAX_CONNECTION_MAP];
	long int spike_messages_in;
	double energy, time;
	double energy_spike_message, latency_spike_message;
	int map_count;
};

struct synapse_processor
{
	char name[MAX_FIELD_LEN];
	int model, spikes_buffer, weight_bits;
	long int spikes_processed;
	double energy, time;
	double energy_spike_op, energy_memory_access;
	double latency_spike_op, latency_memory_access;
};

struct dendrite_processor
{
	char name[MAX_FIELD_LEN];
	double energy, time;
};

struct soma_processor
{
	FILE *noise_stream;
	char name[MAX_FIELD_LEN];
	int model, leak_towards_zero, reset_mode, reverse_reset_mode;
	int noise_type;
	long int neuron_updates, neurons_fired, neuron_count;
	double energy, time;
	double energy_update_neuron, latency_update_neuron;
	double energy_access_neuron, latency_access_neuron;
	double energy_spiking, latency_spiking;
};

struct axon_output
{
	// The axon output points to a number of axons, stored at the
	//  post-synaptic core. A neuron can point to a number of these
	struct connection_map *map_ptr[ARCH_MAX_CONNECTION_MAP];
	struct tile *t;
	int map_count;
	char name[MAX_FIELD_LEN];

	long int packets_out;
	double energy, time;
	double energy_access, latency_access;
};

struct core
{
	struct tile *t;
	struct neuron **neurons;

	struct axon_input axon_in;
	struct synapse_processor synapse[ARCH_MAX_UNITS];
	struct dendrite_processor dendrite;
	struct soma_processor soma[ARCH_MAX_UNITS];
	struct axon_output axon_out;

	char name[MAX_FIELD_LEN];
	struct message next_message;  // Since last spike
	double energy, blocked_until, latency_after_last_message;
	int id, offset, buffer_pos, is_blocking, soma_count, synapse_count;
	int neuron_count, message_count;
	int curr_axon;
};

struct tile
{
	struct core cores[ARCH_MAX_CORES_PER_TILE];
	struct tile *links[ARCH_MAX_LINKS];
	char name[MAX_FIELD_LEN];
	double energy;
	double energy_east_hop, latency_east_hop;
	double energy_west_hop, latency_west_hop;
	double energy_north_hop, latency_north_hop;
	double energy_south_hop, latency_south_hop;
	double blocked_until;
	long int hops, messages_received, total_neurons_fired;
	long int east_hops, west_hops, north_hops, south_hops;
	int id, x, y, core_count, is_blocking;
	int width; // For now just support 2 dimensions
};

struct architecture
{
	struct tile tiles[ARCH_MAX_TILES];
	char name[MAX_FIELD_LEN];
	int noc_width, noc_height, tile_count, core_count;
	int is_init;
};

#include "description.h"

struct architecture *arch_init(void);
void arch_free(struct architecture *const arch);
int arch_create_noc(struct architecture *const arch, struct attributes *attr, const int attribute_count);
int arch_create_tile(struct architecture *const arch, struct attributes *attr, const int attribute_count);
int arch_create_core(struct architecture *const arch, struct tile *const t, struct attributes *attr, const int attribute_count);
void arch_create_axon_in(struct core *const c, const char *const name, const struct attributes *const attr, const int attribute_count);
void arch_create_synapse(struct core *const c, const char *const name, const struct attributes *const attr, const int attribute_count);
void arch_create_soma(struct core *const c, const char *const name, struct attributes *attr, const int attribute_count);
void arch_create_axon_out(struct core *const c, struct attributes *attr, const int attribute_count);
void arch_create_connection_maps(struct architecture *const arch);
void arch_create_core_connection_map(struct core *const core);
void arch_print_connection_map_summary(struct architecture *const arch);
int arch_map_neuron(struct neuron *const n, struct core *c);
void arch_map_neuron_connections(struct neuron *const n);
void arch_allocate_connection_map(struct neuron *const pre_neuron, struct core *const post_core, const int connection_count);
void arch_add_connection_to_map(struct connection *const con, struct core *const post_core);
int arch_parse_neuron_model(const char *model_str);
int arch_parse_synapse_model(const char *model_str);
void arch_init_message(struct message *m);

#endif
