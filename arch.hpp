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

#ifndef ARCH_HEADER_INCLUDED_
#define ARCH_HEADER_INCLUDED_

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <list>

#include "network.hpp"
#include "description.hpp"

using namespace std;

// Hard define maximum defined h/w sizes
// TODO: better dynamically define or allocate these numbers, so that we can
//  support a range of architectures seamlessly. At the moment, a large amount
//  of memory is needed if we want to support lots of large cores
#define ARCH_MAX_TILES 32
#define ARCH_MAX_CORES_PER_TILE 4
#define ARCH_MAX_CORES 128
//#define ARCH_MAX_X (64)
//#define ARCH_MAX_Y (64)
//#define ARCH_MAX_COMPARTMENTS 256

// Loihi
#define ARCH_MAX_X (8*4)
#define ARCH_MAX_Y (4*2)

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

enum synapse_models
{
	SYNAPSE_CUBA,
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

struct Message
{
	Neuron *src_neuron;
	Message *next;
	double generation_delay, network_delay, receive_delay;
	double blocked_latency;
	double sent_timestamp, received_timestamp, processed_timestamp;
	long int timestep;
	int spikes, hops;
	int src_x, dest_x, src_y, dest_y;
	int dest_tile_id, dest_core_offset, dest_axon_id;
	bool dummy_message, in_noc;

	Message();
};

struct MessageFifo
{
	int count;
	Message *head, *tail;
	struct MessageFifo *next;  // For priority queue of core fifos
};

enum NoiseType
{
	NOISE_NONE = -1,
	NOISE_FILE_STREAM,
	// TODO: implement different random noise generators
};

struct AxonInUnit
{
	std::string name;
	long int spike_messages_in;
	double energy, time;
	double energy_spike_message, latency_spike_message;
	int parent_tile_id;
};

struct SynapseUnit
{
	std::string name;
	int model, weight_bits;
	long int spikes_processed;
	double energy, time;
	double energy_spike_op, energy_memory_access;
	double latency_spike_op, latency_memory_access;
};

struct DendriteUnit
{
	std::string name;
	double energy, time;
};

struct SomaUnit
{
	FILE *noise_stream;
	std::string name;
	long int neuron_updates, neurons_fired, neuron_count;
	double energy, time;
	double energy_update_neuron, latency_update_neuron;
	double energy_access_neuron, latency_access_neuron;
	double energy_spiking, latency_spiking;
	int leak_towards_zero, reset_mode, reverse_reset_mode;
	int noise_type;
};

struct AxonOutUnit
{
	// The axon output points to a number of axons, stored at the
	//  post-synaptic core. A neuron can point to a number of these
	std::string name;
	long int packets_out;
	double energy, time;
	double energy_access, latency_access;
	int parent_tile_id;
};

struct AxonInModel
{
	// List of all neuron connections to send spike to
	AxonInUnit *axon_in_unit;
	std::vector<int> synapse_addresses;

	Message *message;
	long int last_updated;
	int spikes_received, active_synapses;
};

struct AxonOutModel
{
	// List of all neuron connections to send spike to
	AxonOutUnit *axon_out_unit;
	int dest_axon_id, dest_tile_id, dest_core_offset;
	int src_neuron_id;
};

struct Core
{
	std::vector<AxonInUnit> axon_in_hw;
	std::vector<SynapseUnit> synapse;
	std::vector<DendriteUnit> dendrite;
	std::vector <SomaUnit> soma;
	std::vector<AxonOutUnit> axon_out_hw;

	std::vector<Message *> messages_in;
	std::vector<AxonInModel> axons_in;
	std::vector<Neuron *> neurons;
	std::vector<Connection *> synapses;
	std::vector<AxonOutModel> axons_out;

	std::string name;
	Message next_message;  // Since last spike
	size_t max_neurons;
	double energy, latency_after_last_message;
	int parent_tile_id, id, offset, buffer_pos, soma_count, synapse_count;
	int message_count;
};

struct Tile
{
	std::vector<Core> cores;
	std::string name;
	double energy;
	double energy_east_hop, latency_east_hop;
	double energy_west_hop, latency_west_hop;
	double energy_north_hop, latency_north_hop;
	double energy_south_hop, latency_south_hop;
	long int hops, messages_received, total_neurons_fired;
	long int east_hops, west_hops, north_hops, south_hops;
	int id, x, y;
	int width; // For now just support 2 dimensions
};

class Architecture
{
public:
	std::vector<Tile> tiles;
	std::string name;
	int noc_width, noc_height, noc_buffer_size;
	bool noc_init;

	// TODO: get rid of these
	bool spike_vector_on;
	std::vector<std::vector<int>> spike_vector;

	Architecture();
	int get_core_count();
};

#include "description.hpp"

void free_spike_vector(Architecture* arch);
void arch_free(Architecture *const arch);
int arch_create_noc(Architecture &arch, const std::vector<Attribute> &attr);
int arch_create_tile(Architecture &arch, const std::vector<Attribute> &attr);
int arch_create_core(Architecture &arch, Tile &t, const std::vector<Attribute> &attr);
void arch_create_axon_in(Core &c, const std::string &name, const std::vector<Attribute> &attr);
void arch_create_synapse(Core &c, const std::string &name, const std::vector<Attribute> &attr);
void arch_create_soma(Core &c, const std::string &name, const std::vector<Attribute> &attr);
void arch_create_axon_out(Core &c, const std::vector<Attribute> &attr);
void arch_create_axons(Architecture &arch);
void arch_create_core_axon(Core &core);
void arch_print_axon_summary(Architecture &arch);
int arch_map_neuron(Neuron &n, Core &c);
void arch_map_neuron_connections(Neuron &n);
void arch_allocate_axon(Neuron &pre_neuron, Core &post_core);
void arch_add_connection_to_axon(Connection &con, Core &post_core);
int arch_parse_neuron_model(const std::string &model_str);
int arch_parse_synapse_model(const std::string &model_str);
void arch_init_message(Message &m);

#endif
