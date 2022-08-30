// network.h - (spiking) neural network functionality. Spiking neural
//  networks are represented as groups of neurons. A neuron group might have a
//  bunch of neurons all with the same properties (and common hardware).
//  Each neuron has its own state and a set of connections to other neurons.
//  These structures have optional links to hardware for performance simulation.
//  Here we include different neuron, synapse and dendrite models.
#ifndef NETWORK_HEADER_INCLUDED_
#define NETWORK_HEADER_INCLUDED_

#define NETWORK_MAX_NEURON_GROUPS 1024
#define NETWORK_INVALID_NID -1

// This defines the structure of the format for specifying neurons in a
//  network for this simulator. Each row in the file represents data for a
//  unique neuron in the network
//
// The format is a number of neuron fields, followed by a variable number of
//  connections (each with CONNECTION_FIELDS entries).  This means each command
//  is variable length depending on the number of connections for that neuron.
enum neuron_config_format
{
	NEURON_COMMAND = 0,
	NEURON_GROUP_ID,
	NEURON_ID,
	NEURON_RECORD_SPIKES, // Record any spikes to a CSV file
	NEURON_RECORD_VOLTAGE, // Write the voltage to a CSV file every timestep
	NEURON_FORCE_UPDATE,
	NEURON_FIELDS,
};

enum connection_config_format
{
	// This is the structure of the CSV format for specifying synaptic
	//  connections in a network for this simulator.  Each row represents
	//  a unique connection.
	CONNECTION_DEST_GID = 0,
	CONNECTION_DEST_NID,
	CONNECTION_WEIGHT,
	CONNECTION_FIELDS,
};

enum input_types
{
	INPUT_EVENT,
	INPUT_RATE,
	INPUT_POISSON,
};

struct neuron_id
{
	unsigned int group, neuron;
};

struct neuron
{
	struct connection *connections;
	struct neuron_group *group;

	double potential, current,  bias, reset, threshold;
	double potential_decay, potential_time_const;
	double current_decay, current_time_const;
	int id, is_init, fired, post_connection_count, spike_count;
	int log_spikes, log_voltage, update_needed, force_update;
};

struct connection
{
	struct neuron *post_neuron, *pre_neuron;
	double weight;
	// TODO: delay, and any other synapse/connection properties
	int id;
};

struct input
{
	struct connection *connections;
	double rate, spike_val; // rate only applies to poisson and rate-based
	int post_connection_count, type, id;
};

struct neuron_group
{
	// A neuron group is a collection of neurons that share common
	//  parameters (and possibly hardware). All neurons must be based on the
	//  same models. If implemented in hardware, they also must share common
	//  hardware i.e. the same core and processor blocks in the core
	struct neuron *neurons;

	// Mapped hardware, which is optionally set
	struct core *core;
	struct synapse_processor *synapse;
	struct dendrite_processor *dendrite;
	struct soma_processor *soma;
	struct axon_output *axon_out;
	struct axon_input *axon_in;

	int id, neuron_count;
	double default_threshold, default_reset;
};

struct hardware_mapping
{
	struct core *core;
	struct synapse_processor *synapse;
	struct dendrite_processor *dendrite;
	struct soma_processor *soma;
	struct axon_output *axon_out;
	struct axon_input *axon_in;
};

struct network
{
	struct neuron_group groups[NETWORK_MAX_NEURON_GROUPS];
	struct input *external_inputs;
	int neuron_group_count, external_input_count;
};

#include "arch.h"
void network_init(struct network *const net);
int network_create_neuron(struct neuron *const n, const int log_spikes, const int log_voltages, const int force_update, const int connection_count);
int network_create_neuron_group(struct network *net,  const unsigned int neuron_count, const double threshold, const double reset);
struct neuron *network_id_to_neuron_ptr(struct network *const net, const struct neuron_id id);
int network_map_neuron_group(struct neuron_group *const group, const struct hardware_mapping map);
int net_create_inputs(struct network *const net, const int input_count, const int input_type);
int net_create_input_node(struct input *const in, const int connection_count);
void net_set_input(struct network *const net, const int input_id, const double rate);
void network_free(struct network *const net);

#endif
