// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.h - (spiking) neural network functionality. Spiking neural
//  networks are represented as groups of neurons. A neuron group might have a
//  bunch of neurons all with the same properties (and common hardware).
//  Each neuron has its own state and a set of connections to other neurons.
//  These structures have links to hardware for performance simulation.
//  Here we include different neuron, synapse and dendrite models.
#ifndef NETWORK_HEADER_INCLUDED_
#define NETWORK_HEADER_INCLUDED_

#define NETWORK_MAX_NEURON_GROUPS 1024
#define NETWORK_INVALID_NID -1
#define MAX_FIELDs 128
#define MAX_FIELD_LEN 64

#include <stdint.h>
#include "plugins.hpp"

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
	int group, neuron;
};

struct neuron
{
	struct neuron_group *group;
	struct connection *connections_out;
	struct connection_map *maps_in;
	struct connection_map **maps_out;

	// Mapped hardware
	struct core *core, *post_synaptic_cores;
	struct soma_processor *soma_hw;

	class Base_Soma *soma_class;

	char soma_hw_name[MAX_FIELD_LEN];

	// Track the timestep each hardware unit was last updated
	int id, is_init, fired, connection_out_count;
	int max_connections_out, log_spikes, log_potential, update_needed;
	int force_update, spike_count;
	int soma_last_updated, dendrite_last_updated;
	int maps_in_count, maps_out_count;

	double dendritic_current_decay, processing_latency;
	double current, charge;
	Neuron_Status neuron_status;
	int forced_spikes;

	// LIF specific
	// unsigned int random_range_mask;
	// double potential, current, charge, bias;
	// double reset, reverse_reset, threshold, reverse_threshold;
	// double leak_decay, leak_bias, potential_time_const;
	// double dendritic_current_decay, processing_latency;
	// End of LIF specific
};

struct connection
{
	struct neuron *post_neuron, *pre_neuron;
	struct synapse_processor *synapse_hw;
	// TODO: create a table of hw unit names and just index into this
	//  global table (read only)
	char synapse_hw_name[MAX_FIELD_LEN];
	double weight, current, synaptic_current_decay;
	int id, delay;
};

struct input
{
	struct connection *connections;
	double rate, spike_val; // rate only applies to poisson and rate-based
	int post_connection_count, type, id, send_spike;
};

struct neuron_group
{
	// A neuron group is a collection of neurons that share common
	//  parameters (and possibly hardware). All neurons must be based on the
	//  same models. If implemented in hardware, they also must share common
	//  hardware i.e. the same core and processor blocks in the core
	struct neuron *neurons;
	char default_soma_hw_name[MAX_FIELD_LEN];
	char default_synapse_hw_name[MAX_FIELD_LEN];
	int id, neuron_count;
	int default_log_potential, default_log_spikes;
	int default_max_connections_out, default_force_update;

	int reset_mode, reverse_reset_mode;
	
	double default_threshold, default_reset;
	double default_reverse_threshold, default_reverse_reset;
	double default_leak_decay, default_leak_bias;
};

struct network
{
	struct neuron_group groups[NETWORK_MAX_NEURON_GROUPS];
	struct input *external_inputs;
	int neuron_group_count, external_input_count;
};

struct architecture;
struct attributes;
struct core;


void network_init(struct network *const net);
void network_free(struct network *const net);
int network_create_neuron(struct neuron *const n, struct attributes *attr, const int attribute_count);
int network_create_neuron_group(struct network *net, const int neuron_count, struct attributes *attr, const int attribute_count);
struct neuron *network_id_to_neuron_ptr(struct network *const net, const struct neuron_id id);
int network_create_inputs(struct network *const net, const int input_count, const int input_type);
int network_create_input_node(struct input *const in, const int connection_count);
void network_set_input(struct network *const net, const int input_id, const double rate);
int network_parse_reset_mode(const char *str);
int network_connect_neurons(struct connection *const con, struct neuron *const src, struct neuron *const dest, struct attributes *attr, const int attribute_count);
void network_check_mapped(struct network *const net);

#endif
