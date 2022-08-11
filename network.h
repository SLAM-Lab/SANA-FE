#ifndef NETWORK_HEADER_INCLUDED
#define NETWORK_HEADER_INCLUDED

#define NETWORK_MAX_NEURON_GROUPS 1024

// This defines the structure of the csv format for specifying neurons in a
//  network for this simulator.   Each row in the file represents data for a
//  unique neuron in hardware.
//
// The format is a number of neuron fields, followed by a variable number of
//  synapses (each with SYNAPSE_FIELDS entries).  This means each line of the
//  csv file is variable length depending on the number of synapses for that
//  neuron.
enum neuron_config_format
{
	NEURON_ID = 0,
	COMPARTMENT_ID, // Which hardware neuron to map to
	THRESHOLD_VOLTAGE, // Spiking threshold value
	RESET_VOLTAGE,
	RECORD_SPIKES, // Record any spikes to a CSV file
	RECORD_VOLTAGE, // Write the voltage to a CSV file every timestep
	FORCE_UPDATE,
	NEURON_FIELDS,
};

enum connection_config_format
{
	// This is the structure of the CSV format for specifying synaptic
	//  connections in a network for this simulator.  Each row represents
	//  a unique connection.
	CONNECTION_DEST_NID = 0,
	CONNECTION_WEIGHT,
	CONNECTION_FIELDS,
};

struct neuron_id
{
	unsigned int group, neuron;
};

struct neuron
{
	int id, neuron_used, fired, post_connection_count, spike_count;
	int log_spikes, log_voltage, update_needed, force_update;
	double potential, current,  bias, reset, threshold;
	double potential_decay, potential_time_const;
	double current_decay, current_time_const;

	struct connection *connections;
	struct neuron_group *group;
};

struct connection
{
	int id;
	struct neuron *post_neuron, *pre_neuron;
	double weight;
	// TODO: delay, and any other synapse/connection properties
};

struct input
{
	int send_spike, post_connection_count;
	struct connection *connections;
};

struct neuron_group
{
	// A neuron group is a collection of neurons that share common
	//  parameters (and possibly hardware). All neurons must be based on the
	//  same models. If implemented in hardware, they also must share common
	//  hardware i.e. the same core and processor blocks in the core
	struct neuron *neurons;
	int id, neuron_count;
	double default_threshold, default_reset;

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
void network_create_empty(struct architecture *arch);
//void network_read_csv(FILE *fp, struct network *net, struct architecture *arch);
void network_init(struct network *net);

#endif
