// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.c
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "arch.h"
#include "print.h"
#include "network.h"

int total_connection_count = 0;

int network_create_neuron_group(struct network *net,
				const unsigned int neuron_count,
				const double threshold, const double reset)
{
	struct neuron_group *group;
	int id;

	id = net->neuron_group_count;
	assert(id < NETWORK_MAX_NEURON_GROUPS);
	net->neuron_group_count++;

	group = &(net->groups[id]);
	group->neurons = (struct neuron *)
				malloc(sizeof(struct neuron) * neuron_count);
	if (group->neurons == NULL)
	{
		INFO("Error: Couldn't allocate neuron group %d\n", id);
		exit(1);
	}

	group->neuron_count = neuron_count;
	group->default_threshold = threshold;
	group->default_reset = reset;

	for (int i = 0; i < group->neuron_count; i++)
	{
		struct neuron *n = &(group->neurons[i]);

		n->id = i;
		n->group = group;
		n->post_connection_count = 0;
		n->log_spikes = 0;
		n->log_voltage = 0;
		n->force_update = 0;

		n->fired = 0;
		n->potential = 0.0;
		n->current = 0.0;
		n->charge = 0.0;
		n->bias = 0.0;
		n->threshold = group->default_threshold;
		n->reset = group->default_reset;
		n->update_needed = 0;
		n->spike_count = 0;
		n->connections = NULL;

		n->charge_buffer = 0.0;
		n->d_currents_buffer[0] = 0.0;
		n->current_buffer = 0.0;
		n->fired_buffer = 0;

		// Initially the neuron is not mapped to anything
		n->core = NULL;
		n->axon_in = NULL;
		n->synapse_hw = NULL;
		n->dendrite_hw = NULL;
		n->soma_hw = NULL;
		n->axon_out = NULL;

		n->is_init = 0;
	}

	INFO("Added neuron group gid:%d count:%d threshold:%lf reset:%lf\n",
		group->id, group->neuron_count, group->default_threshold,
		group->default_reset);

	return id;
}

int network_create_neuron(struct neuron *const n,
					const double bias,
					const int log_spikes,
					const int log_voltages,
					const int force_update,
					const int connection_count)
{
	// Each hardware timestep corresponds to a simulation of the spiking
	//  network for dt seconds. This relates to the LIF time constant.
	const double dt = 1.0e-3; // Seconds
	assert(connection_count >= 0);
	assert(n != NULL);

	if (n->is_init)
	{
		INFO("Error: Trying to redefine neuron twice.\n");
		return NETWORK_INVALID_NID;
	}

	n->bias = bias;
	n->log_spikes = log_spikes;
	n->log_voltage = log_voltages;
	n->force_update = force_update;
	n->update_needed = n->force_update;
	n->post_connection_count = connection_count;
	total_connection_count += connection_count;
	// TODO: Hard coded LIF / CUBA time constants for now
	n->current_time_const = 10.0e-3;
	n->potential_time_const = 10.0e-3;
	n->current_decay =
		exp(-dt / n->current_time_const);
	n->potential_decay =
		exp(-dt / n->potential_time_const);
	n->current_decay = 1.0;
	n->potential_decay = 1.0;

	assert(n->connections == NULL);
	n->connections = (struct connection *)
			malloc(sizeof(struct connection) * connection_count);
	if (n->connections == NULL)
	{
		INFO("Error: Couldn't allocate connection memory.\n");
		return NETWORK_INVALID_NID;
	}

	n->soma_last_updated = 0;

	// Zero initialize all connections
	for (int i = 0; i < n->post_connection_count; i++)
	{
		struct connection *con = &(n->connections[i]);
		con->pre_neuron = NULL;
		con->post_neuron = NULL;
		con->weight = 0.0;
	}

	INFO("Created neuron: gid:%d nid:%d force:%d thresh:%lf con:%d\n",
					n->group->id, n->id, n->force_update,
					n->threshold, n->post_connection_count);
	n->is_init = 1;
	return n->id;
}

int network_map_neuron(struct neuron *const n,
					const struct hardware_mapping map)
{
	// Map the neuron to hardware units
	assert(map.core != NULL);
	assert(map.core->neurons != NULL);
	assert(n->core == NULL);
	n->core = map.core;
	INFO("mapping core %d to neuron:%d\n",
		map.core->id, map.core->neuron_count);
	map.core->neurons[map.core->neuron_count] = n;
	map.core->neuron_count++;

	n->axon_in = map.axon_in;
	n->synapse_hw = map.synapse_hw;
	n->dendrite_hw = map.dendrite_hw;
	n->soma_hw = map.soma_hw;
	n->axon_out = map.axon_out;

	return 1;
}

int net_create_inputs(struct network *const net, const int input_count,
							const int input_type)
{
	assert(input_count > 0);
	assert(net->external_inputs == NULL);

	net->external_inputs =
		(struct input *) malloc(sizeof(struct input) * input_count);
	if (net->external_inputs == NULL)
	{
		INFO("Error: Couldn't allocate memory for network inputs.\n");
		return -1;
	}
	net->external_input_count = input_count;

	for (int i = 0; i < input_count; i++)
	{
		struct input *in = &(net->external_inputs[i]);

		in->id = i;
		in->connections = NULL;
		in->spike_val = 0.0;
		in->post_connection_count = 0;
		in->type = input_type;
	}

	return 0;
}

int net_create_input_node(struct input *const in, const int connection_count)
{
	TRACE("Creating input node with %d connections.\n", connection_count);

	assert(in->connections == NULL);
	in->connections = (struct connection *)
		malloc(sizeof(struct connection) * connection_count);
	if (in->connections == NULL)
	{
		INFO("Error: Couldn't allocate connection memory.\n");
		return NETWORK_INVALID_NID;
	}

	// Zero initialize all connections
	in->post_connection_count = connection_count;
	for (int i = 0; i < in->post_connection_count; i++)
	{
		struct connection *con = &(in->connections[i]);
		con->pre_neuron = NULL;
		con->post_neuron = NULL;
		con->weight = 0.0;
	}

	return in->id;
}

void network_init(struct network *const net)
{
	net->neuron_group_count = 0;
	net->external_input_count = 0;
	net->external_inputs = NULL;

	for (int i = 0; i < NETWORK_MAX_NEURON_GROUPS; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		group->id = i;
		group->neuron_count = 0;
	}
}

void network_free(struct network *const net)
{
	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		// First free all the allocated connections in each neuron
		for (int j = 0; j < group->neuron_count; j++)
		{
			free(group->neurons[j].connections);
		}
		// Finally free the neurons allocated in the group
		free(net->groups[i].neurons);
	}

	for (int i = 0; i < net->external_input_count; i++)
	{
		struct input *in = &(net->external_inputs[i]);
		free(in->connections);
	}
	free(net->external_inputs);

	return;
}

struct neuron *network_id_to_neuron_ptr(struct network *const net,
					const struct neuron_id id)
{
	struct neuron_group *group;
	struct neuron *neuron;

	if (id.group < net->neuron_group_count)
	{
		INFO("ERROR: Group %d > max %d.\n",
					id.group, net->neuron_group_count);
		exit(1);
	}

	group = &(net->groups[id.group]);
	if (id.neuron < group->neuron_count)
	{
		INFO("ERROR: Neuron %d > max %d.\n",
						id.neuron, group->neuron_count);
		exit(1);
	}
	neuron = &(group->neurons[id.neuron]);

	return neuron;
}

void net_set_input(struct network *const net, const int input_id,
							const double rate)
{
	struct input *in = &(net->external_inputs[input_id]);
	in->rate = rate;
}
