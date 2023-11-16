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

int network_create_neuron_group(struct network *net, const int neuron_count,
	struct attributes *attr, const int attribute_count)
{
	struct neuron_group *group;
	int id, ret;

	id = net->neuron_group_count;
	assert(id < NETWORK_MAX_NEURON_GROUPS);
	net->neuron_group_count++;

	group = &(net->groups[id]);
	group->neurons =
		(struct neuron *) malloc(sizeof(struct neuron) * neuron_count);
	if (group->neurons == NULL)
	{
		INFO("Error: Couldn't allocate neuron group %d\n", id);
		exit(1);
	}

	group->neuron_count = neuron_count;

	group->default_max_connections_out = 0;
	group->default_log_potential = 0; // Disabled by default
	group->default_log_spikes = 0;
	group->default_threshold = 1.0;
	group->default_reverse_threshold = -1.0;
	group->default_reset = 0.0;
	group->default_reverse_reset = 0.0;
	group->reset_mode = NEURON_RESET_HARD;
	group->reverse_reset_mode = NEURON_NO_RESET; // Disabled by default
	group->default_force_update = 0;
	// Default is no leak (potential decay), i.e., the potential for the
	//  next timestep is 100% of the previous timestep's
	group->default_leak_decay = 1.0;
	group->default_leak_bias = 0.0;
	for (int i = 0; i < attribute_count; i++)
	{
		struct attributes *a = &(attr[i]);

		ret = -1;
		if (strncmp("threshold", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(
				a->value_str, "%lf", &group->default_threshold);
		}
		else if (strncmp("reverse_threshold", a->key, MAX_FIELD_LEN) ==
			0)
		{
			ret = sscanf(a->value_str, "%lf",
				&group->default_reverse_threshold);
		}
		else if (strncmp("reset", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(
				a->value_str, "%lf", &group->default_reset);
		}
		else if (strncmp("reverse_reset", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%lf",
				&group->default_reverse_reset);
		}
		else if (strncmp("reset_mode", a->key, MAX_FIELD_LEN) == 0)
		{
			group->reset_mode =
				network_parse_reset_mode(a->value_str);
			// Was parsed successfully if we got here
			ret = 1;
		}
		else if (strncmp("reverse_reset_mode", a->key, MAX_FIELD_LEN) ==
			0)
		{
			group->reverse_reset_mode =
				network_parse_reset_mode(a->value_str);
			// Was parsed successfully if we got here
			ret = 1;
		}
		else if (strncmp("leak_decay", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%lf",
				&group->default_leak_decay);
		}
		else if (strncmp("leak_bias", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(
				a->value_str, "%lf", &group->default_leak_bias);
		}
		else if (strncmp("log_v", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%d",
				&group->default_log_potential);
		}
		else if (strncmp("log_spikes", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(
				a->value_str, "%d", &group->default_log_spikes);
		}
		else if (strncmp("connections_out", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%d",
				&group->default_max_connections_out);
		}
		else if (strncmp("force_update", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%d",
				&group->default_force_update);
		}

		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a->key,
				a->value_str);
			exit(1);
		}
	}

	// Initialize all neurons in this group
	for (int i = 0; i < group->neuron_count; i++)
	{
		struct neuron *n = &(group->neurons[i]);

		n->id = i;
		n->group = group;
		n->connection_out_count = 0;

		// Initialize neuron using group attributes
		n->log_spikes = group->default_log_spikes;
		n->log_potential = group->default_log_potential;
		n->force_update = group->default_force_update;
		n->max_connections_out = group->default_max_connections_out;
		// Default connections
		n->reset = group->default_reset;
		n->reverse_reset = group->default_reset;
		n->threshold = group->default_threshold;
		n->reverse_threshold = group->default_reverse_threshold;

		n->fired = 0;
		n->potential = 0.0;
		n->current = 0.0;
		n->bias = 0.0;

		n->leak_decay = group->default_leak_decay;
		n->leak_bias = group->default_leak_bias;

		n->update_needed = 0;
		n->spike_count = 0;
		n->connections_out = NULL;

		// Initially the neuron is not mapped to anything
		n->core = NULL;
		n->axon_in = NULL;
		n->synapse_hw = NULL;
		n->dendrite_hw = NULL;
		n->soma_hw = NULL;
		n->axon_out = NULL;

		n->maps_in = NULL;
		n->maps_out = NULL;
		n->maps_in_count = 0;
		n->maps_out_count = 0;

		n->is_init = 0;
	}

	INFO("Created neuron group gid:%d count:%d "
	     "threshold:%lf neg threshold:%lf "
	     "pos reset:%lf reset mode:%d "
	     "neg reset:%lf neg reset mode:%d\n",
		group->id, group->neuron_count, group->default_threshold,
		group->default_reverse_threshold, group->default_reset,
		group->reset_mode, group->default_reverse_reset,
		group->reverse_reset_mode);

	return id;
}

int network_create_neuron(struct neuron *const n, struct attributes *attr,
	const int attribute_count)
{
	// Each hardware timestep corresponds to a simulation of the spiking
	//  network for dt seconds. This relates to the LIF time constant.
	assert(n != NULL);
	if (n->is_init)
	{
		INFO("Error: Trying to redefine neuron %d.\n", n->id);
		return NETWORK_INVALID_NID;
	}

	n->potential = 0.0;

	/*** Set attributes ***/
	n->bias = 0.0;
	n->random_range_mask = 0;
	n->dendrite_count = 1;
	for (int i = 0; i < attribute_count; i++)
	{
		struct attributes *a = &(attr[i]);
		int ret = -1;

		if (strncmp("bias", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%lf", &(n->bias));
		}
		else if (strncmp("reset", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%lf", &(n->reset));
		}
		else if (strncmp("reverse_reset", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%lf", &(n->reverse_reset));
		}
		else if (strncmp("threshold", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%lf", &(n->threshold));
		}
		else if (strncmp("reverse_threshold", a->key, MAX_FIELD_LEN) ==
			0)
		{
			ret = sscanf(
				a->value_str, "%lf", &(n->reverse_threshold));
		}
		else if (strncmp("connections_out", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(
				a->value_str, "%d", &(n->max_connections_out));
		}
		else if (strncmp("log_spikes", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%d", &(n->log_spikes));
		}
		else if (strncmp("log_v", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%d", &(n->log_potential));
		}
		else if (strncmp("force_update", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%d", &(n->force_update));
		}
		else if (strncmp("dendrites", a->key, MAX_FIELD_LEN) == 0)
		{
			ret = sscanf(a->value_str, "%d", &(n->dendrite_count));
		}
		else
		{
			INFO("Attribute %s not supported.\n", a->key);
			exit(1);
		}

		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a->key,
				a->value_str);
			exit(1);
		}
	}

	// Set the initial update state, no spikes can arrive before the first
	//  time-step but we can force the neuron to update, or bias it
	n->update_needed = (n->force_update || (fabs(n->bias) > 0.0));

	n->soma_last_updated = 0;
	n->dendrite_last_updated = 0;

	assert(n->connections_out == NULL);
	TRACE1("Allocating memory (%lu b) for connections\n",
		sizeof(struct connection) * n->max_connections_out);
	if (n->max_connections_out > 0)
	{
		n->connections_out = (struct connection *) malloc(
			sizeof(struct connection) * n->max_connections_out);
		if (n->connections_out == NULL)
		{
			INFO("Error: Couldn't allocate connection memory.\n");
			exit(1);
		}
	}

	// Zero initialize all connections
	for (int i = 0; i < n->max_connections_out; i++)
	{
		struct connection *con = &(n->connections_out[i]);
		con->id = i;
		con->current = 0.0;
		con->pre_neuron = NULL;
		con->post_neuron = NULL;
		con->weight = 0.0;
		con->delay = 0.0;
		con->synaptic_current_decay = 0.0;
	}

	if (n->dendrite_count < 1)
	{
		INFO("Error: Number of dendrites must be >= 1 (%d)\n",
			n->dendrite_count);
		exit(1);
	}
	INFO("Allocating %d dendritic compartments.\n", n->dendrite_count);
	n->dendrite_potentials = (double *) malloc(
		sizeof(double) * n->dendrite_count);
	n->next_dendrite_potentials = (double *) malloc(
		sizeof(double) * n->dendrite_count);
	n->dendrite_weights = (double *) malloc(
		sizeof(double) * n->dendrite_count * n->dendrite_count);
	if ((n->dendrite_potentials == NULL) ||
		(n->next_dendrite_potentials == NULL) ||
		(n->dendrite_weights == NULL))
	{
		INFO("Error: Couldn't allocate dendrite memory.\n");
		exit(1);
	}
	for (int i = 0; i < n->dendrite_count; i++)
	{
		n->dendrite_potentials[i] = 0.0;
		n->next_dendrite_potentials[i] = 0.0;
		for (int j = 0; j < n->dendrite_count; j++)
		{
			int idx = (i*n->dendrite_count) + j;
			n->dendrite_weights[idx] = 0.0;
		}
	}

	TRACE1("Created neuron: gid:%d nid:%d force:%d thresh:%lf\n",
		n->group->id, n->id, n->force_update, n->threshold);
	n->is_init = 1;
	return n->id;
}

int network_create_dendrite(struct neuron *const n, const int dendrite_id,
	struct attributes *attr, const int attribute_count)
{
	// Each hardware timestep corresponds to a simulation of the spiking
	//  network for dt seconds. This relates to the LIF time constant.
	assert(n != NULL);
	assert(dendrite_id >= 0);
	assert(dendrite_id < n->dendrite_count);

	/*** Set attributes ***/
	n->dendrite_potentials[dendrite_id] = 0.0;

	for (int i = 0; i < attribute_count; i++)
	{
		struct attributes *a = &(attr[i]);
		int ret = -1;

		if (strncmp("leak", a->key, MAX_FIELD_LEN) == 0)
		{
			// Index into the adjacency matrix of weights for the
			//  dendrite of this neuron. The leak will set the
			//  recurrent edge connected this tap to itself
			int idx = (dendrite_id*n->dendrite_count) + dendrite_id;
			ret = sscanf(a->value_str, "%lf",
				&(n->dendrite_weights[idx]));
		}
		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a->key,
				a->value_str);
			exit(1);
		}
	}
	return RET_OK;
}

int network_connect_neurons(struct connection *const con,
	struct neuron *const src, struct neuron *const dest,
	const int dendrite_id, struct attributes *attr,
	const int attribute_count)
{
	assert(con != NULL);
	con->pre_neuron = src;
	con->post_neuron = dest;
	con->weight = 1.0;
	con->post_dendrite_id = dendrite_id;

	total_connection_count++;
	for (int i = 0; i < attribute_count; i++)
	{
		struct attributes *a = &(attr[i]);
		int ret = -1;

		if ((strncmp("w", a->key, MAX_FIELD_LEN) == 0) ||
			(strncmp("weight", a->key, MAX_FIELD_LEN) == 0))
		{
			ret = sscanf(a->value_str, "%lf", &con->weight);
		}
		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a->key,
				a->value_str);
			exit(1);
		}
	}

	TRACE1("\tAdded con %d.%d->%d.%d (w:%lf)\n", con->pre_neuron->group->id,
		con->pre_neuron->id, con->post_neuron->group->id,
		con->post_neuron->id, con->weight);
	return RET_OK;
}

int network_connect_dendrites(struct neuron *const n, const int src_dendrite_id,
	const int dest_dendrite_id, struct attributes *attr,
	const int attribute_count)
{
	// Index into adjacency matrix of weights
	assert(n != NULL);
	int idx = (src_dendrite_id*n->dendrite_count) + dest_dendrite_id;

	INFO("src:%d dest:%d idx:%d\n", src_dendrite_id, dest_dendrite_id, idx);
	assert(src_dendrite_id != dest_dendrite_id);
	assert(src_dendrite_id < n->dendrite_count);
	assert(dest_dendrite_id < n->dendrite_count);

	for (int i = 0; i < attribute_count; i++)
	{
		struct attributes *a = &(attr[i]);
		int ret = -1;

		if ((strncmp("w", a->key, MAX_FIELD_LEN) == 0) ||
			(strncmp("weight", a->key, MAX_FIELD_LEN) == 0))
		{
			int src_identity = (src_dendrite_id*n->dendrite_count) +
				src_dendrite_id;
			double weight;
			ret = sscanf(a->value_str, "%lf", &weight);
			// Update the src compartment, subtracting some of its
			//  potential which is added to the destination
			//  compartment. Formulate this for the matrix mul
			n->dendrite_weights[idx] = weight;
			n->dendrite_weights[src_identity] -= weight;
		}
		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a->key,
				a->value_str);
			exit(1);
		}
	}
	return RET_OK;
}


/*
int network_create_inputs(struct network *const net, const int input_count,
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

int network_create_input_node(struct input *const in, const int connection_count)
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
		con->current = 0.0;
		con->synaptic_current_decay = 0.0;
		con->delay = 0.0;
	}

	return in->id;
}
*/

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
			free(group->neurons[j].connections_out);
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

struct neuron *network_id_to_neuron_ptr(
	struct network *const net, const struct neuron_id id)
{
	struct neuron_group *group;
	struct neuron *neuron;

	if (id.group < net->neuron_group_count)
	{
		INFO("ERROR: Group %d > max %d.\n", id.group,
			net->neuron_group_count);
		exit(1);
	}

	group = &(net->groups[id.group]);
	if (id.neuron < group->neuron_count)
	{
		INFO("ERROR: Neuron %d > max %d.\n", id.neuron,
			group->neuron_count);
		exit(1);
	}
	neuron = &(group->neurons[id.neuron]);

	return neuron;
}

void network_set_input(
	struct network *const net, const int input_id, const double rate)
{
	struct input *in = &(net->external_inputs[input_id]);
	in->rate = rate;
}

void network_check_mapped(struct network *const net)
{
	// Check that all network neurons are mapped to a physical core
	//  If a neuron is not, print an error message and stop the simulation
	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);

			if (n->core == NULL)
			{
				INFO("Error: Neuron %d.%d not mapped to H/W.\n",
					group->id, n->id);
				exit(1);
			}
		}
	}
}

int network_map_hardware(struct neuron *n, struct core *c)
{
	// Map the neuron to hardware units
	assert(n != NULL);
	assert(c != NULL);
	assert(c->neurons != NULL);
	assert(n->core == NULL);

	n->core = c;
	TRACE1("Mapping neuron %d to core %d\n", n->id, c->id);
	c->neurons[c->neuron_count] = n;
	c->neuron_count++;

	n->axon_in = &c->axon_in;
	n->synapse_hw = &c->synapse;
	n->dendrite_hw = &c->dendrite;
	n->soma_hw = &c->soma;
	n->axon_out = &c->axon_out;

	return RET_OK;
}

int network_parse_reset_mode(const char *str)
{
	int reset_mode = -1;

	if (strcmp(str, "none") == 0)
	{
		reset_mode = NEURON_NO_RESET;
	}
	else if (strcmp(str, "soft") == 0)
	{
		reset_mode = NEURON_RESET_SOFT;
	}
	else if (strcmp(str, "hard") == 0)
	{
		reset_mode = NEURON_RESET_HARD;
	}
	else if (strcmp(str, "saturate") == 0)
	{
		reset_mode = NEURON_RESET_SATURATE;
	}
	else
	{
		INFO("Error: reset mode not recognized.");
		exit(1);
	}

	return reset_mode;
}
