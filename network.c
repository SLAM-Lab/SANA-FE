// network.c
// Utility functions for creating user defined spiking networks
//#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "network.h"
#include "sim.h"

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

	// Initially the group of neurons is not mapped to anything
	group->core = NULL;
	group->axon_in = NULL;
	group->synapse = NULL;
	group->dendrite = NULL;
	group->soma = NULL;
	group->axon_out = NULL;

	return id;
}

void network_add_connections(struct neuron *const src,
				const struct connection connections[],
				const int connection_count)
{	
	// TODO: we need to pass neuron connections all in one go
	//  so that we can allocate the right amount of space 
	// So need to get the connections just for this 
	src->connections = (struct connection *)
			malloc(sizeof(struct connection) * connection_count);
	if (src->connections == NULL)
	{
		INFO("Error: Couldn't allocate connections for "
			"neuron (%d.%d).\n", src->group->id, src->id);
		exit(1);
	}

	// Copy connections across
	for (int i = 0; i < connection_count; i++)
	{
		src->connections[i] = connections[i];
	}

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

/*
// Mapping should happen outside of this file, because it's a combination
//  of both architecture and network
void network_map_neuron_group(struct network *net, struct architecture *arch,
				const int group_id, const int core_id, const int axon_in, const int synapse, const int  )
{
	assert(id < net->neuron_group_count);

	group = &(net->groups[group_id]);
	
	// First find the core
	for (int i = 0; i < arch->tile_count)

	group->core = arch->tiles
}
*/

void network_init(struct network *net)
{
	net->neuron_group_count = 0;
	net->external_input_count = 0;

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);

			n->id = i;
			n->post_connection_count = 0;
			n->log_spikes = 0;
			n->log_voltage = 0;
			n->spike_count = 0;

			n->fired = 0;
			n->potential = 0.0;
			n->current = 0.0;
			n->bias = 0.0;
			n->threshold = group->default_threshold;
			n->reset = group->default_reset;
			n->update_needed = 0;
			n->neuron_used = 0;
		}
	}
}

void network_free(struct network *net)
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

	return;
}
