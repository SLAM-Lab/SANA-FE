// arch.c: Create a neuromorphic design based on an architecture description
//  In this simulator and architecture description is a list of different
//  blocks.
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#include "arch.h"
#include "sim.h"

static struct core *arch_get_core(struct architecture *const arch, const int tile_id, const int core_id);

void arch_init(struct architecture *const arch)
{
	arch->tile_count = 0;
}

int arch_create_axon_in(struct architecture *const arch,
						const int tile_id,
						const int core_id)
{
	struct core *c;
	struct axon_input *in;
	int count;

	c = arch_get_core(arch, tile_id, core_id);
	count = c->axon_in_count;
	in = &(c->axon_in[count]);

	in->id = count;
	in->energy = 0.0;
	in->time = 0.0;
	in->packet_size = 0;
	// We already know a valid tile was given at this point
	in->t = &(arch->tiles[tile_id]);

	c->axon_in_count++;

	return in->id;
}

int arch_create_axon_out(struct architecture *const arch,
						const int tile_id,
						const int core_id)
{
	struct core *c;
	struct axon_output *out;
	int count;

	c = arch_get_core(arch, tile_id, core_id);
	count = c->axon_out_count;
	out = &(c->axon_out[count]);

	out->id = count;
	out->total_packets_sent = 0;
	out->energy = 0.0;
	out->time = 0.0;

	TRACE("Created axon output %d (c:%d t:%d)\n", out->id, core_id,
								tile_id);
	c->axon_out_count++;

	return out->id;
}

int arch_create_core(struct architecture *const arch,
						const int tile_id)
{
	struct tile *t;
	struct core *c;
	unsigned int core_id;
	
	if (tile_id >= arch->tile_count)
	{
		INFO("Error: Invalid tile id: %d.\n", tile_id);
		exit(1);
	}
	t = &(arch->tiles[tile_id]);
	core_id = t->core_count;
	t->core_count++;

	c = &(t->cores[core_id]);
	c->id = core_id;

	for (int i = 0; i < ARCH_MAX_PROCESSORS; i++)
	{
		c->axon_in[i].id = i;
		c->axon_in[i].t = t;
		c->synapse[i].id = i;
		c->dendrite[i].id = i;
		c->soma[i].id = i;
		c->axon_out[i].id = i;
	}
	c->energy = 0.0;
	c->time = 0.0;

	c->axon_in_count = 0;
	c->synapse_count = 0;
	c->dendrite_count = 0;
	c->soma_count = 0;
	c->axon_out_count = 0;

	TRACE("Core created id:%d (tile:%u).\n", c->id, t->id);
	return c->id;
}

int arch_create_tile(struct architecture *const arch)
{
	struct tile *t;
	int id;
	
	if (arch->tile_count >= ARCH_MAX_TILES)
	{
		INFO("Error: Only %d tiles supported.\n", ARCH_MAX_TILES);
		exit(1);
	}

	id = arch->tile_count;
	arch->tile_count++;
	t = &(arch->tiles[id]);

	t->id = id;
	t->energy = 0.0;
	t->time = 0.0;
	t->energy_east_west_hop = 0.0;
	t->time_east_west_hop = 0.0;
	t->energy_north_south_hop = 0.0;
	t->time_north_south_hop = 0.0;
	t->x = 0;
	t->y = 0;
	t->core_count = 0;

	TRACE("Tile created id:%u.\n", t->id);
	return t->id;
}

static struct core *arch_get_core(struct architecture *const arch,
					const int tile_id, const int core_id)
{
	struct tile *t;
	struct core *c;

	if (tile_id >= arch->tile_count)
	{
		INFO("Error: Accessing invalid tile: %d (max:%d)\n",
						tile_id, arch->tile_count);
		exit(1);
	}
	t = &(arch->tiles[tile_id]);

	if (core_id >= t->core_count)
	{
		INFO("Error: Accessing invalid core: %d (max:%d)\n",
						core_id, t->core_count);
		exit(1);
	}
	c = &(t->cores[core_id]);

	return c;
}

/*
void arch_free(struct architecture *arch)
{
	// Free the architecture structure, deallocate all memory and mark
	//  the structure as uninitialized
	arch->initialized = 0;

	// Free any neuron data
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);

		assert(n != NULL);
		free(n->synapses);
		n->synapses = NULL;
	}

	for (int i = 0; i < arch->max_axon_outputs; i++)
	{
		struct axon_output *axon_out = &(arch->axon_outputs[i]);

		assert(axon_out != NULL);
		free(axon_out->packets_sent);
		axon_out->packets_sent = NULL;
	}

	for (int i = 0; i < arch->max_external_inputs; i++)
	{
		struct input *in = &(arch->external_inputs[i]);

		assert(in != NULL);
		free(in->synapses);
		in->synapses = NULL;
	}

	// Free all memory
	free(arch->neurons);
	free(arch->mem_blocks);
	free(arch->routers);
	free(arch->axon_inputs);
	free(arch->axon_outputs);
	free(arch->external_inputs);
	free(arch->timers);

	// Reset all pointers and counters
	arch->neurons = NULL;
	arch->mem_blocks = NULL;
	arch->routers = NULL;
	arch->axon_inputs = NULL;
	arch->axon_outputs = NULL;
	arch->external_inputs = NULL;
	arch->timers = NULL;

	arch->max_neurons = 0;
	arch->max_mem_blocks = 0;
	arch->max_routers = 0;
	arch->max_timers = 0;
	arch->max_axon_inputs = 0;
	arch->max_axon_outputs = 0;
	arch->max_external_inputs = 0;
}
*/

// TODO: delete, neurons aren't modelled in arch
/*
static void arch_parse_neuron(struct architecture *arch,
					char fields[][ARCH_MAX_FIELD_LEN])
{
	struct range field_list[ARCH_MAX_VALUES];
	struct axon_input *axon_in;
	struct axon_output *axon_out;
	double *timer;
	int ret, list_len, mem_id, axon_in_id, axon_out_id, timer_id;
	char block_type;

	assert(fields && fields[0]);
	block_type = fields[0][0];
	assert(block_type == 'n');

	ret = sscanf(fields[2], "%d", &mem_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't mem block (%s).\n", fields[2]);
		exit(1);
	}
	assert(mem_id < arch->max_mem_blocks);

	ret = sscanf(fields[3], "%d", &axon_in_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon in (%s).\n", fields[3]);
		exit(1);
	}
	assert(axon_in_id < arch->max_axon_inputs);
	axon_in = &(arch->axon_inputs[axon_in_id]);

	ret = sscanf(fields[4], "%d", &axon_out_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon out (%s).\n", fields[4]);
		exit(1);
	}
	assert(axon_out_id < arch->max_axon_outputs);
	axon_out = &(arch->axon_outputs[axon_out_id]);

	ret = sscanf(fields[5], "%d", &timer_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse timer (%s).\n", fields[5]);
		exit(1);
	}
	assert(timer_id < arch->max_timers);
	timer = &(arch->timers[timer_id]);
	assert(timer != NULL);

	// Copy the neuron a number of times
	list_len = arch_parse_list(fields[1], &(field_list[0]));
	for (int i = 0; i < list_len; i++)
	{
		struct range *r = &(field_list[i]);

		TRACE("min:%u max%u\n", r->min, r->max);
		for (int id = r->min; id <= r->max; id++)
		{
			struct neuron *n;

			assert(id <= arch->max_neurons);
			n = &(arch->neurons[id]);
			// Copy any common neuron variables
			//n->mem_block = mem_block;
			n->axon_in = axon_in;
			n->axon_out = axon_out;
			n->time = timer;
		}
		TRACE("neuron parsed n:%d-%d mem(%d) i(%d) o(%d).\n",
			r->min, r->max, mem_id, axon_in_id, axon_out_id);
	}
}
*/
