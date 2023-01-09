// arch.c
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#include "print.h"
#include "arch.h"

void arch_init(struct architecture *const arch)
{
	arch->tile_count = 0;
	arch->time_barrier = 0.0;
}

void arch_free(struct architecture *const arch)
{
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);

			free(c->neurons);
			free(c->spikes_sent_per_core);
			free(c->axon_map);
			c->neurons = NULL;
			c->spikes_sent_per_core = NULL;
			c->axon_map = NULL;
		}
	}
}

int arch_create_noc(struct architecture *const arch, const int width,
							const int height)
{
	int tile_id = 0;
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			struct tile *t = &(arch->tiles[tile_id]);
			int north_x, north_y, east_x, east_y, south_x, south_y;
			int west_x, west_y, link_count;

			tile_id++;
			t->x = x;
			t->y = y;

			north_x = t->x;
			north_y = t->y - 1;
			east_x = t->x + 1;
			east_y = t->y;
			south_x = t->x;
			south_y = t->y + 1;
			west_x = t->x - 1;
			west_y = t->y;

			link_count = 0;
			TRACE("tid:%d (x:%d,y:%d)\n", t->id, t->x, t->y);
			if (north_y >= 0)
			{
				int lid = (north_y * width) + north_x;
				t->links[link_count] = &(arch->tiles[lid]);
				link_count++;
			}
			if (east_x < width)
			{
				int lid = (east_y * width) + east_x;
				t->links[link_count] = &(arch->tiles[lid]);
				link_count++;
			}
			if (south_y < height)
			{
				int lid = (south_y * width) + south_x;
				t->links[link_count] = &(arch->tiles[lid]);
				link_count++;
			}
			if (west_x >= 0)
			{
				int lid = (west_y * width) + west_x;
				t->links[link_count] = &(arch->tiles[lid]);
				link_count++;
			}
			assert(link_count > 0);
			assert(link_count <= 4);
			for (int i = 0; i < link_count; i++)
			{
				assert(t->links[i] != NULL);
				TRACE("\tlink[%d]->%d\n", i, (t->links[i])->id);
			}
		}
	}

	INFO("NoC created, mesh, width:%d height:%d.\n", width, height);
	return 0;
}

int arch_create_tile(struct architecture *const arch,
			const double energy_east_west_hop,
			const double energy_north_south_hop,
			const double time_east_west_hop,
			const double time_north_south_hop)
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

	// TODO: generalize this so each link has an associated energy and
	//  latency, set one level up when defining the noc mesh
	t->energy_east_west_hop = energy_east_west_hop;
	t->time_east_west_hop = time_east_west_hop;
	t->energy_north_south_hop = energy_north_south_hop;
	t->time_north_south_hop = time_north_south_hop;
	t->x = 0;
	t->y = 0;
	t->core_count = 0;

	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		struct core *c = &(t->cores[i]);

		c->id = i;
		c->t = t;
	}

	TRACE("Tile created id:%d.\n", t->id);
	return t->id;
}

int arch_create_core(struct architecture *const arch, struct tile *const t)
{
	struct core *c;
	unsigned int core_id;

	assert(t != NULL);
	core_id = t->core_count;
	t->core_count++;

	c = &(t->cores[core_id]);
	c->id = core_id;
	c->t = t;

	for (int i = 0; i < ARCH_MAX_PROCESSORS; i++)
	{
		c->axon_in[i].id = i;
		c->axon_in[i].energy = 0.0;
		c->axon_in[i].time = 0.0;
		c->axon_in[i].t = t;

		c->synapse[i].id = i;
		c->synapse[i].energy = 0.0;
		c->synapse[i].time = 0.0;

		c->dendrite[i].id = i;
		c->dendrite[i].energy = 0.0;
		c->dendrite[i].time = 0.0;

		c->soma[i].id = i;
		c->soma[i].energy = 0.0;
		c->soma[i].time = 0.0;

		c->axon_out[i].id = i;
		c->axon_out[i].energy = 0.0;
		c->axon_out[i].time = 0.0;
		c->axon_out[i].t = t;
	}

	c->neuron_count = 0;
	c->curr_neuron = 0;
	c->neurons = (struct neuron **) malloc(sizeof(struct neuron *) * 1024);
	c->spikes_sent_per_core = (int *) malloc(sizeof(int) * 128);
	c->axon_map = (struct core **) malloc(sizeof(struct core *) * 128);
	//c->message_processing_time =
	//			(struct double **) malloc(sizeof(double) * 128);
	for (int i = 0; i < 1024; i++)
	{
		c->neurons[i] = NULL;
	}

	c->energy = 0.0;
	c->time = 0.0;

	c->axon_in_count = 0;
	c->synapse_count = 0;
	c->dendrite_count = 0;
	c->soma_count = 0;
	c->axon_out_count = 0;

	TRACE("Core created id:%d (tile:%d).\n", c->id, t->id);
	return c->id;
}

int arch_create_axon_in(struct architecture *const arch, struct core *const c)
{
	struct axon_input *in;
	int count = c->axon_in_count;

	if (c->axon_in_count >= ARCH_MAX_PROCESSORS)
	{
		INFO("Error: Max %d axon inputs\n", ARCH_MAX_PROCESSORS);
		return ARCH_INVALID_ID;
	}

	in = &(c->axon_in[count]);
	in->energy = 0.0;
	in->time = 0.0;
	in->packet_size = 0;
	// We already know a valid tile was given at this point
	in->t = c->t;

	c->axon_in_count++;

	TRACE("Created axon input %d (c:%d.%d)\n",
						in->id, c->t->id, c->id);

	return in->id;
}

int arch_create_synapse(struct architecture *const arch, struct core *const c,
					const double energy_spike_op,
					const double time_spike_op)
{
	struct synapse_processor *s;
	int count = c->synapse_count;

	if (c->synapse_count >= ARCH_MAX_PROCESSORS)
	{
		INFO("Error: Max %d synapse processors\n", ARCH_MAX_PROCESSORS);
		return ARCH_INVALID_ID;
	}

	s = &(c->synapse[count]);
	s->energy = 0.0;
	s->time = 0.0;

	s->energy_spike_op = energy_spike_op;
	s->time_spike_op = time_spike_op;

	c->synapse_count++;
	s->buffer_in = 0;
	TRACE("Created synapse processor %d (c:%d.%d)\n",
						s->id, c->t->id, c->id);

	return s->id;
}

int arch_create_soma(struct architecture *const arch, struct core *const c,
			int model,
			double energy_active_neuron_update,
			double time_active_neuron_update,
			double energy_inactive_neuron_update,
			double time_inactive_neuron_update,
			double energy_spiking,
			double time_spiking)
{
	struct soma_processor *s;
	int count = c->soma_count;

	if (c->synapse_count >= ARCH_MAX_PROCESSORS)
	{
		INFO("Error: Max %d synapse processors\n", ARCH_MAX_PROCESSORS);
		return ARCH_INVALID_ID;
	}

	s = &(c->soma[count]);
	s->energy = 0.0;
	s->time = 0.0;

	s->model = model;
	s->energy_active_neuron_update = energy_active_neuron_update;
	s->time_active_neuron_update = time_active_neuron_update;
	s->energy_inactive_neuron_update = energy_inactive_neuron_update;
	s->time_inactive_neuron_update = time_inactive_neuron_update;
	s->energy_spiking = energy_spiking;
	s->time_spiking = time_spiking;

	s->buffer_in = 1;

	c->soma_count++;
	TRACE("Created soma processor %d (c:%d.%d)\n",
						s->id, c->t->id, c->id);

	return s->id;
}

int arch_create_axon_out(struct architecture *const arch, struct core *const c,
			const double spike_energy, const double spike_time)
{
	struct axon_output *out;
	int count = c->axon_out_count;

	if (count >= ARCH_MAX_PROCESSORS)
	{
		INFO("Error: Max %d axon outputs\n", ARCH_MAX_PROCESSORS);
		return ARCH_INVALID_ID;
	}

	out = &(c->axon_out[count]);
	out->packets_out = 0;
	out->energy = 0.0;
	out->time = 0.0;
	out->energy_spike_within_tile = spike_energy;
	out->time_spike_within_tile = spike_time;
	out->buffer_in = 0;

	// Track the tile the axon interfaces with
	out->t = c->t;

	TRACE("Created axon output %d (c:%d.%d)\n", out->id, c->t->id, c->id);

	return out->id;
}

/*
static struct core *arch_get_core(struct architecture *const arch,
					const int tile_id, const int core_id)
{
	struct tile *t;
	struct core *c;

	// Make some extra sanity checks on the user input, there's no
	//  guarantees the ids point to a valid tile or core
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
*/
