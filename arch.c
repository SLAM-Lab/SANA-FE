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
	t->is_blocking = 0;

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
	c->is_blocking = 1;

	c->axon_in.energy = 0.0;
	c->axon_in.time = 0.0;
	c->axon_in.t = t;

	c->synapse.energy = 0.0;
	c->synapse.time = 0.0;

	c->dendrite.energy = 0.0;
	c->dendrite.time = 0.0;

	c->soma.energy = 0.0;
	c->soma.time = 0.0;

	c->axon_out.energy = 0.0;
	c->axon_out.time = 0.0;
	c->axon_out.t = t;

	c->neuron_count = 0;
	c->curr_neuron = 0;
	c->neurons = (struct neuron **) malloc(sizeof(struct neuron *) * 1024);
	c->spikes_sent_per_core = (int *) malloc(sizeof(int) * 128);
	c->axon_map = (struct core **) malloc(sizeof(struct core *) * 128);
	for (int i = 0; i < 1024; i++)
	{
		c->neurons[i] = NULL;
	}

	c->energy = 0.0;
	c->time = 0.0;

	c->buffer_pos = BUFFER_SOMA;

	TRACE("Core created id:%d (tile:%d).\n", c->id, t->id);
	return c->id;
}

void arch_create_axon_in(struct architecture *const arch, struct core *const c)
{
	struct axon_input *in;

	in = &(c->axon_in);
	in->energy = 0.0;
	in->time = 0.0;
	in->packet_size = 0;
	// We already know a valid tile was given at this point
	in->t = c->t;

	TRACE("Created axon input %d (c:%d.%d)\n",
						in->id, c->t->id, c->id);

	return;
}

void arch_create_synapse(struct architecture *const arch, struct core *const c,
						const int weight_bits,
						const int word_bits,
						const double energy_spike_op,
						const double time_spike_op,
						const double energy_memory_read,
						const double time_memory_read)
{
	struct synapse_processor *s;

	s = &(c->synapse);
	s->energy = 0.0;
	s->time = 0.0;

	s->energy_spike_op = energy_spike_op;
	s->time_spike_op = time_spike_op;

	// The word size is the number of bits accessed with each memory read.
	//  The weight size is the number of bits for a single synaptic weight.
	//  A single memory read might return multiple weights
	s->weight_bits = weight_bits;
	s->word_bits = word_bits;
	// Round up to the nearest word
	s->weights_per_word = (word_bits + (weight_bits - 1)) / weight_bits;

	TRACE("Created synapse processor %d (c:%d.%d)\n",
						s->id, c->t->id, c->id);

	return;
}

void arch_create_soma(struct architecture *const arch, struct core *const c,
			int model,
			double energy_active_neuron_update,
			double time_active_neuron_update,
			double energy_inactive_neuron_update,
			double time_inactive_neuron_update,
			double energy_spiking,
			double time_spiking)
{
	struct soma_processor *s;

	s = &(c->soma);
	s->energy = 0.0;
	s->time = 0.0;

	s->model = model;
	s->energy_active_neuron_update = energy_active_neuron_update;
	s->time_active_neuron_update = time_active_neuron_update;
	s->energy_inactive_neuron_update = energy_inactive_neuron_update;
	s->time_inactive_neuron_update = time_inactive_neuron_update;
	s->energy_spiking = energy_spiking;
	s->time_spiking = time_spiking;

	TRACE("Created soma processor %d (c:%d.%d)\n",
						s->id, c->t->id, c->id);

	return;
}

void arch_create_axon_out(struct architecture *const arch, struct core *const c,
			const double access_energy, const double access_time)
{
	struct axon_output *out;

	out = &(c->axon_out);
	out->packets_out = 0;
	out->energy = 0.0;
	out->time = 0.0;
	out->energy_access = access_energy;
	out->time_access = access_time;

	// Track the tile the axon interfaces with
	out->t = c->t;

	TRACE("Created axon output %d (c:%d.%d)\n", out->id, c->t->id, c->id);

	return;
}
