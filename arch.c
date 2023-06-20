// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.c
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#include "print.h"
#include "arch.h"
#include "network.h"

/*** All attributes ***/
//  If implementing a new feature, it most likely should be supported in the
//   relevant struct and parsed here from the description
int arch_build_block(struct description_block *block)
{
	/*
	switch (block->type)
	{
	case ARCH_DESCRIPTION_ARCH:
		if (arch_is_field("topology:", d->content[line]))
		{
			arch_parse_field(d->content[line], block->arch->name);
		}
		else if (arch_is_field("dimensions:", d->content[line]))
		{
			block->arch->noc_dimensions =
				arch_parse_field_int(d->content[line]);
		}
		else if (arch_is_field("width:", d->content[line]))
		{
			block->arch->noc_width =
				arch_parse_field_int(d->content[line]);
		}
		else if (arch_is_field("height:", d->content[line]))
		{
			block->arch->noc_height =
				arch_parse_field_int(d->content[line]);
		}
		else
		{
			// Skip line
		}
		break;
	default:
		break;
	}
	return line+1;
	*/
	INFO("not implemented.\n");
	return 1;
}

/*** All architecture init and handling code ***/
struct architecture *arch_init(void)
{
	struct architecture *arch;

	arch = (struct architecture *) malloc(sizeof(struct architecture));
	if (arch == NULL)
	{
		INFO("Error: Architecture couldn't be created.\n");
		exit(1);
	}

	arch->tile_count = 0;
	arch->time_barrier = 0.0;
	arch->initialized = 0;

	for (int i = 0; i < ARCH_MAX_TILES; i++)
	{
		struct tile *t = &(arch->tiles[i]);

		t->energy = 0.0;
		t->time = 0.0;
		t->energy_east_west_hop = 0.0;
		t->time_east_west_hop = 0.0;
		t->energy_north_south_hop = 0.0;
		t->time_north_south_hop = 0.0;
		t->energy_spike_within_tile = 0.0;
		t->time_spike_within_tile = 0.0;
		t->blocked_until = 0.0;
		t->id = -1;
		t->x = -1;
		t->y = -1;
		t->core_count = 0;
		t->is_blocking = 0;
		t->max_dimensions = 0;
		t->width = 0;

		for (int j = 0; j < ARCH_MAX_CORES; j++)
		{
			struct core *c = &(t->cores[j]);
			c->t = NULL;
			c->next_timing = NULL;
			c->energy = 0.0;
			c->time = 0.0;
			c->blocked_until = 0.0;
			c->id = -1;
			c->buffer_pos = 0;
			c->is_blocking = 0;
			c->neuron_count = 0;
			c->curr_neuron = 0;
			c->neurons_left = 0;
			c->curr_axon = 0;

			c->axon_in.energy = 0.0;
			c->axon_in.time = 0.0;
			c->axon_in.packets_in = 0;
			c->axon_in.packet_size = 0;
			c->axon_in.packets_buffer = 0;
			c->axon_in.spikes_buffer = 0;
			c->axon_in.map_count = 0;

			c->synapse.spikes_buffer = 0;
			c->synapse.weights_per_word = 0;
			c->synapse.word_bits = 0;
			c->synapse.weight_bits = 0;
			c->synapse.total_spikes = 0;
			c->synapse.memory_reads = 0;
			c->synapse.energy = 0.0;
			c->synapse.time = 0.0;
			c->synapse.energy_spike_op = 0.0;
			c->synapse.time_spike_op = 0.0;
			c->synapse.energy_memory_access = 0.0;
			c->synapse.time_memory_access = 0.0;

			c->dendrite.energy = 0.0;
			c->dendrite.time = 0.0;

			c->soma.energy = 0.0;
			c->soma.time = 0.0;
			c->axon_out.energy = 0.0;
			c->axon_out.time = 0.0;
		}
	}

	return arch;
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
			c->neurons = NULL;

			for (int k = 0; k < c->axon_in.map_count; k++)
			{
				struct axon_map *a = &(c->axon_in.map[k]);
				free(a->connections);
				a->connections = NULL;
			}
		}
	}
	free(arch);
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
	c->neurons = (struct neuron **) malloc(sizeof(struct neuron *) *
							ARCH_MAX_COMPARTMENTS);
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
	in->map_count = 0;
	// We already know a valid tile was given at this point
	in->t = c->t;

	TRACE("Created axon input (c:%d.%d)\n", c->t->id, c->id);

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

	s->energy_memory_access = energy_memory_read;
	s->time_memory_access = time_memory_read;

	// The word size is the number of bits accessed with each memory read.
	//  The weight size is the number of bits for a single synaptic weight.
	//  A single memory read might return multiple weights
	s->weight_bits = weight_bits;
	s->word_bits = word_bits;
	// Round up to the nearest word
	s->weights_per_word = (word_bits + (weight_bits - 1)) / weight_bits;

	TRACE("Created synapse processor (c:%d.%d)\n", c->t->id, c->id);

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
	// TODO: figure this parameter for setting the leak
	s->leak_towards_zero = 1;

	TRACE("Created soma processor (c:%d.%d)\n", c->t->id, c->id);

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
	out->map_count = 0;

	// Track the tile the axon interfaces with
	out->t = c->t;

	TRACE("Created axon output (c:%d.%d)\n", c->t->id, c->id);

	return;
}

int arch_map_neuron(struct neuron *const n, const struct hardware_mapping map)
{
	// Map the neuron to hardware units
	assert(map.core != NULL);
	assert(map.core->neurons != NULL);
	assert(n->core == NULL);
	n->core = map.core;
	TRACE("mapping neuron %d to core %d\n", n->id, map.core->id);
	map.core->neurons[map.core->neuron_count] = n;
	map.core->neuron_count++;

	n->axon_in = map.axon_in;
	n->synapse_hw = map.synapse_hw;
	n->dendrite_hw = map.dendrite_hw;
	n->soma_hw = map.soma_hw;
	n->axon_out = map.axon_out;

	return 1;
}

void arch_create_axon_maps(struct architecture *const arch)
{
	INFO("Creating axon map.\n");
	// Initialize the axon structures
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			for (int i = 0; i < ARCH_MAX_AXON_MAP; i++)
			{
				c->axon_in.map[i].connection_count = 0;
				c->axon_in.map[i].spikes_received = 0;
				c->axon_in.map[i].receive_latency = 0.0;
				c->axon_out.map_ptr[i] = NULL;
				c->axon_in.map[i].last_updated = 0;
				c->axon_in.map[i].pre_neuron = NULL;
			}
		}
	}

	// TODO: think of a meaningful way of refactoring this monolithic
	//  block
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			for (int k = 0; k < c->neuron_count; k++)
			{
				struct neuron *pre_neuron = c->neurons[k];
				// Track the connections
				int connection_count[ARCH_MAX_TILES*ARCH_MAX_CORES];
				struct core *cores[ARCH_MAX_TILES*ARCH_MAX_CORES];

				for (int x = 0; x < ARCH_MAX_TILES * ARCH_MAX_CORES; x++)
				{
					connection_count[x] = 0;
					cores[x] = NULL;
				}

				TRACE("Creating axons for neuron %d\n", k);
				// Count how many connections go to each core
				//  from this neuron. This will be used to
				//  allocate the axon maps and the connections
				//  in each map
				for (int conn = 0;
					conn < pre_neuron->connection_out_count;
					conn++)
				{
					TRACE("Looking at connection id: %d\n",
									conn);
					struct connection *curr =
						&(pre_neuron->connections_out[conn]);
					struct core *dest_core =
							curr->post_neuron->core;
					// TODO: HACK want a unique core identifier beyond the
					//  number within the tile. Need to index the correct
					//  core number here
					int core_number = (dest_core->t->id * 4) + dest_core->id;
					connection_count[core_number]++;
					cores[core_number] = dest_core;
					TRACE("Connected to dest core: %d\n",
								core_number);
				}

				int axon_count = 0;
				for (int x = 0; x < ARCH_MAX_TILES * ARCH_MAX_CORES; x++)
				{
					// Now for each connected core, create
					//  a new axon map at the destination
					//  core. Then link this axon to output
					//  of the source core. Finally update
					//  the presynaptic neuron and
					//  postsynaptic neuron to account for
					//  this
					if (connection_count[x] > 0)
					{
						// Create the axon map, and add
						//  it to the map at the dest
						//  core.
						//struct neuron *post_neuron = curr->post_neuron;
						struct core *dest_core = cores[x];
						struct axon_input *axon_in = &(dest_core->axon_in);
						int map_count = axon_in->map_count++;
						TRACE("axon in map count:%d for core:%d.%d, adding %d connections\n",
							map_count, dest_core->id, dest_core->t->id, connection_count[x] );
						struct axon_map *a = &(axon_in->map[map_count]);

						//INFO("Adding connection to core.\n");
						// Allocate the entry and its connections
						//INFO("Axon has %d connections; allocating %lu bytes\n",
						//	connection_count[x], connection_count[x] * sizeof(struct connection));
						a->connections = malloc(connection_count[x] * sizeof(struct connection));
						if (a->connections == NULL)
						{
							INFO("Error: Couldn't allocate axon memory.\n");
							exit(1);
						}

						// Now create the link to this map in the
						//  pre-synaptic core
						map_count = c->axon_out.map_count++;
						c->axon_out.map_ptr[map_count] = a;
						if (pre_neuron->maps_out == NULL)
						{
							TRACE("Setting neuron nid:%d axon out.\n", pre_neuron->id);
							pre_neuron->maps_out =
								&(c->axon_out.map_ptr[map_count]);
							assert(pre_neuron->maps_out != NULL);
							assert(pre_neuron->maps_out[0] != NULL);
						}
						pre_neuron->maps_out_count++;
						TRACE("nid:%d.%d cid:%d.%d added one output axon, axon out map_count:%d, neuron out map count:%d.\n",
							pre_neuron->group->id, pre_neuron->id, c->t->id, c->id, c->axon_out.map_count, pre_neuron->maps_out_count);
						axon_count++;
					}
				}
				//INFO("nid:%d axon count: %d\n", pre_neuron->id, axon_count);
				assert(axon_count < ARCH_MAX_AXON_MAP);

				for (int conn = 0; conn < pre_neuron->connection_out_count; conn++)
				{
					// For each connection, add the connection to the axon.
					//  Also track the axon in the post synaptic neuron
					struct connection *curr =
						&(pre_neuron->connections_out[conn]);
					struct core *p = curr->post_neuron->core;
					//INFO("adding connection:%d\n", conn);
					// just add to the current axon
					int map_count = p->axon_in.map_count;
					if (map_count <= 0 || map_count > ARCH_MAX_AXON_MAP)
					{
						TRACE("map_count:%d\n", map_count);
					}
					assert(map_count > 0);
					assert(map_count <= ARCH_MAX_AXON_MAP);
					//INFO("adding to connection to axon:%d\n",
					//			map_count - 1);
					// Access the most recently created axon
					//  for the core
					struct axon_map *a =
						&(p->axon_in.map[map_count-1]);
					a->connections[a->connection_count++] =
									curr;
					a->pre_neuron = pre_neuron;

					// Update the post synaptic neuron to
					//  track
					if (curr->post_neuron->maps_in == NULL)
					{
						// Point to the first mapping
						curr->post_neuron->maps_in = a;
					}
					// We might add a bunch of connections
					//  from another core coming into this
					//  one, then we need to update and
					//  track
					curr->post_neuron->maps_in_count++;
				}
			}
		}
	}

	TRACE("Created all axons\n");

	// Print summary about axons created
	// TODO: refactor?
	INFO("** Mapping summary **\n");
	int in_count, out_count, core_count, core_used;
	in_count = out_count = core_count = 0;
	for (int i = 0; i < arch->tile_count; i++)
	{
		// For debug only, print the axon maps
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);

			core_used = 0;
			for (int k = 0; k < c->neuron_count; k++)
			{
				//struct neuron *n = c->neurons[k];
				//TRACE("\tnid:%d.%d ", n->group->id, n->id);
				//TRACE("i:%d o:%d\n",
				//	n->maps_in_count, n->maps_out_count);
				core_used = 1;
			}

			if (core_used)
			{
				INFO("cid:%d.%d n:%d i:%d o:%d\n", t->id, c->id,
					c->neuron_count, c->axon_in.map_count,
					c->axon_out.map_count);
				in_count += c->axon_in.map_count;
				out_count += c->axon_out.map_count;
				core_count++;
			}
		}
	}
	INFO("Total cores: %d\n", core_count);
	INFO("Average in map count: %lf\n", (double) in_count / core_count);
	INFO("Average out map count: %lf\n", (double) out_count / core_count);
}

/* *** Architecture Description File Parsing *** */

static struct description_block ARCH_BLOCKS[ARCH_MAX_BLOCK_COUNT];
static int ARCH_BLOCK_COUNT;

void arch_print_description(struct description_block *arch_description,
								const int level)
{
	if (level == 0)
	{
		INFO("Printing arch description.\n");
	}
	if (arch_description->child != NULL)
	{
		arch_print_description(arch_description->child, level+1);
	}

	if (level > 0)
	{
		for (int t = 0; t < level; t++)
		{
			printf("\t");
		}
		printf("%s:[", arch_description->name);
		for (int i = 0; i < arch_description->attribute_count; i++)
		{
			struct attributes *attr = &(arch_description->attributes[i]);
			printf("%s:%s,", attr->key, attr->value_str);
		}
		printf("]\n");
	}

	if (arch_description->next_sibling != NULL)
	{
		arch_print_description(arch_description->next_sibling, level);
	}

	return;
}

int arch_parse_file(FILE *fp,struct description_block *arch_description)
{
	char line_buf[ARCH_MAX_DESCRIPTION_LINE+1];

	assert(fp != NULL);
	memset(line_buf, 0, sizeof(line_buf));
	// Use a dummy block to track the entire design
	arch_description->child = NULL;
	arch_description->parent = NULL;
	arch_description->instances = 1;
	arch_description->name[0] = '\0';
	arch_description->type = ARCH_DESCRIPTION_TOP;
	arch_description->attribute_count = 0;

	// Reset global counter
	ARCH_BLOCK_COUNT = 0;
	// Ignore everything preceding or following the 'architecture' block
	while (arch_get_description_line(fp, line_buf) >= 0)
	{
		INFO("line:%s", line_buf);
		if (arch_is_field("architecture:", line_buf))
		{
			arch_description->child = arch_parse_list(fp,
							ARCH_DESCRIPTION_ARCH,
							arch_description);
			return arch_build_arch(arch_description);
		}
		// else ignore everything preceding the main architecture block
	}

	// The block should be defined exactly once. Currently I don't bother
	//  checking for multiple definitions (the first one will be used)
	INFO("Error: No architecture block defined.\n");
	return 1;
}

int arch_get_description_line(FILE *fp, char *line)
{
	char line_buf[ARCH_MAX_DESCRIPTION_LINE];
	int valid_line, indent;
	memset(line_buf, 0, sizeof(line_buf));

	indent = 0;
	// Keep reading in lines until the file stops or we find a valid line
	valid_line = 0;
	while ((!valid_line) && fgets(line_buf, ARCH_MAX_DESCRIPTION_LINE, fp))
	{
		char *str = line_buf;
		valid_line = 0;
		INFO("Read line:%s\n", line_buf);

		// Get the first non-space character after indents
		for (int i = 0; i < ARCH_MAX_DESCRIPTION_LINE; i++)
		{
			if (*str != ' ')
			{
				// Only spaces can be used to indent
				break;
			}
			else
			{
				indent++;
				str++;
			}
		}
		// Check to see if this is a valid line - skip spaces, comments
		if (((*str) == '\0') || ((*str) == '#') || isspace(*str))
		{
			INFO("Ignoring line:%s", line_buf);
		}
		else
		{
			valid_line = 1;
			strncpy(line, str, ARCH_MAX_DESCRIPTION_LINE-indent);
		}
	}

	if (!valid_line)
	{
		// Negative indent means no line was found
		indent = -1;
	}
	// else don't touch the block struct at all

	return indent;
}

struct description_block *arch_parse_list(FILE *fp, const int list_type,
				struct description_block *const parent)
{
	char line_buf[ARCH_MAX_DESCRIPTION_LINE];
	struct description_block *head, *sibling, *block;

	// If the parent block already has children, we will add onto this list
	sibling = parent->child;
	head = sibling;
	memset(line_buf, 0, sizeof(line_buf));

	// Parse all tiles in the list
	while(1)
	{
		int line_len, indent, block_finished;

		indent = arch_get_description_line(fp, line_buf);
		TRACE("Parsing:%s (indent:%d)\n", line_buf, indent);
		line_len = strnlen(line_buf, ARCH_MAX_DESCRIPTION_LINE) +indent;
		// Reset the fp to point to the beginning of the line
		fseek(fp, 0-line_len, SEEK_CUR);

		block_finished = (indent <= parent->indent);
		if (block_finished)
		{
			TRACE("Block finished.\n");
			return head;
		}
		else if (line_buf[indent] != '-')
		{
			// Not dealing with a list entry, so just parse block
			TRACE("Parsing single block.\n");
			block = arch_parse_block(fp, list_type, parent,sibling);
			return block;
		}
		else
		{
			TRACE("Parsing block list.\n");
			// Store a linked list of blocks
			block = arch_parse_block(fp, list_type, parent,sibling);
			sibling = block;
			head = block;
		}
	}

	return head;
}

struct description_block *arch_parse_block(FILE *fp,
	const int block_type, struct description_block *const parent,
	struct description_block *const sibling)
{
	char line_buf[ARCH_MAX_DESCRIPTION_LINE];
	int block_finished;
	struct description_block *block;

	// Initialize a new block
	block = &(ARCH_BLOCKS[ARCH_BLOCK_COUNT]);
	ARCH_BLOCK_COUNT++;
	block->parent = parent;
	block->next_sibling = sibling;
	block->type = block_type;
	block->instances = 1;
	block->attribute_count = 0;

	block->indent = arch_get_description_line(fp, line_buf);
	block_finished = (block->indent <= parent->indent);
	TRACE("Created block:%d indent:%d, parent->indent:%d block_finished:%d\n",
		ARCH_BLOCK_COUNT-1, block->indent, parent->indent,
		block_finished);
	while (!block_finished)
	{
		char *str = line_buf;

		TRACE("Parsing line:%s (indent:%d finished:%d)\n",
			line_buf, block->indent, block_finished);
		if (str[0] == '-')
		{
			// Skip any list characters
			str++;
			while (*str == ' ') str++;
		}

		// Parse name field
		if (arch_is_field("name:", str))
		{
			block->instances = arch_parse_name(str, block->name);
		}
		// Parse attributes field
		else if (arch_is_field("attributes:", str))
		{
			arch_parse_attributes(fp, block);
		}
		// Parse any sub-blocks
		else if (arch_is_field("tile:", str))
		{
			block->child = arch_parse_list(fp,
				ARCH_DESCRIPTION_TILE, block);
		}
		else if (arch_is_field("core:", str))
		{
			block->child = arch_parse_list(fp,
				ARCH_DESCRIPTION_CORE, block);
		}
		// These need to be siblings
		else if (arch_is_field("axon_in:", str))
		{
			block->child = arch_parse_list(fp,
				ARCH_DESCRIPTION_AXON_IN, block);
		}
		else if (arch_is_field("synapse:", str))
		{
			block->child = arch_parse_list(fp,
				ARCH_DESCRIPTION_SYNAPSE, block);
		}
		else if (arch_is_field("dendrite:", str))
		{
			block->child = arch_parse_list(fp,
				ARCH_DESCRIPTION_DENDRITE, block);
		}
		else if (arch_is_field("soma:", str))
		{
			block->child = arch_parse_list(fp,
				ARCH_DESCRIPTION_SOMA, block);
		}
		else if (arch_is_field("axon_out:", str))
		{
			block->child = arch_parse_list(fp,
				ARCH_DESCRIPTION_AXON_OUT, block);
		}
		block->indent = arch_get_description_line(fp, line_buf);
		block_finished = (block->indent <= parent->indent);
	}

	// Put the line back
	int line_len = strnlen(line_buf,
				(ARCH_MAX_DESCRIPTION_LINE + block->indent));
	fseek(fp, 0-line_len, SEEK_CUR);

	return block;
}

int arch_parse_attributes(FILE *fp, struct description_block *block)
{
	char line_buf[ARCH_MAX_DESCRIPTION_LINE];
	int indent, attributes_finished, line_len;
	memset(line_buf, 0, sizeof(line_buf));

	indent = arch_get_description_line(fp, line_buf);
	attributes_finished = (indent <= block->indent);
	while (!attributes_finished)
	{
		int a = block->attribute_count;
		char *curr = line_buf;
		struct attributes *attr = &(block->attributes[a]);
		// Store a key:value pair
		// key = ;
		while ((*curr != '\0') && (*curr != ':')) curr++;
		if (*curr == ':')
		{
			*curr = '\0';
			curr++;
		}
		else
		{
			INFO("Error: Badly formed attribute (%s).\n", line_buf);
			exit(1);
		}
		strncpy(attr->key, line_buf, ARCH_MAX_FIELD_LEN);
		// Copy, stripping trailing whitespace
		while ((*curr != '\0') && isspace(*curr)) curr++;
		for (int i=0; i < ARCH_MAX_FIELD_LEN; i++)
		{
			if (((curr[i]) == '\0') || isspace(curr[i]))
			{
				break;
			}
			else
			{
				attr->value_str[i] = curr[i];
			}
		}

		TRACE("Parsed pair: %s,%s\n", attr->key, attr->value_str);

		block->attribute_count++;
		indent = arch_get_description_line(fp, line_buf);
		TRACE("Line: %s (indent:%d block->indent:%d)\n",
			line_buf, indent, block->indent);
		attributes_finished = (indent <= block->indent);
	}

	line_len = strnlen(line_buf, ARCH_MAX_DESCRIPTION_LINE) + indent;
	fseek(fp, 0-line_len, SEEK_CUR);
	TRACE("Finished parsing attributes.\n");

	return block->attribute_count;
}

int arch_is_field(const char *fieldname, char *str)
{
	// Returns 1 (true) if the string matches the given fieldname, 0 (false)
	//  otherwise
	return (strncmp(str, fieldname,
			strnlen(fieldname, ARCH_MAX_DESCRIPTION_LINE)) ==0);
}

int arch_parse_name(char *str, char *dest)
{
	char name_field[ARCH_MAX_FIELD_LEN];
	char *curr;
	int range_start, range_end, instances;

	arch_parse_field(str, name_field);
	// Check to see if a range has been specified, otherwise assume just
	//  one instance
	instances = 1;
	for (curr = name_field; *curr != '\0'; curr++)
	{
		if (curr[0] == '[')
		{
			if (sscanf(curr, "[%d..%d]",
					&range_start, &range_end) == 2)
			{
				instances = (range_end - range_start) + 1;
				INFO("Block has %d instances,start:%d,end:%d\n",
					instances, range_start, range_end);
				// Replace the opening bracket with null; chop
				//  this off
				while (*curr != '[')
				{
					curr++;
				}
				*curr = '\0';
			}
			else
			{
				INFO("Error: invalid range (%s)\n", curr);
				exit(1);
			}
		}
	}

	// TODO: refactor
	// Copy the name across, discarding any trailing whitespace
	for (int i=0; i < ARCH_MAX_FIELD_LEN; i++)
	{
		if (((name_field[i]) == '\0') || isspace(name_field[i]))
		{
			break;
		}
		else
		{
			dest[i] = name_field[i];
		}
	}

	INFO("Parsed name:'%s' (instances:%d).\n", dest, instances);
	return instances;
}


int arch_parse_field_int(char *str)
{
	char field_str[ARCH_MAX_FIELD_LEN];

	arch_parse_field(str, field_str);
	return strtol(str, NULL, 10);
}

void arch_parse_field(char *str, char *field)
{
	// Extract a value from a <field>:<value> pair, passed in a string.
	//  Store in the field string.
	char *curr = str;

	// Read up to the colon operator
	while (*curr != '\0' && *curr != ':')
	{
		curr++;
	}

	if (*curr == ':')
	{
		curr++;
	}
	else
	{
		INFO("Error: Not a valid field: %s\n", str);
		exit(1);
	}

	// Ignore whitespace before the field value
	while (*curr != '\0' && isspace(*curr))
	{
		curr++;
	}

	if (*curr == '\0')
	{
		INFO("Error: invalid field (%s)\n", str);
		exit(1);
	}

	// Copy the value across
	while (*curr != '\0' && !isspace(*curr))
	{
		*(field++) = *(curr++);
	}

	// End the string
	*field = '\0';
	TRACE("Parsed field %s from: %s\n", field, str);

	return;
}

int arch_build_arch(struct description_block *arch_description)
{
	return 0;
}
