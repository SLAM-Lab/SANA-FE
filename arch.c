// arch.c: Create a neuromorphic design based on an architecture description
//  In this simulator and architecture description is a list of different
//  blocks.

// TODO: in the future we might be able to skip this whole step, and use a
//  Python interface to directly initialize the C structures
// TODO: get rid of the tokenization code, I think its safe just to assume that
//  each entry is separated by a single space (it's meant to be a dumb machine
//  readable format anyway)
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#include "arch.h"
#include "sim.h"

static void arch_read_line(struct architecture *arch, char *line);
static struct architecture arch_init(FILE *fp);

static void arch_parse_neuron(struct architecture *arch, char fields[][ARCH_MAX_FIELD_LEN]);
//static void arch_parse_synapse(struct architecture *arch, char fields[][ARCH_MAX_FIELD_LEN]);
static void arch_parse_axon_input(struct architecture *arch, char fields[][ARCH_MAX_FIELD_LEN]);
static void arch_parse_axon_output(struct architecture *arch, char fields[][ARCH_MAX_FIELD_LEN]);
static void arch_parse_router(struct architecture *arch, char fields[][ARCH_MAX_FIELD_LEN]);

struct architecture arch_read_file(FILE *fp)
{
	// Read an architecture description and return an initialized structure
	//  representing the neuromorphic architecture
	struct architecture arch = arch_init(fp);
	char line[ARCH_LINE_LEN];

	while (fgets(line, ARCH_LINE_LEN, fp))
	{
		arch_read_line(&arch, line);
	}

	return arch;
}

// TODO: I don't know if this is a wasted effort. We could get away with only
//  supporting ranges, and the only place ranges are really useful is for
//  neuron compartments (of which there are many). This is a lot of complicated
//  code and structures for debatable gain. Maybe get rid of this
static int arch_parse_list(const char *field, struct range *list)
{
	char value[ARCH_MAX_VALUE_DIGITS];
	int is_range, list_entries, v, curr_char;

	list_entries = 0;
	curr_char = 0;
	is_range = 0;
	v = 0;
	for (curr_char = 0; curr_char < ARCH_MAX_FIELD_LEN; curr_char++)
	{
		// Parse a field, which is a comma-separated list containing
		// 1) values e.g. (3)
		// 2) value ranges e.g. (3..7)
		if ((field[curr_char] == ',') || (field[curr_char] == '\0'))
		{
			struct range *curr_range = &(list[list_entries]);

			value[v] = '\0';
			// Parse the previous number
			if (is_range)
			{
				int ret = sscanf(value, "%u..%u",
					&(curr_range->min), &(curr_range->max));
				if (ret < 2)
				{
					INFO("Error: Parsing range (%s).\n",
									field);
					exit(1);
				}
				TRACE("Parsed range of values (%u to %u)\n",
					curr_range->min, curr_range->max);
			}
			else
			{
				int ret = sscanf(value, "%u",
							&curr_range->min);
				if (ret < 1)
				{
					INFO("Error: Parsing value (%s).\n",
									field);
					exit(1);
				}
				// Just one value is read i.e. min == max
				curr_range->max = curr_range->min;
				TRACE("Parsed value (%u)\n", curr_range->min);
			}

			v = 0;
			is_range = 0;
			list_entries++;
			if (field[curr_char] == '\0')
			{
				// End of string, stop processing
				break;
			}
		}
		else // copy string into temporary buffer for parsing
		{
			value[v] = field[curr_char];
			if ((field[curr_char] == '.') &&
						(field[curr_char+1] == '.'))
			{
				// This entry is a range of values e.g. 1..3
				is_range = 1;
			}
			v++;
		}
	}

	return list_entries;
}

static void arch_read_line(struct architecture *arch, char *line)
{
	char fields[ARCH_MAX_FIELDS][ARCH_MAX_FIELD_LEN];
	char *curr_field;
	char block_type;

	block_type = line[0];
	// Sanity check input
	if ((line == NULL) || (block_type == '\0') || (block_type == '#'))
	{
		TRACE("Warning: Empty line - skipping.\n");
		return;
	}

	// Read all space separated, machine readable fields
	curr_field = strtok(line, " ");
	for (int i = 0; i < ARCH_MAX_FIELDS; i++)
	{
		if (curr_field == NULL)
		{
			TRACE("Parsed %d fields.\n", i);
			break;
		}
		strncpy(fields[i], curr_field, ARCH_MAX_FIELD_LEN);
		curr_field = strtok(NULL, " ");
	}

	// Process list of identifiers
	switch (block_type)
	{
	case '\0':
	case '#':
	case '\n':
		// Line is a comment
		break;
	case 'c':
		arch_parse_neuron(arch, fields);
		break;
	case 's':
		//arch_parse_synapse(arch, fields);
		break;
	case 'r':
		arch_parse_router(arch, fields);
		break;
	case 'i':
		arch_parse_axon_input(arch, fields);
		break;
	case 'o':
		arch_parse_axon_output(arch, fields);
		break;
	case 'e':
		//arch_parse_external_input(arch, id_list, line);
		break;
	default:
		TRACE("Warning: unrecognized unit (%c) - skipping.\n",
							block_type);
		break;
	}
}

static unsigned int arch_get_count(const char *field)
{
	const char *curr = field;
	unsigned int values;
	// Count how many values the field has

	// A field can either be a number, a comma separated list,
	//  or a range of values. For simplicity, this format doesn't allow you
	//  to combine lists and ranges
	values = 0;
	while (values < ARCH_MAX_VALUES)
	{
		// Either count the number of commas for lists, but if a range
		//  parse the left and right hand values. Otherwise it's just a
		//  single value.
		if (curr == NULL)
		{
			break;
		}
		if ((curr[0] == '\n') || (curr[0] == '\0'))
		{
			values++;
			break;
		}
		else if (isdigit(curr[0]))
		{
			curr++;
		}
		else if (curr[0] == ',')
		{
			values++;
			curr++;
		}
		else if ((curr[0] == '.') && (curr[1] == '.'))
		{
			int ret;
			unsigned int lval, rval;

			ret = sscanf(field, "%d..%d", &lval, &rval);
			if (ret)
			{
				// Calculate the range (inclusive), e.g.
				// 1..3 -> 1,2,3 which is 3 values
				values = (rval - lval) + 1;
				TRACE("Parsed range (%d..%d).\n", lval, rval);
			}
			else
			{
				values = 0;
				TRACE("Warning: invalid range (%s).\n", curr);
			}
			break;
		}
		else
		{
			TRACE("Warning: invalid field (%s).\n", field);
			break;
		}
	}

	TRACE("Parsed %d values.\n", values);
	return values;
}

static struct architecture arch_init(FILE *fp)
{
	// Parse the architecture file, counting the number of units needed
	//  and then allocating the memory
	struct architecture arch;
	char line[ARCH_LINE_LEN];
	char *field;

	arch.compartments = NULL;
	arch.mem_blocks = NULL;
	arch.routers = NULL;
	arch.axon_inputs = NULL;
	arch.axon_outputs = NULL;
	arch.external_inputs = NULL;

	arch.max_compartments = 0;
	arch.max_mem_blocks = 0;
	arch.max_routers = 0;
	arch.max_timers = 0;
	arch.max_axon_inputs = 0;
	arch.max_axon_outputs = 0;
	arch.max_external_inputs = 131072; // HACK

	// TODO: refactor - arch_count_units(fp) ?
	while (fgets(line, ARCH_LINE_LEN, fp))
	{
		unsigned int unit_count;
		char block_type = line[0];

		field = strtok(line, " ");
		field = strtok(NULL, " ");
		unit_count = arch_get_count(field);
		TRACE("unit count:%d.\n", unit_count);

		// TODO: can refactor so the switch sets an enum?
		//  That doesn't really add anything though since we'd need
		//  another switch
		switch (block_type)
		{
		case '\0':
		case '#':
		case '\n':
			// Line is a comment
			continue;
		case 'c':
			arch.max_compartments += unit_count;
			break;
		case 'm':
			arch.max_mem_blocks += unit_count;
			break;
		case 'r':
			arch.max_routers += unit_count;
			break;
		case 'i':
			arch.max_axon_inputs += unit_count;
			break;
		case 'o':
			arch.max_axon_outputs += unit_count;
			break;
		case 'e':
			arch.max_external_inputs += unit_count;
			break;
		case 't':
			arch.max_timers += unit_count;
			break;
		case 'd':
			break;
		case 's':
			break;
		default:
			TRACE("Warning: unrecognized unit (%c) - skipping.\n",
								block_type);
			break;
		}
	}

	INFO("Parsed %d compartments.\n", arch.max_compartments);

	// Reset file pointer to start of file
	fseek(fp, 0, SEEK_SET);

	// Based on the number of different units, allocate enough memory to
	//  simulate this design
	INFO("Allocating memory for %d compartments.\n", arch.max_compartments);
	arch.compartments = (struct compartment *)
		malloc(arch.max_compartments * sizeof(struct compartment));
	arch.mem_blocks = (struct mem *)
		malloc(arch.max_mem_blocks * sizeof(struct mem));
	arch.routers = (struct router *)
		malloc(arch.max_routers * sizeof(struct router));
	arch.axon_inputs = (struct axon_input *)
		malloc(arch.max_axon_inputs * sizeof(struct axon_input));
	arch.axon_outputs = (struct axon_output *)
		malloc(arch.max_axon_outputs * sizeof(struct axon_output));
	arch.timers = (double *) malloc(arch.max_timers * sizeof(double));

	if ((arch.compartments == NULL) || (arch.mem_blocks == NULL) ||
		(arch.routers == NULL) || (arch.timers == NULL) ||
		(arch.axon_inputs == NULL) || (arch.axon_outputs == NULL))
	{
		INFO("Error: Failed to allocate compartment.\n");
		exit(1);
	}

	for (int i = 0; i < arch.max_compartments; i++)
	{
		struct compartment *c = &(arch.compartments[i]);

		c->id = i;
		c->synapses = NULL;
		c->axon_in = NULL;
		c->axon_out = NULL;

		c->fired = 0;
		c->update_needed = 0;
		c->compartment_used = 0;
	}

	for (int i = 0; i < arch.max_mem_blocks; i++)
	{
		struct mem *mem = &(arch.mem_blocks[i]);

		mem->id = i;
		TRACE("Allocating synapse memory for block %d.\n", mem->id);
	}

	for (int i = 0; i < arch.max_axon_outputs; i++)
	{
		struct axon_output *axon_out = &(arch.axon_outputs[i]);

		axon_out->packets_sent = (unsigned int *)
			malloc(arch.max_axon_inputs * sizeof(unsigned int));
		if (axon_out->packets_sent == NULL)
		{
			INFO("Error: Failed to allocate axon memory.\n");
			exit(1);
		}
	}

	INFO("Allocating memory for %d external inputs.\n",
						arch.max_external_inputs);
	arch.external_inputs = (struct input *)
		malloc(arch.max_external_inputs * sizeof(struct input));
	if (arch.external_inputs == NULL)
	{
		INFO("Error: Failed to allocate input memory.\n");
		exit(1);
	}

	// Zero initialize the input nodes
	for (int i = 0; i < arch.max_external_inputs; i++)
	{
		struct input *in = &(arch.external_inputs[i]);

		in->synapses = NULL;
		in->post_connection_count = 0;
		in->send_spike = 0;
	}

	for (int i = 0; i < arch.max_timers; i++)
	{
		arch.timers[i] = 0.0;
	}

	arch.initialized = 1;

	return arch;
}

void arch_free(struct architecture *arch)
{
	// Free the architecture structure, deallocate all memory and mark
	//  the structure as uninitialized
	arch->initialized = 0;

	// Free any neuron data
	for (int i = 0; i < arch->max_compartments; i++)
	{
		struct compartment *c = &(arch->compartments[i]);

		assert(c != NULL);
		free(c->synapses);
		c->synapses = NULL;
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
	free(arch->compartments);
	free(arch->mem_blocks);
	free(arch->routers);
	free(arch->axon_inputs);
	free(arch->axon_outputs);
	free(arch->external_inputs);
	free(arch->timers);

	// Reset all pointers and counters
	arch->compartments = NULL;
	arch->mem_blocks = NULL;
	arch->routers = NULL;
	arch->axon_inputs = NULL;
	arch->axon_outputs = NULL;
	arch->external_inputs = NULL;
	arch->timers = NULL;

	arch->max_compartments = 0;
	arch->max_mem_blocks = 0;
	arch->max_routers = 0;
	arch->max_timers = 0;
	arch->max_axon_inputs = 0;
	arch->max_axon_outputs = 0;
	arch->max_external_inputs = 0;
}

static void arch_parse_neuron(struct architecture *arch,
					char fields[][ARCH_MAX_FIELD_LEN])
{
	struct range field_list[ARCH_MAX_VALUES];
	//struct mem *mem_block;
	struct axon_input *axon_in;
	struct axon_output *axon_out;
	double *timer;
	int ret, list_len, mem_id, axon_in_id, axon_out_id, timer_id;
	char block_type;

	assert(fields && fields[0]);
	block_type = fields[0][0];
	assert(block_type == 'c');

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

	// Copy the neuron a number of times
	list_len = arch_parse_list(fields[1], &(field_list[0]));
	for (int i = 0; i < list_len; i++)
	{
		struct range *r = &(field_list[i]);

		TRACE("min:%u max%u\n", r->min, r->max);
		for (int id = r->min; id <= r->max; id++)
		{
			struct compartment *c;

			assert(id <= arch->max_compartments);
			c = &(arch->compartments[id]);
			// Copy any common compartment variables
			//c->mem_block = mem_block;
			c->axon_in = axon_in;
			c->axon_out = axon_out;
			c->time = timer;
		}
		TRACE("Compartment parsed n:%d-%d mem(%d) i(%d) o(%d).\n",
			r->min, r->max, mem_id, axon_in_id, axon_out_id);
	}
}

static void arch_parse_axon_input(struct architecture *arch,
					char fields[][ARCH_MAX_FIELD_LEN])
{
	struct axon_input *axon_in;
	int ret, id, router_id;
	char type;

	assert(fields && fields[0]);
	type = fields[0][0];
	assert(type == 'i');


	ret = sscanf(fields[1], "%u", &id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon input id (%s).\n", fields[1]);
		exit(1);
	}
	assert(id < arch->max_axon_inputs);
	axon_in = &(arch->axon_inputs[id]);
	axon_in->id = id;

	ret = sscanf(fields[2], "%u", &router_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon router (%s).\n", fields[2]);
		exit(1);
	}
	assert(router_id < arch->max_routers);
	axon_in->r = &(arch->routers[router_id]);

	ret = sscanf(fields[3], "%d", &axon_in->fan_in);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon input fan in (%s).\n",
								fields[3]);
		exit(1);
	}

	TRACE("Axon input parsed i:%d r(%d) fan_in(%d)\n", axon_in->id,
					axon_in->r->id, axon_in->fan_in);
}

static void arch_parse_axon_output(struct architecture *arch,
					char fields[][ARCH_MAX_FIELD_LEN])
{
	struct axon_output *axon_out;
	int ret, id, router_id;
	char type;

	assert(fields && fields[0]);
	type = fields[0][0];
	assert(type == 'o');

	ret = sscanf(fields[1], "%d", &id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output id (%s).\n", fields[1]);
		exit(1);
	}
	assert(id < arch->max_axon_outputs);
	axon_out = &(arch->axon_outputs[id]);
	axon_out->id = id;

	ret = sscanf(fields[2], "%d", &router_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon router (%s).\n", fields[2]);
		exit(1);
	}
	assert(router_id < arch->max_routers);
	axon_out->r = &(arch->routers[router_id]);

	ret = sscanf(fields[3], "%d", &axon_out->fan_out);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output fan out (%s).\n",
								fields[3]);
		exit(1);
	}

	TRACE("Axon output parsed o:%d r(%d) fan_out(%d)\n", axon_out->id,
					axon_out->r->id, axon_out->fan_out);
}

static void arch_parse_router(struct architecture *arch,
					char fields[][ARCH_MAX_FIELD_LEN])
{
	struct router *r;
	int id, ret;

	ret = sscanf(fields[1], "%u", &id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output id (%s).\n", fields[1]);
		exit(1);
	}
	assert(id < arch->max_routers);
	r = &(arch->routers[id]);
	r->id = id;

	// TODO: change to single scanf (everywhere, not just in this function)
	ret = sscanf(fields[2], "%u", &(r->x));
	ret += sscanf(fields[3], "%u", &(r->y));
	if (ret < 2)
	{
		INFO("Error: couldn't parse router r:%d.\n", id);
		exit(1);
	}
	TRACE("Router parsed r:%u x(%d) y(%d).\n", r->id, r->x, r->y);
}

/*
static void arch_parse_external_inputs(struct architecture *arch, char *line)
{

}

void arch_parse_synapse(struct architecture *arch, char *line)
{

}

*/
