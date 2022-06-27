// arch.c: Create a neuromorphic design based on an architecture description
//  In this simulator and architecture description is a list of different
//  blocks.
// TODO: in the future we might be able to skip this whole step, and use a
//  Python interface to directly initialize the C structures
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

enum block_type
{
	NEURON = 0,
	synapse_mem,
	DENDRITE,
	AXON_IN,
	AXON_OUT,
	ROUTER,
	TIMER,
	EXTERNAL_INPUT,
	N_BLOCKS,
};

/*
const char block_letters[] =
{
	'n', // neuron
	's', // synapse block
	'd', // dendrite
	'i', // axon inputs
	'o', // axon outputs
	'r', // routers
	't', // timing
	'\0',
};
*/

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
			INFO("Parsed %d fields.\n", i);
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
	case 'n':
		arch_parse_neuron(arch, fields);
		break;
	case 's':
		//arch_parse_synapse(arch, fields);
		break;
	case 'r':
		//arch_parse_router(arch, id_list, line);
		break;
	case 'i':
		//arch_parse_axon_input(arch, id_list, line);
		break;
	case 'o':
		//arch_parse_axon_output(arch, id_list, line);
		break;
	case 't':
		//arch_parse_timer(arch, id_list, line);
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
				INFO("Parsed range (%d..%d).\n", lval, rval);
			}
			else
			{
				values = 0;
				INFO("Error: invalid range (%s).\n", curr);
			}
			break;
		}
		else
		{
			INFO("Error: invalid field (%s).\n", field);
			break;
		}
	}

	INFO("Parsed %d values.\n", values);
	return values;
}

static struct architecture arch_init(FILE *fp)
{
	// Parse the architecture file, counting the number of units needed
	//  and then allocating the memory
	struct architecture arch;
	char line[ARCH_LINE_LEN];
	char *field;

	arch.neurons = NULL;
	arch.mem_blocks = NULL;
	arch.routers = NULL;
	arch.timers = NULL;
	arch.axon_inputs = NULL;
	arch.axon_outputs = NULL;
	arch.external_inputs = NULL;

	arch.max_neurons = 0;
	arch.max_mem_blocks = 0;
	arch.max_routers = 0;
	arch.max_timers = 0;
	arch.max_axon_inputs = 0;
	arch.max_axon_outputs = 0;
	arch.max_external_inputs = 0;

	// TODO: refactor - arch_count_units(fp) ?
	while (fgets(line, ARCH_LINE_LEN, fp))
	{
		unsigned int unit_count;
		char block_type = line[0];

		field = strtok(line, " ");
		field = strtok(NULL, " ");
		INFO("debug: field: %s.\n", field);
		unit_count = arch_get_count(field);
		INFO("unit count:%d.\n", unit_count);

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
		case 'n':
			arch.max_neurons += unit_count;
			break;
		case 's':
			arch.max_mem_blocks += unit_count;
			break;
		case 'r':
			arch.max_routers += unit_count;
			break;
		case 't':
			arch.max_timers += unit_count;
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
		default:
			INFO("Warning: unrecognized unit (%c) - skipping.\n",
								block_type);
			break;
		}
	}

	INFO("Parsed %d neurons.\n", arch.max_neurons);

	// Reset file pointer to start of file
	fseek(fp, 0, SEEK_SET);

	// Based on the number of different units, allocate enough memory to
	//  simulate this design
	INFO("Allocating memory for %d neurons.\n", arch.max_neurons);
	arch.neurons = (struct neuron *)
		malloc(arch.max_neurons * sizeof(struct neuron));
	arch.mem_blocks = (struct synapse_mem *)
		malloc(arch.max_mem_blocks * sizeof(struct synapse_mem));
	arch.routers = (struct router *)
		malloc(arch.max_routers * sizeof(struct router));
	arch.timers = (struct timer *)
		malloc(arch.max_timers * sizeof(struct timer));
	arch.axon_inputs = (struct axon_input *)
		malloc(arch.max_axon_inputs * sizeof(struct axon_input));
	arch.axon_outputs = (struct axon_output *)
		malloc(arch.max_axon_outputs * sizeof(struct axon_output));

	if (arch.neurons == NULL || arch.mem_blocks == NULL ||
		arch.routers == NULL || arch.timers == NULL ||
		arch.axon_inputs == NULL || arch.axon_outputs == NULL)
	{
		INFO("Error: failed to allocate neuron.\n");
		exit(1);
	}

	for (int i = 0; i < arch.max_neurons; i++)
	{
		struct neuron *n = &(arch.neurons[i]);
		n->id = i;
		n->synapses = NULL;
		n->axon_in = NULL;
		n->axon_out = NULL;

		n->fired = 0;
		n->active = 0;
	}

	// Allocate memory for the synapse blocks
	for (int i = 0; i < arch.max_mem_blocks; i++)
	{
		struct synapse_mem *mem = &(arch.mem_blocks[i]);

		mem->id = i;
		INFO("Allocating synapse memory for block %d.\n", mem->id);
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

	arch.initialized = 1;

	return arch;
}

void arch_free(struct architecture *arch)
{
	// Free the architecture structure, deallocate all memory and mark
	//  the structure as uninitialized
	arch->initialized = 0;

	// Free any data inside the neuron
	for (int i = 0; i < arch->max_neurons; i++)
	{
		struct neuron *n = &(arch->neurons[i]);

		free(n->synapses);
		n->synapses = NULL;
	}

	// Free all memory
	free(arch->neurons);
	free(arch->mem_blocks);
	free(arch->routers);
	free(arch->timers);
	free(arch->axon_inputs);
	free(arch->axon_outputs);
	free(arch->external_inputs);

	// Reset all pointers and counters
	arch->neurons = NULL;
	arch->mem_blocks = NULL;
	arch->routers = NULL;
	arch->timers = NULL;
	arch->axon_inputs = NULL;
	arch->axon_outputs = NULL;
	arch->external_inputs = NULL;

	arch->max_neurons = 0;
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
	struct synapse_mem *mem_block;
	int ret, list_len, mem_id;
	char block_type;

	assert(fields && fields[0]);
	block_type = fields[0][0];
	assert(block_type == 'n');

	ret = sscanf(fields[2], "%u", &mem_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse value (%s).\n", fields[2]);
		exit(1);
	}
	assert(mem_id < arch->max_mem_blocks);
	mem_block = &(arch->mem_blocks[mem_id]);

	// Copy the neuron a number of times
	list_len = arch_parse_list(fields[1], &(field_list[0]));
	for (int i = 0; i < list_len; i++)
	{
		struct range *r = &(field_list[i]);

		INFO("min:%u max%u\n", r->min, r->max);
		for (unsigned int id = r->min; id <= r->max; id++)
		{
			struct neuron *n;

			assert(id <= arch->max_neurons);
			n = &(arch->neurons[id]);
			// Copy any common compartment variables
			n->mem_block = mem_block;

			TRACE("Neuron parsed n:%d\n", id);
		}
	}
}

/*
static void arch_parse_external_inputs(struct architecture *arch, char *line)
{
	
}

void arch_parse_synapse(struct architecture *arch, char *line)
{

}

void arch_parse_router(struct architecture *arch, char *line)
{

}

void arch_parse_timer(struct architecture *arch, char *line)
{

}
*/