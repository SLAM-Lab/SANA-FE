// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// description.c - Parse architecture description.
//
// An alternative, native architecture parser that does not generate an
//  intermediate format. Python is preferred if possible (more robust
//  and more flexible YAML parsing). This only parses the arch description.
//  I haven't supported the SNN description format.
//
// Note: this code is not used currently, but was saved in case it might be.
// The Python way of generating an intermediate format is preferred. This has
//  been kept in case there is a need to build the simulator with no Python
// dependencies. It will probably need rework to be integrated again

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "description.hpp"
#include "arch.hpp"
#include ".parse.hpp"
#include "print.hpp"

static struct description_block ARCH_BLOCKS[ARCH_MAX_BLOCK_COUNT];
static int ARCH_BLOCK_COUNT;

int arch_parse_topology(char *str, void *topology)
{
	int *val = (int *) topology;

	*val = 0;
	return 0;
}

int arch_build_arch(struct architecture *arch,
			struct description_block *arch_description)
{
	arch_description->arch = arch;
	if (arch_description->type == ARCH_DESCRIPTION_TOP)
	{
		arch_build_arch(arch, arch_description->child);
	}
	else if (arch_description->type == ARCH_DESCRIPTION_ARCH)
	{
		arch_build_arch(arch, arch_description->child);
		// Build the NoC interconnect last after building tiles etc
		arch_build_noc(arch_description);
	}
	else
	{
		for (int i = 0; i < arch_description->instances; i++)
		{
			// Build the current block and then all of its
			//  siblings (same hierarchy level) and then all of its
			//  children (lower hierarchy level)
			arch_build_block(arch_description);
			if (arch_description->next_sibling != NULL)
			{
				arch_description->next_sibling->t =
							arch_description->t;
				arch_description->next_sibling->c =
							arch_description->c;
				arch_build_arch(arch,
						arch_description->next_sibling);
			}
			if (arch_description->child != NULL)
			{
				arch_description->child->t =
							arch_description->t;
				arch_description->child->c =
							arch_description->c;
				arch_build_arch(arch, arch_description->child);
			}
		}
	}

	return 0;
}

int arch_parse_file(FILE *fp,struct architecture *arch,
			struct description_block *arch_description)
{
	char line_buf[ARCH_MAX_DESCRIPTION_LINE+1];

	assert(fp != NULL);
	memset(line_buf, 0, sizeof(line_buf));
	// Use a dummy block to track the entire design
	arch_description->child = NULL;
	arch_description->parent = NULL;
	arch_description->instances = 1;
	arch_description->indent = -1;
	arch_description->name[0] = '\0';
	arch_description->type = ARCH_DESCRIPTION_TOP;
	arch_description->attribute_count = 0;\

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
			return arch_build_arch(arch, arch_description);
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
		TRACE("Read line:%s\n", line_buf);

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
			TRACE("Ignoring line:%s", line_buf);
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
	sibling = NULL;
	head = NULL;
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
	char line_buf[ARCH_MAX_DESCRIPTION_LINE], *str;
	int block_finished, indent;
	struct description_block *block;

	// Initialize a new block
	block = &(ARCH_BLOCKS[ARCH_BLOCK_COUNT]);
	ARCH_BLOCK_COUNT++;
	block->parent = parent;
	block->next_sibling = sibling;
	block->type = block_type;
	block->instances = 1;
	block->attribute_count = 0;

	block->arch = NULL;
	block->t = NULL;
	block->c = NULL;
	block->indent = arch_get_description_line(fp, line_buf);

	str = line_buf;
	if (str[0] == '-')
	{
		// Skip any list characters when parsing this block
		str++;
		block->indent++;
		while (*str == ' ')
		{
			str++;
			block->indent++;
		}
	}

	indent = block->indent;
	block_finished = (indent <= parent->indent);
	TRACE("Created block:%d indent:%d, parent->indent:%d block_finished:%d\n",
		ARCH_BLOCK_COUNT-1, block->indent, parent->indent,
		block_finished);
	while (!block_finished)
	{
		TRACE("Parsing line:%s, field:%s (indent:%d parent->indent:%d,"
							"finished:%d)\n",
			line_buf, str, block->indent, parent->indent,
			block_finished);


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
			block->child = arch_parse_block(fp,
				ARCH_DESCRIPTION_AXON_IN, block, block->child);
		}
		else if (arch_is_field("synapse:", str))
		{
			block->child = arch_parse_block(fp,
				ARCH_DESCRIPTION_SYNAPSE, block, block->child);
		}
		else if (arch_is_field("dendrite:", str))
		{
			block->child = arch_parse_block(fp,
				ARCH_DESCRIPTION_DENDRITE, block, block->child);
		}
		else if (arch_is_field("soma:", str))
		{
			block->child = arch_parse_block(fp,
				ARCH_DESCRIPTION_SOMA, block, block->child);
		}
		else if (arch_is_field("axon_out:", str))
		{
			block->child = arch_parse_block(fp,
				ARCH_DESCRIPTION_AXON_OUT, block, block->child);
		}
		indent = arch_get_description_line(fp, line_buf);
		block_finished = (indent <= parent->indent);
		str = line_buf;
		TRACE("Parsed line:%s (indent:%d, finished:%d)\n", str,
							indent, block_finished);
	}

	// Put the line back
	int line_len = strnlen(line_buf, ARCH_MAX_DESCRIPTION_LINE) + indent;
	fseek(fp, 0-line_len, SEEK_CUR);

	TRACE("Finished parsing block %s\n", block->name);
	return block;
}

int arch_parse_attributes(FILE *fp, struct description_block *block)
{
	char line_buf[ARCH_MAX_DESCRIPTION_LINE];
	int indent, attributes_finished, line_len;
	memset(line_buf, 0, sizeof(line_buf));

	INFO("Parsing attributes.\n");
	indent = arch_get_description_line(fp, line_buf);
	attributes_finished = (indent <= block->indent);
	INFO("Line: %s (indent:%d block->indent:%d)\n",
		line_buf, indent, block->indent);
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
		strncpy(attr->key, line_buf, MAX_FIELD_LEN);
		// Copy, stripping trailing whitespace
		while ((*curr != '\0') && isspace(*curr)) curr++;
		for (int i=0; i < MAX_FIELD_LEN; i++)
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
	char name_field[MAX_FIELD_LEN];
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
				TRACE("Block has %d instances,start:%d,end:%d\n",
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
	for (int i=0; i < MAX_FIELD_LEN; i++)
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

	TRACE("Parsed name:'%s' (instances:%d).\n", dest, instances);
	return instances;
}

int arch_parse_field_int(char *str, int *val)
{
	char val_str[MAX_FIELD_LEN];

	arch_parse_field(str, val_str);
	*val = strtol(str, NULL, 10);

	// TODO: add error checking
	return 1;
}

int arch_parse_field(char *str, char *field)
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

	return 1;
}


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

/*** Build routines ***/

int arch_build_block(struct description_block *block)
{
	int tid, cid;
	// Create the hw block
	switch (block->type)
	{
	case ARCH_DESCRIPTION_TILE:
		tid = arch_build_tile_block(block);
		block->t = &(block->arch->tiles[tid]);
		break;
	case ARCH_DESCRIPTION_CORE:
		cid = arch_build_core_block(block);
		block->c = &(block->t->cores[cid]);
		break;
	case ARCH_DESCRIPTION_AXON_IN:
		arch_build_axon_in_block(block);
		break;
	case ARCH_DESCRIPTION_SYNAPSE:
		arch_build_synapse_block(block);
		break;
	case ARCH_DESCRIPTION_DENDRITE:
		arch_build_dendrite_block(block);
		break;
	case ARCH_DESCRIPTION_SOMA:
		arch_build_soma_block(block);
		break;
	case ARCH_DESCRIPTION_AXON_OUT:
		arch_build_axon_out_block(block);
		break;
	default:
		break;
	}
	return 1;
}

int arch_build_tile_block(struct attributes *attributes)
{
	int blocking;
	double energy_broadcast, latency_broadcast;
	double energy_east_west_hop, latency_east_west_hop;
	double energy_north_south_hop, latency_north_south_hop;
	double energy_spike_within_tile, latency_spike_within_tile;

	blocking = 0;
	energy_broadcast = 0.0;
	latency_broadcast = 0.0;
	energy_east_west_hop = 0.0;
	latency_east_west_hop = 0.0;
	energy_north_south_hop = 0.0;
	latency_north_south_hop = 0.0;
	energy_spike_within_tile = 0.0;
	latency_spike_within_tile = 0.0;

	for (int i = 0; i < block->attribute_count; i++)
	{
		struct attributes *attr = &(block->attributes[i]);

		if (strncmp("blocking", attr->key,
						MAX_FIELD_LEN) == 0)
		{
			blocking = (strncmp(attr->value_str,
					"True", MAX_FIELD_LEN) == 0);
		}
		else if (strncmp("energy_broadcast", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &energy_broadcast);
		}
		else if (strncmp("latency_broadcast", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &latency_broadcast);
		}
		else if (strncmp("energy_east_west", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &energy_east_west_hop);
		}
		else if (strncmp("latency_east_west", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &latency_east_west_hop);
		}
		else if (strncmp("energy_north_south", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &energy_north_south_hop);
		}
		else if (strncmp("latency_north_south", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &latency_north_south_hop);
		}
		else if (strncmp("energy_spike_within_tile", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &energy_spike_within_tile);
		}
		else if (strncmp("latency_spike_within_tile", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &latency_spike_within_tile);
		}
	}

	return arch_create_tile(block->arch, blocking,
			energy_spike_within_tile, energy_east_west_hop,
			energy_north_south_hop, latency_spike_within_tile,
			latency_east_west_hop, latency_north_south_hop);
}

int arch_build_core_block(struct description_block *block)
{
	int blocking = 0;

	for (int i = 0; i < block->attribute_count; i++)
	{
		struct attributes *attr = &(block->attributes[i]);

		if (strncmp("blocking", attr->key, MAX_FIELD_LEN) == 0)
		{
			blocking = (strncmp(attr->value_str, "True",
						MAX_FIELD_LEN) == 0);
		}
	}

	return arch_create_core(block->arch, block->t, blocking);
}

int arch_build_axon_in_block(struct description_block *block)
{
	arch_create_axon_in(block->arch, block->c);
	return 0;
}

int arch_build_synapse_block(struct description_block *block)
{
	int weight_bits, word_bits;
	double energy_memory_read, latency_memory_read;
	double energy_spike_op, latency_spike_op;

	weight_bits = 8;
	word_bits = 64;
	energy_memory_read = 0.0;
	latency_memory_read = 0.0;
	energy_spike_op = 0.0;
	latency_spike_op = 0.0;

	for (int i = 0; i < block->attribute_count; i++)
	{
		struct attributes *attr = &(block->attributes[i]);

		if (strncmp("weight_bits", attr->key, MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%d", &weight_bits);
		}
		if (strncmp("word_bits", attr->key, MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%d", &word_bits);
		}
		else if (strncmp("energy_memory", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &energy_memory_read);
		}
		else if (strncmp("latency_memory", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &latency_memory_read);
		}
		else if (strncmp("energy_spike", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &energy_spike_op);
		}
		else if (strncmp("latency_spike", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &latency_spike_op);
		}
	}

	arch_create_synapse(block->arch, block->c, weight_bits,
		word_bits, energy_spike_op, latency_spike_op,
		energy_memory_read, latency_memory_read);
	return 0;
}

int arch_build_dendrite_block(struct description_block *block)
{
	//arch_create_dendrite(block->arch, block->c);
	return 0;
}

int arch_build_soma_block(struct attributes *attr)
{
	int model;
	double energy_active_neuron_update, latency_active_neuron_update;
	double energy_inactive_neuron_update, latency_inactive_neuron_update;
	double energy_spiking, latency_spiking;

	model = NEURON_LIF;
	energy_active_neuron_update = 0.0;
	latency_active_neuron_update = 0.0;
	energy_inactive_neuron_update = 0.0;
	latency_inactive_neuron_update = 0.0;
	energy_spiking = 0.0;
	latency_spiking = 0.0;

	for (int i = 0; i < block->attribute_count; i++)
	{
		struct attributes *attr = &(block->attributes[i]);

		if (strncmp("model", attr->key, MAX_FIELD_LEN) == 0)
		{
			model = arch_parse_neuron_model(attr->value_str);
		}
		else if (strncmp("energy_active", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf",
						&energy_active_neuron_update);
		}
		else if (strncmp("latency_active", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf",
					&latency_active_neuron_update);
		}
		else if (strncmp("energy_inactive", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf",
					&energy_inactive_neuron_update);
		}
		else if (strncmp("latency_inactive", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf",
					&latency_inactive_neuron_update);
		}
		else if (strncmp("energy_spiking", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &energy_spiking);
		}
		else if (strncmp("latency_spiking", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &latency_spiking);
		}
	}

	arch_create_soma(block->arch, block->c, model,
		energy_active_neuron_update, latency_active_neuron_update,
		energy_inactive_neuron_update, latency_inactive_neuron_update,
		energy_spiking, latency_spiking);
	return 0;
}

int arch_build_axon_out_block(struct description_block *block)
{
	double access_energy, access_latency;

	access_energy = 0.0;
	access_latency = 0.0;

	for (int i = 0; i < block->attribute_count; i++)
	{
		struct attributes *attr = &(block->attributes[i]);

		if (strncmp("energy", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &access_energy);
		}
		else if (strncmp("latency", attr->key,
			MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%lf", &access_latency);
		}
	}

	arch_create_axon_out(block->arch, block->c,
				access_energy, access_latency);
	return 0;
}

int arch_build_noc(struct description_block *block)
{
	assert(block->type == ARCH_DESCRIPTION_ARCH);
	int dimensions, width, height;

	// Default values
	dimensions = 2;
	width = -1;
	height = -1;

	if (block->arch->tile_count <= 0)
	{
		// The NoC interconnect is built after tiles are all defined
		//  This is because we link the tiles together in the NoC
		//  mesh
		INFO("Error: NoC must be built after tiles defined.\n");
		exit(1);
	}

	for (int i = 0; i < block->attribute_count; i++)
	{
		struct attributes *attr = &(block->attributes[i]);

		if (strncmp("dimensions", attr->key, MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%d", &dimensions);
		}
		else if (strncmp("width", attr->key, MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%d", &width);
		}
		else if (strncmp("height", attr->key, MAX_FIELD_LEN) == 0)
		{
			sscanf(attr->value_str, "%d", &height);
		}
	}
	return arch_create_noc(block->arch, width, height);
}
