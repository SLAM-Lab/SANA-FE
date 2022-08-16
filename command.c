#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "print.h"
#include "arch.h"
#include "network.h"
#include "command.h"

int command_parse_command(char fields[][MAX_FIELD_LEN],
				const int field_count, struct network *net,
				struct architecture *arch)
{
	int ret;
	char command_type;

	command_type = fields[0][0];
	// Sanity check input
	if ((command_type == '\0') || (command_type == '\n') ||
							(command_type == '#'))
	{
		INFO("Warning: No command, skipping.\n");
		return COMMAND_OK;
	}

	// Process the command based on the type
	switch (command_type)
	{
	case 'g': // Add neuron group
		ret = command_parse_neuron_group(net, fields, field_count);
		break;
	case 'n': // Add neuron
		ret = command_parse_neuron(net, fields, field_count);
		break;
	case '@': // Add network-on-chip (noc)
		ret = command_parse_noc(arch, fields, field_count);
		break;
	case 't':
		ret = command_parse_tile(arch, fields, field_count);
		break;
	case 'c':
		ret = command_parse_core(arch, fields, field_count);
		break;
	case '&': // Map neuron to H/W
		ret = command_map_hardware(net, arch, fields, field_count);
		break;
	case 'i':
		ret = command_parse_axon_input(arch, fields, field_count);
		break;
	case 'o':
		ret = command_parse_axon_output(arch, fields, field_count);
		break;
	/*
	case 'e':
		net_parse_external_input(net, fields);
		break;
	case 's':
		arch_parse_synapse(arch, fields);
		break;

	case '+':
		arch_parse_soma(arch, fields);
		break;


	case '*':
		// TODO: step simulation
		break;
	*/
	default:
		TRACE("Warning: unrecognized unit (%c) - skipping.\n",
							command_type);
		ret = COMMAND_OK;
		break;
	}

	return ret;
}

int command_parse_line(char *line, char fields[][MAX_FIELD_LEN],
			struct network *net, struct architecture *arch)
{
	char *token;
	int field_count;

	assert(net != NULL);
	assert(arch != NULL);
	// Zero initialize all fields, each entry is variable length
	//  The neuron fields are fixed, but this is follow by a
	//  variable number of connections (0-MAX_FANOUT).  The synaptic
	//  fields are null delimited
	for (int i = 0; i < MAX_FIELDS; i++)
	{
		fields[i][0] = '\0';
	}


	field_count = 0;
	token = strtok(line, TOKEN_SEPERATORS);
	while (token != NULL)
	{
		// Only parse the neuron fields, the rest of the fields
		//  contain the synaptic data which we'll read next pass
		strncpy(fields[field_count], token, MAX_FIELD_LEN);
		token = strtok(NULL, TOKEN_SEPERATORS);
		field_count++;
	}

	// TODO: max_fields is specific to the thing we're decoding
	// i.e. move this to neuron / SNN specific files
	/*
	if (field_count < MAX_FIELDS)
	{
		TRACE("nid:%d Number of fields read < %d, "
			"ignoring line.\n", neuron_count,
							NEURON_FIELDS);
		continue;
	}
	*/

	return command_parse_command(fields, field_count, net, arch);
}

int command_read_file(FILE *fp, struct network *net,
						struct architecture *arch)
{
	char *line;
	char (*fields)[MAX_FIELD_LEN];
	int ret;

	assert(net != NULL);
	assert(arch != NULL);
	
	fields = (char (*)[MAX_FIELD_LEN])
			malloc(sizeof(char[MAX_FIELD_LEN]) * MAX_FIELDS);
	line = (char *) malloc(sizeof(char) * MAX_CSV_LINE);
	ret = COMMAND_OK;

	if ((line == NULL) || (fields == NULL))
	{
		INFO("Error: Couldn't allocate memory for text input.\n");
		ret = COMMAND_FAIL;
	}
	else
	{
		while (fgets(line, MAX_CSV_LINE, fp) && (ret != COMMAND_FAIL))
		{
			ret = command_parse_line(line, fields, net, arch);
		}
	}
	
	free(line);
	free(fields);

	return ret;
}

int command_parse_noc(struct architecture *arch, char fields[][MAX_FIELD_LEN],
							const int field_count)
{
	int dimensions, width, height, tile_count, ret;

	// Here we can define the links between tiles, all the stuff like whether
	//  it's a mesh, tree can be done externally. Here we just apply raw
	//  links. Maybe we can even remove NoC links dynamically. So we parse
	//  whether to create (1) or remove (0) a link between two tiles, or
	//  something like this
	INFO("Parsing NoC.\n");
	if (strncmp("mesh", fields[1], strlen("mesh")) != 0)
	{
		INFO("Error: Only 'mesh' NoC type supported (%s).\n",
								fields[1]);
		return COMMAND_FAIL;
	}

	ret = sscanf(fields[2], "%d", &dimensions);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output id (%s).\n", fields[1]);
		exit(1);
	}

	if (dimensions != 2)
	{
		INFO("Error: Only 2 dimensions supported (%d).\n", dimensions);
	}

	ret = sscanf(fields[3], "%d", &width);
	ret += sscanf(fields[4], "%d", &height);
	if (ret < 2)
	{
		INFO("Error: Couldn't parse dimensions (%s).\n", fields[1]);
		exit(1);
	}

	tile_count = width * height;
	if (tile_count != arch->tile_count)
	{
		INFO("Error: width*height (%d) != tile count (%d)\n",
						tile_count, arch->tile_count);
	}

	return arch_create_noc(arch, width, height);
}

int command_parse_tile(struct architecture *const arch,
			char fields[][MAX_FIELD_LEN], const int field_count)
{
	TRACE("Parsing tile.\n");
	return arch_create_tile(arch);
}

int command_parse_core(struct architecture *arch,
			char fields[][MAX_FIELD_LEN],
			const int field_count)
{
	struct tile *t;
	int tile_id;
	int ret;

	TRACE("Parsing core.\n");
	if (field_count < 2)
	{
		INFO("Error: Invalid <add core> command; not enough fields.\n");
		return COMMAND_FAIL;
	}
	ret = sscanf(fields[1], "%d", &tile_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse tile id (%s).\n", fields[1]);
		return COMMAND_FAIL;
	}

	if (tile_id >= arch->tile_count)
	{
		INFO("Error: Invalid tile id: %d.\n", tile_id);
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);

	return arch_create_core(arch, t);
}

int command_parse_neuron_group(struct network *const net,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	double threshold, reset;
	int ret, neuron_count;

	if (field_count < 3)
	{
		INFO("Error: Invalid <add neuron group> command; "
							"not enough fields.\n");
		exit(1);
	}
	ret = sscanf(fields[1], "%d", &neuron_count);
	ret += sscanf(fields[2], "%lf", &threshold);
	ret += sscanf(fields[3], "%lf", &reset);
	
	if (ret < 3)
	{
		INFO("Error: Couldn't parse command.\n");
		return COMMAND_FAIL;
	}
	else
	{
		INFO("Creating neuron group with %d neurons.\n", neuron_count);
		return network_create_neuron_group(net, neuron_count, threshold,
									reset);
	}
}

int command_parse_neuron(struct network *const net,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	struct neuron_group *group;
	struct neuron *n;
	int neuron_group_id, neuron_id, connection_count, log_spikes;
	int log_voltage, force_update, ret;

	TRACE("Parsing neuron.\n");
	if (field_count < NEURON_FIELDS)
	{
		INFO("Error: Invalid <neuron> command; (%d/%d) fields.\n",
						field_count, NEURON_FIELDS);
		return COMMAND_FAIL;
	}
	ret = sscanf(fields[NEURON_GROUP_ID], "%d", &neuron_group_id);
	ret += sscanf(fields[NEURON_ID], "%d", &neuron_id);
	ret += sscanf(fields[NEURON_RECORD_SPIKES], "%d", &log_spikes);
	ret += sscanf(fields[NEURON_RECORD_VOLTAGE], "%d", &log_voltage);
	ret += sscanf(fields[NEURON_FORCE_UPDATE], "%d", &force_update);

	if (ret < (NEURON_FIELDS-1))
	{
		INFO("Error: Couldn't parse neuron command.\n");
		return COMMAND_FAIL;	
	}

	connection_count = (field_count - NEURON_FIELDS) / CONNECTION_FIELDS;
	INFO("Parsed neuron gid:%d nid:%d log s:%d log v:%d force:%d "
		"connections:%d\n", neuron_group_id, neuron_id, log_spikes,
				log_voltage, force_update, connection_count);

	if (neuron_group_id >= net->neuron_group_count)
	{
		INFO("Error: group (%d) > group count (%d)",
				neuron_group_id, net->neuron_group_count);
		return COMMAND_FAIL;
	}
	group = &(net->groups[neuron_group_id]);


	if (neuron_id >= group->neuron_count)
	{
		INFO("Error: neuron (%d) > neuron count (%d)",
					neuron_id, group->neuron_count);
		return COMMAND_FAIL;
	}
	n = &(group->neurons[neuron_id]);
	ret = network_create_neuron(n, log_spikes, log_voltage, force_update,
							connection_count);
	if (ret == NETWORK_INVALID_NID)
	{
		return COMMAND_FAIL;
	}
	
	TRACE("nid:%d creating %d connections.\n", n->id, connection_count);
	// Parse each synapse between neurons, those neurons must already exit
	for (int i = 0; i < connection_count; i++)
	{
		// Parse each connection
		struct neuron_group *dest_group;
		struct neuron *src, *dest;
		struct connection *con;
		double weight;
		int dest_gid, dest_nid;
		const int curr_field = NEURON_FIELDS + (i*CONNECTION_FIELDS);
		assert(curr_field < field_count);
		con = &(n->connections[i]);

		ret = sscanf(fields[curr_field+CONNECTION_DEST_GID],
							"%d", &dest_gid);
		ret += sscanf(fields[curr_field+CONNECTION_DEST_NID],
							"%d", &dest_nid);
		TRACE("weight field:%s\n", fields[curr_field+CONNECTION_WEIGHT]);
		ret += sscanf(fields[curr_field+CONNECTION_WEIGHT],
							"%lf", &weight);
		if (ret < 3)
		{
			INFO("Error: Couldn't parse synapse.\n");
			return COMMAND_FAIL;
		}

		src = n;
		dest_group = &(net->groups[dest_gid]);
		dest = &(dest_group->neurons[dest_nid]);

		con->pre_neuron = src;
		con->post_neuron = dest;
		con->weight = weight;
		TRACE("\tAdded con %d.%d->%d.%d (w:%lf)\n",
			con->pre_neuron->group->id, con->pre_neuron->id,
			con->post_neuron->group->id, con->post_neuron->id,
			con->weight);
	}

	return COMMAND_OK;
}

int command_map_hardware(struct network *const net, struct architecture *arch,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	struct neuron_group *group;
	struct tile *t;
	struct core *c;
	struct hardware_mapping map;
	int neuron_group_id, tile_id, core_id, axon_in_id, synapse_id;
	int dendrite_id, soma_id, axon_out_id, ret;
	
	TRACE("Parsing mapping.\n");
	if (field_count < 9)
	{
		INFO("Error: Invalid <map hw> command; (%d/9) fields.\n",
								field_count);
		return COMMAND_FAIL;
	}
	ret = sscanf(fields[1], "%d", &neuron_group_id);
	ret += sscanf(fields[2], "%d", &tile_id);
	ret += sscanf(fields[3], "%d", &core_id);
	ret += sscanf(fields[4], "%d", &axon_in_id);
	ret += sscanf(fields[5], "%d", &synapse_id);
	ret += sscanf(fields[6], "%d", &dendrite_id);
	ret += sscanf(fields[7], "%d", &soma_id);
	ret += sscanf(fields[8], "%d", &axon_out_id);
	if (ret < 8)
	{
		INFO("Error: Couldn't parse command.\n");
		return COMMAND_FAIL;	
	}
	group = &(net->groups[neuron_group_id]);

	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);
	map.core = c;
	map.axon_in = &(c->axon_in[axon_in_id]);
	map.synapse = &(c->synapse[synapse_id]);
	map.dendrite = &(c->dendrite[dendrite_id]);
	map.soma = &(c->soma[soma_id]);
	map.axon_out = &(c->axon_out[axon_out_id]);

	return network_map_neuron_group(group, map);
}

int command_parse_axon_input(struct architecture *const arch,
						char fields[][MAX_FIELD_LEN],
							const int field_count)
{
	struct tile *t;
	struct core *c;
	int tile_id, core_id, ret;

	ret = sscanf(fields[1], "%d", &tile_id);
	ret += sscanf(fields[2], "%d", &core_id);
	if (ret < 2)
	{
		INFO("Error: Couldn't parse axon input.\n");
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);

	return arch_create_axon_in(arch, c);
}

int command_parse_axon_output(struct architecture *arch,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	struct tile *t;
	struct core *c;
	double spike_energy, spike_time;
	int tile_id, core_id, ret;

	ret = sscanf(fields[1], "%d", &tile_id);
	ret += sscanf(fields[2], "%d", &core_id);
	ret += sscanf(fields[3], "%lf", &spike_energy);
	ret += sscanf(fields[4], "%lf", &spike_time);
	if (ret < 4)
	{
		INFO("Error: Couldn't parse axon input.\n");
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);

	INFO("Parsed axon output: spike energy:%e time:%e\n",
						spike_energy, spike_time);

	return arch_create_axon_out(arch, c, spike_energy, spike_time);
}
