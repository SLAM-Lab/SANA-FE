#include <assert.h>
#include <stdlib.h>

#include "description.h"
#include "arch.h"
#include "print.h"
#include "command.h"

int description_parse_file(
	FILE *fp, struct network *net, struct architecture *arch)
{
	char *line;
	char(*fields)[MAX_FIELD_LEN];
	int ret;

	assert((arch != NULL) || (net != NULL));

	fields = (char(*)[MAX_FIELD_LEN]) malloc(
		sizeof(char[MAX_FIELD_LEN]) * MAX_FIELDS);
	line = (char *) malloc(sizeof(char) * MAX_LINE);
	ret = RET_OK;

	if ((line == NULL) || (fields == NULL))
	{
		INFO("Error: Couldn't allocate memory for text input.\n");
		ret = RET_FAIL;
	}
	else
	{
		while (fgets(line, MAX_LINE, fp) && (ret != RET_FAIL))
		{
			ret = description_read_line(line, fields, net, arch);
		}
	}

	free(line);
	free(fields);

	return ret;
}

int description_read_line(char *line, char fields[][MAX_FIELD_LEN],
	struct network *net, struct architecture *arch)
{
	char *token;
	int field_count, last_char_idx;

	assert((net != NULL) || (arch != NULL));
	// Zero initialize all fields, each entry is variable length
	//  The neuron fields are fixed, but this is follow by a
	//  variable number of connections (0-MAX_FANOUT).  The synaptic
	//  fields are null delimited
	for (int i = 0; i < MAX_FIELDS; i++)
	{
		fields[i][0] = '\0';
	}

	field_count = 0;
	last_char_idx = strlen(line) - 1;
	assert(line[last_char_idx] == '\n');
	line[last_char_idx] = '\0';
	token = strtok(line, TOKEN_SEPERATORS);
	while (token != NULL)
	{
		strncpy(fields[field_count], token, MAX_FIELD_LEN);
		token = strtok(NULL, TOKEN_SEPERATORS);
		field_count++;
	}

	if (arch != NULL)
	{
		if (net != NULL)
		{
			assert(arch->is_init);
			// The network file description requires an initialized arch
			//  and a pointer to the network struct
			return description_read_network_entry(
				fields, field_count, arch, net);
		}
		return description_read_arch_entry(fields, field_count, arch);
	}
	else
	{
		INFO("File must be either SNN or arch description.\n");
		exit(1);
	}
}

int description_read_arch_entry(char fields[][MAX_FIELD_LEN],
	const int field_count, struct architecture *arch)
{
	struct attributes attributes[128];
	char name[MAX_FIELD_LEN];
	int attribute_count = 0;
	struct tile *t;
	struct core *c;
	int ret, tile_id, core_offset, first_field;
	char entry_type;

	entry_type = fields[0][0];
	// Sanity check input
	if ((entry_type == '\0') || (entry_type == '\n') || (entry_type == '#'))
	{
		TRACE1("Warning: No entry, skipping\n");
		return RET_OK;
	}

	t = NULL;
	c = NULL;
	first_field = 1;
	if (entry_type != '@')
	{
		strncpy(name, fields[1], MAX_FIELD_LEN);
		first_field++;
	}
	if (entry_type != '@' && entry_type != 't')
	{
		ret = sscanf(fields[2], "%d", &tile_id);
		if (ret < 1)
		{
			INFO("Error: Couldn't parse tile ID (%s)\n", fields[2]);
			exit(1);
		}
		t = &(arch->tiles[tile_id]);
		first_field++;
	}
	if ((entry_type != '@') && (entry_type != 't') && (entry_type != 'c'))
	{
		ret = sscanf(fields[3], "%d", &core_offset);
		if (ret < 1)
		{
			INFO("Error: Couldn't parse core ID (%s)\n", fields[3]);
			exit(1);
		}

		assert(t != NULL);
		c = &(t->cores[core_offset]);
		first_field++;
	}

	for (int i = first_field; i < field_count; i++)
	{
		char *key, *value_str;
		struct attributes *attr = &(attributes[attribute_count]);

		TRACE2("Parsing field:%s\n", fields[i]);
		key = strtok(fields[i], "=");
		value_str = strtok(NULL, "=");
		if ((key == NULL) || (value_str == NULL))
		{
			INFO("Invalid attribute: %s\n", fields[i]);
			exit(1);
		}
		strncpy(attr->key, key, MAX_FIELD_LEN);
		strncpy(attr->value_str, value_str, MAX_FIELD_LEN);
		TRACE2("Parsed attribute: %s:%s\n", attr->key, attr->value_str);
		attribute_count++;
	}

	// Process the command and create the unit
	switch (entry_type)
	{
	case '@':
		ret = arch_create_noc(arch, attributes, attribute_count);
		break;
	case 't':
		ret = arch_create_tile(arch, attributes, attribute_count);
		break;
	case 'c':
		ret = arch_create_core(arch, t, attributes, attribute_count);
		break;
	case 'i':
		arch_create_axon_in(c, name, attributes, attribute_count);
		ret = RET_OK;
		break;
	case 's':
		arch_create_synapse(c, name, attributes, attribute_count);
		ret = RET_OK;
		break;
	case 'd':
		// TODO: support dendritic ops
		ret = RET_OK;
		break;
	case '+':
		arch_create_soma(c, name, attributes, attribute_count);
		ret = RET_OK;
		break;
	case 'o':
		arch_create_axon_out(c, attributes, attribute_count);
		ret = RET_OK;
		break;
	default:
		INFO("Warning: unrecognized unit (%c) - skipping.\n",
			entry_type);
		ret = RET_OK;
		break;
	}

	return ret;
}

int description_read_network_entry(char fields[][MAX_FIELD_LEN],
	const int field_count, struct architecture *arch, struct network *net)
{
	struct attributes attributes[DESCRIPTION_MAX_ATTRIBUTES];
	int attribute_count = 0;
	struct neuron_group *group, *dest_group;
	struct neuron *n, *dest;
	struct tile *t;
	struct core *c;
	int ret, neuron_group_id, neuron_id, neuron_count;
	int dest_group_id, dest_neuron_id, tile_id, core_offset;
	char entry_type;

	entry_type = fields[0][0];
	// Sanity check input
	if ((entry_type == '\0') || (entry_type == '\n') || (entry_type == '#'))
	{
		TRACE1("Warning: No entry, skipping\n");
		return RET_OK;
	}

	neuron_count = 0;
	neuron_group_id = -1;
	neuron_id = -1;
	tile_id = -1;
	core_offset = -1;
	dest_group_id = -1;
	dest_neuron_id = -1;
	group = NULL;
	n = NULL;
	c = NULL;
	dest_group = NULL;
	dest = NULL;
	if (entry_type == 'g' || entry_type == 'x')
	{
		ret = sscanf(fields[1], "%d", &neuron_count);
		if (ret < 1)
		{
			INFO("Error: Couldn't parse count (%s)\n", fields[0]);
			exit(1);
		}
	}
	else if (entry_type == '&')
	{
		ret = sscanf(fields[1], "%d.%d@%d.%d", &neuron_group_id,
			&neuron_id, &tile_id, &core_offset);
		if (ret < 4)
		{
			INFO("Error couldn't parse mapping.\n");
			exit(1);
		}
		if (tile_id >= arch->tile_count)
		{
			INFO("Error: Tile (%d) >= tile count (%d)\n", tile_id,
				arch->tile_count);
			exit(1);
		}
		t = &(arch->tiles[tile_id]);

		if (core_offset >= t->core_count)
		{
			INFO("Error: Core (%d) >= core count (%d)\n",
				core_offset, t->core_count);
			exit(1);
		}
		c = &(t->cores[core_offset]);
	}
	else if (entry_type == 'e')
	{
		// Edge on SNN graph (e.g., connection between neurons)
		ret = sscanf(fields[1], "%d.%d->%d.%d", &neuron_group_id,
			&neuron_id, &dest_group_id, &dest_neuron_id);
		if (ret < 4)
		{
			INFO("Error couldn't parse connection / edge.\n");
			exit(1);
		}
		if (dest_group_id > -1)
		{
			if (dest_group_id >= net->neuron_group_count)
			{
				INFO("Error: Group (%d) >= group count (%d).\n",
					dest_group_id, net->neuron_group_count);
				return RET_FAIL;
			}
			dest_group = &(net->groups[dest_group_id]);

			TRACE3("Parsed neuron gid:%d nid:%d\n", dest_group_id,
				neuron_id);
			if (dest_neuron_id > -1)
			{
				if (dest_neuron_id >= dest_group->neuron_count)
				{
					INFO("Error: Neuron (%d) >= "
					     "group neurons (%d).\n",
						dest_neuron_id,
						dest_group->neuron_count);
					return RET_FAIL;
				}
				dest = &(dest_group->neurons[dest_neuron_id]);
			}
		}
	}
	else // parse neuron or input node
	{
		ret = sscanf(fields[1], "%d.%d", &neuron_group_id, &neuron_id);
		if (ret < 2)
		{
			INFO("Error: Couldn't parse neuron (%s)\n", fields[0]);
			exit(1);
		}
	}

	if (neuron_group_id > -1)
	{
		if (neuron_group_id >= net->neuron_group_count)
		{
			INFO("Error: Group (%d) >= group count (%d).\n",
				neuron_group_id, net->neuron_group_count);
			return RET_FAIL;
		}
		group = &(net->groups[neuron_group_id]);
		TRACE3("Parsed neuron gid:%d nid:%d\n", neuron_group_id,
			neuron_id);
		if (neuron_id > -1)
		{
			if (neuron_id >= group->neuron_count)
			{
				INFO("Error: Neuron (%d) >= "
				     "group neurons (%d).\n",
					neuron_id, group->neuron_count);
				return RET_FAIL;
			}
			n = &(group->neurons[neuron_id]);
		}
	}

	// Parse attributes from fields
	for (int i = 2; i < field_count; i++)
	{
		char *key, *value_str;
		struct attributes *attr = &(attributes[attribute_count]);

		TRACE3("Parsing field:%s\n", fields[i]);
		key = strtok(fields[i], "=");
		value_str = strtok(NULL, "=");
		if ((key == NULL) || (value_str == NULL))
		{
			INFO("Invalid attribute: %s\n", fields[i]);
			exit(1);
		}
		strncpy(attr->key, key, MAX_FIELD_LEN);
		strncpy(attr->value_str, value_str, MAX_FIELD_LEN);
		TRACE3("Parsed attribute: %s:%s\n", attr->key, attr->value_str);
		attribute_count++;
	}

	// Process the entry
	int src_is_input = 0;
	switch (fields[0][0])
	{
	case 'g': // Add neuron group
		ret = network_create_neuron_group(
			net, neuron_count, attributes, attribute_count);
		break;
	case 'n': // Add neuron
		ret = network_create_neuron(n, attributes, attribute_count);
		break;
	case 'e':
		if (src_is_input)
		{
			//ret = network_connect_neurons()
			INFO("not implemented.\n");
			ret = RET_OK;
		}
		else
		{
			struct connection *con;

			assert(n != NULL);
			int edge_id = (n->connection_out_count)++;
			if (edge_id >= n->max_connections_out)
			{
				INFO("Edge (%d) >= neuron connections (%d)\n",
					edge_id, n->max_connections_out);
				exit(1);
			}
			con = &(n->connections_out[edge_id]);
			ret = network_connect_neurons(
				con, n, dest, attributes, attribute_count);
		}
		break;
	case 'x': // Add group of external inputs
		INFO("not implemented\n");
		//ret = network_create_extern_input_group(net, fields, field_count);
		break;
	case '<': // Add single input node
		INFO("not implemented");
		//ret = network_create_extern_input_node(net, fields, field_count);
		break;
	case '&': // Map neuron to hardware
		ret = arch_map_neuron(n, c);
		break;
	}

	return ret;
}

/*
int command_parse_neuron_group(struct network *const net,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	double threshold, reset, reverse_threshold, reverse_reset, leak;
	int ret, neuron_count, reset_mode, reverse_reset_mode;

	if (field_count < 8)
	{
		INFO("Error: Invalid <add neuron group> command; "
							"not enough fields.\n");
		exit(1);
	}
	ret = sscanf(fields[1], "%d", &neuron_count);
	ret += sscanf(fields[2], "%lf", &threshold);
	ret += sscanf(fields[3], "%lf", &reset);
	ret += sscanf(fields[4], "%lf", &reverse_threshold);
	ret += sscanf(fields[5], "%lf", &reverse_reset);
	ret += sscanf(fields[6], "%lf", &leak);
	reset_mode = command_parse_reset_mode(fields[7]);
	reverse_reset_mode = command_parse_reset_mode(fields[8]);

	if (ret < 6)
	{
		INFO("Error: Couldn't parse command.\n");
		return RET_FAIL;
	}
	else
	{
		TRACE("Creating neuron group with %d neurons.\n", neuron_count);
		return network_create_neuron_group(net, neuron_count, threshold,
			reset, reverse_threshold, reverse_reset, leak,
						reset_mode, reverse_reset_mode);
	}
}

int description_parse_neuron(struct network *const net, struct architecture *arch,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	struct neuron_group *group;
	struct neuron *n;
	double bias;
	int neuron_group_id, neuron_id, connection_count;
	int log_spikes, log_potential, force_update, ret;

	TRACE("Parsing neuron.\n");
	if (field_count < NEURON_FIELDS)
	{
		INFO("Error: Invalid <neuron> command; (%d/%d) fields.\n",
						field_count, NEURON_FIELDS);
		return RET_FAIL;
	}
	ret = sscanf(fields[NEURON_GROUP_ID], "%d", &neuron_group_id);
	ret += sscanf(fields[NEURON_ID], "%d", &neuron_id);
	ret += sscanf(fields[NEURON_BIAS], "%lf", &bias);
	ret += sscanf(fields[NEURON_RECORD_SPIKES], "%d", &log_spikes);
	ret += sscanf(fields[NEURON_RECORD_potential], "%d", &log_potential);
	ret += sscanf(fields[NEURON_FORCE_UPDATE], "%d", &force_update);

	if (ret < (NEURON_FIELDS-1))
	{
		INFO("Error: Couldn't parse neuron command.\n");
		return RET_FAIL;
	}

	if (neuron_group_id >= net->neuron_group_count)
	{
		INFO("Error: Group (%d) >= group count (%d).\n",
				neuron_group_id, net->neuron_group_count);
		return RET_FAIL;
	}
	group = &(net->groups[neuron_group_id]);

	connection_count = (field_count - NEURON_FIELDS) / CONNECTION_FIELDS;
	TRACE("Parsed neuron gid:%d nid:%d log s:%d log v:%d force:%d "
		"connections:%d\n", neuron_group_id, neuron_id, log_spikes,
				log_potential, force_update, connection_count);

	if (neuron_id >= group->neuron_count)
	{
		INFO("Error: Neuron (%d) >= group neurons (%d).\n",
					neuron_id, group->neuron_count);
		return RET_FAIL;
	}
	n = &(group->neurons[neuron_id]);
	ret = network_create_neuron(n, bias, log_spikes, log_potential,
				force_update, connection_count);
	if (ret == NETWORK_INVALID_NID)
	{
		return RET_FAIL;
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
		con = &(n->connections_out[i]);

		ret = sscanf(fields[curr_field+CONNECTION_DEST_GID],
							"%d", &dest_gid);
		ret += sscanf(fields[curr_field+CONNECTION_DEST_NID],
							"%d", &dest_nid);
		TRACE("group id:%s\n", fields[curr_field+CONNECTION_DEST_GID]);
		TRACE("neuron id:%s\n", fields[curr_field+CONNECTION_DEST_NID]);
		TRACE("weight field:%s\n", fields[curr_field+CONNECTION_WEIGHT]);
		ret += sscanf(fields[curr_field+CONNECTION_WEIGHT],
							"%lf", &weight);
		if (ret < 3)
		{
			INFO("group id:%s\n",
				fields[curr_field+CONNECTION_DEST_GID]);
			INFO("neuron id:%s\n",
				fields[curr_field+CONNECTION_DEST_NID]);
			INFO("weight field:%s\n",
				fields[curr_field+CONNECTION_WEIGHT]);

			INFO("Error: Couldn't parse synapse.\n");
			return RET_FAIL;
		}

		src = n;
		assert(dest_gid < net->neuron_group_count);
		dest_group = &(net->groups[dest_gid]);
		assert(dest_nid < dest_group->neuron_count);
		dest = &(dest_group->neurons[dest_nid]);

		con->pre_neuron = src;
		con->post_neuron = dest;
		con->weight = weight;
		TRACE("\tAdded con %d.%d->%d.%d (w:%lf)\n",
			con->pre_neuron->group->id, con->pre_neuron->id,
			con->post_neuron->group->id, con->post_neuron->id,
			con->weight);
	}

	return RET_OK;
}


int command_parse_extern_input_group(struct network *const net,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	int input_type, input_count, ret;

	TRACE("Parsing group of inputs.\n");
	if (net->external_inputs != NULL)
	{
		INFO("Error: Inputs already created.\n");
		return RET_FAIL;
	}

	if (field_count < 3)
	{
		INFO("Error: Invalid <input group> command; (%d/3) fields.\n",
								field_count);
		return RET_FAIL;
	}

	ret = sscanf(fields[1], "%d", &input_count);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse input count.\n");
		return RET_FAIL;
	}

	// Parse the type of inputs to use. For now just allow one type of
	//  input. If, in future, we wanted to mix and match, I could move this
	//  stuff to the input node
	if (strstr("event", fields[2]))
	{
		input_type = INPUT_EVENT;
	}
	else if (strstr("poisson", fields[2]) == 0)
	{
		input_type = INPUT_RATE;
	}
	else if (strstr("rate", fields[2]) == 0)
	{
		input_type = INPUT_POISSON;
	}
	else
	{
		INFO("Error: Invalid neuron type (%s).\n", fields[2]);
		return RET_FAIL;
	}

 	return net_create_inputs(net, input_count, input_type);
}

int command_parse_extern_input_node(struct network *const net,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	struct input *in;
	int input_id, connection_count, ret;

	TRACE("Parsing input node.\n");
	if (field_count < 2)
	{
		INFO("Error: Invalid <input> command; (%d/2) fields.\n",
							field_count);
		return RET_FAIL;
	}
	ret = sscanf(fields[1], "%d", &input_id);
	if ((ret < 1) || (input_id < 0) ||
					(input_id >= net->external_input_count))
	{
		INFO("Error: Invalid input ID (%s).\n", fields[1]);
		return RET_FAIL;
	}
	in = &(net->external_inputs[input_id]);

	connection_count = (field_count - 2) / CONNECTION_FIELDS;
 	ret = net_create_input_node(in, connection_count);
	if (ret == NETWORK_INVALID_NID)
	{
		return RET_FAIL;
	}
	TRACE("Parsed input iid:%d connections:%d\n",
					input_id, connection_count);
	if (connection_count <= 0)
	{
		INFO("Warning: input %d with no outgoing connections.\n",
								input_id);
		return RET_OK;
	}

	TRACE("iid:%d creating %d connections.\n",
					in->id, in->post_connection_count);
	// Parse each synapse between neurons, those neurons must already exit
	for (int i = 0; i < connection_count; i++)
	{
		// Parse each connection
		struct neuron_group *dest_group;
		struct neuron *src, *dest;
		struct connection *con;
		double weight;
		int dest_gid, dest_nid;
		const int curr_field = 2 + (i*CONNECTION_FIELDS);
		assert(curr_field < field_count);
		con = &(in->connections[i]);

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
			return RET_FAIL;
		}

		src = NULL; // Is not valid when using an input node to feed
		dest_group = &(net->groups[dest_gid]);
		dest = &(dest_group->neurons[dest_nid]);

		con->pre_neuron = src;
		con->post_neuron = dest;

		// TODO: clean up hacks, quantization effects on Loihi
		con->weight = (int) (weight * 64*127.0);
		con->weight = con->weight / (64*127.0);
		TRACE("\tAdded con in[%d]->%d.%d (w:%lf)\n", in->id,
			con->post_neuron->group->id, con->post_neuron->id,
			con->weight);
	}

	return RET_OK;
}
*/

int description_parse_command(char fields[][MAX_FIELD_LEN],
	const int field_count, struct network *net, struct architecture *arch,
	struct simulation *sim)
{
	int ret;
	char command_type;

	command_type = fields[0][0];
	// Sanity check input
	if ((command_type == '\0') || (command_type == '\n') ||
		(command_type == '#'))
	{
		INFO("Warning: No command, skipping.\n");
		return RET_OK;
	}

	// Process the command based on the type
	switch (command_type)
	{
	case '$': // Input vector (either rates or individual spikes)
		ret = command_parse_input_spikes(net, fields, field_count);
		break;
	case '*': // Step simulator
		ret = command_parse_step_sim(net, arch, sim);
		break;
	default:
		INFO("Warning: unrecognized unit (%c) - skipping.\n",
			command_type);
		ret = RET_OK;
		break;
	}

	return ret;
}
