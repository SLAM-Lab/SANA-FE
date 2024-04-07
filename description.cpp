#include <cassert>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>

#include <string>
#include <regex>

#include <vector>
#include <list>

#include "description.hpp"
#include "arch.hpp"
#include "print.hpp"
#include "command.hpp"

int description_parse_file(
	std::fstream &fp, struct architecture &arch)
{
	std::string line;
	int ret = RET_OK;

	line.reserve(DEFAULT_LINE_LEN);
	while (std::getline(fp, line) && (ret == RET_OK))
	{
		INFO("Parsing line: %s\n", line.c_str());
		std::vector<std::string> fields = description_get_fields(line);
		for (auto f: fields)
		{
			INFO("\tField:%s\n", f.c_str());
		}

		description_read_arch_entry(fields, arch);
	}

	return ret;
}

int description_parse_file(
	std::fstream &fp, struct network &net, struct architecture &arch)
{
	std::string line;
	int ret = RET_OK;

	line.reserve(DEFAULT_LINE_LEN);
	while (std::getline(fp, line) && (ret == RET_OK))
	{
		INFO("Parsing line: %s\n", line.c_str());
		std::vector<std::string> fields = description_get_fields(line);
		INFO("%ld fields.\n", fields.size());
		for (auto f: fields)
		{
			INFO("\tField:%s\n", f.c_str());
		}
		description_read_network_entry(fields, arch, net);
	}

	return ret;
}

std::vector<std::string> description_get_fields(const std::string &line)
{
	// Get all the fields from a line of text. Every field is separated by
	//  whitespace and has the format <attribute>=<value>
	// Returns a vector of field strings
	const auto field_delimiters = R"( |\t|\n)";
	const std::regex re(field_delimiters);

	// Find the positions of any delimiters in the line
	std::sregex_token_iterator it { line.begin(), line.end(), re, -1 };
	// Use this list to find all seperate fields
	std::vector<std::string> fields { it, {} };
	// Finally, remove any empty fields
	fields.erase(std::remove_if(fields.begin(), fields.end(),
		[](std::string const &s) { return s.size() == 0; }));

	return fields;
}

int description_read_arch_entry(
	const std::vector<std::string> &fields, struct architecture &arch)
{
	int ret= RET_OK;
	/*
	std::vector<attributes> attributes;
	std::string name;
	struct tile *t;
	struct core *c;
	int ret, tile_id, core_offset, first_field;
	char entry_type;

	auto it = fields.front();
	entry_type = it[0];
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
		name = (*it);
		it = it.next();
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
		it = it.next();
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
		it = it.next();
	}

	while (it != fields.end())
	{
		char *key, *value_str;
		std::string key, value_str;

		TRACE2("Parsing field:%s\n", fields[i]);
		key = strtok(*it, "=");
		value_str = strtok(NULL, "=");
		if ((key == NULL) || (value_str == NULL))
		{
			INFO("Invalid attribute: %s\n", fields[i]);
			exit(1);
		}
		strncpy(attr->key, key, MAX_FIELD_LEN);
		strncpy(attr->value_str, value_str, MAX_FIELD_LEN);
		TRACE2("Parsed attribute: %s:%s\n", attr->key, attr->value_str);
		attributes.push_back(attribute);
		it = it.next();
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
	*/
	return ret;
}

int description_read_network_entry(
	const std::vector<std::string> &fields, struct architecture &arch,
	struct network &net)
{
	std::list<attribute> attributes;
	struct neuron_group *group, *dest_group;
	struct neuron *n, *dest;
	struct tile *t;
	struct core *c;
	int ret, neuron_group_id, neuron_id, neuron_count;
	int dest_group_id, dest_neuron_id, tile_id, core_offset;

	const char entry_type = fields[0][0];
	// Sanity check input
	if ((entry_type == '\0') || (entry_type == '\n') || (entry_type == '#'))
	{
		TRACE1("Warning: No entry, skipping\n");
		return RET_OK;
	}

	if (fields.size() < 2)
	{
		INFO("Error: fields < 2 (%ld)", fields.size());
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
	if (entry_type == 'g')
	{
		std::istringstream ss(fields[1]);
		ss >> neuron_count;
		if (ss.fail())
		{
			INFO("Error: Couldn't parse count (%s)\n",
				fields[1].c_str());
			exit(1);
		}
	}
	else if (entry_type == '&')
	{
		ret = sscanf(fields[1].c_str(), "%d.%d@%d.%d", &neuron_group_id,
			&neuron_id, &tile_id, &core_offset);
		if (ret < 4)
		{
			INFO("Error couldn't parse mapping.\n");
			exit(1);
		}
		if (tile_id >= arch.tile_count)
		{
			INFO("Error: Tile (%d) >= tile count (%d)\n", tile_id,
				arch.tile_count);
			exit(1);
		}
		t = &(arch.tiles[tile_id]);

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
		ret = sscanf(fields[1].c_str(), "%d.%d->%d.%d",
			&neuron_group_id, &neuron_id, &dest_group_id,
			&dest_neuron_id);
		if (ret < 4)
		{
			INFO("Error couldn't parse connection / edge.\n");
			exit(1);
		}
		if (dest_group_id > -1)
		{
			if (dest_group_id >= net.neuron_group_count)
			{
				INFO("Error: Group (%d) >= group count (%d).\n",
					dest_group_id, net.neuron_group_count);
				return RET_FAIL;
			}
			dest_group = &(net.groups[dest_group_id]);

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
		ret = sscanf(fields[1].c_str(), "%d.%d", &neuron_group_id,
			&neuron_id);
		if (ret < 2)
		{
			INFO("Error: Couldn't parse neuron (%s)\n",
				fields[0].c_str());
			exit(1);
		}
	}

	if (neuron_group_id > -1)
	{
		if (neuron_group_id >= net.neuron_group_count)
		{
			INFO("Error: Group (%d) >= group count (%d).\n",
				neuron_group_id, net.neuron_group_count);
			return RET_FAIL;
		}
		group = &(net.groups[neuron_group_id]);
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
	for (std::vector<std::string>::size_type i = 2; i < fields.size(); i++)
	{
		INFO("Parsing field:%s\n", fields[i].c_str());

		if ((fields[i].length() < 3))
		{
			INFO("Error: Invalid field: %s\n", fields[i].c_str());
			continue;
		}

		int pos = fields[i].find_first_of('=');
		std::string key = fields[i].substr(0, pos);
		std::string value_str = fields[i].substr(pos+1);

		if ((key.length() == 0) || (value_str.length() == 0))
		{
			INFO("Invalid attribute: %s\n", fields[i].c_str());
			continue;
		}

		attribute attr = { key, value_str };
		INFO("Parsed attribute: %s:%s\n", attr.key.c_str(),
			attr.value_str.c_str());
		attributes.push_back(attr);
	}

	// Process the entry
	switch (entry_type)
	{

	case 'g': // Add neuron group
		ret = network_create_neuron_group(
			net, neuron_count, attributes);
		break;
	case 'n': // Add neuron
		ret = network_create_neuron(n, attributes);
		break;
	case 'e':
	{
		connection con;
		assert(n != NULL);
		// Zero initialize all connections TODO: put in constructors
		con.id = n->connections_out.size();
		con.current = 0.0;
		con.pre_neuron = NULL;
		con.post_neuron = NULL;
		con.weight = 0.0;
		con.delay = 0.0;
		con.synaptic_current_decay = 0.0;
		con.synapse_hw = NULL;
		n->connections_out.push_back(con);
		ret = network_connect_neurons(
			n->connections_out[n->connections_out.size()-1], n, dest,
			attributes);
		break;
	}
	case '&': // Map neuron to hardware
		ret = arch_map_neuron(n, c);
		break;
	default:
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
*/
