#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <string_view>
#include <charconv>

#include "arch.hpp"
#include "network.hpp"
#include "description.hpp"
#include "print.hpp"

using namespace sanafe;

int sanafe::description_parse_arch_file(
	std::ifstream &fp, Architecture &arch)
{
	std::vector<std::string_view> fields;
	fields.reserve(32);
	std::string line;
	line.reserve(DEFAULT_LINE_LEN);
	int line_number = 1;

	while (std::getline(fp, line))
	{
		TRACE1("Parsing line: %s\n", line.c_str());
		description_get_fields(fields, line);
#ifdef DEBUG
		for (auto f: fields)
		{
			TRACE1("\tField:%s\n", f.c_str());
		}
#endif
		if (fields.size() > 0)
		{
			description_read_arch_entry(fields, arch, line_number);
		}
		line_number++;
	}
	INFO("File parsed.\n");
	return RET_OK;
}

int sanafe::description_parse_net_file(
	std::ifstream &fp, struct Network &net, Architecture &arch)
{
	std::vector<std::string_view> fields;
	fields.reserve(32);
	std::string line;
	line.reserve(DEFAULT_LINE_LEN);
	int line_number = 1;
	while (std::getline(fp, line))
	{
		TRACE1("Parsing line: %s\n", line.c_str());
		description_get_fields(fields, line);

		TRACE1("%ld fields.\n", fields.size());
#ifdef DEBUG
		for (auto f: fields)
		{
			TRACE1("\tField:%s\n", f.c_str());
		}
#endif
		if (fields.size() > 0)
		{
			description_read_network_entry(
				fields, arch, net, line_number);
		}
		line_number++;
	}

	return RET_OK;
}

void sanafe::description_get_fields(
	std::vector<std::string_view> &fields, const std::string &line)
{
	// Get all the fields from a line of text. Every field is separated by
	//  whitespace and has the format <Attribute>=<value>
	// Returns a vector of field strings
	const char *delim = " \t\r\n";

	fields.clear();
	std::string_view line_buffer(line);
	auto field_start = line_buffer.find_first_not_of(delim);
 	while (field_start != std::string_view::npos)
	{
		auto field_end = line_buffer.find_first_of(delim, field_start);

		const std::string_view new_field = line_buffer.substr(
			field_start, field_end - field_start);
		if (field_end != field_start)
		{
			fields.push_back(new_field);
		}
		field_start = line_buffer.find_first_not_of(delim, field_end);
	}

	return;
}

size_t sanafe::field_to_int(const std::string_view &field)
{
	size_t val = 0;
	auto [ptr, error_code] = std::from_chars(field.data(),
		field.data() + field.size(),
		val);
	if (error_code != std::errc())
	{
		std::string error_str =
			"Error: Couldn't parse integer val for field:" +
			std::string(field);
		throw std::runtime_error(std::string(error_str));
	}

	return val;
}

int sanafe::description_read_arch_entry(
	const std::vector<std::string_view> &fields, Architecture &arch,
	const int line_number)
{
	int ret = RET_OK;

	std::vector<Attribute> attributes;
	std::string name;
	Tile *t;
	Core *c;
	int tile_id, core_offset, first_field;

	assert(fields.size() > 0);
	const char entry_type = fields[0][0];
	// Sanity check input
	if ((entry_type == '\0') || (entry_type == '\n') ||
		(entry_type == '#') || (entry_type == '\r'))
	{
		TRACE1("Warning: No entry, skipping\n");
		return RET_OK;
	}

	t = NULL;
	c = NULL;
	first_field = 1;
	if (entry_type != '@')
	{
		name = fields[1];
		first_field++;
	}
	if (entry_type != '@' && entry_type != 't')
	{
		tile_id = field_to_int(fields[2]);
		t = &(arch.tiles[tile_id]);
		first_field++;
	}
	if ((entry_type != '@') && (entry_type != 't') && (entry_type != 'c'))
	{
		core_offset = field_to_int(fields[3]);
		assert(t != nullptr);
		c = &(t->cores[core_offset]);
		first_field++;
	}

	// Parse attributes from fields
	for (std::vector<std::string>::size_type i = first_field;
		i < fields.size(); i++)
	{
		TRACE1("Parsing field:%s\n", std::string(fields[i]).c_str());

		if ((fields[i].length() < 3))
		{
			INFO("Error: Invalid field: %s\n",
				std::string(fields[i]).c_str());
			continue;
		}

		int pos = fields[i].find_first_of('=');
		std::string key = std::string(fields[i].substr(0, pos));
		std::string value_str = std::string(fields[i].substr(pos+1));

		if ((key.length() == 0) || (value_str.length() == 0))
		{
			INFO("Invalid attribute: %s\n",
				std::string(fields[i]).c_str());
			continue;
		}

		Attribute attr = { key, value_str };
		TRACE1("Parsed attribute: %s:%s\n", attr.key.c_str(),
			attr.value_str.c_str());
		attributes.push_back(attr);
	}

	// Process the command and create the unit
	switch (entry_type)
	{
	case '@':
		ret = arch_create_noc(arch, attributes);
		break;
	case 't':
		ret = arch_create_tile(arch, attributes);
		break;
	case 'c':
		ret = arch_create_core(arch, *t, attributes);
		break;
	case 'i':
		arch_create_axon_in(*c, name.c_str(), attributes);
		ret = RET_OK;
		break;
	case 's':
		arch_create_synapse(*c, name.c_str(), attributes);
		ret = RET_OK;
		break;
	case 'd':
		// TODO: support dendritic ops
		ret = RET_OK;
		break;
	case '+':
		arch_create_soma(*c, name.c_str(), attributes);
		ret = RET_OK;
		break;
	case 'o':
		arch_create_axon_out(*c, attributes);
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

void sanafe::parse_neuron_field(
	const std::string_view &neuron_field,
	std::vector<NeuronGroup>::size_type &group_id,
	std::vector<Neuron>::size_type &neuron_id)
{
	const auto pos = neuron_field.find('.');
	if (pos == std::string_view::npos)
	{
		throw std::runtime_error("Invalid neuron format");
	}

	const auto group_str = neuron_field.substr(0, pos);
	group_id = field_to_int(group_str);
	const auto neuron_str = neuron_field.substr(pos + 1,
		neuron_field.size());
	neuron_id = field_to_int(neuron_str);

	return;
}

void sanafe::parse_core_field(
	const std::string_view &core_field,
	std::vector<NeuronGroup>::size_type &tile_id,
	std::vector<Neuron>::size_type &core_offset)
{
	const auto pos = core_field.find('.');
	if (pos == std::string_view::npos)
	{
		throw std::runtime_error("Invalid neuron format");
	}

	const auto tile_str = core_field.substr(0, pos);
	tile_id = field_to_int(tile_str);
	const auto core_str = core_field.substr(pos + 1,
		core_field.size());
	core_offset = field_to_int(core_str);

	return;
}

void sanafe::parse_edge_field(
	const std::string_view &edge_field,
	std::vector<NeuronGroup>::size_type &group_id,
	std::vector<Neuron>::size_type &neuron_id,
	std::vector<NeuronGroup>::size_type &dest_group_id,
	std::vector<Neuron>::size_type &dest_neuron_id)
{
	const auto pos = edge_field.find("->");
	if (pos == std::string_view::npos)
	{
		throw std::runtime_error("Invalid edge format");
	}
	else if ((pos + 1) >= edge_field.size())
	{
		throw std::runtime_error("Invalid edge format");
	}

	const auto src_neuron_address = edge_field.substr(0, pos);
	parse_neuron_field(src_neuron_address, group_id, neuron_id);

	const auto dest_neuron_address = edge_field.substr(
		pos+2, edge_field.size());
	parse_neuron_field(dest_neuron_address, dest_group_id, dest_neuron_id);

	return;
}

void sanafe::parse_mapping_field(
	const std::string_view &mapping_field,
	std::vector<NeuronGroup>::size_type &group_id,
	std::vector<Neuron>::size_type &neuron_id,
	std::vector<Tile>::size_type &tile_id,
	std::vector<Core>::size_type &core_offset)
{
	const auto pos = mapping_field.find("@");
	if (pos == std::string_view::npos)
	{
		throw std::runtime_error("Invalid mapping format");
	}
	else if (pos >= mapping_field.size())
	{
		throw std::runtime_error("Invalid mapping format");
	}

	const auto neuron_address = mapping_field.substr(0, pos);
	parse_neuron_field(neuron_address, group_id, neuron_id);

	const auto core_address = mapping_field.substr(
		pos+1, mapping_field.size());
	parse_neuron_field(core_address, tile_id, core_offset);

	return;
}

int sanafe::description_read_network_entry(
	const std::vector<std::string_view> &fields, Architecture &arch,
	Network &net, const int line_number)
{
	std::vector<Attribute> attributes;
	attributes.reserve(16);
	NeuronGroup *dest_group;
	Neuron *n, *dest;
	Tile *t;
	Core *c;
	int ret, neuron_count;
	std::vector<Tile>::size_type tile_id;
	std::vector<Core>::size_type core_offset;
	std::vector<NeuronGroup>::size_type neuron_group_id, dest_group_id;
	std::vector<Neuron>::size_type neuron_id, dest_neuron_id;
	bool group_set, neuron_set;

	assert(fields.size() > 0);
	const char entry_type = fields[0][0];
	// Sanity check input
	if ((entry_type == '\0') || (entry_type == '\n') ||
		(entry_type == '#') || (entry_type == '\r'))
	{
		TRACE1("Warning: No entry, skipping line %d\n", line_number);
		return RET_OK;
	}

	if (fields.size() < 2)
	{
		INFO("Error: fields < 2 (%ld)", fields.size());
	}

	neuron_count = 0;
	n = nullptr;
	c = nullptr;
	dest_group = nullptr;
	dest = nullptr;

	group_set = false;
	neuron_set = false;

	if (entry_type == 'g')
	{
		neuron_count = field_to_int(fields[1]);
	}
	else if (entry_type == '&')
	{
		parse_mapping_field(fields[1],
			neuron_group_id, neuron_id, tile_id, core_offset);
		/*
		if (!success)
		{
			std::ostringstream error_str;
			error_str << "Error: Line ";
			error_str << line_number;
			error_str << "Couldn't parse mapping. (Format should ";
			error_str << "be src_gid.src_nid@dst_gid.dst_nic)";
			throw std::runtime_error(error_str.str());
		}
		*/
		if (tile_id >= arch.tiles.size())
		{
			INFO("Error: Line %d: Tile (%lu) >= tile count (%lu)\n",
				line_number, tile_id, arch.tiles.size());
			exit(1);
		}
		t = &(arch.tiles[tile_id]);

		if (core_offset >= t->cores.size())
		{
			INFO("Error: Line %d: Core (%lu) >= core count (%lu)\n",
				line_number, core_offset, t->cores.size());
			exit(1);
		}
		c = &(t->cores[core_offset]);
		group_set = true;
		neuron_set = true;
	}
	else if (entry_type == 'e')
	{
		// Edge on SNN graph (e.g., connection between neurons)
		parse_edge_field(fields[1], neuron_group_id, neuron_id,
			dest_group_id, dest_neuron_id);
		/*
		if (!success)
		{
			INFO("Error: Line %d: Couldn't parse connection / edge.\n",
				line_number);
			exit(1);
		}
		*/
		if (dest_group_id >= net.groups.size())
		{
			INFO("Error: Line %d: Group (%lu) >= group count (%lu).\n",
				line_number, dest_group_id, net.groups.size());
			return RET_FAIL;
		}
		dest_group = &(net.groups[dest_group_id]);

		TRACE1("Parsed neuron gid:%lu nid:%lu\n", dest_group_id,
			neuron_id);
		if (dest_neuron_id >= dest_group->neurons.size())
		{
			INFO("Error: Line %d: Trying to access neuron "
				"(%d.%lu) but group %d only "
				"allocates %lu neurons.\n",
				line_number, dest_group->id, dest_neuron_id,
				dest_group->id,
				dest_group->neurons.size());
			return RET_FAIL;
		}
		dest = &(dest_group->neurons[dest_neuron_id]);
		group_set = true;
		neuron_set = true;
	}
	else if (entry_type == 'n') // parse neuron
	{
		//ret = sscanf(std::string(fields[1]).c_str(), "%lu.%lu", &neuron_group_id,
		//	&neuron_id);
		parse_neuron_field(fields[1], neuron_group_id, neuron_id);
		/*
		if (!success)
		{
			INFO("Error: Line %d: Couldn't parse neuron (%s)\n",
				line_number, std::string(fields[0]).c_str());
			exit(1);
		}
		*/
		group_set = true;
		neuron_set = true;
	}
	else
	{
		INFO("Error: Line %d: Invalid entry type (%s)",
			line_number, std::string(fields[0]).c_str());
	}

	if (group_set)
	{
		if (neuron_group_id >= net.groups.size())
		{
			INFO("Error: Line %d: Group (%lu) >= group count (%lu).\n",
				line_number, neuron_group_id, net.groups.size());
			return RET_FAIL;
		}
		NeuronGroup &group = net.groups[neuron_group_id];
		TRACE1("Parsed neuron gid:%lu nid:%lu\n", neuron_group_id,
			neuron_id);
		if (neuron_set)
		{
			if (neuron_id >= group.neurons.size())
			{
                               INFO("Error: Line %d: Trying to access neuron "
                                       "(%d.%lu) but group %d only "
                                       "allocates %lu neuron(s).\n",
					line_number,
					group.id, neuron_id, group.id,
					group.neurons.size());
				return RET_FAIL;
			}
			n = &(group.neurons[neuron_id]);
		}
	}

	// Parse attributes from fields
	for (size_t i = 2; i < fields.size(); i++)
	{
		TRACE1("Parsing field:%s\n", fields[i].c_str());

		if ((fields[i].length() < 3))
		{
			INFO("Error: Line %d: Invalid field: %s\n",
				line_number, std::string(fields[i]).c_str());
			continue;
		}

		const auto pos = fields[i].find_first_of('=');
		std::string key = std::string(fields[i].substr(0, pos));
		std::string value_str = std::string(fields[i].substr(pos+1));

		if ((key.length() == 0) || (value_str.length() == 0))
		{
			INFO("Error: Line %d: Invalid attribute: %s\n",
				line_number, std::string(fields[i]).c_str());
			continue;
		}

		Attribute attr = { key, value_str };
		TRACE1("Parsed attribute: %s:%s\n", attr.key.c_str(),
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
		ret = network_create_neuron(net, *n, attributes);
		break;
	case 'e':
	{
		assert(n != nullptr);
		// Zero initialize all connections TODO: put in constructors
		n->connections_out.push_back(
			Connection(n->connections_out.size()));
		ret = network_connect_neurons(
			n->connections_out[n->connections_out.size()-1],
			*n, *dest, attributes);
		break;
	}
	case '&': // Map neuron to hardware
		ret = arch_map_neuron(net, *n, *c);
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

int description_parse_neuron(struct network *const net, Architecture *arch,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	struct neuron_group *group;
	Neuron *n;
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
		Neuron *src, *dest;
		Connection *con;
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
	const int field_count, struct network *net, Architecture *arch,
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
