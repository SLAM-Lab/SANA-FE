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
#include <unordered_map>

#include "arch.hpp"
#include "network.hpp"
#include "description.hpp"
#include "print.hpp"

using namespace sanafe;

const int default_line_len = 4096;
const int default_fields = 32;

int sanafe::description_parse_arch_file(
	std::ifstream &fp, Architecture &arch)
{
	std::vector<std::string_view> fields;
	fields.reserve(default_fields);
	std::string line;
	line.reserve(default_line_len);
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
	fields.reserve(default_fields);
	std::string line;
	line.reserve(default_line_len);
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

void sanafe::description_read_arch_entry(
	const std::vector<std::string_view> &fields, Architecture &arch,
	const int line_number)
{
	std::unordered_map<std::string, std::string> attributes;
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
		return;
	}

	t = NULL;
	c = NULL;
	first_field = 1;
	tile_id = -1;
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
	for (size_t i = first_field; i < fields.size(); i++)
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
		attributes.insert({key, value_str});
	}

	// Process the command and create the unit
	switch (entry_type)
	{
	case '@':
		arch.set_noc_attributes(attributes);
		break;
	case 't':
		arch.create_tile(attributes);
		break;
	case 'c':
		arch.create_core(tile_id, attributes);
		break;
	case 'i':
		c->create_axon_in(name, attributes);
		break;
	case 's':
		c->create_synapse(name, attributes);
		break;
	case 'd':
		// TODO: support dendritic ops
		break;
	case '+':
		c->create_soma(name, attributes);
		break;
	case 'o':
		c->create_axon_out(name, attributes);
		break;
	default:
		INFO("Warning: unrecognized unit (%c) - skipping.\n",
			entry_type);
		break;
	}

	return;
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

void sanafe::description_read_network_entry(
	const std::vector<std::string_view> &fields, Architecture &arch,
	Network &net, const int line_number)
{
	std::unordered_map<std::string, std::string> attributes;
	attributes.reserve(16);
	NeuronGroup *group, *dest_group;
	Neuron *n, *dest;
	Tile *t;
	Core *c;
	int neuron_count;
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
		return;
	}

	if (fields.size() < 2)
	{
		INFO("Error: fields < 2 (%ld)", fields.size());
	}

	neuron_count = 0;
	n = nullptr;
	group = nullptr;
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
		if (dest_group_id >= net.groups.size())
		{
			INFO("Error: Line %d: Group (%lu) >= group count (%lu).\n",
				line_number, dest_group_id, net.groups.size());
			throw std::invalid_argument("Invalid group id");
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
			throw std::invalid_argument("Invalid nid");
		}
		dest = &(dest_group->neurons[dest_neuron_id]);
		group_set = true;
		neuron_set = true;
	}
	else if (entry_type == 'n') // parse neuron
	{
		parse_neuron_field(fields[1], neuron_group_id, neuron_id);
		group_set = true;
		neuron_set = true;
	}
	else
	{
		INFO("Error: Line %d: Invalid entry type (%s)",
			line_number, std::string(fields[0]).c_str());
		throw std::invalid_argument("Invalid description entry type");
	}

	if (group_set)
	{
		if (neuron_group_id >= net.groups.size())
		{
			INFO("Error: Line %d: Group (%lu) >= group count (%lu).\n",
				line_number, neuron_group_id, net.groups.size());
			throw std::invalid_argument("Invalid group id");
		}
		group = &(net.groups[neuron_group_id]);
		TRACE1("Parsed neuron gid:%lu nid:%lu\n", neuron_group_id,
			neuron_id);
		if (neuron_set)
		{
			if (neuron_id >= group->neurons.size())
			{
                               INFO("Error: Line %d: Trying to access neuron "
                                       "(%d.%lu) but group %d only "
                                       "allocates %lu neuron(s).\n",
					line_number,
					group->id, neuron_id, group->id,
					group->neurons.size());
				throw std::invalid_argument(
					"Invalid neuron id");

			}
			n = &(group->neurons[neuron_id]);
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

		attributes.insert({key, value_str});
	}

	// Process the entry
	switch (entry_type)
	{

	case 'g': // Add neuron group
		net.create_neuron_group(neuron_count, attributes);
		break;
	case 'n': // Add neuron
		n->set_attributes(attributes);
		break;
	case 'e':
	{
		assert(n != nullptr);
		// Zero initialize all connections
		n->connect_to_neuron(*dest, attributes);
		break;
	}
	case '&': // Map neuron to hardware
		c->map_neuron(*n);
		break;
	default:
		break;
	}

	return;
}
