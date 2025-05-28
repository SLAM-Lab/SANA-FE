// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// netlist.cpp
#include <charconv>
#include <cstddef>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include <ryml.hpp> // NOLINT(misc-include-cleaner)

#include "arch.hpp"
#include "attribute.hpp"
#include "netlist.hpp"
#include "network.hpp"
#include "print.hpp"
#include "yaml_common.hpp"
#include "yaml_snn.hpp"

// Netlist (v1) SNN description format. Supported for back-compatability.
//  This format is useful for extremely large network files, as parsing this
//  simpler format requires less memory than YAML parsing.
constexpr int default_line_len = 4096;
constexpr int default_fields = 32;

sanafe::SpikingNetwork sanafe::netlist_parse_file(
        std::ifstream &fp, Architecture &arch)
{
    SpikingNetwork net("");

    std::vector<std::string_view> fields;
    fields.reserve(default_fields);
    std::string line;
    line.reserve(default_line_len);
    int line_number = 1;
    while (std::getline(fp, line))
    {
        TRACE1(DESCRIPTION, "Parsing line: %s\n", line.c_str());
        netlist_get_fields(fields, line);

        TRACE1(DESCRIPTION, "%ld fields.\n", fields.size());

#if (DEBUG_LEVEL_DESCRIPTION > 0)
        for (auto &f : fields)
        {
            std::cout << "\tField:" << std::string(f) << "\n";
        }
#endif
        if (!fields.empty())
        {
            netlist_read_network_entry(fields, arch, net, line_number);
        }
        line_number++;
    }

    return net;
}

void sanafe::netlist_get_fields(
        std::vector<std::string_view> &fields, const std::string &line)
{
    // Get all the fields from a line of text. Every field is separated by
    //  whitespace and has the format <Attribute>=<value>
    // Returns a vector of field strings
    const char *delim = " \t\r\n";

    fields.clear();
    const std::string_view line_buffer(line);
    auto field_start = line_buffer.find_first_not_of(delim);
    while (field_start != std::string_view::npos)
    {
        auto field_end = line_buffer.find_first_of(delim, field_start);

        const std::string_view new_field =
                line_buffer.substr(field_start, field_end - field_start);
        if (field_end != field_start)
        {
            fields.push_back(new_field);
        }
        field_start = line_buffer.find_first_not_of(delim, field_end);
    }
}

std::pair<std::string, size_t> sanafe::netlist_parse_neuron_field(
        const std::string_view &neuron_field)
{
    const auto pos = neuron_field.find('.');
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Error: Invalid neuron format");
    }

    const std::string group_id(neuron_field.substr(0, pos));
    const auto neuron_str = neuron_field.substr(pos + 1, neuron_field.size());
    const size_t neuron_id = field_to_int(neuron_str);

    return {group_id, neuron_id};
}

std::pair<size_t, size_t> sanafe::netlist_parse_core_field(
        const std::string_view &core_field)
{
    const auto pos = core_field.find('.');
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Error: Invalid neuron format");
    }

    auto tile_str = core_field.substr(0, pos);
    const size_t tile_id = field_to_int(tile_str);
    auto core_str = core_field.substr(pos + 1, core_field.size());
    const size_t core_offset = field_to_int(core_str);

    return {tile_id, core_offset};
}

std::tuple<std::string, size_t, std::string, size_t>
sanafe::netlist_parse_edge_field(const std::string_view &edge_field)
{
    // Edge description entries support two formats, to represent
    //  neuron-neuron connections
    // i.e. Connection: e group.neuron->group.neuron:compartment <attributes>
    //   Note that the destination compartment is optional (default=0)
    // Split the source and destination neuron addresses
    TRACE1(DESCRIPTION, "Parsing edge.\n");
    const auto pos = edge_field.find("->");
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Invalid edge format");
    }

    // Parse the source group, neuron and optional compartment identifiers
    //  from the source neuron substring (before the "->")
    const auto src_neuron_address = edge_field.substr(0, pos);
    const auto [group_id, neuron_id] =
            netlist_parse_neuron_field(src_neuron_address);
    // Parse the destination group, neuron and compartment
    //  identifiers from the substring after "->"
    const auto dest_neuron_address =
            edge_field.substr(pos + 2, edge_field.size());
    const auto [dest_group_id, dest_neuron_id] =
            netlist_parse_neuron_field(dest_neuron_address);

    return {group_id, neuron_id, dest_group_id, dest_neuron_id};
}

std::tuple<std::string, size_t, size_t, size_t>
sanafe::netlist_parse_mapping_field(const std::string_view &mapping_field)
{
    const auto pos = mapping_field.find('@');
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Invalid mapping format");
    }
    if (pos >= mapping_field.size())
    {
        throw std::runtime_error("Invalid mapping format");
    }

    const auto neuron_address = mapping_field.substr(0, pos);
    auto [group_id, neuron_id] = netlist_parse_neuron_field(neuron_address);

    const auto core_address =
            mapping_field.substr(pos + 1, mapping_field.size());
    auto [tile_id, core_offset] = netlist_parse_core_field(core_address);

    return {group_id, neuron_id, tile_id, core_offset};
}

void sanafe::netlist_read_network_entry(
        const std::vector<std::string_view> &fields, Architecture &arch,
        SpikingNetwork &net, const int line_number)
{
    const char entry_type = fields[0][0];
    // Sanity check input
    if ((entry_type == '\0') || (entry_type == '\n') || (entry_type == '#') ||
            (entry_type == '\r'))
    {
        TRACE1(DESCRIPTION, "Warning: No entry, skipping line %d\n",
                line_number);
        return;
    }

    if (fields.size() < 2)
    {
        INFO("Error: fields < 2 (%ld)", fields.size());
    }

    if (entry_type == 'g')
    {
        netlist_read_group(fields, net, line_number);
    }
    else if (entry_type == '&')
    {
        netlist_read_mapping(fields, arch, net, line_number);
    }
    else if (entry_type == 'e')
    {
        netlist_read_edge(fields, net, line_number);
    }
    else if (entry_type == 'n') // parse neuron
    {
        netlist_read_neuron(fields, net, line_number);
    }
    else
    {
        INFO("Error: Line %d: Invalid entry type (%s)", line_number,
                std::string(fields[0]).c_str());
        throw std::invalid_argument("Invalid description entry type");
    }
}

std::map<std::string, sanafe::ModelAttribute> sanafe::netlist_parse_attributes(
        const std::vector<std::string_view> &attribute_fields,
        const int line_number)
{
    std::map<std::string, ModelAttribute> attributes{};
    // Parse attributes from remaining fields
    if (attribute_fields.empty() || attribute_fields.at(0).empty())
    {
        return attributes;
    }

    const std::string_view first_field = attribute_fields.at(0);
    const char first_char = first_field.at(0);
    if (first_char == '[' || first_char == '{')
    {
        return netlist_parse_embedded_json(attribute_fields, line_number);
    }

    for (const auto &field : attribute_fields)
    {
        netlist_parse_attribute_field(field, attributes, line_number);
    }

    return attributes;
}

void sanafe::netlist_parse_attribute_field(const std::string_view &field,
        std::map<std::string, ModelAttribute> &attributes,
        const int line_number)
{
    TRACE1(DESCRIPTION, "Parsing field:%s\n", std::string(field).c_str());
    if ((field.length() < 3))
    {
        INFO("Error: Line %d: Invalid field: %s\n", line_number,
                std::string(field).c_str());
        return;
    }
    const auto pos = field.find_first_of('=');
    if (pos == std::string::npos)
    {
        INFO("Error: Line %d: Missing '=' in field: %s\n", line_number,
                std::string(field).c_str());
        return;
    }
    ModelAttribute attribute;
    const std::string key(field.substr(0, pos));
    const std::string value_str(field.substr(pos + 1));

    attribute.name = key;
    attribute.value = netlist_parse_attribute_value(value_str);

    if ((key.empty()) || (value_str.empty()))
    {
        INFO("Error: Line %d: Invalid attribute: %s\n", line_number,
                std::string(field).c_str());
        return;
    }

    attributes.insert({key, attribute});
}

std::variant<bool, int, double, std::string, std::vector<sanafe::ModelAttribute>>
sanafe::netlist_parse_attribute_value(std::string value_str)
{
    std::stringstream int_ss(value_str);
    std::stringstream float_ss(value_str);
    std::stringstream bool_ss(value_str);

    int decoded_int = 0;
    // Use the extractor operator to attempt to parse each type in turn,
    //  returning at the first successful parsed value. Note that we check
    //  for end of file (eof) to make sure we didn't only partially parse the
    //  value string
    if ((int_ss >> decoded_int) && int_ss.eof())
    {
        TRACE1(DESCRIPTION, "Parsed integer: %d.\n", decoded_int);
        return decoded_int;
    }

    double decoded_double = std::numeric_limits<double>::quiet_NaN();
    if ((float_ss >> decoded_double) && float_ss.eof())
    {
        TRACE1(DESCRIPTION, "Parsed float: %e.\n", decoded_double);
        return decoded_double;
    }

    bool decoded_bool = false;
    if ((bool_ss >> decoded_bool) && bool_ss.eof())
    {
        TRACE1(DESCRIPTION, "Parsed bool: %d.\n", decoded_bool);
        return decoded_bool;
    }

    TRACE1(DESCRIPTION, "Parsed string: %s\n", value_str.c_str());
    return std::move(value_str);
}

char sanafe::netlist_get_closing_char(const char opening_char)
{
    switch (opening_char)
    {
    case '[':
        return ']';
    case '{':
        return '}';
    default:
        INFO("Error: Invalid opening character (%c)\n", opening_char);
        throw std::runtime_error("Invalid opening character\n");
    }
}

size_t sanafe::netlist_embedded_json_end_pos(const char opening_char,
        const std::string &all_fields, const int line_number)
{
    const char closing_char = netlist_get_closing_char(opening_char);
    int nested_level = 0;

    size_t end_pos = 0;
    while (end_pos < all_fields.length())
    {
        const char ch = all_fields[end_pos];
        if (ch == opening_char)
        {
            ++nested_level;
        }
        else if (ch == closing_char)
        {
            --nested_level;
        }

        if (nested_level < 1)
        {
            break;
        }
        ++end_pos;
    }

    if (nested_level > 0)
    {
        INFO("Error: JSON attributes weren't terminated (%s) on line:%d.\n",
                all_fields.c_str(), line_number);
        throw std::invalid_argument("JSON attributes weren't terminated.");
    }
    if (all_fields.length() > (end_pos + 2))
    {
        INFO("Warning: Text found after delimiter, end of line:%d ignored\n",
                line_number);
    }

    return end_pos;
}

std::map<std::string, sanafe::ModelAttribute>
sanafe::netlist_parse_embedded_json(
        const std::vector<std::string_view> &attribute_fields,
        const int line_number)
{
    std::map<std::string, ModelAttribute> model_attributes{};

    // Concatenate fields again, and search for the terminating character
    std::string all_fields;
    for (auto field : attribute_fields)
    {
        all_fields += field;
        all_fields += ' ';
    }
    if (all_fields.empty())
    {
        return model_attributes;
    }

    const char opening_char = all_fields[0];
    const size_t end_pos = netlist_embedded_json_end_pos(
            opening_char, all_fields, line_number);
    all_fields.resize(end_pos + 1);

    // Create YAML parser, ignore linter include warnings on RapidYAML classes
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    const ryml::Parser parser(
            &event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = ryml::parse_in_place(all_fields.data());
    const ryml::ConstNodeRef yaml_node = tree.rootref();
    // NOLINTEND(misc-include-cleaner)

    model_attributes =
            description_parse_model_attributes_yaml(parser, yaml_node);

    return model_attributes;
}

void sanafe::netlist_read_group(const std::vector<std::string_view> &fields,
        sanafe::SpikingNetwork &net, const int line_number)
{
    const size_t neuron_count = field_to_int(fields[1]);
    std::string neuron_group_id = std::to_string(net.groups.size());

    TRACE1(DESCRIPTION, "Parsed neuron gid:%s\n", neuron_group_id.c_str());

    NeuronConfiguration neuron_config{};
    const auto attribute_fields =
            std::vector<std::string_view>(fields.begin() + 2, fields.end());
    auto attributes = netlist_parse_attributes(attribute_fields, line_number);
    neuron_config.model_attributes = attributes;
    if (attributes.find("synapse_hw_name") != attributes.end())
    {
        neuron_config.default_synapse_hw_name =
                static_cast<std::string>(attributes["synapse_hw_name"]);
    }
    if (attributes.find("dendrite_hw_name") != attributes.end())
    {
        neuron_config.dendrite_hw_name =
                static_cast<std::string>(attributes["dendrite_hw_name"]);
    }
    if (attributes.find("soma_hw_name") != attributes.end())
    {
        neuron_config.soma_hw_name =
                static_cast<std::string>(attributes["soma_hw_name"]);
    }
    if (attributes.find("force_update") != attributes.end())
    {
        neuron_config.force_soma_update =
                static_cast<bool>(attributes["force_update"]);
    }
    if (attributes.find("log_spikes") != attributes.end())
    {
        neuron_config.log_spikes = static_cast<bool>(attributes["log_spikes"]);
    }
    if (attributes.find("log_v") != attributes.end())
    {
        neuron_config.log_potential = static_cast<bool>(attributes["log_v"]);
    }

    TRACE1(DESCRIPTION, "Creating neuron group:%s with count:%zu\n",
            neuron_group_id.c_str(), neuron_count);
    net.create_neuron_group(
            std::move(neuron_group_id), neuron_count, neuron_config);
}

void sanafe::netlist_read_neuron(const std::vector<std::string_view> &fields,
        SpikingNetwork &net, const int line_number)
{
    std::string neuron_group_id;
    size_t neuron_id = 0;

    std::tie(neuron_group_id, neuron_id) =
            netlist_parse_neuron_field(fields[1]);
    NeuronGroup &group = net.groups.at(neuron_group_id);
    if (neuron_id >= group.neurons.size())
    {
        INFO("Error: Line %d: Trying to access neuron "
             "(%s.%lu) but group %s only "
             "allocates %lu neuron(s).\n",
                line_number, group.name.c_str(), neuron_id, group.name.c_str(),
                group.neurons.size());
        throw std::invalid_argument("Invalid neuron id");
    }

    NeuronConfiguration neuron_config{};
    const auto attribute_fields =
            std::vector<std::string_view>(fields.begin() + 2, fields.end());
    auto attributes = netlist_parse_attributes(attribute_fields, line_number);
    neuron_config.model_attributes = attributes;
    if (attributes.find("synapse_hw_name") != attributes.end())
    {
        neuron_config.default_synapse_hw_name =
                static_cast<std::string>(attributes["synapse_hw_name"]);
    }
    if (attributes.find("dendrite_hw_name") != attributes.end())
    {
        neuron_config.dendrite_hw_name =
                static_cast<std::string>(attributes["dendrite_hw_name"]);
    }
    if (attributes.find("soma_hw_name") != attributes.end())
    {
        neuron_config.soma_hw_name =
                static_cast<std::string>(attributes["soma_hw_name"]);
    }
    if (attributes.find("force_update") != attributes.end())
    {
        neuron_config.force_soma_update =
                static_cast<bool>(attributes["force_update"]);
    }
    if (attributes.find("log_spikes") != attributes.end())
    {
        neuron_config.log_spikes = static_cast<bool>(attributes["log_spikes"]);
    }
    if (attributes.find("log_v") != attributes.end())
    {
        neuron_config.log_potential = static_cast<bool>(attributes["log_v"]);
    }

    Neuron &neuron = group.neurons.at(neuron_id);
    neuron.set_attributes(neuron_config);
}

void sanafe::netlist_read_edge(const std::vector<std::string_view> &fields,
        sanafe::SpikingNetwork &net, const int line_number)
{
    std::string neuron_group_id;
    std::string dest_group_id{};
    size_t neuron_id = 0;
    size_t dest_neuron_id = 0;

    // Edge on SNN graph (e.g., connection between neurons)
    std::tie(neuron_group_id, neuron_id, dest_group_id, dest_neuron_id) =
            netlist_parse_edge_field(fields[1]);
    NeuronGroup &group = net.groups.at(neuron_group_id);
    if (neuron_id >= group.neurons.size())
    {
        INFO("Error: Line %d: Trying to access neuron "
             "(%s.%lu) but group %s only "
             "allocates %lu neuron(s).\n",
                line_number, group.name.c_str(), neuron_id, group.name.c_str(),
                group.neurons.size());
        throw std::invalid_argument("Invalid neuron id");
    }
    Neuron &neuron = group.neurons.at(neuron_id);

    if (net.groups.find(dest_group_id) == net.groups.end())
    {
        INFO("Error: Line %d: Group (%s) not in groups.\n", line_number,
                dest_group_id.c_str());
        throw std::invalid_argument("Invalid group id");
    }
    NeuronGroup &dest_group = net.groups.at(dest_group_id);

    TRACE1(DESCRIPTION, "Parsed neuron gid:%s nid:%lu\n", dest_group_id.c_str(),
            neuron_id);
    if (dest_neuron_id >= dest_group.neurons.size())
    {
        INFO("Error: Line %d: Trying to access neuron "
             "(%s.%lu) but group %s only "
             "allocates %lu neurons.\n",
                line_number, dest_group.name.c_str(), dest_neuron_id,
                dest_group.name.c_str(), dest_group.neurons.size());
        throw std::invalid_argument("Invalid nid");
    }
    Neuron &dest_neuron = dest_group.neurons[dest_neuron_id];

    const auto attribute_fields =
            std::vector<std::string_view>(fields.begin() + 2, fields.end());
    auto attributes = netlist_parse_attributes(attribute_fields, line_number);

    const size_t idx = neuron.connect_to_neuron(dest_neuron);
    Connection &con = neuron.edges_out[idx];
    con.synapse_attributes = attributes;
    con.dendrite_attributes = attributes;
}

void sanafe::netlist_read_mapping(const std::vector<std::string_view> &fields,
        sanafe::Architecture &arch, sanafe::SpikingNetwork &net,
        const int line_number)
{
    std::string neuron_group_id{};
    size_t neuron_id{0};
    size_t tile_id{0};
    size_t core_offset{0};

    std::tie(neuron_group_id, neuron_id, tile_id, core_offset) =
            netlist_parse_mapping_field(fields[1]);
    if (tile_id >= arch.tiles.size())
    {
        INFO("Error: Line %d: Tile (%lu) >= tile count (%lu)\n", line_number,
                tile_id, arch.tiles.size());
        throw std::runtime_error("Error: Couldn't parse mapping.");
    }
    TileConfiguration &tile = arch.tiles[tile_id];
    auto *tile_ptr = &tile;

    if (core_offset >= tile_ptr->cores.size())
    {
        INFO("Error: Line %d: Core (%lu) >= core count (%lu)\n", line_number,
                core_offset, tile_ptr->cores.size());
        throw std::runtime_error("Error: Couldn't parse mapping.");
    }
    const CoreConfiguration &core = tile_ptr->cores[core_offset];

    NeuronGroup &group = net.groups.at(neuron_group_id);
    if (neuron_id >= group.neurons.size())
    {
        INFO("Error: Line %d: Trying to access neuron "
             "(%s.%lu) but group %s only "
             "allocates %lu neuron(s).\n",
                line_number, group.name.c_str(), neuron_id, group.name.c_str(),
                group.neurons.size());
        throw std::invalid_argument("Invalid neuron id");
    }

    Neuron &neuron = group.neurons.at(neuron_id);
    const auto attribute_fields =
            std::vector<std::string_view>(fields.begin() + 2, fields.end());
    auto attributes = netlist_parse_attributes(attribute_fields, line_number);
    neuron.map_to_core(core);
}

std::string sanafe::netlist_group_to_netlist(const NeuronGroup &group)
{
    // TODO: support attribute specific to certain h/w units
    std::string entry = "g " + std::to_string(group.neurons.size());

    if (group.default_neuron_config.default_synapse_hw_name.has_value() &&
            !group.default_neuron_config.default_synapse_hw_name.value().empty())
    {
        entry += " synapse_hw_name=" +
                group.default_neuron_config.default_synapse_hw_name.value();
    }
    if (group.default_neuron_config.dendrite_hw_name.has_value() &&
            !group.default_neuron_config.dendrite_hw_name.value().empty())
    {
        entry += " dendrite_hw_name=" +
                group.default_neuron_config.dendrite_hw_name.value();
    }
    if (group.default_neuron_config.force_dendrite_update.has_value() &&
            group.default_neuron_config.force_dendrite_update.value())
    {
        entry += " force_dendrite_update=" +
                std::to_string(static_cast<int>(group.default_neuron_config
                                .force_dendrite_update.value()));
    }
    if (group.default_neuron_config.force_soma_update.has_value() &&
            group.default_neuron_config.force_soma_update.value())
    {
        entry += " force_soma_update=" +
                std::to_string(static_cast<int>(
                        group.default_neuron_config.force_soma_update.value()));
    }
    if (group.default_neuron_config.force_synapse_update.has_value() &&
            group.default_neuron_config.force_synapse_update.value())
    {
        entry += " force_synapse_update=" +
                std::to_string(static_cast<int>(group.default_neuron_config
                                .force_synapse_update.value()));
    }
    if (group.default_neuron_config.log_potential.has_value() &&
            group.default_neuron_config.log_potential.value())
    {
        entry += " log_potential=" +
                std::to_string(static_cast<int>(
                        group.default_neuron_config.log_potential.value()));
    }
    if (group.default_neuron_config.log_spikes.has_value() &&
            group.default_neuron_config.log_spikes.value())
    {
        entry += " log_spikes=" +
                std::to_string(static_cast<int>(
                        group.default_neuron_config.log_spikes.value()));
    }
    if (group.default_neuron_config.soma_hw_name.has_value() &&
            !group.default_neuron_config.soma_hw_name.value().empty())
    {
        entry += " soma_hw_name=" +
                group.default_neuron_config.soma_hw_name.value();
    }

    TRACE2(NET, "saving attributes\n");
    const std::map<std::string, ModelAttribute> no_default_attributes{};
    entry += netlist_attributes_to_netlist(
            group.default_neuron_config.model_attributes,
            no_default_attributes);

    return entry;
}

void sanafe::add_string_attribute_if_unique(std::string &entry,
        const std::string &attr_name, const std::string &neuron_value,
        const std::optional<std::string> &group_default)
{
    if (!neuron_value.empty() &&
            is_unique_attribute(group_default, neuron_value))
    {
        entry += " " + attr_name + "=" + neuron_value;
    }
}

void sanafe::add_bool_attribute_if_unique(std::string &entry,
        const std::string &attr_name, bool neuron_value,
        const std::optional<bool> &group_default)
{
    if (neuron_value && is_unique_attribute(group_default, neuron_value))
    {
        entry += " " + attr_name + "=" +
                std::to_string(static_cast<int>(neuron_value));
    }
}

std::string sanafe::netlist_neuron_to_netlist(const Neuron &neuron,
        const SpikingNetwork &net,
        const std::map<std::string, size_t> &group_name_to_id)
{
    // TODO: support attributes specific to certain h/w units
    TRACE1(DESCRIPTION, "Saving neuron nid:%s.%zu to netlist\n",
            neuron.parent_group_name.c_str(), neuron.get_id());
    const NeuronGroup &parent_group = net.groups.at(neuron.parent_group_name);

    const size_t group_id = group_name_to_id.at(parent_group.name);
    const auto &default_config = parent_group.default_neuron_config;
    std::string entry = "n " + std::to_string(group_id) + "." +
            std::to_string(neuron.offset);

    // Add hardware name attributes if they differ from group defaults
    add_string_attribute_if_unique(entry, "soma_hw_name", neuron.soma_hw_name,
            default_config.soma_hw_name);
    add_string_attribute_if_unique(entry, "synapse_hw_name",
            neuron.default_synapse_hw_name,
            default_config.default_synapse_hw_name);
    add_string_attribute_if_unique(entry, "dendrite_hw_name",
            neuron.dendrite_hw_name, default_config.dendrite_hw_name);

    // Add force update flags if they differ from group defaults
    add_bool_attribute_if_unique(entry, "force_synapse_update",
            neuron.force_synapse_update, default_config.force_synapse_update);
    add_bool_attribute_if_unique(entry, "force_dendrite_update",
            neuron.force_dendrite_update, default_config.force_dendrite_update);
    add_bool_attribute_if_unique(entry, "force_soma_update",
            neuron.force_soma_update, default_config.force_soma_update);

    // Add logging flags if they differ from group defaults
    add_bool_attribute_if_unique(
            entry, "log_spikes", neuron.log_spikes, default_config.log_spikes);
    add_bool_attribute_if_unique(entry, "log_potential", neuron.log_potential,
            default_config.log_potential);

    entry += netlist_attributes_to_netlist(neuron.model_attributes,
            parent_group.default_neuron_config.model_attributes);

    return entry;
}

std::string sanafe::netlist_mapping_to_netlist(const Neuron &neuron,
        const std::map<std::string, size_t> &group_name_to_id)
{
    std::string entry{};

    if (neuron.core_address.has_value())
    {
        const CoreAddress &address = neuron.core_address.value();
        const size_t group_id = group_name_to_id.at(neuron.parent_group_name);
        entry = "& " + std::to_string(group_id) + "." +
                std::to_string(neuron.offset) + "@" +
                std::to_string(address.parent_tile_id) + "." +
                std::to_string(address.offset_within_tile);
    }
    else
    {
        TRACE1(NET, "No mapping defined\n");
    }
    return entry;
}

std::string sanafe::netlist_connection_to_netlist(const Connection &con,
        const std::map<std::string, size_t> &group_name_to_id)
{
    // TODO: support attributes specific to only synapse or dendrite h/w
    std::string pre_neuron_offset_str;
    if (con.pre_neuron.neuron_offset.has_value())
    {
        pre_neuron_offset_str = ".";
        pre_neuron_offset_str +=
                std::to_string(con.pre_neuron.neuron_offset.value());
    }
    std::string post_neuron_offset_str;
    if (con.post_neuron.neuron_offset.has_value())
    {
        post_neuron_offset_str = ".";
        post_neuron_offset_str +=
                std::to_string(con.post_neuron.neuron_offset.value());
    }

    const size_t src_group_id = group_name_to_id.at(con.pre_neuron.group_name);
    const size_t dest_group_id =
            group_name_to_id.at(con.post_neuron.group_name);
    std::string entry = "e " + std::to_string(src_group_id) +
            pre_neuron_offset_str + "->" + std::to_string(dest_group_id) +
            post_neuron_offset_str;

    const std::map<std::string, ModelAttribute> no_default_attributes{};
    entry += netlist_attributes_to_netlist(
            con.synapse_attributes, no_default_attributes);

    return entry;
}

std::string sanafe::netlist_attributes_to_netlist(
        const std::map<std::string, sanafe::ModelAttribute> &model_attributes,
        const std::map<std::string, sanafe::ModelAttribute> &default_attributes)
{
    std::string attribute_str{};

    TRACE1(DESCRIPTION, "Parsing attributes\n");
    bool nested_attributes = false;
    for (const auto &[key, attribute] : model_attributes)
    {
        if (std::holds_alternative<std::vector<ModelAttribute>>(
                    attribute.value))
        {
            // One or more attributes is nested
            nested_attributes = true;
            break;
        }
    }

    if (nested_attributes)
    {
        TRACE2(DESCRIPTION, "Parsing nested attributes\n");
        ryml::Tree tree; // NOLINT(misc-include-cleaner)
        ryml::NodeRef root = tree.rootref(); // NOLINT(misc-include-cleaner)
        root |= ryml::MAP; // NOLINT(misc-include-cleaner)
        root |= ryml::FLOW_SL; // NOLINT(misc-include-cleaner)

        yaml_serialize_model_attributes(
                default_attributes, root, model_attributes);
        std::ostringstream ss;
        ss << tree;
        attribute_str = " " + ss.str();
        INFO("attribute str:%s\n", attribute_str.c_str());
    }
    else
    {
        TRACE2(DESCRIPTION, "Parsing attributes using normal format\n");
        for (const auto &attribute : model_attributes)
        {
            attribute_str += " ";
            attribute_str += attribute.first;
            attribute_str += "=";
            attribute_str += attribute.second.print();
        }
    }

    return attribute_str;
}

size_t sanafe::field_to_int(const std::string_view &field)
{
    size_t val = 0;
    auto [ptr, error_code] =
            std::from_chars(field.data(), field.data() + field.size(), val);
    if (error_code != std::errc())
    {
        const std::string error_str =
                "Error: Couldn't parse integer val for field:" +
                std::string(field);
        throw std::runtime_error(error_str);
    }

    return val;
}
