// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric> // For std::accumulate
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <ryml.hpp>
#include <ryml_std.hpp>

#include "arch.hpp"
#include "description.hpp"
#include "network.hpp"
#include "print.hpp"
#include "chip.hpp"

// TODO: create an option where we can directly construct the mapped network on
//  the hardware, without creating an intermediate abstract spiking network
//  This will save time and memory for the larger applications.

sanafe::DescriptionParsingError::DescriptionParsingError(
        const std::string &error, const ryml::Parser &parser,
        const ryml::ConstNodeRef &node)
        : std::invalid_argument(error)
{
    ryml::Location pos = parser.location(node);
    message = "Error: " + error + " (Line " + std::to_string(pos.line + 1) +
            ':' + std::to_string(pos.col + 1) + ").\n";
}

const char *sanafe::DescriptionParsingError::what() const noexcept
{
    return message.c_str();
}

void sanafe::check_key(const ryml::Parser &parser,
        const ryml::ConstNodeRef node, const std::string &key)
{
    if (!node.is_map())
    {
        throw DescriptionParsingError(
                "Node should be a mapping\n. For more info on YAML mappings "
                "refer to the YAML 1.2 specification, 7.4.2 'Flow Mappings' "
                "and 8.2.2 'Block Mappings'",
                parser, node);
    }
    const ryml::ConstNodeRef child = node.find_child(key.c_str());
    if (child.invalid())
    {
        const std::string message = "Value for key '" + key + "' not defined";
        throw DescriptionParsingError(message, parser, node);
    }
}

template <typename T>
T sanafe::description_required_field(const ryml::Parser &parser,
        const ryml::ConstNodeRef node, const std::string &key)
{
    // Wrapper around YAML library for field=map[key], adding more error prints
    if (node.invalid())
    {
        const std::string message = "Invalid node when looking up key: " + key;
        throw std::runtime_error(message);
    }
    const ryml::ConstNodeRef field_node = node.find_child(key.c_str());
    if (field_node.invalid())
    {
        const std::string message = "Key '" + key + "' does not exist";
        throw DescriptionParsingError(message, parser, node);
    }
    if (!field_node.has_val())
    {
        const std::string message = "'" + key + "' value should be a scalar";
        throw DescriptionParsingError(message, parser, field_node);
    }

    T field;
    // Efficiently convert to type T by trying the RapidYAML reader.
    //  If read() fails, it returns false and execution falls through
    if (c4::yml::read(field_node, &field))
    {
        return field; // type T
    }

    const std::string message = "Could not cast field '" + key +
            "' to type: " + description_get_type_string(typeid(field));
    throw DescriptionParsingError(message, parser, field_node);
}

std::string sanafe::description_get_type_string(const std::type_info &type)
{
    if (type == typeid(bool))
    {
        return "bool";
    }
    if (type == typeid(int))
    {
        return "int";
    }
    if (type == typeid(size_t))
    {
        return "size_t";
    }
    if (type == typeid(double))
    {
        return "double";
    }
    if (type == typeid(std::string))
    {
        return "string";
    }

    // Not a scalar type; fall back to default name which may be mangled
    return type.name();
}

void sanafe::description_parse_axon_in_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef axon_in_node, CoreConfiguration &parent_core)
{
    auto name = description_required_field<std::string>(
            parser, axon_in_node, "name");
    const ryml::ConstNodeRef &attributes =
            axon_in_node.find_child("attributes");
    if (attributes.invalid())
    {
        throw DescriptionParsingError(
                "No attributes section defined", parser, axon_in_node);
    }
    const AxonInPowerMetrics in_metrics =
            description_parse_axon_in_attributes_yaml(parser, attributes);
    parent_core.create_axon_in(name, in_metrics);
}

sanafe::AxonInPowerMetrics sanafe::description_parse_axon_in_attributes_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    AxonInPowerMetrics axon_in_metrics;
    axon_in_metrics.energy_message_in = description_required_field<double>(
            parser, attributes, "energy_message_in");
    axon_in_metrics.latency_message_in = description_required_field<double>(
            parser, attributes, "latency_message_in");

    return axon_in_metrics;
}

void sanafe::description_parse_synapse_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef synapse_node, CoreConfiguration &parent_core)
{
    auto name = description_required_field<std::string>(
            parser, synapse_node, "name");
    auto model = description_parse_synapse_attributes_yaml(
            parser, synapse_node.find_child("attributes"));
    parent_core.create_synapse(name, model);
}

sanafe::ModelInfo sanafe::description_parse_synapse_attributes_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    ModelInfo model_details;
    model_details.name = description_required_field<std::string>(
            parser, attributes, "model");
    model_details.model_parameters =
            description_parse_model_parameters_yaml(parser, attributes);

    return model_details;
}

void sanafe::description_parse_dendrite_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef dendrite_node, CoreConfiguration &parent_core)
{
    auto dendrite_name = description_required_field<std::string>(
            parser, dendrite_node, "name");
    std::pair<int, int> dendrite_range = {0, 0};
    if (dendrite_name.find("..") != std::string::npos)
    {
        dendrite_range = description_parse_range_yaml(dendrite_name);
    }

    for (int d = dendrite_range.first; d <= dendrite_range.second; ++d)
    {
        std::string name(dendrite_name);
        if (dendrite_name.find("..") != std::string::npos)
        {
            name = dendrite_name.substr(0, dendrite_name.find('[')) + '[' +
                    std::to_string(d) + ']';
        }
        const ryml::ConstNodeRef attributes =
                dendrite_node.find_child("attributes");
        if (attributes.invalid())
        {
            throw DescriptionParsingError(
                    "No attributes section defined", parser, dendrite_node);
        }

        ModelInfo model_details;
        model_details.name = description_required_field<std::string>(
                parser, attributes, "model");
        const ryml::ConstNodeRef plugin_path_node =
                attributes.find_child("plugin");
        if (!plugin_path_node.invalid())
        {
            if (plugin_path_node.has_val())
            {
                std::string plugin_path;
                plugin_path_node >> plugin_path;
                model_details.plugin_library_path = plugin_path;
            }
            else
            {
                throw DescriptionParsingError(
                        "Expected plugin path to be string", parser,
                        plugin_path_node);
            }
        }
        model_details.model_parameters =
                description_parse_model_parameters_yaml(parser, attributes);
        parent_core.create_dendrite(name, model_details);
    }
}

void sanafe::description_parse_soma_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef soma_node, CoreConfiguration &parent_core)
{
    auto soma_name =
            description_required_field<std::string>(parser, soma_node, "name");
    std::pair<int, int> soma_range = {0, 0};
    if (soma_name.find("..") != std::string::npos)
    {
        soma_range = description_parse_range_yaml(soma_name);
    }

    for (int s = soma_range.first; s <= soma_range.second; ++s)
    {
        std::string name(soma_name);
        if (soma_name.find("..") != std::string::npos)
        {
            name = soma_name.substr(0, soma_name.find('[')) + '[' +
                    std::to_string(s) + ']';
        }
        const ryml::ConstNodeRef attributes =
                soma_node.find_child("attributes");
        if (attributes.invalid())
        {
            throw DescriptionParsingError(
                    "No attributes section defined", parser, soma_node);
        }
        std::string model_str;

        ModelInfo model_details;
        model_details.name = description_required_field<std::string>(
                parser, attributes, "model");
        model_details.model_parameters =
                description_parse_model_parameters_yaml(parser, attributes);
        const ryml::ConstNodeRef plugin_path_node =
                attributes.find_child("plugin");
        if (!plugin_path_node.invalid())
        {
            if (plugin_path_node.has_val())
            {
                std::string plugin_path;
                plugin_path_node >> plugin_path;
                model_details.plugin_library_path = plugin_path;
            }
            else
            {
                throw DescriptionParsingError(
                        "Expected plugin path to be string", parser,
                        plugin_path_node);
            }
        }

        if (!attributes.find_child("noise").invalid())
        {
            // TODO: support optional noise arg again alongside the plugin mechanism
            /*
            s.noise_type = NOISE_FILE_STREAM;
            s.noise_stream = fopen(value_str.c_str(), "r");
            TRACE1(DESCRIPTION, "Opening noise str: %s\n", value_str.c_str());
            if (s.noise_stream == NULL)
            {
                INFO("Error: Failed to open noise stream: %s.\n",
                    value_str.c_str());
                exit(1);
            }
            */
        }
        parent_core.create_soma(name, model_details);
    }
}

void sanafe::description_parse_axon_out_section(const ryml::Parser &parser,
        const ryml::ConstNodeRef axon_out_node, CoreConfiguration &parent_core)
{
    auto axon_out_name = description_required_field<std::string>(
            parser, axon_out_node, "name");

    const auto &attributes = axon_out_node.find_child("attributes");
    if (attributes.invalid())
    {
        throw DescriptionParsingError(
                "No attributes section defined", parser, axon_out_node);
    }
    AxonOutPowerMetrics power_metrics;
    power_metrics.energy_message_out = description_required_field<double>(
            parser, attributes, "energy_message_out");
    power_metrics.latency_message_out = description_required_field<double>(
            parser, attributes, "latency_message_out");

    parent_core.create_axon_out(axon_out_name, power_metrics);
}

void sanafe::description_parse_core_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef core_node, const size_t parent_tile_id,
        Architecture &arch)
{
    auto core_name =
            description_required_field<std::string>(parser, core_node, "name");
    std::pair<int, int> core_range = {0, 0};
    if (core_name.find("..") != std::string::npos)
    {
        core_range = description_parse_range_yaml(core_name);
    }

    for (int c = core_range.first; c <= core_range.second; c++)
    {
        const std::string name = core_name.substr(0, core_name.find('[')) +
                '[' + std::to_string(c) + ']';
        const CorePipelineConfiguration pipeline_config =
                description_parse_core_pipeline_yaml(
                        parser, core_node["attributes"]);
        CoreConfiguration &core =
                arch.create_core(name, parent_tile_id, pipeline_config);

        if (!core_node.find_child("axon_in").invalid())
        {
            const ryml::ConstNodeRef axon_in_node = core_node["axon_in"];
            if (axon_in_node.is_seq())
            {
                for (const auto &axon : axon_in_node)
                {
                    description_parse_axon_in_section_yaml(parser, axon, core);
                }
            }
            else
            {
                description_parse_axon_in_section_yaml(
                        parser, axon_in_node, core);
            }
        }
        else
        {
            const std::string error = "No axon in section defined";
            throw DescriptionParsingError(error, parser, core_node);
        }

        if (!core_node.find_child("synapse").invalid())
        {
            const ryml::ConstNodeRef synapses = core_node["synapse"];
            if (synapses.is_seq())
            {
                for (const auto &syn : synapses)
                {
                    description_parse_synapse_section_yaml(parser, syn, core);
                }
            }
            else
            {
                description_parse_synapse_section_yaml(parser, synapses, core);
            }
        }
        else
        {
            const std::string error = "No synapse section defined";
            throw DescriptionParsingError(error, parser, core_node);
        }
        if (!core_node.find_child("dendrite").invalid())
        {
            const ryml::ConstNodeRef dendrite_node = core_node["dendrite"];
            if (dendrite_node.is_seq())
            {
                for (const auto &dendrite : dendrite_node)
                {
                    description_parse_dendrite_section_yaml(
                            parser, dendrite, core);
                }
            }
            else
            {
                description_parse_dendrite_section_yaml(
                        parser, dendrite_node, core);
            }
        }
        else
        {
            const std::string error = "No dendrite section defined";
            throw DescriptionParsingError(error, parser, core_node);
        }

        if (!core_node.find_child("soma").invalid())
        {
            const ryml::ConstNodeRef soma_node = core_node["soma"];
            if (soma_node.is_seq())
            {
                for (const auto &soma : soma_node)
                {
                    description_parse_soma_section_yaml(parser, soma, core);
                }
            }
            else
            {
                description_parse_soma_section_yaml(parser, soma_node, core);
            }
        }
        else
        {
            const std::string error = "No soma section defined";
            throw DescriptionParsingError(error, parser, core_node);
        }

        if (!core_node.find_child("axon_out").invalid())
        {
            const ryml::ConstNodeRef axon_out_node = core_node["axon_out"];
            if (axon_out_node.is_seq())
            {
                for (const auto &axon : axon_out_node)
                {
                    description_parse_axon_out_section(parser, axon, core);
                }
            }
            else
            {
                description_parse_axon_out_section(parser, axon_out_node, core);
            }
        }
        else
        {
            const std::string error = "No axon out seciont defined";
            throw DescriptionParsingError(error, parser, core_node);
        }
    }
}

sanafe::CorePipelineConfiguration sanafe::description_parse_core_pipeline_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    CorePipelineConfiguration pipeline_config{};

    pipeline_config.buffer_position = pipeline_parse_buffer_pos_str(
            description_required_field<std::string>(
                    parser, attributes, "buffer_position"));
    pipeline_config.max_neurons_supported = description_required_field<int>(
            parser, attributes, "max_neurons_supported");

    return pipeline_config;
}

sanafe::TilePowerMetrics sanafe::description_parse_tile_metrics_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    TilePowerMetrics tile_metrics;

    tile_metrics.energy_north_hop = description_required_field<double>(
            parser, attributes, "energy_north_hop");
    tile_metrics.latency_north_hop = description_required_field<double>(
            parser, attributes, "latency_north_hop");

    tile_metrics.energy_east_hop = description_required_field<double>(
            parser, attributes, "energy_east_hop");
    tile_metrics.latency_east_hop = description_required_field<double>(
            parser, attributes, "latency_east_hop");

    tile_metrics.energy_south_hop = description_required_field<double>(
            parser, attributes, "energy_south_hop");
    tile_metrics.latency_south_hop = description_required_field<double>(
            parser, attributes, "latency_south_hop");

    tile_metrics.energy_west_hop = description_required_field<double>(
            parser, attributes, "energy_west_hop");
    tile_metrics.latency_west_hop = description_required_field<double>(
            parser, attributes, "latency_west_hop");

    return tile_metrics;
}

void sanafe::description_parse_tile_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef tile_node, Architecture &arch)
{
    std::string tile_name;
    tile_node["name"] >> tile_name;
    std::pair<int, int> range = {0, 0};

    if (tile_name.find("..") != std::string::npos)
    {
        range = description_parse_range_yaml(tile_name);
    }

    for (int t = range.first; t <= range.second; t++)
    {
        const std::string name = tile_name.substr(0, tile_name.find('[')) +
                "[" + std::to_string(t) + "]";
        const TilePowerMetrics power_metrics =
                description_parse_tile_metrics_yaml(
                        parser, tile_node["attributes"]);
        TileConfiguration &new_tile = arch.create_tile(name, power_metrics);
        if (tile_node.find_child("core").invalid())
        {
            const std::string error = "No core section defined";
            throw DescriptionParsingError(error, parser, tile_node);
        }
        const ryml::ConstNodeRef core_section = tile_node["core"];
        if (core_section.is_seq())
        {
            for (const auto &core : core_section)
            {
                description_parse_core_section_yaml(
                        parser, core, new_tile.id, arch);
            }
        }
        else // Is a single core
        {
            description_parse_core_section_yaml(
                    parser, core_section, new_tile.id, arch);
        }
    }
}

sanafe::NetworkOnChipConfiguration
sanafe::description_parse_noc_configuration_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef noc_attributes)
{
    NetworkOnChipConfiguration noc;
    noc.width_in_tiles =
            description_required_field<int>(parser, noc_attributes, "width");
    noc.height_in_tiles =
            description_required_field<int>(parser, noc_attributes, "height");
    noc.link_buffer_size = description_required_field<int>(
            parser, noc_attributes, "link_buffer_size");

    return noc;
}

sanafe::Architecture sanafe::description_parse_arch_section_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef arch_node)
{
    if (arch_node.invalid())
    {
        throw DescriptionParsingError(
                "No top-level architecture section defined.\n", parser,
                arch_node);
    }
    std::string arch_name;
    arch_node["name"] >> arch_name;
    if (arch_name.find('[') != std::string::npos)
    {
        throw DescriptionParsingError(
                "Multiple architectures not supported", parser, arch_node);
    }
    NetworkOnChipConfiguration noc = description_parse_noc_configuration_yaml(
            parser, arch_node["attributes"]);
    Architecture new_arch(arch_name, noc);
    if (!arch_node.find_child("tile").invalid())
    {
        const ryml::ConstNodeRef tiles = arch_node["tile"];
        const bool is_list_of_tiles = tiles.is_seq();
        if (is_list_of_tiles)
        {
            for (const auto &tile : tiles)
            {
                description_parse_tile_section_yaml(parser, tile, new_arch);
            }
        }
        else // Only one tile defined
        {
            description_parse_tile_section_yaml(parser, tiles, new_arch);
        }
    }
    else
    {
        throw DescriptionParsingError(
                "No tile section defined", parser, arch_node);
    }

    return new_arch;
}

sanafe::Architecture sanafe::description_parse_arch_file_yaml(std::ifstream &fp)
{
    if (!fp.is_open())
    {
        throw std::runtime_error("Error opening file\n");
    }
    // Get file size
    fp.seekg(0, std::ios::end);
    std::streampos file_size = fp.tellg();
    fp.seekg(0, std::ios::beg);

    // Allocate memory
    std::string file_content;
    file_content.reserve(file_size);

    // Read the file
    file_content.assign((std::istreambuf_iterator<char>(fp)),
            std::istreambuf_iterator<char>());
    fp.close();

    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    INFO("Loading YAML information from file.\n");
    ryml::Tree top_level_yaml =
            ryml::parse_in_place(&parser, file_content.data());
    INFO("YAML information loaded from file.\n");

    if (top_level_yaml["architecture"].invalid())
    {
        throw DescriptionParsingError(
                "No architecture section defined", parser, top_level_yaml);
    }
    return description_parse_arch_section_yaml(
            parser, top_level_yaml["architecture"]);
}

sanafe::SpikingNetwork sanafe::description_parse_network_file_yaml(
        std::ifstream &fp, Architecture &arch)
{
    if (!fp.is_open())
    {
        throw std::runtime_error("Error opening file\n");
    }

    // Get file size
    fp.seekg(0, std::ios::end);
    const std::streampos file_size = fp.tellg();
    fp.seekg(0, std::ios::beg);

    // Allocate memory
    std::string file_content;
    file_content.reserve(file_size);

    // Read the file
    file_content.assign((std::istreambuf_iterator<char>(fp)),
            std::istreambuf_iterator<char>());
    fp.close();
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree top_level_yaml =
            ryml::parse_in_place(&parser, file_content.data());
    INFO("Loading network YAML information from file.\n");
    ryml::Tree tree = ryml::parse_in_place(file_content.data());
    INFO("Network YAML information loaded from file.\n");

    ryml::ConstNodeRef yaml_node = tree.rootref();
    if (yaml_node.is_map())
    {
        if (yaml_node.find_child("network").invalid())
        {
            throw DescriptionParsingError(
                    "No top-level 'network' section defined", parser,
                    yaml_node);
        }
        SpikingNetwork net = description_parse_network_section_yaml(
                parser, yaml_node["network"]);
        if (yaml_node.find_child("mappings").invalid())
        {
            throw DescriptionParsingError(
                    "No 'mappings' section defined", parser, yaml_node);
        }
        description_parse_mapping_section_yaml(
                parser, yaml_node["mappings"], arch, net);

        return net;
    }
    throw DescriptionParsingError(
            "Mapped network file has invalid format", parser, yaml_node);
}

sanafe::SpikingNetwork sanafe::description_parse_network_section_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef net_node)
{
    std::string net_name;
    if (!net_node.find_child("name").invalid())
    {
        net_node["name"] >> net_name;
        if (net_name.find('[') != std::string::npos)
        {
            throw DescriptionParsingError(
                    "Multiple networks not supported", parser, net_node);
        }
    }
    else
    {
        INFO("Warning: No network name given; leaving name empty.\n");
    }

    INFO("Parsing network: %s\n", net_name.c_str());

    SpikingNetwork new_net(std::move(net_name));
    description_parse_neuron_group_section_yaml(
            parser, net_node["groups"], new_net);
    description_parse_edges_section_yaml(parser, net_node["edges"], new_net);

    return new_net;
}

void sanafe::description_parse_neuron_group_section_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef groups_node,
        SpikingNetwork &net)
{
    INFO("Parsing neuron groups.\n");
    if (groups_node.is_seq())
    {
        for (const auto &group : groups_node)
        {
            description_parse_group(parser, group, net);
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Neuron group section does not define a list of groups", parser,
                groups_node);
    }
}

void sanafe::description_parse_edges_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef edges_node, SpikingNetwork &net)
{
    INFO("Parsing edges section.\n");
    if (edges_node.is_seq())
    {
        for (const auto list_entry : edges_node)
        {
            for (const auto edge_node : list_entry)
            {
                std::string edge_description;
                edge_node >> ryml::key(edge_description);
                description_parse_edge(
                        edge_description, parser, edge_node, net);
            }
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Edges section does not define a list of edges", parser,
                edges_node);
    }
}

void sanafe::description_parse_group(const ryml::Parser &parser,
        const ryml::ConstNodeRef neuron_group_node, SpikingNetwork &net)
{
    const auto group_name = description_required_field<std::string>(
            parser, neuron_group_node, "name");
    INFO("Parsing neuron group: %s\n", group_name.c_str());

    if (neuron_group_node.find_child("neurons").invalid())
    {
        throw DescriptionParsingError(
                "No neurons section defined.", parser, neuron_group_node);
    }
    const auto &neurons_node = neuron_group_node["neurons"];
    const size_t neuron_count = description_count_neurons(parser, neurons_node);

    NeuronConfiguration default_neuron_config{};
    if (!neuron_group_node.find_child("attributes").invalid())
    {
        TRACE1(DESCRIPTION, "Parsing neuron group attributes\n");
        default_neuron_config = description_parse_neuron_attributes_yaml(
                parser, neuron_group_node["attributes"]);
    }
    NeuronGroup &group = net.create_neuron_group(
            group_name, neuron_count, default_neuron_config);
    TRACE1(DESCRIPTION, "Parsing neuron section\n");
    description_parse_neuron_section_yaml(parser, neurons_node, group);
}

size_t sanafe::description_count_neurons(
        const ryml::Parser &parser, const ryml::ConstNodeRef neuron_node)
{
    size_t neuron_count{0UL};

    if (neuron_node.is_seq())
    {
        for (const auto &neuron_entry : neuron_node)
        {
            if (neuron_entry.is_map() || neuron_entry.is_seq())
            {
                for (const auto neuron_description : neuron_entry)
                {
                    std::string id;
                    neuron_description >> ryml::key(id);
                    const bool is_range = (id.find("..") != std::string::npos);
                    if (is_range)
                    {
                        const auto range = description_parse_range_yaml(id);
                        neuron_count += (range.second - range.first) + 1;
                    }
                    else // if entry defines a single neuron
                    {
                        ++neuron_count;
                    }
                }
            }
            else
            {
                std::string id;
                neuron_entry >> id;
                const bool is_range = (id.find("..") != std::string::npos);
                if (is_range)
                {
                    const auto range = description_parse_range_yaml(id);
                    neuron_count += (range.second - range.first) + 1;
                }
                else // if entry defines a single neuron
                {
                    ++neuron_count;
                }
            }
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Invalid neuron format, should be list", parser, neuron_node);
    }

    INFO("Counted %zu neurons\n", neuron_count);
    return neuron_count;
}

void sanafe::description_parse_neuron_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef neuron_node, NeuronGroup &neuron_group)
{
    if (neuron_node.is_seq())
    {
        for (const auto &list_entry : neuron_node)
        {
            for (const auto &neuron_description : list_entry)
            {
                // Iterate, but there should only be one mapping per list entry
                std::string id;
                neuron_description >> ryml::key(id);
                description_parse_neuron(
                        id, parser, neuron_description, neuron_group);
            }
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Invalid neuron format, should be list", parser, neuron_node);
    }
}

void sanafe::description_parse_neuron(const std::string &id,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes,
        NeuronGroup &neuron_group)
{
    std::pair<size_t, size_t> range;
    TRACE1(DESCRIPTION, "Parsing neuron(s): %s\n", id.c_str());
    const NeuronConfiguration config = description_parse_neuron_attributes_yaml(
            parser, attributes, neuron_group.default_neuron_config);
    const bool is_range = (id.find("..") != std::string::npos);
    if (is_range)
    {
        range = description_parse_range_yaml(id);
        for (size_t instance = range.first; instance <= range.second;
                ++instance)
        {
            Neuron &n = neuron_group.neurons.at(instance);
            n.configure(config);
        }
    }
    else
    {
        const size_t nid = std::stoull(id);
        Neuron &n = neuron_group.neurons.at(nid);
        n.configure(config);
    }
}

sanafe::NeuronConfiguration sanafe::description_parse_neuron_attributes_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes,
        const NeuronConfiguration &default_template)
{
    NeuronConfiguration neuron_template = default_template;

    if (attributes.is_seq())
    {
        // Ordered list format, recursively parse attributes for each element
        for (const auto &attribute : attributes)
        {
            neuron_template = description_parse_neuron_attributes_yaml(
                    parser, attribute, neuron_template);
        }
        return neuron_template;
    }

    if (!attributes.find_child("log_potential").invalid())
    {
        bool log_potential;
        attributes["log_potential"] >> log_potential;
        neuron_template.log_potential = log_potential;
    }
    if (!attributes.find_child("log_spikes").invalid())
    {
        bool log_spikes;
        attributes["log_spikes"] >> log_spikes;
        neuron_template.log_spikes = log_spikes;
    }
    if (!attributes.find_child("force_synapse_update").invalid())
    {
        bool force_synapse_update;
        attributes["force_synapse_update"] >> force_synapse_update;
        neuron_template.force_synapse_update = force_synapse_update;
    }
    if (!attributes.find_child("force_dendrite_update").invalid())
    {
        bool force_dendrite_update;
        attributes["force_dendrite_update"] >> force_dendrite_update;
        neuron_template.force_dendrite_update = force_dendrite_update;
    }
    if (!attributes.find_child("force_soma_update").invalid())
    {
        bool force_soma_update;
        attributes["force_soma_update"] >> force_soma_update;
        neuron_template.force_soma_update = force_soma_update;
    }
    if (!attributes.find_child("synapse_hw_name").invalid())
    {
        std::string synapse_hw_name;
        attributes["synapse_hw_name"] >> synapse_hw_name;
        neuron_template.default_synapse_hw_name = synapse_hw_name;
    }
    if (!attributes.find_child("dendrite_hw_name").invalid())
    {
        std::string dendrite_hw_name;
        attributes["dendrite_hw_name"] >> dendrite_hw_name;
        neuron_template.dendrite_hw_name = dendrite_hw_name;
    }
    if (!attributes.find_child("soma_hw_name").invalid())
    {
        std::string soma_hw_name;
        attributes["soma_hw_name"] >> soma_hw_name;
        neuron_template.soma_hw_name = soma_hw_name;
    }

    // Parse and add shared parameters, which are defined alongside attributes
    auto model_params =
            description_parse_model_parameters_yaml(parser, attributes);
    for (auto &[key, parameter] : model_params)
    {
        parameter.forward_to_dendrite = true;
        parameter.forward_to_soma = true;
        neuron_template.model_parameters[key] = parameter;
    }
    // Parse h/w unit specific model parameters
    if (!attributes.find_child("dendrite").invalid())
    {
        auto dendrite_parameters = description_parse_model_parameters_yaml(
                parser, attributes["dendrite"]);
        for (auto &[key, param] : dendrite_parameters)
        {
            param.forward_to_synapse = false;
            param.forward_to_soma = false;
            neuron_template.model_parameters[key] = param;
        }
    }
    if (!attributes.find_child("soma").invalid())
    {
        auto soma_parameters = description_parse_model_parameters_yaml(
                parser, attributes["soma"]);
        for (auto &[key, parameter] : soma_parameters)
        {
            parameter.forward_to_synapse = false;
            parameter.forward_to_dendrite = false;
            neuron_template.model_parameters[key] = parameter;
        }
    }

    return neuron_template;
}

std::string_view description_trim_whitespace(const std::string_view input)
{
    constexpr auto whitespace = " \t\n\r";
    auto start = input.find_first_not_of(whitespace);
    auto end = input.find_last_not_of(whitespace);
    if ((start == std::string::npos) || (end == std::string::npos))
    {
        return "";
    }
    return input.substr(start, end - start + 1);
}

std::tuple<sanafe::NeuronAddress, sanafe::NeuronAddress>
sanafe::description_parse_edge_description(const std::string_view &description)
{
    auto arrow_pos = description.find("->");
    if (arrow_pos == std::string::npos)
    {
        throw std::runtime_error(
                "Edge is not formatted correctly: " + std::string(description));
    }

    const std::string_view source_part =
            description_trim_whitespace(description.substr(0, arrow_pos));
    const std::string_view target_part =
            description_trim_whitespace(description.substr(arrow_pos + 2));

    const auto source_dot_pos = source_part.find('.');
    const auto target_dot_pos = target_part.find('.');

    const bool source_neuron_defined = (source_dot_pos != std::string::npos);
    const bool target_neuron_defined = (target_dot_pos != std::string::npos);
    if (source_neuron_defined && !target_neuron_defined)
    {
        throw std::runtime_error(
                "No target neuron defined in edge:" + std::string(description));
    }
    else if (target_neuron_defined && !source_neuron_defined)
    {
        throw std::runtime_error(
                "No target neuron defined in edge:" + std::string(description));
    }

    NeuronAddress source;
    source.group_name = source_part.substr(0, source_dot_pos);
    if (source_neuron_defined)
    {
        source.neuron_id = std::stoull(
                std::string(source_part.substr(source_dot_pos + 1)));
    }

    NeuronAddress target;
    target.group_name = target_part.substr(0, target_dot_pos);
    if (target_neuron_defined)
    {
        target.neuron_id = std::stoull(
                std::string(target_part.substr(target_dot_pos + 1)));
    }

    return std::make_tuple(source, target);
}

void sanafe::description_parse_edge(const std::string &description,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node,
        SpikingNetwork &net)
{
    // Description has format src_group.src_neuron -> tgt_group.tgt_neuron
    const auto [source_address, target_address] =
            description_parse_edge_description(description);

    const bool is_hyper_edge = !source_address.neuron_id.has_value();
    if (is_hyper_edge)
    {
        description_parse_hyperedge(
                source_address, target_address, parser, attributes_node, net);
    }
    else
    {
        description_parse_neuron_connection(
                source_address, target_address, parser, attributes_node, net);
    }
}

void sanafe::description_parse_neuron_connection(
        const NeuronAddress &source_address,
        const NeuronAddress &target_address, const ryml::Parser &parser,
        const ryml::ConstNodeRef attributes_node, SpikingNetwork &net)
{
    if (net.groups.find(source_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid source neuron group:" + source_address.group_name;
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    NeuronGroup &source_group = net.groups.at(source_address.group_name);
    if (source_address.neuron_id >= source_group.neurons.size())
    {
        const std::string error =
                "Invalid source neuron id: " + source_address.group_name + "." +
                std::to_string(source_address.neuron_id.value());
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    Neuron &source_neuron =
            source_group.neurons[source_address.neuron_id.value()];

    if (net.groups.find(target_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid target neuron group:" + target_address.group_name;
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    NeuronGroup &target_group = net.groups.at(target_address.group_name);

    if (target_address.neuron_id >= target_group.neurons.size())
    {
        const std::string error =
                "Invalid target neuron id: " + target_address.group_name + "." +
                std::to_string(target_address.neuron_id.value());
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    Neuron &target_neuron =
            target_group.neurons.at(target_address.neuron_id.value());

    const size_t edge_idx = source_neuron.connect_to_neuron(target_neuron);
    Connection &edge = source_neuron.edges_out[edge_idx];
    description_parse_edge_attributes(edge, parser, attributes_node);
}

void sanafe::description_parse_hyperedge(const NeuronAddress &source_address,
        const NeuronAddress &target_address, const ryml::Parser &parser,
        const ryml::ConstNodeRef hyperedge_node, SpikingNetwork &net)
{
    if (net.groups.find(source_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid source neuron group:" + source_address.group_name;
        throw DescriptionParsingError(error, parser, hyperedge_node);
    }
    NeuronGroup &source_group = net.groups.at(source_address.group_name);

    if (net.groups.find(target_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid target neuron group:" + target_address.group_name;
        throw DescriptionParsingError(error, parser, hyperedge_node);
    }

    NeuronGroup &target_group = net.groups.at(target_address.group_name);

    // First parse the lists of attributes, this is common for all hyperedge
    //  connectivity
    const auto type = description_required_field<std::string>(
            parser, hyperedge_node, "type");

    std::map<std::string, std::vector<ModelParam>> attribute_lists{};
    if (type == "conv2d")
    {
        Conv2DParameters convolution{};
        for (const auto &attribute : hyperedge_node.children())
        {
            if (attribute.key() == "input_height")
            {
                attribute >> convolution.input_height;
            }
            else if (attribute.key() == "input_width")
            {
                attribute >> convolution.input_width;
            }
            else if (attribute.key() == "input_channels")
            {
                attribute >> convolution.input_channels;
            }
            else if (attribute.key() == "kernel_width")
            {
                attribute >> convolution.kernel_width;
            }
            else if (attribute.key() == "kernel_height")
            {
                attribute >> convolution.kernel_height;
            }
            else if (attribute.key() == "kernel_count")
            {
                attribute >> convolution.kernel_count;
            }
            else if (attribute.key() == "stride_width")
            {
                attribute >> convolution.stride_width;
            }
            else if (attribute.key() == "stride_height")
            {
                attribute >> convolution.stride_height;
            }
            else if (attribute.key() == "synapse")
            {
                for (const auto &synapse_param_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelParam> attribute_list;
                    for (const auto &model_param_node : synapse_param_node)
                    {
                        ModelParam value = description_parse_parameter_yaml(
                                parser, model_param_node);
                        value.forward_to_dendrite = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(value);
                    }
                    std::string attribute_name;
                    synapse_param_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = attribute_list;
                }
            }
            else if (attribute.key() == "dendrite")
            {
                for (const auto &dendrite_param_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelParam> attribute_list;
                    for (const auto &model_param_node : dendrite_param_node)
                    {
                        ModelParam value = description_parse_parameter_yaml(
                                parser, model_param_node);
                        value.forward_to_synapse = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(value);
                    }
                    std::string attribute_name;
                    dendrite_param_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = attribute_list;
                }
            }
            else if (attribute.key() != "type")
            {
                // TODO: refactor
                std::vector<ModelParam> attribute_list;
                for (const auto &model_param_node : attribute)
                {
                    ModelParam value = description_parse_parameter_yaml(
                            parser, model_param_node);
                    attribute_list.push_back(value);
                }
                std::string attribute_name;
                attribute >> ryml::key(attribute_name);
                attribute_lists[attribute_name] = attribute_list;
            }
        }

        source_group.connect_neurons_conv2d(
                target_group, attribute_lists, convolution);
    }
    else if (type == "dense")
    {
        for (const auto &attribute : hyperedge_node.children())
        {
            // TODO: refactor
            if (attribute.key() == "synapse")
            {
                for (const auto &synapse_param_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelParam> attribute_list;
                    for (const auto &model_param_node : synapse_param_node)
                    {
                        ModelParam value = description_parse_parameter_yaml(
                                parser, model_param_node);
                        value.forward_to_dendrite = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(value);
                    }
                    std::string attribute_name;
                    synapse_param_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = attribute_list;
                }
            }
            else if (attribute.key() == "dendrite")
            {
                for (const auto &dendrite_param_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelParam> attribute_list;
                    for (const auto &model_param_node : dendrite_param_node)
                    {
                        ModelParam value = description_parse_parameter_yaml(
                                parser, model_param_node);
                        value.forward_to_synapse = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(value);
                    }
                    std::string attribute_name;
                    dendrite_param_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = attribute_list;
                }
            }
            else if (attribute.key() != "type")
            {
                // TODO: refactor
                std::vector<ModelParam> attribute_list;
                for (const auto &model_param_node : attribute)
                {
                    ModelParam value = description_parse_parameter_yaml(
                            parser, model_param_node);
                    attribute_list.push_back(value);
                }
                std::string attribute_name;
                attribute >> ryml::key(attribute_name);
                attribute_lists[attribute_name] = attribute_list;
            }
        }
        source_group.connect_neurons_dense(target_group, attribute_lists);
    }
    else if (type == "sparse")
    {
    }
    else
    {
        const std::string error = "Invalid hyperedge type: " + type;
        throw DescriptionParsingError(error, parser, hyperedge_node);
    }
}

void sanafe::description_parse_edge_attributes(Connection &edge,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node)
{
    if (attributes_node.is_seq())
    {
        for (const auto &attribute : attributes_node)
        {
            description_parse_edge_attributes(edge, parser, attribute);
        }

        return;
    }

    if (!attributes_node.find_child("synapse").invalid())
    {
        auto synapse_parameters = description_parse_model_parameters_yaml(
                parser, attributes_node["synapse"]);
        for (auto &[key, param] : synapse_parameters)
        {
            param.forward_to_dendrite = false;
            param.forward_to_soma = false;
            edge.synapse_params[key] = param;
        }
    }
    if (!attributes_node.find_child("dendrite").invalid())
    {
        auto dendrite_params = description_parse_model_parameters_yaml(
                parser, attributes_node["dendrite"]);
        for (auto &[key, param] : dendrite_params)
        {
            param.forward_to_synapse = false;
            param.forward_to_soma = false;
            edge.dendrite_params[key] = param;
        }
    }

    const auto shared_model_params =
            description_parse_model_parameters_yaml(parser, attributes_node);
    for (const auto &[key, parameter] : shared_model_params)
    {
        if ((key != "synapse") && (key != "dendrite") && (key != "soma"))
        {
            TRACE2(DESCRIPTION, "Adding con parameter:%s\n", key.c_str());
            edge.synapse_params[key] = parameter;
            edge.dendrite_params[key] = parameter;
        }
    }
}

void sanafe::description_parse_mapping_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef mappings_node, Architecture &arch,
        SpikingNetwork &net)
{
    INFO("Parsing mapping section.\n");
    if (!mappings_node.is_seq())
    {
        throw DescriptionParsingError(
                "Mappings must be given as a sequence / list. Each list entry "
                "must first use the neuron address as a key, followed by the "
                "core address.\nE.g.: "
                "- G.n: {core: 1.1}\nThis maps group G's neuron 'n' "
                "to Tile 1's Core 1, specifically using soma unit 'foo'.)",
                parser, mappings_node);
    }

    for (const auto &mapping : mappings_node)
    {
        if (!mapping.is_map())
        {
            throw DescriptionParsingError(
                    "Expected mapping to be defined in the format: "
                    "<group>.<neuron>: [<attributes>]",
                    parser, mapping);
        }

        int entries = 0;
        // Ordered mappings, so each list entry contains a mapping with a single
        //  key. Only look up the first entry.
        for (const auto &mapping_entry : mapping)
        {
            description_parse_mapping(parser, mapping_entry, arch, net);
            ++entries;
            if (entries > 1)
            {
                throw DescriptionParsingError("Should be one entry per mapping",
                        parser, mapping_entry);
            }
        }
    }
}

void sanafe::description_parse_mapping(const ryml::Parser &parser,
        const ryml::ConstNodeRef mapping_info, Architecture &arch,
        SpikingNetwork &net)
{
    std::string neuron_address;
    mapping_info >> ryml::key(neuron_address);
    auto dot_pos = neuron_address.find_first_of('.');

    std::string group_name = neuron_address.substr(0, dot_pos);
    std::string neuron_str = neuron_address.substr(dot_pos + 1);
    if (net.groups.find(group_name) == net.groups.end())
    {
        std::string error = "Invalid neuron group:" + group_name;
        throw DescriptionParsingError(error, parser, mapping_info);
    }
    //INFO("Mapping neuron: %s.%zu\n", group_name.c_str(), neuron_id);
    NeuronGroup &group = net.groups.at(group_name);

    // TODO: support mapping neuron groups (i.e. no neuron field given)
    // Parse range of neurons to map
    size_t start_id;
    size_t end_id;
    if (neuron_str.find("..") != std::string::npos)
    {
        std::tie(start_id, end_id) = description_parse_range_yaml(neuron_str);
    }
    else
    {
        start_id = std::stoull(neuron_str);
        end_id = start_id;
    }
    for (size_t neuron_id = start_id; neuron_id <= end_id; ++neuron_id)
    {
        if (neuron_id >= group.neurons.size())
        {
            std::string error = "Invalid neuron id: ";
            error += group_name;
            error += '.';
            error += neuron_id;
            throw DescriptionParsingError(error, parser, mapping_info);
        }
        Neuron &neuron = group.neurons[neuron_id];

        std::string core_address;
        description_parse_mapping_info(
                parser, mapping_info, neuron, core_address);
        const auto dot_pos = core_address.find('.');

        const size_t tile_id = std::stoull(core_address.substr(0, dot_pos));
        const size_t core_offset_within_tile =
                std::stoull(core_address.substr(dot_pos + 1));

        if (tile_id >= arch.tiles.size())
        {
            throw DescriptionParsingError(
                    "Tile ID >= tile count", parser, mapping_info);
        }
        TileConfiguration &tile = arch.tiles[tile_id];
        if (core_offset_within_tile >= tile.cores.size())
        {
            throw DescriptionParsingError(
                    "Core ID >= core count", parser, mapping_info);
        }
        CoreConfiguration &core = tile.cores[core_offset_within_tile];
        neuron.map_to_core(core);
    }
}

void sanafe::description_parse_mapping_info(const ryml::Parser &parser,
        const ryml::ConstNodeRef info, Neuron &n, std::string &core_name)
{
    if (info.is_seq())
    {
        for (const auto &field : info)
        {
            description_parse_mapping_info(parser, field, n, core_name);
        }
    }
    else if (!info.is_map())
    {
        throw DescriptionParsingError(
                "Expected attributes to be map", parser, info);
    }
    else
    {
        if (!info.find_child("synapse").invalid())
        {
            info["synapse"] >> n.default_synapse_hw_name;
        }
        if (!info.find_child("dendrite").invalid())
        {
            info["dendrite"] >> n.dendrite_hw_name;
        }
        if (!info.find_child("soma").invalid())
        {
            info["soma"] >> n.soma_hw_name;
        }
        if (!info.find_child("core").invalid())
        {
            info["core"] >> core_name;
        }
    }

    return;
}

std::map<std::string, sanafe::ModelParam>
sanafe::description_parse_model_parameters_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef parameters_node)
{
    std::map<std::string, ModelParam> model_parameters;

    if (parameters_node.is_seq())
    {
        for (const auto &mapping_list_node : parameters_node)
        {
            // Recursive call to flatten list of mappings
            auto new_parameters = description_parse_model_parameters_yaml(
                    parser, mapping_list_node);
            model_parameters.insert(
                    new_parameters.begin(), new_parameters.end());
        }
    }
    else if (parameters_node.is_map())
    {
        for (const auto &node : parameters_node)
        {
            std::string key_str;
            node >> ryml::key(key_str);
            std::set<std::string> unit_specific_keys = {
                    "synapse", "dendrite", "soma"};
            if (unit_specific_keys.find(key_str) == unit_specific_keys.end())
            {
                //INFO("Parsing attribute: %s\n", key.c_str());
                model_parameters[key_str] =
                        description_parse_parameter_yaml(parser, node);
            }
        }
    }
    else
    {
        throw std::invalid_argument(
                "Error: Model parameters must be an ordered map or mapping of"
                "named parameters.\n");
    }

    return model_parameters;
}

sanafe::ModelParam sanafe::description_parse_parameter_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node)
{
    ModelParam attribute;

    if (attribute_node.is_seq())
    {
        // Create an list of unnamed attributes
        std::vector<ModelParam> attribute_list;
        for (const auto &node : attribute_node)
        {
            INFO("Parsing sub-parameter in list.\n");
            ModelParam curr = description_parse_parameter_yaml(parser, node);
            attribute_list.push_back(curr);
        }
        INFO("Setting parameter to an list of %zu unnamed parameters\n",
                attribute_list.size());
        attribute.value = attribute_list;
    }
    else if (attribute_node.is_map())
    {
        // Create a list of named attributes
        std::vector<ModelParam> attribute_list;
        for (const auto &node : attribute_node)
        {
            INFO("Parsing mapping of parameters.\n");
            ModelParam curr = description_parse_parameter_yaml(parser, node);
            std::string key;
            node >> ryml::key(key);
            curr.name = key;
            INFO("Saving to key: %s\n", key.c_str());
            attribute_list.push_back(curr);
        }
        INFO("Setting parameter to a list of %zu named attributes\n",
                attribute_list.size());
        attribute.value = attribute_list;
    }
    else
    {
        if (attribute_node.invalid())
        {
            throw std::invalid_argument("Invalid node.\n");
        }
        int decoded_int;
        double decoded_double;
        bool decoded_bool;
        std::string decoded_str;
        if (c4::yml::read(attribute_node, &decoded_int))
        {
            //INFO("Parsed int: %d.\n", decoded_int);
            attribute.value = decoded_int;
        }
        else if (c4::yml::read(attribute_node, &decoded_double))
        {
            //INFO("Parsed float: %lf.\n", decoded_double);
            attribute.value = decoded_double;
        }
        else if (c4::yml::read(attribute_node, &decoded_bool))
        {
            //INFO("Parsed bool: %d.\n", decoded_bool);
            attribute.value = decoded_bool;
        }
        else if (c4::yml::read(attribute_node, &decoded_str))
        {
            //INFO("Parsed string: %s.\n", decoded_str.c_str());
            attribute.value = std::move(decoded_str);
        }
    }

    return attribute;
}

std::pair<size_t, size_t> sanafe::description_parse_range_yaml(
        const std::string &range_str)
{
    const size_t delimiter_pos = range_str.find("..");

    size_t start_pos;
    size_t end_pos;
    if (range_str.find('[') == std::string::npos)
    {
        start_pos = 0;
    }
    else
    {
        start_pos = range_str.find('[') + 1;
    }

    if (range_str.find(']') == std::string::npos)
    {
        end_pos = range_str.length();
    }
    else
    {
        end_pos = range_str.find(']');
    }

    if ((end_pos <= start_pos) || (delimiter_pos == std::string::npos) ||
            (delimiter_pos <= start_pos) || (delimiter_pos >= end_pos))
    {
        throw std::runtime_error("Invalid range format");
    }

    const std::string first_str =
            range_str.substr(start_pos, delimiter_pos - start_pos);
    const std::string last_str =
            range_str.substr(delimiter_pos + 2, (end_pos - delimiter_pos) - 2);

    size_t first;
    size_t last;
    try
    {
        first = std::stoull(first_str);
        last = std::stoull(last_str);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Invalid range string");
    }
    TRACE1(DESCRIPTION, "Range: %zu to %zu\n", first, last);
    if (first > last)
    {
        throw std::runtime_error("Invalid range; first > last");
    }

    return std::make_pair(first, last);
}

// Markdown (v1) network description format. Supported for back-compatability.
//  This format is useful for extremely large network files, as parsing this
//  simpler format requires less memory than YAML parsing.
constexpr int default_line_len = 4096;
constexpr int default_fields = 32;

sanafe::SpikingNetwork sanafe::description_parse_network_file_netlist(
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
        description_get_fields(fields, line);

        TRACE1(DESCRIPTION, "%ld fields.\n", fields.size());

        if (DEBUG_LEVEL_DESCRIPTION > 0)
        {
            for (auto f : fields)
            {
                fprintf(stdout, "\tField:%s\n", std::string(f).c_str());
            }
        }
        if (fields.size() > 0)
        {
            description_read_network_entry(fields, arch, net, line_number);
        }
        line_number++;
    }

    return net;
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

        const std::string_view new_field =
                line_buffer.substr(field_start, field_end - field_start);
        if (field_end != field_start)
        {
            fields.push_back(new_field);
        }
        field_start = line_buffer.find_first_not_of(delim, field_end);
    }

    return;
}

std::pair<std::string, size_t> sanafe::parse_neuron_field(
        const std::string_view &neuron_field)
{
    const auto pos = neuron_field.find('.');
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Error: Invalid neuron format");
    }

    std::string group_id(neuron_field.substr(0, pos));
    const auto neuron_str = neuron_field.substr(pos + 1, neuron_field.size());
    size_t neuron_id = field_to_int(neuron_str);

    return {group_id, neuron_id};
}

std::pair<size_t, size_t> sanafe::parse_core_field(
        const std::string_view &core_field)
{
    const auto pos = core_field.find('.');
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Error: Invalid neuron format");
    }

    auto tile_str = core_field.substr(0, pos);
    size_t tile_id = field_to_int(tile_str);
    auto core_str = core_field.substr(pos + 1, core_field.size());
    size_t core_offset = field_to_int(core_str);

    return {tile_id, core_offset};
}

std::tuple<std::string, size_t, std::string, size_t> sanafe::parse_edge_field(
        const std::string_view &edge_field)
{
    // Edge description entries support two formats, to represent
    //  neuron-neuron connections and compartment-compartment branches
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
    const auto [group_id, neuron_id] = parse_neuron_field(src_neuron_address);
    // Parse the destination group, neuron and compartment
    //  identifiers from the substring after "->"
    const auto dest_neuron_address =
            edge_field.substr(pos + 2, edge_field.size());
    const auto [dest_group_id, dest_neuron_id] =
            parse_neuron_field(dest_neuron_address);

    return {group_id, neuron_id, dest_group_id, dest_neuron_id};
}

std::tuple<std::string, size_t, size_t, size_t> sanafe::parse_mapping_field(
        const std::string_view &mapping_field)
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
    auto [group_id, neuron_id] = parse_neuron_field(neuron_address);

    const auto core_address =
            mapping_field.substr(pos + 1, mapping_field.size());
    auto [tile_id, core_offset] = parse_core_field(core_address);

    return {group_id, neuron_id, tile_id, core_offset};
}

void sanafe::description_read_network_entry(
        const std::vector<std::string_view> &fields, Architecture &arch,
        SpikingNetwork &net, const int line_number)
{
    std::map<std::string, ModelParam> attributes;
    Neuron *neuron_ptr;
    Neuron *dest_ptr;
    TileConfiguration *tile_ptr;
    CoreConfiguration *core_ptr;
    std::string neuron_group_id;
    std::string dest_group_id;
    size_t tile_id;
    size_t core_offset;
    size_t neuron_id;
    size_t dest_neuron_id;
    int neuron_count;
    bool group_set;
    bool neuron_set;

    assert(fields.size() > 0);
    const char entry_type = fields[0][0];
    // Sanity check input
    if ((entry_type == '\0') || (entry_type == '\n') || (entry_type == '#') ||
            (entry_type == '\r'))
    {
        TRACE1(DESCRIPTION, "Warning: No entry, skipping line %d\n", line_number);
        return;
    }

    if (fields.size() < 2)
    {
        INFO("Error: fields < 2 (%ld)", fields.size());
    }

    neuron_count = 0;
    neuron_ptr = nullptr;
    core_ptr = nullptr;
    dest_ptr = nullptr;

    group_set = false;
    neuron_set = false;

    if (entry_type == 'g')
    {
        neuron_count = field_to_int(fields[1]);
        neuron_group_id = std::to_string(net.groups.size());
    }
    else if (entry_type == '&')
    {
        std::tie(neuron_group_id, neuron_id, tile_id, core_offset) =
                parse_mapping_field(fields[1]);
        if (tile_id >= arch.tiles.size())
        {
            INFO("Error: Line %d: Tile (%lu) >= tile count (%lu)\n",
                    line_number, tile_id, arch.tiles.size());
            throw std::runtime_error("Error: Couldn't parse mapping.");
        }
        TileConfiguration &tile = arch.tiles[tile_id];
        tile_ptr = &tile;

        if (core_offset >= tile_ptr->cores.size())
        {
            INFO("Error: Line %d: Core (%lu) >= core count (%lu)\n",
                    line_number, core_offset, tile_ptr->cores.size());
            throw std::runtime_error("Error: Couldn't parse mapping.");
        }
        CoreConfiguration &core = tile_ptr->cores[core_offset];
        core_ptr = &core;
        group_set = true;
        neuron_set = true;
    }
    else if (entry_type == 'e')
    {
        // Edge on SNN graph (e.g., connection between neurons)
        std::tie(neuron_group_id, neuron_id, dest_group_id, dest_neuron_id) =
                parse_edge_field(fields[1]);
        if (net.groups.find(dest_group_id) == net.groups.end())
        {
            INFO("Error: Line %d: Group (%s) not in groups.\n", line_number,
                    dest_group_id.c_str());
            throw std::invalid_argument("Invalid group id");
        }
        NeuronGroup &dest_group = net.groups.at(dest_group_id);

        TRACE1(DESCRIPTION, "Parsed neuron gid:%s nid:%lu\n",
                dest_group_id.c_str(), neuron_id);
        if (dest_neuron_id >= dest_group.neurons.size())
        {
            INFO("Error: Line %d: Trying to access neuron "
                 "(%s.%lu) but group %s only "
                 "allocates %lu neurons.\n",
                    line_number, dest_group.name.c_str(), dest_neuron_id,
                    dest_group.name.c_str(), dest_group.neurons.size());
            throw std::invalid_argument("Invalid nid");
        }
        dest_ptr = &(dest_group.neurons[dest_neuron_id]);
        group_set = true;
        neuron_set = true;
    }
    else if (entry_type == 'n') // parse neuron
    {
        std::tie(neuron_group_id, neuron_id) = parse_neuron_field(fields[1]);
        group_set = true;
        neuron_set = true;
    }
    else
    {
        INFO("Error: Line %d: Invalid entry type (%s)", line_number,
                std::string(fields[0]).c_str());
        throw std::invalid_argument("Invalid description entry type");
    }

    if (group_set)
    {
        if (net.groups.find(neuron_group_id) == net.groups.end())
        {
            INFO("Error: Line %d: Group (%s) not in groups.\n", line_number,
                    neuron_group_id.c_str());
            throw std::invalid_argument("Invalid group id");
        }
        NeuronGroup &group = net.groups.at(neuron_group_id);
        TRACE1(DESCRIPTION, "Parsed neuron gid:%s nid:%lu\n",
                neuron_group_id.c_str(), neuron_id);
        if (neuron_set)
        {
            if (neuron_id >= group.neurons.size())
            {
                INFO("Error: Line %d: Trying to access neuron "
                     "(%s.%lu) but group %s only "
                     "allocates %lu neuron(s).\n",
                        line_number, group.name.c_str(), neuron_id,
                        group.name.c_str(), group.neurons.size());
                throw std::invalid_argument("Invalid neuron id");
            }
            neuron_ptr = &(group.neurons[neuron_id]);
        }
    }

    // Parse attributes from remaining fields

    // TODO:
    // if the first character is [ or { read using a YAML parser, otherwise
    //  parse as list of space separated fields using <key>=<value> as below
    std::map<std::string, ModelParam> params{};
    for (size_t i = 2; i < fields.size(); i++)
    {
        TRACE1(DESCRIPTION, "Parsing field:%s\n",
                std::string(fields[i]).c_str());

        if ((fields[i].length() < 3))
        {
            INFO("Error: Line %d: Invalid field: %s\n", line_number,
                    std::string(fields[i]).c_str());
            continue;
        }

        const auto pos = fields[i].find_first_of('=');
        std::string key(fields[i].substr(0, pos));
        std::string value_str(fields[i].substr(pos + 1));

        ModelParam parameter;
        int decoded_int;
        double decoded_double;
        bool decoded_bool;
        parameter.name = key;

        std::stringstream int_ss(value_str);
        std::stringstream float_ss(value_str);
        std::stringstream bool_ss(value_str);
        if ((int_ss >> decoded_int) && int_ss.eof())
        {
            TRACE1(DESCRIPTION, "Parsed integer: %d (%s).\n",
                    decoded_int, key.c_str());
            parameter.value = decoded_int;
        }
        else if ((float_ss >> decoded_double) && float_ss.eof())
        {
            TRACE1(DESCRIPTION, "Parsed float: %e (%s).\n",
                    decoded_double, key.c_str());
            //INFO("Parsed float: %e (%s).\n",
            //        decoded_double, key.c_str());
            parameter.value = decoded_double;
        }
        else if ((bool_ss >> decoded_bool) && bool_ss.eof())
        {
            TRACE1(DESCRIPTION, "Parsed bool: %d.\n", decoded_bool);
            parameter.value = decoded_bool;
        }
        else
        {
            // Parsed string
            TRACE1(DESCRIPTION, "Parsed string: %s\n", value_str.c_str());
            parameter.value = std::string(value_str);
        }

        if ((key.length() == 0) || (value_str.length() == 0))
        {
            INFO("Error: Line %d: Invalid attribute: %s\n", line_number,
                    std::string(fields[i]).c_str());
            continue;
        }

        params.insert({key, parameter});
    }

    NeuronConfiguration neuron_config{};
    if (group_set)
    {
        NeuronGroup &group = net.groups.at(neuron_group_id);
        neuron_config = group.default_neuron_config;
    }
    // Process simulator specific keys
    if (params.find("synapse_hw_name") != params.end())
    {
        neuron_config.default_synapse_hw_name =
                static_cast<std::string>(params["synapse_hw_name"]);
    }
    if (params.find("dendrite_hw_name") != params.end())
    {
        neuron_config.dendrite_hw_name =
                static_cast<std::string>(params["dendrite_hw_name"]);
    }
    if (params.find("soma_hw_name") != params.end())
    {
        neuron_config.soma_hw_name =
                static_cast<std::string>(params["soma_hw_name"]);
    }
    if (params.find("force_update") != params.end())
    {
        neuron_config.force_soma_update =
                static_cast<bool>(params["force_update"]);
    }
    if (params.find("log_spikes") != params.end())
    {
        neuron_config.log_spikes = static_cast<bool>(params["log_spikes"]);
    }
    if (params.find("log_v") != params.end())
    {
        neuron_config.log_potential = static_cast<bool>(params["log_v"]);
    }

    // Process the entry
    neuron_config.model_parameters = params;
    switch (entry_type)
    {
    case 'g': // Add neuron group
        net.create_neuron_group(neuron_group_id, neuron_count, neuron_config);
        break;
    case 'n': // Add neuron
        neuron_ptr->configure(neuron_config);
        break;
    case 'e': {
        assert(neuron_ptr != nullptr);
        const size_t idx = neuron_ptr->connect_to_neuron(*dest_ptr);
        Connection &con = neuron_ptr->edges_out[idx];
        con.synapse_params = params;
        con.dendrite_params = params;
    }
    break;
    case '&': // Map neuron to hardware
        neuron_ptr->map_to_core(*core_ptr);
        break;
    default:
        break;
    }

    return;
}

size_t sanafe::field_to_int(const std::string_view &field)
{
    size_t val = 0;
    auto [ptr, error_code] =
            std::from_chars(field.data(), field.data() + field.size(), val);
    if (error_code != std::errc())
    {
        std::string error_str = "Error: Couldn't parse integer val for field:" +
                std::string(field);
        throw std::runtime_error(std::string(error_str));
    }

    return val;
}
