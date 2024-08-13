#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric> // For std::accumulate
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include <ryml.hpp>
#include <ryml_std.hpp>

#include "arch.hpp"
#include "description.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"

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
        const ryml::ConstNodeRef axon_in_node, Core &parent_core)
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
        const ryml::ConstNodeRef synapse_node, Core &parent_core)
{
    auto name = description_required_field<std::string>(
            parser, synapse_node, "name");
    auto [power_metrics, model] = description_parse_synapse_attributes_yaml(
            parser, synapse_node.find_child("attributes"));
    parent_core.create_synapse(name, power_metrics, model);
}

std::pair<sanafe::SynapsePowerMetrics, sanafe::ModelInfo>
sanafe::description_parse_synapse_attributes_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    ModelInfo model;
    model.name = description_required_field<std::string>(
            parser, attributes, "model");
    const ryml::ConstNodeRef plugin_path_node = attributes.find_child("plugin");
    if (!plugin_path_node.invalid())
    {
        if (plugin_path_node.is_val())
        {
            std::string plugin_path;
            plugin_path_node >> plugin_path;
            model.plugin_library_path = plugin_path;
        }
        else
        {
            throw DescriptionParsingError("Expected plugin path to be string",
                    parser, plugin_path_node);
        }
    }

    SynapsePowerMetrics power_metrics;
    power_metrics.energy_process_spike = description_required_field<double>(
            parser, attributes, "energy_process_spike");
    power_metrics.latency_process_spike = description_required_field<double>(
            parser, attributes, "latency_process_spike");

    return {power_metrics, model};
}

void sanafe::description_parse_dendrite_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef dendrite_node, Core &parent_core)
{
    auto dendrite_name = description_required_field<std::string>(
            parser, dendrite_node, "name");
    const ryml::ConstNodeRef attributes =
            dendrite_node.find_child("attributes");
    if (attributes.invalid())
    {
        throw DescriptionParsingError(
                "No attributes section defined", parser, dendrite_node);
    }

    ModelInfo model;
    model.name = description_required_field<std::string>(
            parser, attributes, "model");
    const ryml::ConstNodeRef plugin_path_node = attributes.find_child("plugin");
    if (!plugin_path_node.invalid())
    {
        if (plugin_path_node.is_val())
        {
            std::string plugin_path;
            plugin_path_node >> plugin_path;
            model.plugin_library_path = plugin_path;
        }
        else
        {
            throw DescriptionParsingError("Expected plugin path to be string",
                    parser, plugin_path_node);
        }
    }

    DendritePowerMetrics power_metrics;
    power_metrics.energy_access = description_required_field<double>(
            parser, attributes, "energy_access");
    power_metrics.latency_access = description_required_field<double>(
            parser, attributes, "latency_access");
    parent_core.create_dendrite(dendrite_name, power_metrics, model);
}

void sanafe::description_parse_soma_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef soma_node, Core &parent_core)
{
    auto soma_name =
            description_required_field<std::string>(parser, soma_node, "name");
    const ryml::ConstNodeRef attributes = soma_node.find_child("attributes");
    if (attributes.invalid())
    {
        throw DescriptionParsingError(
                "No attributes section defined", parser, soma_node);
    }
    std::string model_str;

    ModelInfo model;
    model.name = description_required_field<std::string>(
            parser, attributes, "model");
    const ryml::ConstNodeRef plugin_path_node = attributes.find_child("plugin");
    if (!plugin_path_node.invalid())
    {
        if (plugin_path_node.is_val())
        {
            std::string plugin_path;
            plugin_path_node >> plugin_path;
            model.plugin_library_path = plugin_path;
        }
        else
        {
            throw DescriptionParsingError("Expected plugin path to be string",
                    parser, plugin_path_node);
        }
    }

    SomaPowerMetrics power_metrics;
    power_metrics.energy_update_neuron = description_required_field<double>(
            parser, attributes, "energy_update_neuron");
    power_metrics.latency_update_neuron = description_required_field<double>(
            parser, attributes, "latency_update_neuron");
    power_metrics.energy_access_neuron = description_required_field<double>(
            parser, attributes, "energy_access_neuron");
    power_metrics.latency_access_neuron = description_required_field<double>(
            parser, attributes, "latency_access_neuron");
    power_metrics.energy_spike_out = description_required_field<double>(
            parser, attributes, "energy_spike_out");
    power_metrics.latency_spike_out = description_required_field<double>(
            parser, attributes, "latency_spike_out");

    if (!attributes.find_child("noise").invalid())
    {
        // TODO: support optioanl noise arg again alongside the plugin mechanism
        /*
        s.noise_type = NOISE_FILE_STREAM;
        s.noise_stream = fopen(value_str.c_str(), "r");
        TRACE1("Opening noise str: %s\n", value_str.c_str());
        if (s.noise_stream == NULL)
        {
            INFO("Error: Failed to open noise stream: %s.\n",
                value_str.c_str());
            exit(1);
        }
        */
    }
    parent_core.create_soma(std::move(soma_name), power_metrics, model);
}

void sanafe::description_parse_axon_out_section(const ryml::Parser &parser,
        const ryml::ConstNodeRef axon_out_node, Core &parent_core)
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
        Core &core = arch.create_core(name, parent_tile_id, pipeline_config);

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
    CorePipelineConfiguration pipeline_config;

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
        Tile &new_tile = arch.create_tile(name, power_metrics);
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

sanafe::Network sanafe::description_parse_network_file_yaml(
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
        Network net = description_parse_network_section_yaml(
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

sanafe::Network sanafe::description_parse_network_section_yaml(
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

    Network new_net(std::move(net_name));
    description_parse_neuron_group_section_yaml(
            parser, net_node["groups"], new_net);
    description_parse_edges_section_yaml(parser, net_node["edges"], new_net);

    return new_net;
}

void sanafe::description_parse_neuron_group_section_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef groups_node,
        Network &net)
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
        const ryml::ConstNodeRef edges_node, Network &net)
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
        const ryml::ConstNodeRef neuron_group_node, Network &net)
{
    const auto group_name = description_required_field<std::string>(
            parser, neuron_group_node, "name");
    INFO("Parsing neuron group: %s\n", group_name.c_str());

    const auto &neurons_node = neuron_group_node["neurons"];
    if (neuron_group_node.find_child("size").invalid())
    {
        throw DescriptionParsingError(
                "Neuron group 'size' (aka neuron count) not given", parser,
                neuron_group_node);
    }
    size_t neuron_count;
    neuron_group_node["size"] >> neuron_count;

    NeuronTemplate default_neuron_config{};
    if (!neuron_group_node.find_child("attributes").invalid())
    {
        TRACE1("Parsing neuron group attributes\n");
        default_neuron_config = description_parse_neuron_attributes_yaml(
                parser, neuron_group_node["attributes"]);
    }
    NeuronGroup &group = net.create_neuron_group(
            group_name, neuron_count, default_neuron_config);
    TRACE1("Parsing neuron section\n");
    description_parse_neuron_section_yaml(parser, neurons_node, group);
}

void sanafe::description_parse_neuron_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef neuron_node, NeuronGroup &neuron_group)
{
    if (neuron_node.is_seq())
    {
        for (const auto list_entry : neuron_node)
        {
            for (const auto neuron_description : list_entry)
            {
                // There should only be one mapping per list entry
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
                "Invalid neuron format", parser, neuron_node);
    }
}

void sanafe::description_parse_neuron(const std::string &id,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes,
        NeuronGroup &neuron_group)
{
    std::pair<size_t, size_t> range;
    TRACE1("Parsing neuron(s): %s\n", id.c_str());
    const NeuronTemplate config = description_parse_neuron_attributes_yaml(
            parser, attributes, neuron_group.default_neuron_config);
    const bool is_range = (id.find("..") != std::string::npos);
    if (is_range)
    {
        range = description_parse_range_yaml(id);
        for (size_t instance = range.first; instance <= range.second;
                ++instance)
        {
            Neuron &n = neuron_group.neurons[instance];
            n.set_attributes(config);
        }
    }
    else
    {
        const size_t nid = std::stoull(id);
        Neuron &n = neuron_group.neurons[nid];
        n.set_attributes(config);
    }
}

sanafe::NeuronTemplate sanafe::description_parse_neuron_attributes_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes,
        const NeuronTemplate &default_template)
{
    NeuronTemplate neuron_template = default_template;

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
        attributes["log_potential"] >> neuron_template.log_potential;
    }
    if (!attributes.find_child("log_spikes").invalid())
    {
        attributes["log_spikes"] >> neuron_template.log_spikes;
    }
    if (!attributes.find_child("force_update").invalid())
    {
        attributes["force_update"] >> neuron_template.force_update;
    }
    if (!attributes.find_child("synapse_hw_name").invalid())
    {
        attributes["synapse_hw_name"] >>
                neuron_template.default_synapse_hw_name;
    }
    if (!attributes.find_child("soma_hw_name").invalid())
    {
        attributes["soma_hw_name"] >> neuron_template.soma_hw_name;
    }

    // Parse and add shared parameters, which are defined alongside attributes
    auto model_params =
            description_parse_model_parameters_yaml(parser, attributes);
    for (auto &[key, parameter] : model_params)
    {
        neuron_template.dendrite_model_params.insert({key, parameter});
        neuron_template.soma_model_params.insert({key, parameter});
    }
    // Parse and add unit specific model parameters defined under 'dendrite' or
    //  'soma' keys
    if (!attributes.find_child("dendrite").invalid())
    {
        auto dendrite_params = description_parse_model_parameters_yaml(
                parser, attributes["dendrite"]);
        neuron_template.dendrite_model_params.insert(
                dendrite_params.begin(), dendrite_params.end());
    }
    if (!attributes.find_child("soma").invalid())
    {
        auto soma_params = description_parse_model_parameters_yaml(
                parser, attributes["soma"]);
        neuron_template.soma_model_params.insert(
                soma_params.begin(), soma_params.end());
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
    if ((source_dot_pos == std::string::npos) ||
            (target_dot_pos == std::string::npos))
    {
        throw std::runtime_error(
                "Edge is not formatted correctly: " + std::string(description));
    }

    NeuronAddress source;
    source.group_name = source_part.substr(0, source_dot_pos);
    source.neuron_id =
            std::stoull(std::string(source_part.substr(source_dot_pos + 1)));

    NeuronAddress target;
    target.group_name = target_part.substr(0, target_dot_pos);
    target.neuron_id =
            std::stoull(std::string(target_part.substr(target_dot_pos + 1)));

    return std::make_tuple(source, target);
}

void sanafe::description_parse_edge(const std::string &description,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node,
        Network &net)
{
    // Description has format src_group.src_neuron -> tgt_group.tgt_neuron
    auto [source_address, target_address] =
            description_parse_edge_description(description);

    if (net.groups.find(source_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid source neuron group:" + source_address.group_name;
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    NeuronGroup &src_group = net.groups.at(source_address.group_name);
    if (source_address.neuron_id > src_group.neurons.size())
    {
        const std::string error =
                "Invalid source neuron id: " + source_address.group_name + "." +
                std::to_string(source_address.neuron_id);
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    Neuron &src_neuron = src_group.neurons[source_address.neuron_id];

    if (net.groups.find(target_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid target neuron group:" + target_address.group_name;
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    NeuronGroup &dst_group = net.groups.at(target_address.group_name);

    if (target_address.neuron_id > dst_group.neurons.size())
    {
        const std::string error =
                "Invalid target neuron id: " + target_address.group_name + "." +
                std::to_string(target_address.neuron_id);
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    Neuron &dst_neuron = dst_group.neurons.at(target_address.neuron_id);

    Connection &con = src_neuron.connect_to_neuron(dst_neuron);
    // TODO: create a ConnectionTemplate? and pass this when connecting two
    //  neurons instead. The problem with using a reference is when the vector
    //  of connections expands, any references may become invalidated...
    description_parse_edge_attributes(con, parser, attributes_node);
}

void sanafe::description_parse_edge_attributes(Connection &con,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node)
{
    if (attributes_node.is_seq())
    {
        for (const auto &attribute : attributes_node)
        {
            description_parse_edge_attributes(con, parser, attribute);
        }

        return;
    }

    if (!attributes_node.find_child("synapse").invalid())
    {
        con.synapse_params = description_parse_model_parameters_yaml(
                parser, attributes_node["synapse"]);
    }
    if (!attributes_node.find_child("dendrite").invalid())
    {
        con.dendrite_params = description_parse_model_parameters_yaml(
                parser, attributes_node["dendrite"]);
    }

    const auto shared_model_params =
            description_parse_model_parameters_yaml(parser, attributes_node);
    for (const auto &[key, parameter] : shared_model_params)
    {
        if ((key != "synapse") && (key != "dendrite") && (key != "soma"))
        {
            con.synapse_params.insert({key, parameter});
            con.dendrite_params.insert({key, parameter});
        }
    }
}

void sanafe::description_parse_mapping_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef mappings_node, Architecture &arch,
        Network &net)
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
                    "<group>.<neuron>: {core: <tile>.<core>}",
                    parser, mapping);
        }

        for (const auto &entry : mapping)
        {
            std::string neuron_address;
            entry >> ryml::key(neuron_address);
            const auto dot_pos = neuron_address.find('.');

            const std::string group_name = neuron_address.substr(0, dot_pos);
            const size_t neuron_id =
                    std::stoull(neuron_address.substr(dot_pos + 1));
            if (net.groups.find(group_name) == net.groups.end())
            {
                const std::string error = "Invalid neuron group:" + group_name;
                throw DescriptionParsingError(error, parser, mapping);
            }
            //INFO("Mapping neuron: %s.%zu\n", group_name.c_str(), neuron_id);
            NeuronGroup &group = net.groups.at(group_name);
            if (neuron_id >= group.neurons.size())
            {
                std::string error = "Invalid neuron id: ";
                error += group_name;
                error += '.';
                error += neuron_id;
                throw DescriptionParsingError(error, parser, mapping);
            }
            Neuron &neuron = group.neurons[neuron_id];

            description_parse_mapping(parser, neuron, entry, arch);
        }
    }
}

void sanafe::description_parse_mapping(const ryml::Parser &parser,
        Neuron &neuron, const ryml::ConstNodeRef mapping_info,
        Architecture &arch)
{
    if (!mapping_info.is_map())
    {
        throw DescriptionParsingError(
                "Invalid mapping: mappings must be defined using a YAML map "
                "<group>.<neuron>: core: <tile>.<core>",
                parser, mapping_info);
    }

    const auto core_address = description_required_field<std::string>(
            parser, mapping_info, "core");
    const auto dot_pos = core_address.find('.');

    const size_t tile_id = std::stoull(core_address.substr(0, dot_pos));
    const size_t core_offset_within_tile =
            std::stoull(core_address.substr(dot_pos + 1));

    if (tile_id >= arch.tiles.size())
    {
        throw DescriptionParsingError(
                "Tile ID >= tile count", parser, mapping_info);
    }
    Tile &tile = arch.tiles[tile_id];
    if (core_offset_within_tile >= tile.cores.size())
    {
        throw DescriptionParsingError(
                "Core ID >= core count", parser, mapping_info);
    }
    Core &core = tile.cores[core_offset_within_tile];
    core.map_neuron(neuron);
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
    TRACE1("Range: %zu to %zu\n", first, last);
    if (first > last)
    {
        throw std::runtime_error("Invalid range; first > last");
    }

    return std::make_pair(first, last);
}
