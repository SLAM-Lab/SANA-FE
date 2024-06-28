#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric> // For std::accumulate
#include <optional>
#include <string>
#include <string_view>
//#include <unordered_set>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "arch.hpp"
#include "description.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"

sanafe::DescriptionParsingError::DescriptionParsingError(
        const std::string &error, const YAML::Mark &pos)
        : std::invalid_argument(error)
        , message("Error: " + error + " (Line " + std::to_string(pos.line + 1) +
                  ':' + std::to_string(pos.column + 1) + ").\n")
{
}

const char *sanafe::DescriptionParsingError::what() const noexcept
{
    return message.c_str();
}

void sanafe::check_key(const YAML::Node &node, const std::string &key)
{
    if (!node.IsMap())
    {
        throw DescriptionParsingError(
                "Node should be a mapping\n. For more info on YAML mappings "
                "refer to the YAML 1.2 specification, 7.4.2 'Flow Mappings' "
                "and 8.2.2 'Block Mappings'",
                node.Mark());
    }
    const YAML::Node &child = node[key];
    if (!child.IsDefined() || child.IsNull())
    {
        const std::string message = "Value for key '" + key + "' not defined";
        throw DescriptionParsingError(message, node.Mark());
    }
}

template <typename T>
T sanafe::description_required_field(
        const YAML::Node &node, const std::string &key)
{
    // Wrapper around YAML-CPP for field=map[key], adding more error prints
    check_key(node, key);
    const YAML::Node &field_node = node[key];
    if (!field_node.IsScalar())
    {
        const std::string message = "'" + key + "' value should be a scalar";
        throw DescriptionParsingError(message, field_node.Mark());
    }

    // Efficiently convert to type T by trying the YAML-CPP decoder.
    //  If decode() fails, it returns false and execution falls through
    T field;
    if (YAML::convert<T>::decode(field_node, field))
    {
        return field; // type T
    }

    const std::string message = "Could not cast field '" +
            YAML::Dump(field_node) + "' (key '" + key +
            "') to type: " + description_get_type_string(typeid(field));
    throw DescriptionParsingError(message, field_node.Mark());
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

void sanafe::description_parse_axon_in_section_yaml(
        const YAML::Node &axon_in_node, Core &parent_core)
{
    auto name = description_required_field<std::string>(axon_in_node, "name");
    check_key(axon_in_node, "attributes");
    const AxonInPowerMetrics in_metrics =
            description_parse_axon_in_attributes_yaml(
                    axon_in_node["attributes"]);
    parent_core.create_axon_in(name, in_metrics);
}

sanafe::AxonInPowerMetrics sanafe::description_parse_axon_in_attributes_yaml(
        const YAML::Node &attributes)
{
    AxonInPowerMetrics axon_in_metrics;
    axon_in_metrics.energy_message_in =
            description_required_field<double>(attributes, "energy_message_in");
    axon_in_metrics.latency_message_in = description_required_field<double>(
            attributes, "latency_message_in");

    return axon_in_metrics;
}

void sanafe::description_parse_synapse_section_yaml(
        const YAML::Node &synapse_node, Core &parent_core)
{
    auto name = description_required_field<std::string>(synapse_node, "name");
    check_key(synapse_node, "attributes");
    auto [power_metrics, model] = description_parse_synapse_attributes_yaml(
            synapse_node["attributes"]);
    parent_core.create_synapse(name, power_metrics, model);
}

std::pair<sanafe::SynapsePowerMetrics, sanafe::ModelInfo>
sanafe::description_parse_synapse_attributes_yaml(const YAML::Node &attributes)
{
    ModelInfo model;
    model.name = description_required_field<std::string>(attributes, "model");
    if (const YAML::Node &plugin_path_node = attributes["plugin"])
    {
        if (plugin_path_node.IsScalar())
        {
            model.plugin_library_path = plugin_path_node.as<std::string>();
        }
        else
        {
            throw DescriptionParsingError("Expected plugin path to be string",
                    plugin_path_node.Mark());
        }
    }

    SynapsePowerMetrics power_metrics;
    power_metrics.energy_process_spike = description_required_field<double>(
            attributes, "energy_process_spike");
    power_metrics.latency_process_spike = description_required_field<double>(
            attributes, "latency_process_spike");

    return {power_metrics, model};
}

void sanafe::description_parse_dendrite_section_yaml(
        const YAML::Node &dendrite_node, Core &parent_core)
{
    auto dendrite_name =
            description_required_field<std::string>(dendrite_node, "name");
    check_key(dendrite_node, "attributes");
    const YAML::Node &attributes = dendrite_node["attributes"];

    ModelInfo model;
    model.name = description_required_field<std::string>(attributes, "model");
    if (const YAML::Node &plugin_path_node = attributes["plugin"])
    {
        if (plugin_path_node.IsScalar())
        {
            model.plugin_library_path = plugin_path_node.as<std::string>();
        }
        else
        {
            throw DescriptionParsingError("Expected plugin path to be string",
                    plugin_path_node.Mark());
        }
    }

    DendritePowerMetrics power_metrics;
    power_metrics.energy_access =
            description_required_field<double>(attributes, "energy_access");
    power_metrics.latency_access =
            description_required_field<double>(attributes, "latency_access");
    parent_core.create_dendrite(dendrite_name, power_metrics, model);
}

void sanafe::description_parse_soma_section_yaml(
        const YAML::Node &soma_node, Core &parent_core)
{
    auto soma_name = description_required_field<std::string>(soma_node, "name");
    const YAML::Node &attributes = soma_node["attributes"];
    std::string model_str;

    ModelInfo model;
    model.name = description_required_field<std::string>(attributes, "model");
    if (const YAML::Node &plugin_path_node = attributes["plugin"])
    {
        if (plugin_path_node.IsScalar())
        {
            model.plugin_library_path = plugin_path_node.as<std::string>();
        }
        else
        {
            throw DescriptionParsingError("Expected plugin path to be string",
                    plugin_path_node.Mark());
        }
    }

    SomaPowerMetrics power_metrics;
    power_metrics.energy_update_neuron = description_required_field<double>(
            attributes, "energy_update_neuron");
    power_metrics.latency_update_neuron = description_required_field<double>(
            attributes, "latency_update_neuron");
    power_metrics.energy_access_neuron = description_required_field<double>(
            attributes, "energy_access_neuron");
    power_metrics.latency_access_neuron = description_required_field<double>(
            attributes, "latency_access_neuron");
    power_metrics.energy_spike_out =
            description_required_field<double>(attributes, "energy_spike_out");
    power_metrics.latency_spike_out =
            description_required_field<double>(attributes, "latency_spike_out");

    if (attributes["noise"].IsDefined())
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

void sanafe::description_parse_axon_out_section(
        const YAML::Node &axon_out_node, Core &parent_core)
{
    auto axon_out_name =
            description_required_field<std::string>(axon_out_node, "name");

    check_key(axon_out_node, "attributes");
    const auto &attributes = axon_out_node["attributes"];
    AxonOutPowerMetrics power_metrics;
    power_metrics.energy_message_out = description_required_field<double>(
            attributes, "energy_message_out");
    power_metrics.latency_message_out = description_required_field<double>(
            attributes, "latency_message_out");

    parent_core.create_axon_out(axon_out_name, power_metrics);
}

void sanafe::description_parse_core_section_yaml(const YAML::Node &core_node,
        const size_t parent_tile_id, Architecture &arch)
{
    auto core_name = description_required_field<std::string>(core_node, "name");
    std::pair<int, int> core_range = {0, 0};
    if (core_name.find("..") != std::string::npos)
    {
        core_range = description_parse_range_yaml(
                core_name, core_node["name"].Mark());
    }

    for (int c = core_range.first; c <= core_range.second; c++)
    {
        const std::string name = core_name.substr(0, core_name.find('[')) +
                '[' + std::to_string(c) + ']';
        const CorePipelineConfiguration pipeline_config =
                description_parse_core_pipeline_yaml(core_node["attributes"]);
        Core &core = arch.create_core(name, parent_tile_id, pipeline_config);

        if (const YAML::Node &axon_in_node = core_node["axon_in"])
        {
            if (axon_in_node.IsSequence())
            {
                for (const auto &axon : axon_in_node)
                {
                    description_parse_axon_in_section_yaml(axon, core);
                }
            }
            else
            {
                description_parse_axon_in_section_yaml(axon_in_node, core);
            }
        }
        else
        {
            const std::string error = "No axon in section defined";
            throw DescriptionParsingError(error, core_node.Mark());
        }

        if (const YAML::Node &synapses = core_node["synapse"])
        {
            if (synapses.IsSequence())
            {
                for (const auto &syn : synapses)
                {
                    description_parse_synapse_section_yaml(syn, core);
                }
            }
            else
            {
                description_parse_synapse_section_yaml(synapses, core);
            }
        }
        else
        {
            const std::string error = "No synapse section defined";
            throw DescriptionParsingError(error, core_node.Mark());
        }
        if (const YAML::Node &dendrite_node = core_node["dendrite"])
        {
            if (dendrite_node.IsSequence())
            {
                for (const auto &dendrite : dendrite_node)
                {
                    description_parse_dendrite_section_yaml(dendrite, core);
                }
            }
            else
            {
                description_parse_dendrite_section_yaml(dendrite_node, core);
            }
        }
        else
        {
            const std::string error = "No dendrite section defined";
            throw DescriptionParsingError(error, core_node.Mark());
        }

        if (const YAML::Node &soma_node = core_node["soma"])
        {
            if (soma_node.IsSequence())
            {
                for (const auto &soma : soma_node)
                {
                    description_parse_soma_section_yaml(soma, core);
                }
            }
            else
            {
                description_parse_soma_section_yaml(soma_node, core);
            }
        }
        else
        {
            const std::string error = "No soma section defined";
            throw DescriptionParsingError(error, core_node.Mark());
        }

        if (const YAML::Node &axon_out_node = core_node["axon_out"])
        {
            if (axon_out_node.IsSequence())
            {
                for (const auto &axon : axon_out_node)
                {
                    description_parse_axon_out_section(axon, core);
                }
            }
            else
            {
                description_parse_axon_out_section(axon_out_node, core);
            }
        }
        else
        {
            const std::string error = "No axon out seciont defined";
            throw DescriptionParsingError(error, core_node.Mark());
        }
    }
}

sanafe::CorePipelineConfiguration sanafe::description_parse_core_pipeline_yaml(
        const YAML::Node &attributes)
{
    CorePipelineConfiguration pipeline_config;

    pipeline_config.buffer_position = pipeline_parse_buffer_pos_str(
            description_required_field<std::string>(
                    attributes, "buffer_position"));
    pipeline_config.max_neurons_supported = description_required_field<int>(
            attributes, "max_neurons_supported");

    return pipeline_config;
}

sanafe::TilePowerMetrics sanafe::description_parse_tile_metrics_yaml(
        const YAML::Node &attributes)
{
    TilePowerMetrics tile_metrics;

    tile_metrics.energy_north_hop =
            description_required_field<double>(attributes, "energy_north_hop");
    tile_metrics.latency_north_hop =
            description_required_field<double>(attributes, "latency_north_hop");

    tile_metrics.energy_east_hop =
            description_required_field<double>(attributes, "energy_east_hop");
    tile_metrics.latency_east_hop =
            description_required_field<double>(attributes, "latency_east_hop");

    tile_metrics.energy_south_hop =
            description_required_field<double>(attributes, "energy_south_hop");
    tile_metrics.latency_south_hop =
            description_required_field<double>(attributes, "latency_south_hop");

    tile_metrics.energy_west_hop =
            description_required_field<double>(attributes, "energy_west_hop");
    tile_metrics.latency_west_hop =
            description_required_field<double>(attributes, "latency_west_hop");

    return tile_metrics;
}

void sanafe::description_parse_tile_section_yaml(
        const YAML::Node &tile_node, Architecture &arch)
{
    auto tile_name = tile_node["name"].as<std::string>();
    std::pair<int, int> range = {0, 0};

    if (tile_name.find("..") != std::string::npos)
    {
        range = description_parse_range_yaml(
                tile_name, tile_node["name"].Mark());
    }

    for (int t = range.first; t <= range.second; t++)
    {
        const std::string name = tile_name.substr(0, tile_name.find('[')) +
                "[" + std::to_string(t) + "]";
        const TilePowerMetrics power_metrics =
                description_parse_tile_metrics_yaml(tile_node["attributes"]);
        Tile &new_tile = arch.create_tile(name, power_metrics);
        if (const YAML::Node &cores = tile_node["core"])
        {
            const bool is_list_of_cores = cores.IsSequence();
            if (is_list_of_cores)
            {
                for (const auto &core : cores)
                {
                    description_parse_core_section_yaml(
                            core, new_tile.id, arch);
                }
            }
            else // Is a single core
            {
                description_parse_core_section_yaml(cores, new_tile.id, arch);
            }
        }
        else
        {
            const std::string error = "No core section defined";
            throw DescriptionParsingError(error, tile_node.Mark());
        }
    }
}

sanafe::NetworkOnChipConfiguration
sanafe::description_parse_noc_configuration_yaml(
        const YAML::Node &noc_attributes)
{
    NetworkOnChipConfiguration noc;
    noc.width_in_tiles =
            description_required_field<int>(noc_attributes, "width");
    noc.height_in_tiles =
            description_required_field<int>(noc_attributes, "height");
    noc.link_buffer_size =
            description_required_field<int>(noc_attributes, "link_buffer_size");

    return noc;
}

sanafe::Architecture sanafe::description_parse_arch_section_yaml(
        const YAML::Node &arch_node)
{
    const auto arch_name = arch_node["name"].as<std::string>();
    if (arch_name.find('[') != std::string::npos)
    {
        throw DescriptionParsingError(
                "Multiple architectures not supported", arch_node.Mark());
    }
    NetworkOnChipConfiguration noc =
            description_parse_noc_configuration_yaml(arch_node["attributes"]);
    Architecture new_arch(arch_name, noc);
    if (const YAML::Node &tiles = arch_node["tile"])
    {
        const bool is_list_of_tiles = tiles.IsSequence();
        if (is_list_of_tiles)
        {
            for (const auto &tile : tiles)
            {
                description_parse_tile_section_yaml(tile, new_arch);
            }
        }
        else // Only one tile defined
        {
            description_parse_tile_section_yaml(tiles, new_arch);
        }
    }
    else
    {
        throw DescriptionParsingError("No tile section defined", tiles.Mark());
    }

    return new_arch;
}

sanafe::Architecture sanafe::description_parse_arch_file_yaml(std::ifstream &fp)
{
    YAML::Node top_level_yaml_node = YAML::Load(fp);
    INFO("YAML information loaded from file.\n");
    if (const YAML::Node &arch_yaml_node = top_level_yaml_node["architecture"])
    {
        return description_parse_arch_section_yaml(arch_yaml_node);
    }
    throw DescriptionParsingError("No top-level architecture section defined",
            top_level_yaml_node.Mark());
}

sanafe::Network sanafe::description_parse_network_file_yaml(std::ifstream &fp)
{
    const YAML::Node yaml_node = YAML::Load(fp);
    INFO("YAML information loaded from file.\n");
    if (yaml_node.IsMap() && (yaml_node["network"].IsDefined()))
    {
        return description_parse_network_section_yaml(yaml_node["network"]);
    }
    throw DescriptionParsingError(
            "No top-level network section defined", yaml_node.Mark());
}

sanafe::Network sanafe::description_parse_network_section_yaml(
        const YAML::Node &net_node)
{
    std::string net_name{};
    // TODO: refactor?
    if (net_node["name"].IsDefined())
    {
        net_name = net_node["name"].as<std::string>();
        if (net_name.find('[') != std::string::npos)
        {
            throw DescriptionParsingError(
                    "Multiple networks not supported", net_node.Mark());
        }
    }
    else
    {
        INFO("Warning: No network name given; leaving name empty.\n");
    }

    Network new_net(std::move(net_name));
    check_key(net_node, "groups");
    check_key(net_node, "edges");
    description_parse_neuron_group_section_yaml(net_node["groups"], new_net);
    description_parse_edges_section_yaml(net_node["edges"], new_net);

    return new_net;
}

void sanafe::description_parse_neuron_group_section_yaml(
        const YAML::Node &groups_node, Network &net)
{
    if (groups_node.IsSequence())
    {
        for (const auto &group : groups_node)
        {
            description_parse_group(group, net);
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Neuron group section does not define a list of groups",
                groups_node.Mark());
    }
}

void sanafe::description_parse_edges_section_yaml(
        const YAML::Node &edges_node, Network &net)
{
    if (edges_node.IsMap())
    {
        for (const auto &edge : edges_node)
        {
            const auto &edge_description = edge.first.as<std::string>();
            //std::cout << "edge:" << edge_description << ':' << edge.second << '\n';
            const YAML::Node &edge_attributes = edge.second;
            description_parse_edge_short_format(
                    edge_description, edge_attributes, net);
        }
    }
    else if (edges_node.IsSequence())
    {
        for (const auto &edge : edges_node)
        {
            description_parse_edge_long_format(edge, net);
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Edges section does not define a list or mapping of edges",
                edges_node.Mark());
    }
}

void sanafe::description_parse_group(
        const YAML::Node &neuron_group_node, Network &net)
{
    const auto group_name =
            description_required_field<std::string>(neuron_group_node, "name");
    check_key(neuron_group_node, "neurons");
    const auto &neurons_node = neuron_group_node["neurons"];

    NeuronTemplate default_neuron_config{};
    if (const YAML::Node &group_attributes = neuron_group_node["attributes"])
    {
        default_neuron_config =
                description_parse_neuron_attributes_yaml(group_attributes);
    }
    NeuronGroup &group =
            net.create_neuron_group(group_name, default_neuron_config);
    description_parse_neuron_section_yaml(neurons_node, group);
}

void sanafe::description_parse_neuron_section_yaml(
        const YAML::Node &neuron_node, NeuronGroup &neuron_group)
{
    if (neuron_node.IsMap())
    {
        for (const auto &it : neuron_node)
        {
            const auto id = it.first.as<std::string>();
            description_parse_neuron_short_format(id, it.second, neuron_group);
        }
    }
    else if (neuron_node.IsSequence())
    {
        for (const auto &node : neuron_node)
        {
            description_parse_neuron_long_format(node, neuron_group);
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Invalid neuron format", neuron_node.Mark());
    }
}

void sanafe::description_parse_neuron_long_format(
        const YAML::Node &neuron_node, NeuronGroup &neuron_group)
{
    const auto &id =
            description_required_field<std::string>(neuron_node, "id");
    std::pair<size_t, size_t> range = {0, 0};
    if (id.find("..") != std::string::npos)
    {
        range = description_parse_range_yaml(
                id, neuron_node["name"].Mark());
    }
    check_key(neuron_node, "attributes");
    const auto &attributes = neuron_node["attributes"];
    // TODO: parse default config
    const NeuronTemplate config = description_parse_neuron_attributes_yaml(
            attributes, neuron_group.default_neuron_config);

    for (size_t instance = range.first; instance <= range.second; ++instance)
    {
        std::string instance_name = id;
        const bool range_given = (range.second > range.first);
        if (range_given)
        {
            const bool bracket_notation =
                    (id.find('[') != std::string::npos);
            if (bracket_notation)  // E.g., n[123] for instance 123
            {
                instance_name = id.substr(0, id.find('[')) +
                        '[' + std::to_string(instance) + ']';
            }
            else  // Special case where id is only integer value in the range
            {
                instance_name = std::to_string(instance);
            }
        }
        neuron_group.create_neuron(instance_name, config);
    }
}

void sanafe::description_parse_neuron_short_format(const std::string &id,
        const YAML::Node &attributes, NeuronGroup &neuron_group)
{
    std::pair<size_t, size_t> range = {0, 0};
    if (id.find("..") != std::string::npos)
    {
        range = description_parse_range_yaml(id, attributes.Mark());
    }
    const NeuronTemplate config = description_parse_neuron_attributes_yaml(
            attributes, neuron_group.default_neuron_config);
    for (size_t instance = range.first; instance <= range.second; ++instance)
    {
        std::string instance_name = id;
        const bool range_given = (range.second > range.first);
        if (range_given)
        {
            const bool bracket_notation =
                    (id.find('[') != std::string::npos);
            if (bracket_notation)  // E.g., n[123] for instance 123
            {
                instance_name = id.substr(0, id.find('[')) +
                        '[' + std::to_string(instance) + ']';
            }
            else  // Special case where id is only integer value in the range
            {
                instance_name = std::to_string(instance);
            }
        }

        neuron_group.create_neuron(instance_name, config);
    }
}

sanafe::NeuronTemplate sanafe::description_parse_neuron_attributes_yaml(
        const YAML::Node &attributes, const NeuronTemplate &default_template)
{
    NeuronTemplate neuron_template = default_template;

    if (attributes.IsSequence())
    {
        // Ordered list format, recursively parse attributes for each element
        return std::accumulate(attributes.begin(), attributes.end(),
                neuron_template,
                [](const NeuronTemplate &config, const YAML::Node &yaml) {
                    return description_parse_neuron_attributes_yaml(
                            yaml, config);
                });
    }

    if (attributes["log_potential"].IsDefined())
    {
        neuron_template.log_potential = attributes["log_potential"].as<bool>();
    }
    if (attributes["log_spikes"].IsDefined())
    {
        neuron_template.log_spikes = attributes["log_spikes"].as<bool>();
    }
    if (attributes["force_update"].IsDefined())
    {
        neuron_template.force_update = attributes["force_update"].as<bool>();
    }
    if (attributes["synapse_hw_name"].IsDefined())
    {
        neuron_template.default_synapse_hw_name =
                attributes["synapse_hw_name"].as<std::string>();
    }
    if (attributes["soma_hw_name"].IsDefined())
    {
        neuron_template.soma_hw_name =
                attributes["soma_hw_name"].as<std::string>();
    }

    // Parse and add shared parameters, which are defined alongside attributes
    auto model_params = description_parse_model_parameters_yaml(attributes);
    for (auto &[key, parameter] : model_params)
    {
        neuron_template.dendrite_model_params.insert({key, parameter});
        neuron_template.soma_model_params.insert({key, parameter});
    }
    // Parse and add unit specific model parameters defined under 'dendrite' or
    //  'soma' keys
    if (attributes["dendrite"].IsDefined())
    {
        auto dendrite_params =
                description_parse_model_parameters_yaml(attributes["dendrite"]);
        neuron_template.dendrite_model_params.insert(
                dendrite_params.begin(), dendrite_params.end());
    }
    if (attributes["soma"].IsDefined())
    {
        auto soma_params =
                description_parse_model_parameters_yaml(attributes["soma"]);
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
        throw DescriptionParsingError(
                "Edge is not formatted correctly.", YAML::Mark());
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
        throw DescriptionParsingError(
                "Edge is not formatted correctly.", YAML::Mark());
    }

    NeuronAddress source;
    source.group_name = source_part.substr(0, source_dot_pos);
    source.neuron_id = source_part.substr(source_dot_pos + 1);

    NeuronAddress target;
    target.group_name = target_part.substr(0, target_dot_pos);
    target.neuron_id = target_part.substr(target_dot_pos + 1);

    return std::make_tuple(source, target);
}

void sanafe::description_parse_edge_short_format(const std::string &description,
        const YAML::Node &attributes_node, Network &net)
{
    if (attributes_node.IsSequence())
    {
        for (const YAML::Node &attr : attributes_node)
        {
            description_parse_edge_short_format(description, attr, net);
        }
    }

    // Description has format src_group.src_neuron -> tgt_group.tgt_neuron
    auto [source_address, target_address] =
            description_parse_edge_description(description);

    // TODO: refactor - all of the below code is the same as the long format
    if (net.groups.find(source_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid source neuron group:" + source_address.group_name;
        throw DescriptionParsingError(error, attributes_node.Mark());
    }
    NeuronGroup &src_group = net.groups.at(source_address.group_name);
    if (src_group.neurons.find(source_address.neuron_id) ==
            src_group.neurons.end())
    {
        const std::string error = "Invalid source neuron id: " +
                source_address.group_name + "." + source_address.neuron_id;
        throw DescriptionParsingError(error, attributes_node.Mark());
    }
    Neuron &src_neuron = src_group.neurons.at(source_address.neuron_id);

    if (net.groups.find(target_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid target neuron group:" + target_address.group_name;
        throw DescriptionParsingError(error, attributes_node.Mark());
    }
    NeuronGroup &dst_group = net.groups.at(target_address.group_name);

    if (dst_group.neurons.find(target_address.neuron_id) ==
            dst_group.neurons.end())
    {
        const std::string error = "Invalid target neuron id: " +
                target_address.group_name + "." + target_address.neuron_id;
        throw DescriptionParsingError(error, attributes_node.Mark());
    }
    Neuron &dst_neuron = dst_group.neurons.at(target_address.neuron_id);

    // TODO: fix this space wastage
    std::map<std::string, ModelParam> synapse_params{};
    std::map<std::string, ModelParam> dendrite_params{};

    if (const YAML::Node &synapse_node = attributes_node["synapse"])
    {
        synapse_params =
                description_parse_model_parameters_yaml(synapse_node);
    }
    if (const YAML::Node &dendrite_node = attributes_node["dendrite"])
    {
        dendrite_params =
                description_parse_model_parameters_yaml(dendrite_node);
    }

    const auto shared_model_params =
            description_parse_model_parameters_yaml(attributes_node);
    for (const auto &[key, parameter] : shared_model_params)
    {
        //static const std::unordered_set<std::string> unit_specific_keys = {
        //    "synapse", "dendrite", "soma"};
        //if (unit_specific_keys.find(key) != unit_specific_keys.end())
        if ((key != "synapse") && (key != "dendrite") && (key != "soma"))
        {
            synapse_params.insert({key, parameter});
            dendrite_params.insert({key, parameter});
        }
    }
    std::optional<std::string> synapse_hw_name;
    src_neuron.connect_to_neuron(
            dst_neuron, synapse_params, dendrite_params, synapse_hw_name);
}

void sanafe::description_parse_edge_long_format(
        const YAML::Node &edge_node, Network &net)
{
    const YAML::Node &src_node = edge_node["src"];
    const YAML::Node &dst_node = edge_node["dst"];

    if (!src_node.IsMap() || !dst_node.IsMap())
    {
        throw DescriptionParsingError(
                "Source or dest node not defined correctly. "
                "Both should be tuples [<group id>, <neuron id>]",
                edge_node.Mark());
    }

    const auto src_group_id =
            description_required_field<std::string>(src_node, "group");
    const auto src_nid =
            description_required_field<std::string>(src_node, "id");

    const auto dst_group_id =
            description_required_field<std::string>(dst_node, "group");
    const auto dst_nid =
            description_required_field<std::string>(dst_node, "id");

    if (net.groups.find(src_group_id) == net.groups.end())
    {
        const std::string error = "Invalid source neuron group:" + src_group_id;
        throw DescriptionParsingError(error, src_node["group"].Mark());
    }
    NeuronGroup &src_group = net.groups.at(src_group_id);
    if (src_group.neurons.find(src_nid) == src_group.neurons.end())
    {
        const std::string error =
                "Invalid source neuron id: " + src_group_id + "." + src_nid;
        throw DescriptionParsingError(error, src_node["id"].Mark());
    }
    Neuron &src_neuron = src_group.neurons.at(src_nid);

    if (net.groups.find(dst_group_id) == net.groups.end())
    {
        const std::string error = "Invalid target neuron group:" + dst_group_id;
        throw DescriptionParsingError(error, dst_node["group"].Mark());
    }
    NeuronGroup &dst_group = net.groups.at(dst_group_id);

    if (dst_group.neurons.find(dst_nid) == dst_group.neurons.end())
    {
        const std::string error =
                "Invalid target neuron id: " + dst_group_id + "." + dst_nid;
        throw DescriptionParsingError(error, dst_node["id"].Mark());
    }
    Neuron &dst_neuron = dst_group.neurons.at(dst_nid);

    std::map<std::string, ModelParam> synapse_params{};
    std::map<std::string, ModelParam> dendrite_params{};
    if (const YAML::Node &attributes_node = edge_node["attributes"])
    {
        if (const YAML::Node &synapse_node = attributes_node["synapse"])
        {
            synapse_params =
                    description_parse_model_parameters_yaml(synapse_node);
        }
        if (const YAML::Node &dendrite_node = attributes_node["dendrite"])
        {
            dendrite_params =
                    description_parse_model_parameters_yaml(dendrite_node);
        }

        std::set<std::string> unit_specific_keys = {
                "synapse", "dendrite", "soma"};
        auto shared_model_params =
                description_parse_model_parameters_yaml(attributes_node);
        for (auto &[key, parameter] : shared_model_params)
        {
            if (unit_specific_keys.find(key) != unit_specific_keys.end())
            {
                synapse_params.insert({key, parameter});
                dendrite_params.insert({key, parameter});
            }
        }
    }
    std::optional<std::string> synapse_hw_name;

    src_neuron.connect_to_neuron(
            dst_neuron, synapse_params, dendrite_params, synapse_hw_name);
}

void sanafe::description_parse_mapping_file_yaml(
        std::ifstream &fp, Architecture &arch, Network &net)
{
    YAML::Node yaml_node = YAML::Load(fp);
    INFO("YAML information loaded from file.\n");
    if (yaml_node.IsMap() && yaml_node["mappings"].IsDefined())
    {
        description_parse_mapping_section_yaml(
                yaml_node["mappings"], arch, net);
        return;
    }
    throw DescriptionParsingError(
            "No top level mapping section defined", yaml_node.Mark());
}

void sanafe::description_parse_mapping_section_yaml(
        const YAML::Node &mappings_node, Architecture &arch, Network &net)
{
    INFO("Parsing mapping section.\n");
    if (!mappings_node.IsSequence())
    {
        throw DescriptionParsingError(
                "Mappings must be given as a sequence / list. Each list entry "
                "must first use the neuron address as a key, followed by the "
                "core address and optionally specific hardware units.\nE.g.: "
                "- G.n: {core: 1.1, soma: foo}\nThis maps group G's neuron 'n' "
                "to Tile 1's Core 1, specifically using soma unit 'foo'.)",
                mappings_node.Mark());
    }

    for (const auto &mapping : mappings_node)
    {
        if (!mapping.IsMap())
        {
            throw DescriptionParsingError(
                    "Expected mapping to be defined with "
                    "the neuron as the key, and the destination hardware as "
                    "the value.",
                    mapping.Mark());
        }

        for (const auto &pair : mapping)
        {
            auto neuron_address = pair.first.as<std::string>();
            const auto dot_pos = neuron_address.find('.');

            const std::string group_name = neuron_address.substr(0, dot_pos);
            const std::string neuron_id = neuron_address.substr(dot_pos + 1);
            if (net.groups.find(group_name) == net.groups.end())
            {
                const std::string error = "Invalid neuron group:" + group_name;
                throw DescriptionParsingError(error, mapping.Mark());
            }
            NeuronGroup &group = net.groups.at(group_name);
            if (group.neurons.find(neuron_id) == group.neurons.end())
            {
                std::string error = "Invalid neuron id: ";
                error += group_name;
                error += '.';
                error += neuron_id;
                throw DescriptionParsingError(error, mapping.Mark());
            }
            Neuron &neuron = group.neurons.at(neuron_id);

            description_parse_mapping(neuron, pair.second, arch);
        }
    }
}

void sanafe::description_parse_mapping(
        Neuron &neuron, const YAML::Node &mapping_info, Architecture &arch)
{
    if (!mapping_info.IsMap())
    {
        throw DescriptionParsingError(
                "Invalid mapping: mappings must be defined using a YAML map "
                "using 'core', 'synapse', 'dendrite' and / or 'soma' keys",
                mapping_info.Mark());
    }

    const auto core_address =
            description_required_field<std::string>(mapping_info, "core");
    const auto dot_pos = core_address.find('.');

    const size_t tile_id = std::stoull(core_address.substr(0, dot_pos));
    const size_t core_offset_within_tile =
            std::stoull(core_address.substr(dot_pos + 1));

    if (mapping_info["soma"].IsDefined())
    {
        neuron.soma_hw_name = mapping_info["soma"].as<std::string>();
    }
    if (mapping_info["dendrite"].IsDefined())
    {
        neuron.dendrite_hw_name = mapping_info["dendrite"].as<std::string>();
    }
    if (mapping_info["synapse"].IsDefined())
    {
        neuron.default_synapse_hw_name =
                mapping_info["synapse"].as<std::string>();
    }

    if (tile_id > arch.tiles.size())
    {
        throw DescriptionParsingError(
                "Tile ID >= tile count", mapping_info.Mark());
    }
    Tile &tile = arch.tiles[tile_id];
    if (core_offset_within_tile > arch.tiles.size())
    {
        throw DescriptionParsingError(
                "Core ID >= core count", mapping_info.Mark());
    }
    Core &core = tile.cores[core_offset_within_tile];
    core.map_neuron(neuron);
}

std::map<std::string, sanafe::ModelParam>
sanafe::description_parse_model_parameters_yaml(
        const YAML::Node &parameters_node)
{
    std::map<std::string, ModelParam> model_parameters;

    if (parameters_node.IsSequence())
    {
        for (const auto &mapping_list_node : parameters_node)
        {
            // Recursive call to flatten list of mappings
            auto new_parameters =
                    description_parse_model_parameters_yaml(mapping_list_node);
            model_parameters.insert(
                    new_parameters.begin(), new_parameters.end());
        }
    }
    else if (parameters_node.IsMap())
    {
        for (const auto &node : parameters_node)
        {
            const std::string &key = node.first.as<std::string>();
            std::set<std::string> unit_specific_keys = {
                    "synapse", "dendrite", "soma"};
            if (unit_specific_keys.find(key) == unit_specific_keys.end())
            {
                //INFO("Parsing attribute: %s\n", key.c_str());
                model_parameters[key] =
                        description_parse_parameter_yaml(node.second);
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
        const YAML::Node &attribute_node)
{
    ModelParam attribute;

    if (attribute_node.IsSequence())
    {
        // Create an list of unnamed attributes
        std::vector<ModelParam> attribute_list;
        for (const auto &node : attribute_node)
        {
            INFO("Parsing sub-parameter in list.\n");
            ModelParam curr = description_parse_parameter_yaml(node);
            attribute_list.push_back(curr);
        }
        INFO("Setting parameter to an list of %zu unnamed parameters\n",
                attribute_list.size());
        attribute.value = attribute_list;
    }
    else if (attribute_node.IsMap())
    {
        // Create a list of named attributes
        std::vector<ModelParam> attribute_list;
        for (const auto &node : attribute_node)
        {
            INFO("Parsing mapping of parameters.\n");
            ModelParam curr = description_parse_parameter_yaml(node.second);
            curr.name = node.first.as<std::string>();
            INFO("Saving to key: %s\n", curr.name.value().c_str());
            attribute_list.push_back(curr);
        }
        INFO("Setting parameter to a list of %zu named attributes\n",
                attribute_list.size());
        attribute.value = attribute_list;
    }
    else
    {
        if (!attribute_node.IsDefined())
        {
            throw std::invalid_argument("Invalid node.\n");
        }
        double decoded_double;
        int decoded_int;
        bool decoded_bool;
        std::string decoded_str;
        if (YAML::convert<int>::decode(attribute_node, decoded_int))
        {
            //INFO("Parsed int: %d.\n", decoded_int);
            attribute.value = decoded_int;
        }
        else if (YAML::convert<double>::decode(attribute_node, decoded_double))
        {
            //INFO("Parsed float: %lf.\n", decoded_double);
            attribute.value = decoded_double;
        }
        else if (YAML::convert<bool>::decode(attribute_node, decoded_bool))
        {
            //INFO("Parsed bool: %d.\n", decoded_bool);
            attribute.value = decoded_bool;
        }
        else if (YAML::convert<std::string>::decode(
                         attribute_node, decoded_str))
        {
            //INFO("Parsed string: %s.\n", decoded_str.c_str());
            attribute.value = std::move(decoded_str);
        }
    }

    return attribute;
}

std::pair<size_t, size_t> sanafe::description_parse_range_yaml(
        const std::string &range_str, const YAML::Mark &entry_mark)
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
        throw DescriptionParsingError("Invalid range format", entry_mark);
    }

    YAML::Mark range_mark(entry_mark);
    range_mark.column = start_pos;
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
        throw DescriptionParsingError("Invalid range string", range_mark);
    }
    TRACE1("Range: %zu to %zu\n", first, last);
    if (first > last)
    {
        throw DescriptionParsingError(
                "Invalid range; first > last", range_mark);
    }

    return std::make_pair(first, last);
}
