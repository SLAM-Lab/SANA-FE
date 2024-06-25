#include <algorithm>
#include <cassert>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
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
                  ':' + std::to_string(pos.column + 1) + ").")
{
}

const char *sanafe::DescriptionParsingError::what() const noexcept
{
    return message.c_str();
}

template <>
YAML::Node sanafe::description_required_field<YAML::Node>(
        const YAML::Node &node, const std::string &key)
{
    // Specialization of YAML-CPP wrapper for node=map[key], see generic
    //  implementation below (for scalar values)
    if (!node.IsMap())
    {
        throw DescriptionParsingError(
                "Node should be a mapping\n. For more info on YAML mappings "
                "refer to the YAML 1.2 specification, 7.4.2 'Flow Mappings' "
                "and 8.2.2 'Block Mappings'",
                node.Mark());
    }
    const YAML::Node &child = node[key];
    if ((child == nullptr) || child.IsNull() || !child.IsDefined())
    {
        const std::string message = "Value for key '" + key + "' not defined";
        throw DescriptionParsingError(message, node.Mark());
    }

    return child;
}

template <typename T>
T sanafe::description_required_field(
        const YAML::Node &node, const std::string &key)
{
    // Wrapper around YAML-CPP for field=map[key], adding more error prints
    const auto &field_node = description_required_field<YAML::Node>(node, key);
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
    const auto &attributes =
            description_required_field<YAML::Node>(axon_in_node, "attributes");
    const AxonInPowerMetrics in_metrics =
            description_parse_axon_in_attributes_yaml(attributes);
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
    const auto &attributes =
            description_required_field<YAML::Node>(synapse_node, "attributes");

    auto [power_metrics, model] =
            description_parse_synapse_attributes_yaml(attributes);
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
    const YAML::Node &attributes =
            description_required_field<YAML::Node>(dendrite_node, "attributes");

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

    if (attributes["noise"] != nullptr)
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

    const auto &attributes =
            description_required_field<YAML::Node>(axon_out_node, "attributes");
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
    if (core_name.find('[') != std::string::npos)
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

        if (const YAML::Node axon_in_node = core_node["axon_in"])
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

        if (const YAML::Node synapses = core_node["synapse"])
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
        if (const YAML::Node dendrite_node = core_node["dendrite"])
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

        if (const YAML::Node soma_node = core_node["soma"])
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

        if (const YAML::Node axon_out_node = core_node["axon_out"])
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

    if (tile_name.find('[') != std::string::npos)
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
        if (const YAML::Node cores = tile_node["core"])
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

std::string sanafe::description_yaml_parsing_error(
        const std::string &detail, const YAML::Mark &location)
{
    return ("Error: " + detail + std::to_string(location.line + 1) + ":" +
            std::to_string(location.column + 1) + ").\n");
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
    if (YAML::Node arch_yaml_node = top_level_yaml_node["architecture"])
    {
        return description_parse_arch_section_yaml(arch_yaml_node);
    }
    throw DescriptionParsingError("No top-level architecture section defined",
            top_level_yaml_node.Mark());
}

sanafe::Network sanafe::description_parse_network_file_yaml(
        std::ifstream &fp)
{
    YAML::Node yaml_node = YAML::Load(fp);
    if (yaml_node.IsMap() && (yaml_node["network"] != nullptr))
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
    if (net_node["name"] != nullptr)
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

    Network new_net(net_name);
    description_parse_neuron_group_section_yaml(
            description_required_field<YAML::Node>(net_node, "groups"), new_net);
    description_parse_edges_section_yaml(
            description_required_field<YAML::Node>(net_node, "edges"), new_net);

    return new_net;
}

void sanafe::description_parse_neuron_group_section_yaml(
        const YAML::Node &groups_node, Network &net)
{
    if (groups_node.IsSequence())
    {
        for (const auto &group: groups_node)
        description_parse_group(group, net);
    }
    else
    {
        throw DescriptionParsingError(
                "Neuron group section does not define a list of groups",
                groups_node.Mark());
    }
}

void sanafe::description_parse_edges_section_yaml(const YAML::Node &edges_node, Network &net)
{
    if (edges_node.IsMap())
    {
        for (const auto &edge: edges_node)
        {
            const auto edge_description = edge.first;
            const YAML::Node edge_attributes = edge.second;
            description_parse_edge_long_format(edges_node, net);
        }
    }
    else if (edges_node.IsSequence())
    {
        for (const auto &edge : edges_node)
        {
            description_parse_edge_short_format(edge, net);
        }
    }
    else
    {
        throw DescriptionParsingError("Edges section does not define a list or mapping of edges");
    }
}

void sanafe::description_parse_group(
        const YAML::Node &neuron_group_node, Network &net)
{
    const auto name =
            description_required_field<std::string>(neuron_group_node, "name");
    const auto &neurons_node = description_required_field<YAML::Node>(
            neuron_group_node, "neurons");
    const size_t group_neuron_count =
            description_count_neurons_yaml(neurons_node);

    NeuronTemplate default_neuron_config{};
    if (const YAML::Node &group_attributes = neuron_group_node["attributes"])
    {
        default_neuron_config = description_parse_neuron_attributes_yaml(
                group_attributes);
    }
    NeuronGroup &group = net.create_neuron_group(
            name, group_neuron_count, default_neuron_config);
    description_parse_neuron_section_yaml(neurons_node, group);
}

void sanafe::description_parse_neuron_section_yaml(
        const YAML::Node &neuron_node,
        NeuronGroup &neuron_group)
{
    if (neuron_node.IsMap())
    {
        for (const auto &it: neuron_node)
        {
            const auto name = it.first.as<std::string>();
            description_parse_neuron_short_format(name, it.second, neuron_group);
        }
    }
    else if (neuron_node.IsSequence())
    {
        for (const auto &neuron: neuron_group)
        {
            description_parse_neuron_long_format(neuron, neuron_group);
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Invalid neuron format", neuron_node.Mark());
    }
}

void sanafe::description_parse_neuron_long_format(const YAML::Node &neuron_node, NeuronGroup &neuron_group)
{
    const auto &neuron_name =
            description_required_field<std::string>(neuron_node, "name");
    std::pair<size_t, size_t> range = {0, 0};
    if (neuron_name.find('[') != std::string::npos)
    {
        range = description_parse_range_yaml(
                neuron_name, neuron_node["name"].Mark());
    }
    const auto &attributes =
            description_required_field<YAML::Node>(neuron_node, "attributes");
    const NeuronTemplate config =
            description_parse_neuron_attributes_yaml(attributes);

    for (size_t id = range.first; id <= range.second; id++)
    {
        const std::string name = neuron_name.substr(0, neuron_name.find('[')) +
                '[' + std::to_string(id) + ']';
        neuron_group.neurons[name].set_attributes(config);
    }
}

void sanafe::description_parse_neuron_short_format(const std::string &neuron_name,
        const YAML::Node &attributes, NeuronGroup &neuron_group)
{
    std::pair<size_t, size_t> range = {0, 0};
    if (neuron_name.find('['))
    {
        range = description_parse_range_yaml(neuron_name, attributes.Mark());
    }
    const NeuronTemplate config = description_parse_neuron_attributes_yaml(attributes);
    for (size_t id = range.first; id <= range.second; id++)
    {
        const std::string name = neuron_name.substr(0, neuron_name.find('[')) +
                '[' + std::to_string(id) + ']';
        neuron_group.neurons[name].set_attributes(config);
    }
}

sanafe::NeuronTemplate sanafe::description_parse_neuron_attributes_yaml(
        const YAML::Node &attributes, const NeuronTemplate &default_template)
{
    NeuronTemplate neuron_template = default_template;

    if (attributes["log_potential"] != nullptr)
    {
        neuron_template.log_potential = attributes["log_potential"].as<bool>();
    }
    if (attributes["log_spikes"] != nullptr)
    {
        neuron_template.log_spikes = attributes["log_spikes"].as<bool>();
    }
    if (attributes["force_update"] != nullptr)
    {
        neuron_template.force_update = attributes["force_update"].as<bool>();
    }
    if (attributes["soma"] != nullptr)
    {
        neuron_template.soma_model_params =
                description_parse_model_parameters_yaml(attributes["soma"]);
    }
    if (attributes["dendrite"] != nullptr)
    {
        neuron_template.dendrite_model_params =
                description_parse_model_parameters_yaml(attributes["dendrite"]);
    }

    return neuron_template;
}

void sanafe::description_parse_edge(
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

    const auto src_group_id = src_node["group"].as<std::string>();
    const auto src_nid = src_node["id"].as<std::string>();
    const auto dst_group_id = dst_node["group"].as<std::string>();
    const auto dst_nid = dst_node["id"].as<std::string>();

    /*
    if (src_group_id >= net.groups.size())
    {
        const std::string error = "Invalid group ID for source neuron: " +
                std::to_string(src_group_id) +
                ">=" + std::to_string(net.groups.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, src_node.Mark()));
    }
    */

    /*
    if (dst_group_id >= net.groups.size())
    {
        const std::string error = "Invalid group ID for source neuron: " +
                std::to_string(dst_group_id) +
                ">=" + std::to_string(net.groups.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, dst_node.Mark()));
    }
    */

    NeuronGroup &src_group = *(net.groups[src_group_id]);
    /*
    if (src_nid >= src_group.neurons.size())
    {
        const std::string error = "Invalid neuron ID for source neuron" +
                std::to_string(src_nid) +
                ">=" + std::to_string(src_group.neurons.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, src_node.Mark()));
    }
    */
    Neuron &src_neuron = src_group.neurons[src_nid];
    NeuronGroup &dst_group = *(net.groups[dst_group_id]);
    /*
    if (dst_nid >= dst_group.neurons.size())
    {
        const std::string error = "Invalid neuron ID for source neuron" +
                std::to_string(dst_nid) +
                ">=" + std::to_string(dst_group.neurons.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, dst_node.Mark()));
    }
    */
    Neuron &dst_neuron = dst_group.neurons[dst_nid];

    std::map<std::string, ModelParam> synapse_params;
    std::map<std::string, ModelParam> dendrite_params;
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
    }
    std::optional<std::string> synapse_hw_name;
    if (edge_node["map_to"] != nullptr)
    {
        synapse_hw_name = edge_node["map_to"].as<std::string>();
    }
    src_neuron.connect_to_neuron(
            dst_neuron, synapse_params, dendrite_params, synapse_hw_name);
}

void sanafe::description_parse_neuron_mapping_subsection_yaml(
        const YAML::Node &mapping_node, Architecture &arch, Neuron &neuron)
{
    const auto tile_id = mapping_node["tile"].as<size_t>();
    const auto core_offset_within_tile = mapping_node["core"].as<size_t>();
    if (mapping_node["soma"] != nullptr)
    {
        neuron.soma_hw_name = mapping_node["soma"].as<std::string>();
    }
    if (mapping_node["dendrite"] != nullptr)
    {
        neuron.dendrite_hw_name = mapping_node["dendrite"].as<std::string>();
    }
    if (mapping_node["synapse"] != nullptr)
    {
        neuron.default_synapse_hw_name =
                mapping_node["synapse"].as<std::string>();
    }

    if (tile_id > arch.tiles.size())
    {
        throw DescriptionParsingError(
                "Tile ID >= tile count", mapping_node.Mark());
    }
    Tile &tile = arch.tiles[tile_id];
    if (core_offset_within_tile > arch.tiles.size())
    {
        throw DescriptionParsingError(
                "Core ID >= core count", mapping_node.Mark());
    }
    Core &core = tile.cores[core_offset_within_tile];
    core.map_neuron(neuron);
}

void sanafe::description_parse_edge_mapping_subsection_yaml(
        const YAML::Node &mapping_node, Connection &edge)
{
    edge.synapse_hw_name = mapping_node.as<std::string>();
}

std::map<std::string, sanafe::ModelParam>
sanafe::description_parse_model_parameters_yaml(
        const YAML::Node &parameters_node)
{
    if (!parameters_node.IsMap())
    {
        throw std::invalid_argument(
                "Error: Model parameters must be a set of named attributes.\n");
    }

    std::map<std::string, ModelParam> model_parameters;
    for (const auto &node : parameters_node)
    {
        const std::string &key = node.first.as<std::string>();
        INFO("Parsing attribute: %s\n", key.c_str());
        model_parameters[key] = description_parse_parameter_yaml(node.second);
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
            INFO("Parsing list of attributes.\n");
            ModelParam curr = description_parse_parameter_yaml(node);
            attribute_list.push_back(curr);
        }
        INFO("Setting attribute to an unnamed list of %zu attributes",
                attribute_list.size());
        attribute.value = attribute_list;
    }
    else if (attribute_node.IsMap())
    {
        // Create a list of named attributes
        std::vector<ModelParam> attribute_list;
        for (const auto &node : attribute_node)
        {
            INFO("Parsing map of attributes.\n");
            ModelParam curr = description_parse_parameter_yaml(node.second);
            curr.name = node.first.as<std::string>();
            INFO("Saving to key: %s\n", curr.name.value().c_str());
            attribute_list.push_back(curr);
        }
        INFO("Setting attribute to a named list of %zu attributes",
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
            INFO("Parsed int: %d.\n", decoded_int);
            attribute.value = decoded_int;
        }
        else if (YAML::convert<double>::decode(attribute_node, decoded_double))
        {
            INFO("Parsed float: %lf.\n", decoded_double);
            attribute.value = decoded_double;
        }
        else if (YAML::convert<bool>::decode(attribute_node, decoded_bool))
        {
            INFO("Parsed bool: %d.\n", decoded_bool);
            attribute.value = decoded_bool;
        }
        else if (YAML::convert<std::string>::decode(
                         attribute_node, decoded_str))
        {
            INFO("Parsed string: %s.\n", decoded_str.c_str());
            attribute.value = std::move(decoded_str);
        }
    }

    return attribute;
}

std::pair<size_t, size_t> sanafe::description_parse_range_yaml(
        const std::string &range_str, const YAML::Mark &entry_mark)
{
    const size_t start_pos = range_str.find('[');
    const size_t end_pos = range_str.find(']');
    const size_t delimiter_pos = range_str.find("..");
    if ((start_pos == std::string::npos) || (end_pos == std::string::npos) ||
            (delimiter_pos == std::string::npos) || (end_pos <= start_pos) ||
            (delimiter_pos <= start_pos) || (delimiter_pos >= end_pos))
    {
        throw DescriptionParsingError("Invalid range format", entry_mark);
    }

    YAML::Mark range_mark(entry_mark);
    range_mark.column = start_pos;
    const std::string first_str =
            range_str.substr(start_pos + 1, delimiter_pos - start_pos - 1);
    const std::string last_str =
            range_str.substr(delimiter_pos + 2, end_pos - 2 - delimiter_pos);

    size_t first;
    size_t last;
    try
    {
        first = std::stoull(first_str);
        last = std::stoull(last_str);
    }
    catch(const std::exception &e)
    {
        throw DescriptionParsingError("Invalid range string", range_mark);
    }

    if (first > last)
    {
        throw DescriptionParsingError(
                "Invalid range; first > last", range_mark);
    }

    return {first, last};

}
