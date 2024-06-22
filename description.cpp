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
#include "print.hpp"

sanafe::DescriptionParsingError::DescriptionParsingError(
        const std::string &error, const YAML::Mark &pos)
        : std::invalid_argument(error)
{
    message = ("Error: " + error + " (Line " + std::to_string(pos.line + 1) +
            ':' + std::to_string(pos.column + 1) + ").");
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
    if (!child.IsDefined())
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

    const std::string error = "Could not cast field '" +
            YAML::Dump(field_node) + "' (key '" + key +
            "') to type: " + description_get_type_string(field);
    throw DescriptionParsingError(error, field_node.Mark());
}

template <typename T>
std::string sanafe::description_get_type_string(const T &value)
{
    if (typeid(value) == typeid(bool))
    {
        return "bool";
    }
    else if (typeid(value) == typeid(int))
    {
        return "int";
    }
    else if (typeid(value) == typeid(size_t))
    {
        return "size_t";
    }
    else if (typeid(value) == typeid(double))
    {
        return "double";
    }
    else if (typeid(value) == typeid(std::string))
    {
        return "string";
    }
    else
    {
        // Not a scalar type; fall back to default name which may be mangled
        return typeid(value).name();
    }
}

void sanafe::description_parse_axon_in_section_yaml(
        const YAML::Node &axon_in_node, Core &parent_core)
{
    const auto name =
            description_required_field<std::string>(axon_in_node, "name");
    const auto &attributes =
            description_required_field<YAML::Node>(axon_in_node, "attributes");
    const AxonInPowerMetrics in_metrics =
            description_parse_axon_in_attributes_yaml(attributes, parent_core);
    parent_core.create_axon_in(name, in_metrics);
}

sanafe::AxonInPowerMetrics sanafe::description_parse_axon_in_attributes_yaml(
        const YAML::Node &attributes, const Core &parent_core)
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
    const auto name =
            description_required_field<std::string>(synapse_node, "name");
    const auto &attributes =
            description_required_field<YAML::Node>(synapse_node, "attributes");

    auto [power_metrics, model] =
            description_parse_synapse_attributes_yaml(attributes, parent_core);
    parent_core.create_synapse(name, power_metrics, model);
}

std::pair<sanafe::SynapsePowerMetrics, sanafe::ModelInfo>
sanafe::description_parse_synapse_attributes_yaml(
        const YAML::Node &attributes, Core &parent_core)
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
            DescriptionParsingError("Expected plugin path to be string",
                    plugin_path_node.Mark());
        }
    }

    SynapsePowerMetrics power_metrics;
    // TODO: possibly remove this code entirely - leave for now
    //power_metrics.energy_memory_access = description_required_field<double>(
    //        attributes, "energy_memory_access");
    //power_metrics.latency_memory_access = description_required_field<double>(
    //    attributes, "latecy_memory_access");
    power_metrics.energy_process_spike = description_required_field<double>(
            attributes, "energy_process_spike");
    power_metrics.latency_process_spike = description_required_field<double>(
            attributes, "latency_process_spike");

    return {power_metrics, model};
}

void sanafe::description_parse_dendrite_section_yaml(
        const YAML::Node &dendrite_node, Core &parent_core)
{
    std::string dendrite_name = dendrite_node["name"].as<std::string>();
    std::replace(dendrite_name.begin(), dendrite_name.end(), ' ', '_');
    std::replace(dendrite_name.begin(), dendrite_name.end(), '\t', '_');

    const YAML::Node &attributes = dendrite_node["attributes"];
    std::string model_str;
    if (attributes["model"])
    {
        model_str = attributes["model"].as<std::string>();
    }
    else
    {
        throw std::invalid_argument("No dendrite model defined.\n");
    }

    double energy_access = 0.0;
    double latency_access = 0.0;
    if (attributes["energy"])
    {
        energy_access = attributes["energy"].as<double>();
    }
    if (attributes["latency"])
    {
        latency_access = attributes["latency"].as<double>();
    }
    parent_core.create_dendrite(
            dendrite_name, model_str, energy_access, latency_access);
}

void sanafe::description_parse_soma_section_yaml(
        const YAML::Node &soma_node, Core &parent_core)
{
    std::string soma_name = soma_node["name"].as<std::string>();
    const YAML::Node &attributes = soma_node["attributes"];
    std::string model_str;
    if (attributes["model"])
    {
        model_str = attributes["model"].as<std::string>();
    }
    else
    {
        throw std::invalid_argument("No dendrite model defined.\n");
    }

    std::optional<std::string> plugin_lib_path;
    if (attributes["plugin"])
    {
        plugin_lib_path =
                std::filesystem::path(attributes["plugin"].as<std::string>());
    }

    // TODO: maybe just replace this with the metric struct rather than
    //  constructing it twice?
    double energy_update_neuron = 0.0;
    double latency_update_neuron = 0.0;
    double energy_access_neuron = 0.0;
    double latency_access_neuron = 0.0;
    double energy_spike_out = 0.0;
    double latency_spike_out = 0.0;
    if (attributes["energy_update_neuron"])
    {
        energy_update_neuron = attributes["energy_update_neuron"].as<double>();
    }
    if (attributes["latency_update_neuron"])
    {
        latency_update_neuron =
                attributes["latency_update_neuron"].as<double>();
    }
    if (attributes["energy_access_neuron"])
    {
        energy_access_neuron = attributes["energy_access_neuron"].as<double>();
    }
    if (attributes["latency_access_neuron"])
    {
        latency_access_neuron =
                attributes["latency_access_neuron"].as<double>();
    }
    if (attributes["energy_spike_out"])
    {
        energy_spike_out = attributes["energy_spike_out"].as<double>();
    }
    if (attributes["latency_spike_out"])
    {
        latency_spike_out = attributes["latency_spike_out"].as<double>();
    }
    if (attributes["noise"])
    {
        // TODO: support noise again alongside the plugin mechanism
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

    SomaPowerMetrics power_metrics(energy_update_neuron, latency_update_neuron,
            energy_access_neuron, latency_access_neuron, energy_spike_out,
            latency_spike_out);
    parent_core.create_soma(std::move(soma_name), std::move(model_str),
            power_metrics, plugin_lib_path);
}

void sanafe::description_parse_axon_out_section(
        const YAML::Node &axon_out_node, Core &parent_core)
{
    std::string axon_out_name = axon_out_node["name"].as<std::string>();
    std::replace(axon_out_name.begin(), axon_out_name.end(), ' ', '_');
    std::replace(axon_out_name.begin(), axon_out_name.end(), '\t', '_');

    // TODO: maybe just replace this with the metric struct
    const YAML::Node &attributes = axon_out_node["attributes"];
    double energy_message = 0.0;
    double latency_message = 0.0;
    if (attributes["energy_message_out"])
    {
        energy_message = attributes["energy_message_out"].as<double>();
    }
    if (attributes["latency_message_out"])
    {
        latency_message = attributes["latency_message_out"].as<double>();
    }

    parent_core.create_axon_out(axon_out_name, energy_message, latency_message);
}

void sanafe::description_parse_core_section_yaml(const YAML::Node &core_node,
        const size_t parent_tile_id, Architecture &arch)
{
    std::string core_name = core_node["name"].as<std::string>();
    std::replace(core_name.begin(), core_name.end(), ' ', '_');
    std::replace(core_name.begin(), core_name.end(), '\t', '_');

    std::pair<int, int> core_range = {0, 0};
    if (core_name.find('[') != std::string::npos)
    {
        core_range = description_parse_range_yaml(core_name);
    }

    for (int c = core_range.first; c <= core_range.second; c++)
    {
        const std::string name = core_name.substr(0, core_name.find("[")) +
                "[" + std::to_string(c) + "]";
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
            const std::string error =
                    "Error: No axon in section defined (line:" +
                    std::to_string(axon_in_node.Mark().line + 1) + ":" +
                    std::to_string(axon_in_node.Mark().column + 1) + ").\n";
            throw std::invalid_argument(error);
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
            const std::string error =
                    "Error: No synapse section defined (line:" +
                    std::to_string(synapses.Mark().line + 1) + ":" +
                    std::to_string(synapses.Mark().column + 1) + ").\n";
        }
        if (const YAML::Node dendrite_node = core_node["dendrite"])
        {
            if (dendrite_node.IsSequence())
            {
                for (auto d : dendrite_node)
                {
                    description_parse_dendrite_section_yaml(d, core);
                }
            }
            else
            {
                description_parse_dendrite_section_yaml(dendrite_node, core);
            }
        }
        else
        {
            const std::string error =
                    "Error: No dendrite section defined (line:" +
                    std::to_string(dendrite_node.Mark().line + 1) + ":" +
                    std::to_string(dendrite_node.Mark().column + 1) + ").\n";
        }

        if (const YAML::Node soma_node = core_node["soma"])
        {
            if (soma_node.IsSequence())
            {
                for (auto s : soma_node)
                {
                    description_parse_soma_section_yaml(s, core);
                }
            }
            else
            {
                description_parse_soma_section_yaml(soma_node, core);
            }
        }
        else
        {
            const std::string error = "Error: No soma section defined (line:" +
                    std::to_string(soma_node.Mark().line + 1) + ":" +
                    std::to_string(soma_node.Mark().column + 1) + ").\n";
        }

        if (const YAML::Node axon_out_node = core_node["axon_out"])
        {
            if (axon_out_node.IsSequence())
            {
                for (auto axon : axon_out_node)
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
            const std::string error =
                    "Error: No axon out section defined (line:" +
                    std::to_string(core_node.Mark().line + 1) + ":" +
                    std::to_string(core_node.Mark().column + 1) + ").\n";
        }
    }
}

sanafe::CorePipelineConfiguration sanafe::description_parse_core_pipeline_yaml(
        const YAML::Node &attributes)
{
    std::string buffer_pos = "soma";
    size_t max_neurons_supported = 1024;
    if (attributes["buffer_position"])
    {
        buffer_pos = attributes["buffer_position"].as<std::string>();
    }
    if (attributes["max_neurons_supported"])
    {
        max_neurons_supported = attributes["buffer_position"].as<size_t>();
    }
    CorePipelineConfiguration pipeline_config(
            buffer_pos, max_neurons_supported);

    return pipeline_config;
}

sanafe::TilePowerMetrics sanafe::description_parse_tile_metrics_yaml(
        const YAML::Node &attributes)
{
    // Per-hop power metrics
    double energy_north = 0.0;
    double latency_north = 0.0;
    double energy_east = 0.0;
    double latency_east = 0.0;
    double energy_west = 0.0;
    double latency_west = 0.0;
    double energy_south = 0.0;
    double latency_south = 0.0;

    if (attributes["energy_north_hop"])
    {
        energy_north = attributes["energy_north_hop"].as<double>();
    }
    if (attributes["latency_north_hop"])
    {
        latency_north = attributes["latency_north_hop"].as<double>();
    }
    if (attributes["energy_east_hop"])
    {
        energy_east = attributes["energy_east_hop"].as<double>();
    }
    if (attributes["latency_east_hop"])
    {
        latency_east = attributes["latency_east_hop"].as<double>();
    }
    if (attributes["energy_south_hop"])
    {
        energy_south = attributes["energy_south_hop"].as<double>();
    }
    if (attributes["latency_south_hop"])
    {
        latency_south = attributes["latency_south_hop"].as<double>();
    }
    if (attributes["energy_west_hop"])
    {
        energy_west = attributes["energy_west_hop"].as<double>();
    }
    if (attributes["latency_west_hop"])
    {
        latency_west = attributes["latency_west_hop"].as<double>();
    }

    TilePowerMetrics tile_metrics(energy_north, latency_north, energy_east,
            latency_east, energy_south, latency_south, energy_west,
            latency_west);
    return tile_metrics;
}

void sanafe::description_parse_tile_section_yaml(
        const YAML::Node &tile_node, Architecture &arch)
{
    auto tile_name = tile_node["name"].as<std::string>();
    std::replace(tile_name.begin(), tile_name.end(), ' ', '_');
    std::replace(tile_name.begin(), tile_name.end(), '\t', '_');

    std::pair<int, int> range = {0, 0};
    if (tile_name.find('[') != std::string::npos)
    {
        range = description_parse_range_yaml(tile_name);
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
            throw std::invalid_argument(description_yaml_parsing_error(
                    "No core section defined", tile_node.Mark()));
        }
    }
}

sanafe::NetworkOnChipConfiguration
sanafe::description_parse_noc_configuration_yaml(
        const YAML::Node &noc_attributes)
{
    size_t width_in_tiles = 1;
    size_t height_in_tiles = 1;
    size_t link_buffer_size = 0;

    if (noc_attributes["width"])
    {
        width_in_tiles = noc_attributes["width"].as<size_t>();
    }
    if (noc_attributes["height"])
    {
        height_in_tiles = noc_attributes["height"].as<size_t>();
    }
    if (noc_attributes["link_buffer_size"])
    {
        link_buffer_size = noc_attributes["link_buffer_size"].as<size_t>();
    }

    const NetworkOnChipConfiguration noc(
            width_in_tiles, height_in_tiles, link_buffer_size);
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
        throw std::invalid_argument(description_yaml_parsing_error(
                "Multiple architectures not supported", arch_node.Mark()));
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
        throw std::invalid_argument(description_yaml_parsing_error(
                "No tile section defined", tiles.Mark()));
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
    else
    {
        throw std::runtime_error(description_yaml_parsing_error(
                "No top-level architecture section defined",
                top_level_yaml_node.Mark()));
    }
}

sanafe::Network sanafe::description_parse_net_file_yaml(
        std::ifstream &fp, Architecture &arch)
{
    YAML::Node yaml_node = YAML::Load(fp);
    if (yaml_node.IsMap() && yaml_node["network"])
    {
        return description_parse_net_section_yaml(yaml_node["network"], arch);
    }
    else
    {
        const std::string error = "Error: No network section defined (line:" +
                std::to_string(yaml_node.Mark().line + 1) + ":" +
                std::to_string(yaml_node.Mark().column + 1) + ").\n";
        throw std::runtime_error(error);
    }
}

sanafe::Network sanafe::description_parse_net_section_yaml(
        const YAML::Node &net_node, Architecture &arch)
{
    std::string net_name = "";
    if (net_node["name"])
    {
        net_name = net_node["name"].as<std::string>();
        if (net_name.find("[") != std::string::npos)
        {
            throw std::runtime_error("Error: Multiple networks not supported");
        }
    }
    else
    {
        INFO("No network name given; leaving name empty.\n");
    }

    Network new_net(net_name);
    if (const YAML::Node group_node = net_node["groups"])
    {
        // Iterate through all tiles
        if (group_node.IsSequence())
        {
            for (const auto &g : group_node)
            {
                description_parse_neuron_group_section_yaml(g, arch, new_net);
            }
        }
        else
        {
            description_parse_neuron_group_section_yaml(
                    group_node, arch, new_net);
        }
    }
    else
    {
        const std::string error = "Error: No group section defined (line:" +
                std::to_string(group_node.Mark().line + 1) + ":" +
                std::to_string(group_node.Mark().column + 1) + ").\n";
        throw std::invalid_argument(error);
    }

    if (const YAML::Node edges_node = net_node["edges"])
    {
        // Iterate through all tiles
        if (edges_node.IsSequence())
        {
            for (const auto &e : edges_node)
            {
                description_parse_edge_section_yaml(e, new_net);
            }
        }
        else
        {
            description_parse_edge_section_yaml(edges_node, new_net);
        }
    }
    else
    {
        const std::string error = "No edge section defined" +
                std::to_string(edges_node.Mark().line + 1) + ":" +
                std::to_string(edges_node.Mark().column + 1) + ").\n";
        throw std::invalid_argument(error);
    }

    return new_net;
}

sanafe::NeuronTemplate sanafe::description_parse_neuron_attributes_yaml(
        const YAML::Node &attributes, const NeuronTemplate &default_template)
{
    NeuronTemplate neuron_template = default_template;

    if (attributes["log_potential"])
    {
        neuron_template.log_potential = attributes["log_potential"].as<bool>();
    }
    if (attributes["log_spikes"])
    {
        neuron_template.log_spikes = attributes["log_spikes"].as<bool>();
    }
    if (attributes["force_update"])
    {
        neuron_template.force_update = attributes["force_update"].as<bool>();
    }
    if (attributes["soma"])
    {
        neuron_template.soma_model_params =
                description_parse_model_parameters_yaml(attributes["soma"]);
    }
    if (attributes["dendrite"])
    {
        neuron_template.dendrite_model_params =
                description_parse_model_parameters_yaml(attributes["dendrite"]);
    }

    return neuron_template;
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

size_t sanafe::description_count_neurons_yaml(const YAML::Node &neurons_node)
{
    // TODO: is this basically just overloading the YAML::Node::size() function..
    //  check to see what this api is doing
    if (neurons_node.IsNull())
    {
        INFO("Warning: No neurons defined for group");
        return 0;
    }
    else if (neurons_node.IsScalar())
    {
        return 1;
    }
    else if (neurons_node.IsSequence())
    {
        return neurons_node.size();
    }
    else
    {
        throw std::invalid_argument("Invalid YAML node type for neurons.\n");
    }
}

void sanafe::description_parse_neuron_group_section_yaml(
        const YAML::Node &neuron_group_node, Architecture &arch, Network &net)
{
    if (!neuron_group_node["name"])
    {
        INFO("Error: No neuron group name given.\n");
        throw std::invalid_argument("Error: No neuron group name given.\n");
    }

    const auto group_name = neuron_group_node["name"].as<std::string>();
    if (!neuron_group_node["neurons"])
    {
        INFO("Error: No neurons defined for neuron group.\n");
        throw std::invalid_argument(
                "Error: No neurons defined for neuron group.\n");
    }
    const YAML::Node neurons_node = neuron_group_node["neurons"];
    const size_t group_neuron_count =
            description_count_neurons_yaml(neurons_node);

    if (!neuron_group_node["attributes"])
    {
        INFO("Error: No neuron group attributes defined.\n");
        throw std::invalid_argument(
                "Error: No neuron group attributes defined.\n");
    }
    const YAML::Node &neuron_group_attributes = neuron_group_node["attributes"];
    const NeuronTemplate default_neuron_config =
            description_parse_neuron_attributes_yaml(neuron_group_attributes);
    NeuronGroup &neuron_group = net.create_neuron_group(
            group_name, group_neuron_count, default_neuron_config);

    if (neurons_node.IsSequence())
    {
        for (const auto &n : neurons_node)
        {
            description_parse_neuron_section_yaml(n, arch, neuron_group);
        }
    }
    else
    {
        description_parse_neuron_section_yaml(neurons_node, arch, neuron_group);
    }

    return;
}

void sanafe::description_parse_edge_section_yaml(
        const YAML::Node &edge_node, Network &net)
{
    const YAML::Node &src_node = edge_node["src"];
    const YAML::Node &dst_node = edge_node["dst"];

    if (!src_node.IsMap() || !dst_node.IsMap())
    {
        throw std::invalid_argument(description_yaml_parsing_error(
                "Source or dest node not defined correctly. "
                "Both should be tuples [<group id>, <neuron id>]",
                edge_node.Mark()));
    }

    const auto src_group_id = src_node["group"].as<size_t>();
    const auto src_nid = src_node["id"].as<size_t>();
    const auto dst_group_id = dst_node["group"].as<size_t>();
    const auto dst_nid = dst_node["id"].as<size_t>();

    if (src_group_id >= net.groups.size())
    {
        const std::string error = "Invalid group ID for source neuron: " +
                std::to_string(src_group_id) +
                ">=" + std::to_string(net.groups.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, src_node.Mark()));
    }

    if (dst_group_id >= net.groups.size())
    {
        const std::string error = "Invalid group ID for source neuron: " +
                std::to_string(dst_group_id) +
                ">=" + std::to_string(net.groups.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, dst_node.Mark()));
    }

    NeuronGroup &src_group = *(net.groups[src_group_id]);
    if (src_nid >= src_group.neurons.size())
    {
        const std::string error = "Invalid neuron ID for source neuron" +
                std::to_string(src_nid) +
                ">=" + std::to_string(src_group.neurons.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, src_node.Mark()));
    }
    Neuron &src_neuron = src_group.neurons[src_nid];
    NeuronGroup &dst_group = *(net.groups[dst_group_id]);
    if (dst_nid >= dst_group.neurons.size())
    {
        const std::string error = "Invalid neuron ID for source neuron" +
                std::to_string(dst_nid) +
                ">=" + std::to_string(dst_group.neurons.size());
        throw std::invalid_argument(
                description_yaml_parsing_error(error, dst_node.Mark()));
    }
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
    if (edge_node["map_to"])
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
    if (mapping_node["soma"])
    {
        neuron.soma_hw_name = mapping_node["soma"].as<std::string>();
    }
    if (mapping_node["dendrite"])
    {
        neuron.dendrite_hw_name = mapping_node["dendrite"].as<std::string>();
    }
    if (mapping_node["synapse"])
    {
        neuron.default_synapse_hw_name =
                mapping_node["synapse"].as<std::string>();
    }

    Tile &tile = arch.tiles[tile_id];
    Core &core = tile.cores[core_offset_within_tile];
    core.map_neuron(neuron);
}

void sanafe::description_parse_edge_mapping_subsection_yaml(
        const YAML::Node &mapping_node, Connection &edge)
{
    edge.synapse_hw_name = mapping_node.as<std::string>();
}

void sanafe::description_parse_neuron_section_yaml(
        const YAML::Node &neuron_node, Architecture &arch,
        NeuronGroup &neuron_group)
{
    size_t neuron_id;
    if (neuron_node["id"])
    {
        neuron_id = neuron_node["id"].as<size_t>();
        if (neuron_id >= neuron_group.neurons.size())
        {
            throw std::invalid_argument("Neuron ID > max neurons in group.\n");
        }
    }
    else
    {
        throw std::invalid_argument("No valid neuron ID given.\n");
    }

    if (neuron_node["attributes"])
    {
        // TODO: we need to have it set to the defaults set by the group, and
        //  only override parameters set in the attributes section
        const NeuronTemplate config = description_parse_neuron_attributes_yaml(
                neuron_node["attributes"]);
        neuron_group.neurons[neuron_id].set_attributes(config);
    }

    if (neuron_node["map_to"])
    {
        description_parse_neuron_mapping_subsection_yaml(
                neuron_node["map_to"], arch, neuron_group.neurons[neuron_id]);
    }
    else
    {
        throw std::invalid_argument("No valid neuron to H/W mapping given.\n");
    }
}

std::pair<size_t, size_t> sanafe::description_parse_range_yaml(
        const std::string &range_str)
{
    const size_t start_pos = range_str.find('[');
    const size_t end_pos = range_str.find(']');
    const size_t delimiter_pos = range_str.find("..");
    if ((start_pos == std::string::npos) || (end_pos == std::string::npos) ||
            (delimiter_pos == std::string::npos) || (end_pos <= start_pos))
    {
        const std::string error = "Error: Invalid range:" + range_str + '\n';
        throw std::invalid_argument(error);
    }
    // Use stringstreams to parse the size_t values from the range substrings
    std::istringstream first_ss(
            range_str.substr(start_pos + 1, delimiter_pos - start_pos - 1));
    std::istringstream last_ss(
            range_str.substr(delimiter_pos + 2, end_pos - 2 - delimiter_pos));
    size_t first;
    first_ss >> first;
    size_t last;
    last_ss >> last;

    return {first, last};
}

/*
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
*/

/*
void sanafe::description_read_arch_entry(
        const std::vector<std::string_view> &fields, Architecture &arch,
        const int line_number)
{
    std::map<std::string, std::string> attributes;
    std::string name;
    Tile *tile_ptr;
    Core *core_ptr;
    int tile_id, core_offset, first_field;

    assert(fields.size() > 0);
    const char entry_type = fields[0][0];
    // Sanity check input
    if ((entry_type == '\0') || (entry_type == '\n') || (entry_type == '#') ||
            (entry_type == '\r'))
    {
        TRACE1("Warning: No entry, skipping\n");
        return;
    }

    tile_ptr = nullptr;
    core_ptr = nullptr;
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
        Tile &t = arch.tiles[tile_id];
        tile_ptr = &t;
        first_field++;
    }
    if ((entry_type != '@') && (entry_type != 't') && (entry_type != 'c'))
    {
        core_offset = field_to_int(fields[3]);
        assert(tile_ptr != nullptr);
        Core &c = tile_ptr->cores_ref[core_offset];
        core_ptr = &c;
        first_field++;
    }

    // Parse attributes from fields
    for (size_t i = first_field; i < fields.size(); i++)
    {
        TRACE1("Parsing field:%s\n", std::string(fields[i]).c_str());

        if ((fields[i].length() < 3))
        {
            INFO("Error: Line: %d Invalid field: %s\n", line_number,
                    std::string(fields[i]).c_str());
            continue;
        }

        int pos = fields[i].find_first_of('=');
        std::string key = std::string(fields[i].substr(0, pos));
        std::string value_str = std::string(fields[i].substr(pos + 1));

        if ((key.length() == 0) || (value_str.length() == 0))
        {
            INFO("Invalid attribute: %s\n", std::string(fields[i]).c_str());
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
        arch.create_tile(name, attributes);
        break;
    case 'c':
        arch.create_core(name, tile_id, attributes);
        break;
    case 'i':
        core_ptr->create_axon_in(name, attributes);
        break;
    case 's':
        core_ptr->create_synapse(name, attributes);
        break;
    case 'd':
        core_ptr->create_dendrite(name, attributes);
        break;
    case '+':
        core_ptr->create_soma(name, attributes);
        break;
    case 'o':
        core_ptr->create_axon_out(name, attributes);
        break;
    default:
        INFO("Warning: unrecognized unit (%c) - skipping.\n", entry_type);
        break;
    }

    return;
}
*/

/*
void sanafe::parse_neuron_with_compartment_field(
        const std::string_view &neuron_field, size_t &group_id,
        size_t &neuron_id, std::optional<size_t> &compartment_id)
{
    parse_neuron_field(neuron_field, group_id, neuron_id);
    const auto pos = neuron_field.find(':');
    if (pos != std::string_view::npos)
    {
        const auto compartment_str = neuron_field.substr(pos + 1);
        compartment_id = field_to_int(compartment_str);
    }
    else
    {
        TRACE1("Compartment not set, leaving default value\n");
    }

    return;
}
*/

/*
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

void sanafe::parse_neuron_field(const std::string_view &neuron_field,
        size_t &group_id, size_t &neuron_id)
{
    const auto pos = neuron_field.find('.');
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Error: Invalid neuron format");
    }

    const auto group_str = neuron_field.substr(0, pos);
    group_id = field_to_int(group_str);
    const auto neuron_str = neuron_field.substr(pos + 1, neuron_field.size());
    neuron_id = field_to_int(neuron_str);

    return;
}
*/

/*
void sanafe::parse_core_field(const std::string_view &core_field,
        size_t &tile_id, size_t &core_offset)
{
    const auto pos = core_field.find('.');
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Error: Invalid neuron format");
    }

    const auto tile_str = core_field.substr(0, pos);
    tile_id = field_to_int(tile_str);
    const auto core_str = core_field.substr(pos + 1, core_field.size());
    core_offset = field_to_int(core_str);

    return;
}
*/

/*
void sanafe::parse_edge_field(const std::string_view &edge_field,
        size_t &group_id, size_t &neuron_id, size_t &dest_group_id,
        size_t &dest_neuron_id)
{
    // Edge description entries support two formats, to represent
    //  neuron-neuron connections and compartment-compartment branches
    // i.e. Connection: e group.neuron->group.neuron:compartment <attributes>
    // i.e. Branch:     e group.neuron:compartment->compartment <attributes>
    //   Note that the destination compartment is optional (default=0)
    // Split the source and destination neuron addresses
    TRACE1("Parsing edge.\n");
    const auto pos = edge_field.find("->");
    if (pos == std::string_view::npos)
    {
        throw std::runtime_error("Invalid edge format");
    }

    // Parse the source group, neuron and optional compartment identifiers
    //  from the source neuron substring (before the "->")
    const auto src_neuron_address = edge_field.substr(0, pos);
    parse_neuron_field(src_neuron_address, group_id, neuron_id);
    // Parse the destination group, neuron and compartment
    //  identifiers from the substring after "->"
    const auto dest_neuron_address =
            edge_field.substr(pos + 2, edge_field.size());
    parse_neuron_field(dest_neuron_address, dest_group_id, dest_neuron_id);

    return;
}
*/

/*
void sanafe::parse_mapping_field(const std::string_view &mapping_field,
        size_t &group_id, size_t &neuron_id, size_t &tile_id,
        size_t &core_offset)
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

    const auto core_address =
            mapping_field.substr(pos + 1, mapping_field.size());
    parse_neuron_field(core_address, tile_id, core_offset);

    return;
}
*/

/*
void sanafe::description_read_network_entry(
        const std::vector<std::string_view> &fields, Architecture &arch,
        Network &net, const int line_number)
{
    std::map<std::string, std::string> attributes;
    NeuronGroup *group_ptr, *dest_group_ptr;
    Neuron *neuron_ptr, *dest_ptr;
    Tile *tile_ptr;
    Core *core_ptr;
    size_t tile_id, core_offset, neuron_group_id, dest_group_id;
    size_t neuron_id, dest_neuron_id;
    int neuron_count;
    bool group_set, neuron_set;

    assert(fields.size() > 0);
    const char entry_type = fields[0][0];
    // Sanity check input
    if ((entry_type == '\0') || (entry_type == '\n') || (entry_type == '#') ||
            (entry_type == '\r'))
    {
        TRACE1("Warning: No entry, skipping line %d\n", line_number);
        return;
    }

    if (fields.size() < 2)
    {
        INFO("Error: fields < 2 (%ld)", fields.size());
    }

    neuron_count = 0;
    neuron_ptr = nullptr;
    group_ptr = nullptr;
    core_ptr = nullptr;
    dest_group_ptr = nullptr;
    dest_ptr = nullptr;

    group_set = false;
    neuron_set = false;

    if (entry_type == 'g')
    {
        neuron_count = field_to_int(fields[1]);
    }
    else if (entry_type == '&')
    {
        parse_mapping_field(
                fields[1], neuron_group_id, neuron_id, tile_id, core_offset);
        if (tile_id >= arch.tiles.size())
        {
            INFO("Error: Line %d: Tile (%lu) >= tile count (%lu)\n",
                    line_number, tile_id, arch.tiles.size());
            exit(1);
        }
        Tile &tile = arch.tiles[tile_id];
        tile_ptr = &tile;

        if (core_offset >= tile_ptr->cores.size())
        {
            INFO("Error: Line %d: Core (%lu) >= core count (%lu)\n",
                    line_number, core_offset, tile_ptr->cores.size());
            exit(1);
        }
        Core &core = tile_ptr->cores[core_offset];
        core_ptr = &core;
        group_set = true;
        neuron_set = true;
    }
    else if (entry_type == 'e')
    {
        // Edge on SNN graph (e.g., connection between neurons)
        parse_edge_field(fields[1], neuron_group_id, neuron_id, dest_group_id,
                dest_neuron_id);
        if (dest_group_id >= net.groups.size())
        {
            INFO("Error: Line %d: Group (%lu) >= group count (%lu).\n",
                    line_number, dest_group_id, net.groups.size());
            throw std::invalid_argument("Invalid group id");
        }
        NeuronGroup &dest_group = *(net.groups[dest_group_id]);
        dest_group_ptr = &dest_group;

        TRACE1("Parsed neuron gid:%lu nid:%lu\n", dest_group_id, neuron_id);
        if (dest_neuron_id >= dest_group_ptr->neurons.size())
        {
            INFO("Error: Line %d: Trying to access neuron "
                 "(%d.%lu) but group %d only "
                 "allocates %lu neurons.\n",
                    line_number, dest_group_ptr->id, dest_neuron_id,
                    dest_group_ptr->id, dest_group_ptr->neurons.size());
            throw std::invalid_argument("Invalid nid");
        }
        dest_ptr = &(dest_group_ptr->neurons[dest_neuron_id]);
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
        INFO("Error: Line %d: Invalid entry type (%s)", line_number,
                std::string(fields[0]).c_str());
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
        group_ptr = net.groups[neuron_group_id].get();
        NeuronGroup &group = *group_ptr;
        TRACE1("Parsed neuron gid:%lu nid:%lu\n", neuron_group_id, neuron_id);
        if (neuron_set)
        {
            if (neuron_id >= group.neurons.size())
            {
                INFO("Error: Line %d: Trying to access neuron "
                     "(%d.%lu) but group %d only "
                     "allocates %lu neuron(s).\n",
                        line_number, group.id, neuron_id, group.id,
                        group.neurons.size());
                throw std::invalid_argument("Invalid neuron id");
            }
            neuron_ptr = &(group.neurons[neuron_id]);
        }
    }

    // Parse attributes from fields
    for (size_t i = 2; i < fields.size(); i++)
    {
        TRACE1("Parsing field:%s\n", fields[i].c_str());

        if ((fields[i].length() < 3))
        {
            INFO("Error: Line %d: Invalid field: %s\n", line_number,
                    std::string(fields[i]).c_str());
            continue;
        }

        const auto pos = fields[i].find_first_of('=');
        std::string key = std::string(fields[i].substr(0, pos));
        std::string value_str = std::string(fields[i].substr(pos + 1));

        if ((key.length() == 0) || (value_str.length() == 0))
        {
            INFO("Error: Line %d: Invalid attribute: %s\n", line_number,
                    std::string(fields[i]).c_str());
            continue;
        }

        attributes.insert({key, value_str});
    }

    // Process the entry
    switch (entry_type)
    {
    case 'g': // Add neuron group
        //net.create_neuron_group("", neuron_count, attributes);
        break;
    case 'n': // Add neuron
        //neuron_ptr->set_attributes(attributes);
        break;
    case 'e':
        assert(neuron_ptr != nullptr);
        // Zero initialize all connections
        neuron_ptr->connect_to_neuron(*dest_ptr, attributes);
        break;
    case '&': // Map neuron to hardware
        core_ptr->map_neuron(*neuron_ptr);
        break;
    default:
        break;
    }

    return;
}
*/