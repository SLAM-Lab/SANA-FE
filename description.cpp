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

std::pair<size_t, size_t> sanafe::description_parse_range(
        const std::string &range_str)
{
    const size_t start_pos = range_str.find("[");
    const size_t end_pos = range_str.find("]");
    const size_t delimiter_pos = range_str.find("..");
    if ((start_pos == std::string::npos) || (end_pos == std::string::npos) ||
            (delimiter_pos == std::string::npos) || (end_pos <= start_pos))
    {
        INFO("Error: Invalid range formatting: %s\n", range_str.c_str());
        throw std::invalid_argument("Error: Invalid range formatting.\n");
    }

    std::istringstream first_ss(
            range_str.substr(start_pos + 1, delimiter_pos - start_pos - 1));
    size_t first;
    first_ss >> first;

    std::istringstream last_ss(
            range_str.substr(delimiter_pos + 2, end_pos - 2 - delimiter_pos));
    size_t last;
    last_ss >> last;

    std::pair<size_t, size_t> range = {first, last};

    return range;
}

void sanafe::description_parse_axon_in_section(
        const YAML::Node &axon_in_node, Core &parent_core)
{
    std::string axon_in_name = axon_in_node["name"].as<std::string>();
    std::replace(axon_in_name.begin(), axon_in_name.end(), ' ', '_');
    std::replace(axon_in_name.begin(), axon_in_name.end(), '\t', '_');
    const YAML::Node &attributes = axon_in_node["attributes"];
    double energy_message = 0.0;
    double latency_message = 0.0;
    if (attributes["energy_message"])
    {
        energy_message = attributes["energy_message"].as<double>();
    }
    if (attributes["latency_message"])
    {
        latency_message = attributes["latency_message"].as<double>();
    }
    parent_core.create_axon_in(axon_in_name, energy_message, latency_message);

    return;
}

void sanafe::description_parse_synapse_section(
        const YAML::Node &synapse_node, Core &parent_core)
{
    std::string synapse_name = synapse_node["name"].as<std::string>();
    std::replace(synapse_name.begin(), synapse_name.end(), ' ', '_');
    std::replace(synapse_name.begin(), synapse_name.end(), '\t', '_');

    const YAML::Node &attributes = synapse_node["attributes"];
    std::string model_str;
    if (attributes["model"])
    {
        model_str = attributes["model"].as<std::string>();
    }
    else
    {
        throw std::invalid_argument("No synapse model defined.\n");
    }
    double energy_memory = 0.0;
    double latency_memory = 0.0;
    double energy_spike = 0.0;
    double latency_spike = 0.0;
    // Parse power metricslatency_message
    if (attributes["energy_memory_access"])
    {
        energy_memory = attributes["energy_memory_access"].as<double>();
    }
    if (attributes["latency_memory_access"])
    {
        latency_memory = attributes["latency_memory_access"].as<double>();
    }
    if (attributes["energy_process_spike"])
    {
        energy_spike = attributes["energy_process_spike"].as<double>();
    }
    if (attributes["latency_process_spike"])
    {
        latency_spike = attributes["latency_process_spike"].as<double>();
    }
    SynapsePowerMetrics power_metrics(
            energy_memory, latency_memory, energy_spike, latency_spike);
    // TODO: I think we need more information on the model.. i.e. possibly an
    //  optional plugin path (indicating its a plugin not a built-in model)
    parent_core.create_synapse(synapse_name, model_str, power_metrics);

    return;
}

void sanafe::description_parse_dendrite_section(
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

void sanafe::description_parse_soma_section(
        const YAML::Node &soma_node, Core &parent_core)
{
    std::string soma_name = soma_node["name"].as<std::string>();
    std::replace(soma_name.begin(), soma_name.end(), ' ', '_');
    std::replace(soma_name.begin(), soma_name.end(), '\t', '_');

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

    // TODO: figure how to define this plugin library parameter
    if (attributes["plugin_lib"])
    {
        //std::string plugin_lib = std::filesystem::path(value_str);
    }

    // TODO: maybe just replace this with the metric struct
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
        latency_update_neuron = attributes["latency_update_neuron"].as<double>();
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
    parent_core.create_soma(soma_name, model_str, power_metrics);
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

void sanafe::description_parse_core_section(const YAML::Node &core_node,
        const size_t parent_tile_id, Architecture &arch)
{
    std::string core_name = core_node["name"].as<std::string>();
    std::replace(core_name.begin(), core_name.end(), ' ', '_');
    std::replace(core_name.begin(), core_name.end(), '\t', '_');

    std::pair<int, int> core_range = {0, 0};
    if (core_name.find("[") != std::string::npos)
    {
        core_range = description_parse_range(core_name);
    }

    for (int c = core_range.first; c <= core_range.second; c++)
    {
        const std::string name = core_name.substr(0, core_name.find("[")) +
                "[" + std::to_string(c) + "]";
        const CorePipelineConfiguration pipeline_config =
                description_parse_core_pipeline(core_node["attributes"]);
        Core &core = arch.create_core(name, parent_tile_id, pipeline_config);

        if (const YAML::Node axon_in_node = core_node["axon_in"])
        {
            if (axon_in_node.IsSequence())
            {
                for (auto a : axon_in_node)
                {
                    description_parse_axon_in_section(a, core);
                }
            }
            else
            {
                description_parse_axon_in_section(axon_in_node, core);
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

        if (const YAML::Node synapse_node = core_node["synapse"])
        {
            if (synapse_node.IsSequence())
            {
                for (auto s : synapse_node)
                {
                    description_parse_synapse_section(s, core);
                }
            }
            else
            {
                description_parse_synapse_section(synapse_node, core);
            }
        }
        else
        {
            const std::string error =
                    "Error: No synapse section defined (line:" +
                    std::to_string(synapse_node.Mark().line + 1) + ":" +
                    std::to_string(synapse_node.Mark().column + 1) + ").\n";
        }
        if (const YAML::Node dendrite_node = core_node["dendrite"])
        {
            if (dendrite_node.IsSequence())
            {
                for (auto d : dendrite_node)
                {
                    description_parse_dendrite_section(d, core);
                }
            }
            else
            {
                description_parse_dendrite_section(dendrite_node, core);
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
                    description_parse_soma_section(s, core);
                }
            }
            else
            {
                description_parse_soma_section(soma_node, core);
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
                for (auto a : axon_out_node)
                {
                    description_parse_axon_out_section(a, core);
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

sanafe::CorePipelineConfiguration sanafe::description_parse_core_pipeline(
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

sanafe::TilePowerMetrics sanafe::description_parse_tile_metrics(
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

void sanafe::description_parse_tile_section(
        const YAML::Node &tile_node, Architecture &arch)
{
    std::string tile_name = tile_node["name"].as<std::string>();
    std::replace(tile_name.begin(), tile_name.end(), ' ', '_');
    std::replace(tile_name.begin(), tile_name.end(), '\t', '_');

    std::pair<int, int> range = {0, 0};
    if (tile_name.find("[") != std::string::npos)
    {
        range = description_parse_range(tile_name);
    }

    for (int t = range.first; t <= range.second; t++)
    {
        const std::string name = tile_name.substr(0, tile_name.find("[")) +
                "[" + std::to_string(t) + "]";
        const TilePowerMetrics power_metrics =
                description_parse_tile_metrics(tile_node["attributes"]);
        Tile &new_tile = arch.create_tile(name, power_metrics);
        if (const YAML::Node core_node = tile_node["core"])
        {
            if (core_node.IsSequence())
            {
                for (auto c : core_node)
                {
                    description_parse_core_section(c, new_tile.id, arch);
                }
            }
            else
            {
                description_parse_core_section(core_node, new_tile.id, arch);
            }
        }
        else
        {
            const std::string error = "Error: No core section defined (line:" +
                    std::to_string(tile_node.Mark().line + 1) + ":" +
                    std::to_string(tile_node.Mark().column + 1) + ").\n";
        }
    }

    return;
}


/*
int sanafe::Architecture::set_noc_attributes(
                          const std::map<std::string, std::string> &attr)
{   noc_init = true;
    TRACE1("NoC created, mesh, width:%d height:%d.\n", noc_width, noc_height);
    return 0;
}
*/

sanafe::NetworkOnChipConfiguration sanafe::description_parse_noc_configuration(
        const YAML::Node &noc_attributes)
{
    int width_in_tiles = 1;
    int height_in_tiles = 1;
    int link_buffer_size = 0;

    if (noc_attributes["width"])
    {
        width_in_tiles = noc_attributes["width"].as<int>();
    }
    if (noc_attributes["height"])
    {
        height_in_tiles = noc_attributes["height"].as<int>();
    }
    if (noc_attributes["link_buffer_size"])
    {
        link_buffer_size = noc_attributes["link_buffer_size"].as<int>();
    }

    const NetworkOnChipConfiguration noc(
            width_in_tiles, height_in_tiles, link_buffer_size);
    return noc;
}

sanafe::Architecture sanafe::description_parse_arch_section(
        const YAML::Node &arch_node)
{
    const std::string arch_name = arch_node["name"].as<std::string>();
    if (arch_name.find("[") != std::string::npos)
    {
        throw std::runtime_error("Error: Multiple architectures not supported");
    }
    NetworkOnChipConfiguration noc =
            description_parse_noc_configuration(arch_node["attributes"]);
    Architecture new_arch(arch_name, noc);
    if (const YAML::Node tile_node = arch_node["tile"])
    {
        // Iterate through all tiles
        if (tile_node.IsSequence())
        {
            for (auto t : tile_node)
            {
                description_parse_tile_section(t, new_arch);
            }
        }
        else
        {
            description_parse_tile_section(tile_node, new_arch);
        }
    }
    else
    {
        const std::string error = "Error: No tile section defined (line:" +
                std::to_string(tile_node.Mark().line + 1) + ":" +
                std::to_string(tile_node.Mark().column + 1) + ").\n";
    }

    return new_arch;
}

sanafe::Architecture sanafe::description_parse_arch_file(std::ifstream &fp)
{
    YAML::Node top_level_yaml_node = YAML::Load(fp);
    if (YAML::Node arch_yaml_node = top_level_yaml_node["architecture"])
    {
        return description_parse_arch_section(arch_yaml_node);
    }
    else
    {
        const std::string error =
                "Error: No architecture section defined (line:" +
                std::to_string(top_level_yaml_node.Mark().line + 1) + ":" +
                std::to_string(top_level_yaml_node.Mark().column + 1) + ").\n";
        throw std::runtime_error(error);
    }
}

/*
int sanafe::description_parse_arch_file_old(std::ifstream &fp, Architecture &arch)
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
        for (auto f : fields)
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
*/

const int default_line_len = 4096;
const int default_fields = 32;

int sanafe::description_parse_net_file(
        std::ifstream &fp, class Network &net, Architecture &arch)
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
        for (auto f : fields)
        {
            TRACE1("\tField:%s\n", f.c_str());
        }
#endif
        if (fields.size() > 0)
        {
            description_read_network_entry(fields, arch, net, line_number);
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
        NeuronGroup &dest_group = net.groups_vec[dest_group_id];
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
        NeuronGroup &group = net.groups_vec[neuron_group_id];
        group_ptr = &group;
        TRACE1("Parsed neuron gid:%lu nid:%lu\n", neuron_group_id, neuron_id);
        if (neuron_set)
        {
            if (neuron_id >= group_ptr->neurons.size())
            {
                INFO("Error: Line %d: Trying to access neuron "
                     "(%d.%lu) but group %d only "
                     "allocates %lu neuron(s).\n",
                        line_number, group_ptr->id, neuron_id, group_ptr->id,
                        group_ptr->neurons.size());
                throw std::invalid_argument("Invalid neuron id");
            }
            neuron_ptr = &(group_ptr->neurons[neuron_id]);
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
        net.create_neuron_group(neuron_count, attributes);
        break;
    case 'n': // Add neuron
        neuron_ptr->set_attributes(attributes);
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
