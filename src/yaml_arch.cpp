
#include <cstddef>
#include <fstream>
#include <ios>
#include <iosfwd>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <c4/yml/event_handler_tree.hpp>
#include <c4/yml/fwd.hpp>
#include <c4/yml/node.hpp>
#include <c4/yml/parse.hpp>
#include <c4/yml/tree.hpp>
#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "arch.hpp"
#include "pipeline.hpp"
#include "print.hpp"
#include "utils.hpp"
#include "yaml_arch.hpp"
#include "yaml_common.hpp"

constexpr std::string_view range_delimiter = "..";

void sanafe::yaml_parse_axon_in(const ryml::Parser &parser,
        const ryml::ConstNodeRef axon_in_node, CoreConfiguration &parent_core,
        const std::string_view & /*type*/, const std::string &name)
{
    const ryml::ConstNodeRef &attributes =
            axon_in_node.find_child("attributes");
    if (attributes.invalid())
    {
        throw YamlDescriptionParsingError(
                "No attributes section defined", parser, axon_in_node);
    }
    const AxonInPowerMetrics in_metrics =
            yaml_parse_axon_in_attributes(parser, attributes);
    parent_core.create_axon_in(name, in_metrics);
}

void sanafe::yaml_parse_axon_out(const ryml::Parser &parser,
        const ryml::ConstNodeRef axon_out_node, CoreConfiguration &parent_core,
        const std::string_view & /*type*/, const std::string &name)
{
    const auto &attributes = axon_out_node.find_child("attributes");
    if (attributes.invalid())
    {
        throw YamlDescriptionParsingError(
                "No attributes section defined", parser, axon_out_node);
    }
    const AxonOutPowerMetrics out_metrics =
            yaml_parse_axon_out_attributes(parser, attributes);

    parent_core.create_axon_out(name, out_metrics);
}

// Generic helper for parsing model info attributes
void sanafe::yaml_parse_processing_unit(const ryml::Parser &parser,
        const ryml::ConstNodeRef node, CoreConfiguration &parent_core,
        const std::string_view &type, const std::string &unit_name)
{
    auto model = yaml_parse_processing_unit_attributes(
            parser, node.find_child("attributes"));

    yaml_merge_or_create_hardware_unit(parent_core, unit_name, model, type);
}

// NOLINTNEXTLINE(misc-include-cleaner)
sanafe::AxonInPowerMetrics sanafe::yaml_parse_axon_in_attributes(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    AxonInPowerMetrics axon_in_metrics;
    axon_in_metrics.energy_message_in = yaml_required_field<double>(
            parser, attributes, "energy_message_in");
    axon_in_metrics.latency_message_in = yaml_required_field<double>(
            parser, attributes, "latency_message_in");

    return axon_in_metrics;
}

sanafe::AxonOutPowerMetrics sanafe::yaml_parse_axon_out_attributes(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    AxonOutPowerMetrics axon_out_metrics;
    axon_out_metrics.energy_message_out = yaml_required_field<double>(
            parser, attributes, "energy_message_out");
    axon_out_metrics.latency_message_out = yaml_required_field<double>(
            parser, attributes, "latency_message_out");

    return axon_out_metrics;
}

sanafe::ModelInfo sanafe::yaml_parse_processing_unit_attributes(
        const ryml::Parser &parser, const ryml::ConstNodeRef &attributes)
{
    ModelInfo model_details;
    model_details.name =
            yaml_required_field<std::string>(parser, attributes, "model");

    // NOLINTNEXTLINE(misc-include-cleaner)
    const ryml::ConstNodeRef energy_node = attributes.find_child("log_energy");
    if (!energy_node.invalid())
    {
        energy_node >> model_details.log_energy;
    }

    const ryml::ConstNodeRef latency_node =
            attributes.find_child("log_latency");
    if (!latency_node.invalid())
    {
        latency_node >> model_details.log_latency;
    }

    // Handle plugin path NOLINTNEXTLINE(misc-include-cleaner)
    const ryml::ConstNodeRef plugin_path_node = attributes.find_child("plugin");
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
            throw YamlDescriptionParsingError(
                    "Expected plugin path to be string", parser,
                    plugin_path_node);
        }
    }

    model_details.model_attributes =
            description_parse_model_attributes_yaml(parser, attributes);
    return model_details;
}

void sanafe::yaml_merge_or_create_hardware_unit(CoreConfiguration &parent_core,
        const std::string &name, ModelInfo &model_details,
        const std::string_view &section)
{
    bool hw_exists = false;
    for (PipelineUnitConfiguration &hw : parent_core.pipeline_hw)
    {
        if (hw.name == name)
        {
            hw_exists = true;
            yaml_set_implements_flag(hw, section);

            // Merge attributes and handle plugin path conflicts
            hw.model_info.model_attributes.merge(
                    model_details.model_attributes);
            if (model_details.plugin_library_path.has_value())
            {
                if (hw.model_info.plugin_library_path.has_value() &&
                        hw.model_info.plugin_library_path !=
                                model_details.plugin_library_path)
                {
                    INFO("Warning: overwriting plugin path:%s\n",
                            model_details.plugin_library_path.value().c_str());
                }
                hw.model_info.plugin_library_path =
                        model_details.plugin_library_path;
            }
            break;
        }
    }

    if (!hw_exists)
    {
        PipelineUnitConfiguration &new_unit =
                parent_core.create_hardware_unit(name, model_details);
        yaml_set_implements_flag(new_unit, section);
    }
}

template <typename ParseFunc>
void sanafe::yaml_parse_pipeline_entry(const ryml::Parser &parser,
        const ryml::ConstNodeRef &unit_node, CoreConfiguration &parent_core,
        const std::string_view &type, ParseFunc parsing_function)
{
    auto name = yaml_required_field<std::string>(parser, unit_node, "name");
    std::pair<int, int> range = {0, 0};

    // Check if name contains range notation (e.g., "foo[0..3]")
    if (name.find(range_delimiter) != std::string::npos)
    {
        range = yaml_parse_range(name);
    }

    // Parse the same entry for each unit in the range
    // Note: We re-parse the YAML node for each iteration rather than
    // copying objects for simplicity, as the iteration count is typically small
    for (int i = range.first; i <= range.second; ++i)
    {
        // Generate unique name for each unit in the range
        std::string unit_name(name);
        if (name.find(range_delimiter) != std::string::npos)
        {
            unit_name = name.substr(0, name.find('[')) + '[' +
                    std::to_string(i) + ']';
        }

        // Delegate to the appropriate parsing function
        parsing_function(parser, unit_node, parent_core, type, unit_name);
    }
}

void sanafe::yaml_set_implements_flag(
        PipelineUnitConfiguration &hw, const std::string_view &section)
{
    // Sets the appropriate implementation flag based on the section type
    //  This determines which type of processing unit (synapse/dendrite/soma) is
    //  being configured. Note that a hardware unit can configured to implement
    //  multiple units' functionality, in which case we would see the same unit
    //  defined in multiple YAML sections
    if (section == "synapse")
    {
        hw.implements_synapse = true;
    }
    else if (section == "dendrite")
    {
        hw.implements_dendrite = true;
    }
    else if (section == "soma")
    {
        hw.implements_soma = true;
    }
    else
    {
        throw std::runtime_error("Section not recognized");
    }
}

void sanafe::description_parse_core_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef core_node, const size_t parent_tile_id,
        Architecture &arch, const std::string_view &name)
{
    // Parses a complete core configuration from YAML
    //  A core contains multiple types of pipeline units (axon_in, synapse,
    //  dendrite, soma, axon_out)
    const CorePipelineConfiguration pipeline_config =
            description_parse_core_pipeline_yaml(
                    parser, core_node["attributes"]);
    CoreConfiguration &core = arch.create_core(
            std::string(name), parent_tile_id, pipeline_config);

    // Define required hardware sections with their parsing functions
    const std::vector<PipelineUnitSectionInfo> required_sections = {
            {"axon_in", yaml_parse_axon_in},
            // Synapse, dendrite and soma units share a common parsing routine
            {"synapse", yaml_parse_processing_unit},
            {"dendrite", yaml_parse_processing_unit},
            {"soma", yaml_parse_processing_unit},
            {"axon_out", yaml_parse_axon_out}};

    for (const auto &section : required_sections)
    {
        const ryml::ConstNodeRef section_node =
                core_node.find_child(std::string(section.name).c_str());
        if (section_node.invalid())
        {
            const std::string error = std::string("No ") +
                    std::string(section.name) + " section defined";
            throw YamlDescriptionParsingError(error, parser, core_node);
        }

        if (section_node.is_seq())
        {
            for (const auto &item_node : section_node)
            {
                yaml_parse_pipeline_entry(parser, item_node, core, section.name,
                        section.parsing_function);
            }
        }
        else
        {
            yaml_parse_pipeline_entry(parser, section_node, core, section.name,
                    section.parsing_function);
        }
    }
}

void sanafe::description_parse_core_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef core_node, const size_t parent_tile_id,
        Architecture &arch)
{
    auto core_name =
            yaml_required_field<std::string>(parser, core_node, "name");
    std::pair<int, int> core_range = {0, 0};

    if (core_name.find(range_delimiter) != std::string::npos)
    {
        core_range = yaml_parse_range(core_name);
    }

    for (int c = core_range.first; c <= core_range.second; c++)
    {
        const std::string name = core_name.substr(0, core_name.find('[')) +
                '[' + std::to_string(c) + ']';
        description_parse_core_yaml(
                parser, core_node, parent_tile_id, arch, name);
    }
}

sanafe::CorePipelineConfiguration sanafe::description_parse_core_pipeline_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    CorePipelineConfiguration pipeline_config{};

    bool buffer_inside_unit = false;

    const ryml::ConstNodeRef buffer_inside_node =
            attributes.find_child("buffer_inside_unit");
    if (!buffer_inside_node.invalid())
    {
        buffer_inside_node >> buffer_inside_unit;
    }
    const ryml::ConstNodeRef energy_node = attributes.find_child("log_energy");
    if (!energy_node.invalid())
    {
        energy_node >> pipeline_config.log_energy;
    }
    const ryml::ConstNodeRef latency_node =
            attributes.find_child("log_latency");
    if (!latency_node.invalid())
    {
        latency_node >> pipeline_config.log_latency;
    }

    pipeline_config.buffer_position = pipeline_parse_buffer_pos_str(
            yaml_required_field<std::string>(
                    parser, attributes, "buffer_position"),
            buffer_inside_unit);
    pipeline_config.max_neurons_supported = yaml_required_field<int>(
            parser, attributes, "max_neurons_supported");

    return pipeline_config;
}

sanafe::TilePowerMetrics sanafe::description_parse_tile_metrics_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    TilePowerMetrics tile_metrics;

    tile_metrics.energy_north_hop =
            yaml_required_field<double>(parser, attributes, "energy_north_hop");
    tile_metrics.latency_north_hop = yaml_required_field<double>(
            parser, attributes, "latency_north_hop");

    tile_metrics.energy_east_hop =
            yaml_required_field<double>(parser, attributes, "energy_east_hop");
    tile_metrics.latency_east_hop =
            yaml_required_field<double>(parser, attributes, "latency_east_hop");

    tile_metrics.energy_south_hop =
            yaml_required_field<double>(parser, attributes, "energy_south_hop");
    tile_metrics.latency_south_hop = yaml_required_field<double>(
            parser, attributes, "latency_south_hop");

    tile_metrics.energy_west_hop =
            yaml_required_field<double>(parser, attributes, "energy_west_hop");
    tile_metrics.latency_west_hop =
            yaml_required_field<double>(parser, attributes, "latency_west_hop");

    if (!attributes.find_child("log_energy").invalid())
    {
        const ryml::ConstNodeRef energy = attributes["log_energy"];
        energy >> tile_metrics.log_energy;
    }
    if (!attributes.find_child("log_latency").invalid())
    {
        const ryml::ConstNodeRef latency = attributes["log_latency"];
        latency >> tile_metrics.log_latency;
    }

    return tile_metrics;
}

void sanafe::description_parse_tile_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef tile_node, Architecture &arch)
{
    std::string tile_name;
    tile_node["name"] >> tile_name;
    std::pair<int, int> range = {0, 0};

    if (tile_name.find(range_delimiter) != std::string::npos)
    {
        range = yaml_parse_range(tile_name);
    }

    for (int t = range.first; t <= range.second; t++)
    {
        std::string name = tile_name.substr(0, tile_name.find('[')) +
                "[" + std::to_string(t) + "]";
        const TilePowerMetrics power_metrics =
                description_parse_tile_metrics_yaml(
                        parser, tile_node["attributes"]);

        const TileConfiguration &new_tile =
                arch.create_tile(std::move(name), power_metrics);

        if (tile_node.find_child("core").invalid())
        {
            const std::string error = "No core section defined";
            throw YamlDescriptionParsingError(error, parser, tile_node);
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
            yaml_required_field<int>(parser, noc_attributes, "width");
    noc.height_in_tiles =
            yaml_required_field<int>(parser, noc_attributes, "height");
    noc.link_buffer_size = yaml_required_field<int>(
            parser, noc_attributes, "link_buffer_size");

    const std::string model_type =
            yaml_optional_field<std::string>(noc_attributes, "sync_model")
                    .value_or("fixed");
    noc.ts_sync_delay_table =
            yaml_parse_sync_delay_table(parser, noc_attributes, model_type);

    return noc;
}

sanafe::LookupTable<double> sanafe::yaml_parse_sync_delay_table(
        const ryml::Parser &parser, const ryml::ConstNodeRef &noc_attributes,
        const std::string &model_type)
{
    LookupTable<double> table;

    if (model_type == "fixed")
    {
        // Single value, defaults to 0.0 if not provided
        const double delay =
                yaml_optional_field<double>(noc_attributes, "latency_sync")
                        .value_or(0.0);
        table.values[0] = delay;
    }
    else if (model_type == "table")
    {
        const ryml::ConstNodeRef delay_node =
                noc_attributes.find_child("latency_sync");
        if (delay_node.invalid())
        {
            throw YamlDescriptionParsingError(
                    "Attribute 'latency_sync' required when "
                    "'table' synchronization model is chosen.",
                    parser, noc_attributes);
        }

        if (delay_node.is_seq()) // YAML list
        {
            // Array of values - tile index is implicit (0, 1, 2, ...)
            size_t tile_index = 0;
            for (const auto &value_node : delay_node)
            {
                double delay_value = 0.0;
                value_node >> delay_value;
                table.values[tile_index++] = delay_value;
            }
        }
        else if (delay_node.is_map()) // YAML mapping
        {
            // Mapping of tile_id -> delay_value
            for (const auto &pair : delay_node)
            {
                size_t tile_id = 0UL;
                double delay_value = 0.0;
                pair >> ryml::key(tile_id);
                pair >> delay_value;
                table.values[tile_id] = delay_value;
            }
        }
        else // if single YAML scalar value given
        {
            double delay_value = 0.0;
            delay_node >> delay_value;
            table.values[0] = delay_value;
        }
        INFO("Setting (%zu) sync latency values \n", table.values.size());
    }
    else
    {
        throw YamlDescriptionParsingError(
                "Unknown sync_model: " + model_type, parser, noc_attributes);
    }

    return table;
}

sanafe::Architecture sanafe::description_parse_arch_section_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef arch_node)
{
    if (arch_node.invalid())
    {
        throw YamlDescriptionParsingError(
                "No top-level architecture section defined.\n", parser,
                arch_node);
    }
    std::string arch_name;
    arch_node["name"] >> arch_name;
    if (arch_name.find('[') != std::string::npos)
    {
        throw YamlDescriptionParsingError(
                "Multiple architectures not supported", parser, arch_node);
    }
    const NetworkOnChipConfiguration noc =
            description_parse_noc_configuration_yaml(
                    parser, arch_node["attributes"]);
    Architecture new_arch(std::move(arch_name), noc);
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
        throw YamlDescriptionParsingError(
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
    const std::streampos file_size = fp.tellg();
    fp.seekg(0, std::ios::beg);

    // Allocate memory
    std::string file_content;
    file_content.reserve(file_size);

    // Read the file
    file_content.assign((std::istreambuf_iterator<char>(fp)),
            std::istreambuf_iterator<char>());
    fp.close();

    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    INFO("Loading YAML information from file.\n");
    ryml::Tree top_level_yaml =
            ryml::parse_in_place(&parser, file_content.data());
    // NOLINTEND(misc-include-cleaner)
    INFO("YAML information loaded from file.\n");

    if (top_level_yaml["architecture"].invalid())
    {
        throw YamlDescriptionParsingError(
                "No architecture section defined", parser, top_level_yaml);
    }
    return description_parse_arch_section_yaml(
            parser, top_level_yaml["architecture"]);
}
