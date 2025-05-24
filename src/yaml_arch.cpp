#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "description.hpp"
#include "pipeline.hpp"
#include "yaml_arch.hpp"

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
    parent_core.create_axon_in(std::move(name), in_metrics);
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

    bool hw_exists = false;
    for (PipelineUnitConfiguration &hw : parent_core.pipeline_hw)
    {
        if (hw.name == name)
        {
            hw_exists = true;
            hw.implements_synapse = true;
            // Merge the attributes from all sections with the same name
            //  and warn if the library is overwritten
            hw.model_info.model_attributes.merge(model.model_attributes);
            if (model.plugin_library_path.has_value())
            {
                if (hw.model_info.plugin_library_path.has_value() &&
                        hw.model_info.plugin_library_path !=
                                model.plugin_library_path)
                {
                    INFO("Warning: overwriting plugin path:%s\n",
                            model.plugin_library_path.value().c_str());
                }
                hw.model_info.plugin_library_path = model.plugin_library_path;
            }
            break;
        }
    }
    if (!hw_exists)
    {
        PipelineUnitConfiguration &synapse =
                parent_core.create_hardware_unit(std::move(name), model);
        synapse.implements_synapse = true;
    }
}

sanafe::ModelInfo sanafe::description_parse_synapse_attributes_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes)
{
    ModelInfo model_details;
    model_details.name = description_required_field<std::string>(
            parser, attributes, "model");
    model_details.model_attributes =
            description_parse_model_attributes_yaml(parser, attributes);

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
        const ryml::ConstNodeRef energy_node =
                attributes.find_child("log_energy");
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
                INFO("Dendrite plugin path found: %s\n", plugin_path.c_str());
                model_details.plugin_library_path = plugin_path;
            }
            else
            {
                throw DescriptionParsingError(
                        "Expected plugin path to be string", parser,
                        plugin_path_node);
            }
        }
        model_details.model_attributes =
                description_parse_model_attributes_yaml(parser, attributes);

        bool hw_exists = false;
        for (PipelineUnitConfiguration &hw : parent_core.pipeline_hw)
        {
            if (hw.name == name)
            {
                hw_exists = true;
                hw.implements_dendrite = true;
                // Merge the attributes from all sections with the same name
                //  and warn if the library is overwritten
                hw.model_info.model_attributes.merge(
                        model_details.model_attributes);
                if (model_details.plugin_library_path.has_value())
                {
                    if (hw.model_info.plugin_library_path.has_value() &&
                            hw.model_info.plugin_library_path !=
                                    model_details.plugin_library_path)
                    {
                        INFO("Warning: overwriting plugin path:%s\n",
                                model_details.plugin_library_path.value()
                                        .c_str());
                    }
                    hw.model_info.plugin_library_path =
                            model_details.plugin_library_path;
                }
            }
        }
        if (!hw_exists)
        {
            PipelineUnitConfiguration &dendrite =
                    parent_core.create_hardware_unit(
                            std::move(name), model_details);
            dendrite.implements_dendrite = true;
        }
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
        const std::string model_str;

        ModelInfo model_details;
        model_details.name = description_required_field<std::string>(
                parser, attributes, "model");
        const ryml::ConstNodeRef energy_node =
                attributes.find_child("log_energy");
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
        model_details.model_attributes =
                description_parse_model_attributes_yaml(parser, attributes);
        const ryml::ConstNodeRef plugin_path_node =
                attributes.find_child("plugin");
        if (!plugin_path_node.invalid())
        {
            if (plugin_path_node.has_val())
            {
                std::string plugin_path;
                plugin_path_node >> plugin_path;
                INFO("Soma plugin path found: %s\n", plugin_path.c_str());
                model_details.plugin_library_path = plugin_path;
            }
            else
            {
                throw DescriptionParsingError(
                        "Expected plugin path to be string", parser,
                        plugin_path_node);
            }
        }

        bool hw_exists = false;
        for (PipelineUnitConfiguration &hw : parent_core.pipeline_hw)
        {
            if (hw.name == name)
            {
                hw_exists = true;
                hw.implements_soma = true;
                // Merge the attributes from all sections with the same name
                //  and warn if the library is overwritten
                hw.model_info.model_attributes.merge(
                        model_details.model_attributes);
                if (model_details.plugin_library_path.has_value())
                {
                    if (hw.model_info.plugin_library_path.has_value() &&
                            hw.model_info.plugin_library_path !=
                                    model_details.plugin_library_path)
                    {
                        INFO("Warning: overwriting plugin path:%s\n",
                                model_details.plugin_library_path.value()
                                        .c_str());
                    }
                    hw.model_info.plugin_library_path =
                            model_details.plugin_library_path;
                }
                break;
            }
        }
        if (!hw_exists)
        {
            PipelineUnitConfiguration &soma = parent_core.create_hardware_unit(
                    std::move(name), model_details);
            soma.implements_soma = true;
        }
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

    parent_core.create_axon_out(std::move(axon_out_name), power_metrics);
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
        CoreConfiguration &core = arch.create_core(
                std::move(name), parent_tile_id, pipeline_config);

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
            description_required_field<std::string>(
                    parser, attributes, "buffer_position"),
            buffer_inside_unit);
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

        const TileConfiguration &new_tile =
                arch.create_tile(std::move(name), power_metrics);

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
