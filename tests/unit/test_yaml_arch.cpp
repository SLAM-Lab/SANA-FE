#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "arch.hpp"
#include "arg_parsing.hpp"
#include "yaml_arch.hpp"
#include "yaml_common.hpp"

// any helper functions go here
namespace
{
ryml::Tree parse_yaml_snippet(
        const std::string &yaml_text, ryml::Parser &parser)
{
    ryml::Tree tree = ryml::parse_in_place(
            &parser, const_cast<char *>(yaml_text.c_str()));
    return tree;
}
}

TEST(YamlArchTest, ParseAxonInAttributes_Valid)
{
    const std::string yaml = R"(
attributes:
  energy_message_in: 7.89
  latency_message_in: 0.12
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = parse_yaml_snippet(yaml, parser);
    auto node = tree["attributes"];
    sanafe::AxonInPowerMetrics axon_in_metrics =
            sanafe::yaml_parse_axon_in_attributes(parser, node);
    EXPECT_DOUBLE_EQ(axon_in_metrics.energy_message_in, 7.89);
    EXPECT_DOUBLE_EQ(axon_in_metrics.latency_message_in, 0.12);
}

TEST(YamlArchTest, ParseAxonInAttributes_Invalid)
{
    const std::string yaml = R"(
attributes:
  energy_message_in: 7.89
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = parse_yaml_snippet(yaml, parser);
    auto node = tree["attributes"];
    EXPECT_THROW(sanafe::yaml_parse_axon_in_attributes(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlArchTest, ParseAxonOutAttributes_Valid)
{
    const std::string yaml = R"(
attributes:
  energy_message_out: 7.89
  latency_message_out: 0.12
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = parse_yaml_snippet(yaml, parser);
    auto node = tree["attributes"];
    sanafe::AxonOutPowerMetrics axon_out_metrics =
            sanafe::yaml_parse_axon_out_attributes(parser, node);
    EXPECT_DOUBLE_EQ(axon_out_metrics.energy_message_out, 7.89);
    EXPECT_DOUBLE_EQ(axon_out_metrics.latency_message_out, 0.12);
}

TEST(YamlArchTest, ParseAxonOutAttributes_Invalid)
{
    const std::string yaml = R"(
attributes:
  energy_message_out: 7.89
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = parse_yaml_snippet(yaml, parser);
    auto node = tree["attributes"];
    EXPECT_THROW(sanafe::yaml_parse_axon_out_attributes(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlArchTest, ParseProcessingUnitAttributesWithPlugin)
{
    const std::string yaml = R"(
attributes:
  model: "testmodel"
  log_energy: true
  log_latency: false
  plugin: "plugin.so"
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = parse_yaml_snippet(yaml, parser);
    auto node = tree["attributes"];

    auto result = sanafe::yaml_parse_processing_unit_attributes(parser, node);
    EXPECT_EQ(result.name, "testmodel");
    EXPECT_TRUE(result.log_energy);
    EXPECT_FALSE(result.log_latency);
    ASSERT_TRUE(result.plugin_library_path.has_value());
    EXPECT_EQ(result.plugin_library_path.value(), "plugin.so");
}

TEST(YamlArchTest, DescriptionParseTileMetricsYaml_Valid)
{
    const std::string yaml = R"(
energy_north_hop: 1.0
latency_north_hop: 2.0
energy_east_hop: 3.0
latency_east_hop: 4.0
energy_south_hop: 5.0
latency_south_hop: 6.0
energy_west_hop: 7.0
latency_west_hop: 8.0
log_energy: true
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    auto result = sanafe::description_parse_tile_metrics_yaml(parser, node);
    EXPECT_DOUBLE_EQ(result.energy_north_hop, 1.0);
    EXPECT_DOUBLE_EQ(result.latency_north_hop, 2.0);
    EXPECT_DOUBLE_EQ(result.energy_east_hop, 3.0);
    EXPECT_DOUBLE_EQ(result.latency_east_hop, 4.0);
    EXPECT_DOUBLE_EQ(result.energy_south_hop, 5.0);
    EXPECT_DOUBLE_EQ(result.latency_south_hop, 6.0);
    EXPECT_DOUBLE_EQ(result.energy_west_hop, 7.0);
    EXPECT_DOUBLE_EQ(result.latency_west_hop, 8.0);
    EXPECT_TRUE(result.log_energy);
}

TEST(YamlArchTest, ParsesBasicArchitecture)
{
    const std::string yaml = R"(
architecture:
  name: minimal_arch
  attributes:
    link_buffer_size: 1
    width: 1
    height: 1
  tile:
    - name: tile0
      attributes:
        energy_north_hop: 1.0
        latency_north_hop: 1.0
        energy_east_hop: 1.0
        latency_east_hop: 1.0
        energy_south_hop: 1.0
        latency_south_hop: 1.0
        energy_west_hop: 1.0
        latency_west_hop: 1.0
      core:
        - name: core0
          attributes:
            buffer_position: soma
            max_neurons_supported: 10
          axon_in:
            - name: axin
              attributes:
                energy_message_in: 0.0
                latency_message_in: 0.0
          synapse:
            - name: syn
              attributes:
                model: current_based
                energy_process_spike: 1.0
                latency_process_spike: 1.0
          dendrite:
            - name: dend
              attributes:
                model: accumulator
                energy_update: 0.0
                latency_update: 0.0
                update_every_timestep: true
          soma:
            - name: soma
              attributes:
                model: leaky_integrate_fire
                energy_access_neuron: 1.0
                latency_access_neuron: 1.0
                energy_update_neuron: 1.0
                latency_update_neuron: 1.0
                energy_spike_out: 1.0
                latency_spike_out: 1.0
          axon_out:
            - name: axout
              attributes:
                energy_message_out: 1.0
                latency_message_out: 1.0
)";

    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    const ryml::ConstNodeRef root = tree.rootref();

    ASSERT_TRUE(root.has_child("architecture"));
    sanafe::Architecture arch = sanafe::description_parse_arch_section_yaml(
            parser, root["architecture"]);

    EXPECT_EQ(arch.tiles.size(), 1);
    EXPECT_EQ(arch.core_count, 1);
    EXPECT_EQ(arch.name, "minimal_arch");
    EXPECT_EQ(arch.noc_width_in_tiles, 1);
    EXPECT_EQ(arch.noc_height_in_tiles, 1);
    EXPECT_EQ(arch.noc_buffer_size, 1);
    EXPECT_EQ(arch.cores().at(0).get().name, "core0[0]");
    // EXPECT_EQ(arch.cores().at(0).get().name, "core0");

    const auto &core = arch.cores().at(0).get();
    EXPECT_EQ(core.name, "core0[0]");

    EXPECT_EQ(core.axon_in.size(), 1);
    EXPECT_EQ(core.axon_in[0].name, "axin");
    EXPECT_DOUBLE_EQ(core.axon_in[0].metrics.energy_message_in, 0.0);
    EXPECT_DOUBLE_EQ(core.axon_in[0].metrics.latency_message_in, 0.0);

    EXPECT_EQ(core.axon_out.size(), 1);
    EXPECT_EQ(core.axon_out[0].name, "axout");
    EXPECT_DOUBLE_EQ(core.axon_out[0].metrics.energy_message_out, 1.0);
    EXPECT_DOUBLE_EQ(core.axon_out[0].metrics.latency_message_out, 1.0);

    const auto &pipeline_hw = core.pipeline_hw;
    ASSERT_EQ(pipeline_hw.size(), 3);

    EXPECT_EQ(pipeline_hw[0].name, "syn");
    EXPECT_EQ(pipeline_hw[0].model_info.name, "current_based");
    EXPECT_TRUE(pipeline_hw[0].implements_synapse);
    EXPECT_DOUBLE_EQ(pipeline_hw[0].model_info.model_attributes.at(
                             "energy_process_spike"),
            1.0);
    EXPECT_DOUBLE_EQ(pipeline_hw[0].model_info.model_attributes.at(
                             "latency_process_spike"),
            1.0);

    EXPECT_EQ(pipeline_hw[1].name, "dend");
    EXPECT_EQ(pipeline_hw[1].model_info.name, "accumulator");
    EXPECT_TRUE(pipeline_hw[1].implements_dendrite);
    EXPECT_DOUBLE_EQ(
            pipeline_hw[1].model_info.model_attributes.at("energy_update"),
            0.0);
    EXPECT_DOUBLE_EQ(
            pipeline_hw[1].model_info.model_attributes.at("latency_update"),
            0.0);
    EXPECT_TRUE(pipeline_hw[1].model_info.model_attributes.at(
            "update_every_timestep"));

    EXPECT_EQ(pipeline_hw[2].name, "soma");
    EXPECT_EQ(pipeline_hw[2].model_info.name, "leaky_integrate_fire");
    EXPECT_TRUE(pipeline_hw[2].implements_soma);
    EXPECT_DOUBLE_EQ(pipeline_hw[2].model_info.model_attributes.at(
                             "energy_access_neuron"),
            1.0);
    EXPECT_DOUBLE_EQ(pipeline_hw[2].model_info.model_attributes.at(
                             "latency_access_neuron"),
            1.0);
    EXPECT_DOUBLE_EQ(pipeline_hw[2].model_info.model_attributes.at(
                             "energy_update_neuron"),
            1.0);
    EXPECT_DOUBLE_EQ(pipeline_hw[2].model_info.model_attributes.at(
                             "latency_update_neuron"),
            1.0);
    EXPECT_DOUBLE_EQ(
            pipeline_hw[2].model_info.model_attributes.at("energy_spike_out"),
            1.0);
    EXPECT_DOUBLE_EQ(
            pipeline_hw[2].model_info.model_attributes.at("latency_spike_out"),
            1.0);
}

TEST(YamlArchTest, ParsesTileRangeNotation)
{
    const std::string yaml = R"(
architecture:
  name: range_test_arch
  attributes:
    link_buffer_size: 1
    width: 3
    height: 1
  tile:
    - name: tile[0..2]
      attributes:
        energy_north_hop: 1.0
        latency_north_hop: 1.0
        energy_east_hop: 1.0
        latency_east_hop: 1.0
        energy_south_hop: 1.0
        latency_south_hop: 1.0
        energy_west_hop: 1.0
        latency_west_hop: 1.0
      core:
        - name: core0
          attributes:
            buffer_position: soma
            max_neurons_supported: 10
          axon_in:
            - name: axin
              attributes:
                energy_message_in: 0.0
                latency_message_in: 0.0
          synapse:
            - name: syn
              attributes:
                model: current_based
                energy_process_spike: 1.0
                latency_process_spike: 1.0
          dendrite:
            - name: dend
              attributes:
                model: accumulator
                energy_update: 0.0
                latency_update: 0.0
          soma:
            - name: soma
              attributes:
                model: leaky_integrate_fire
                energy_access_neuron: 1.0
                latency_access_neuron: 1.0
                energy_update_neuron: 1.0
                latency_update_neuron: 1.0
                energy_spike_out: 1.0
                latency_spike_out: 1.0
          axon_out:
            - name: axout
              attributes:
                energy_message_out: 1.0
                latency_message_out: 1.0
)";

    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    const ryml::ConstNodeRef root = tree.rootref();

    ASSERT_TRUE(root.has_child("architecture"));
    sanafe::Architecture arch = sanafe::description_parse_arch_section_yaml(
            parser, root["architecture"]);

    // tile[0..2] should expand to 3 tiles
    EXPECT_EQ(arch.tiles.size(), 3);
    EXPECT_EQ(arch.name, "range_test_arch");

    EXPECT_EQ(arch.core_count, 3);

    EXPECT_EQ(arch.tiles[0].name, "tile[0]");
    EXPECT_EQ(arch.tiles[1].name, "tile[1]");
    EXPECT_EQ(arch.tiles[2].name, "tile[2]");
}


TEST(YamlArchTest, ParsesCoreRangeNotation)
{
    const std::string yaml = R"(
architecture:
  name: core_range_arch
  attributes:
    link_buffer_size: 1
    width: 1
    height: 1
  tile:
    - name: tile0
      attributes:
        energy_north_hop: 1.0
        latency_north_hop: 1.0
        energy_east_hop: 1.0
        latency_east_hop: 1.0
        energy_south_hop: 1.0
        latency_south_hop: 1.0
        energy_west_hop: 1.0
        latency_west_hop: 1.0
      core:
        - name: core[0..3]
          attributes:
            buffer_position: soma
            max_neurons_supported: 10
          axon_in:
            - name: axin
              attributes:
                energy_message_in: 0.0
                latency_message_in: 0.0
          synapse:
            - name: syn
              attributes:
                model: current_based
                energy_process_spike: 1.0
                latency_process_spike: 1.0
          dendrite:
            - name: dend
              attributes:
                model: accumulator
                energy_update: 0.0
                latency_update: 0.0
          soma:
            - name: soma
              attributes:
                model: leaky_integrate_fire
                energy_access_neuron: 1.0
                latency_access_neuron: 1.0
                energy_update_neuron: 1.0
                latency_update_neuron: 1.0
                energy_spike_out: 1.0
                latency_spike_out: 1.0
          axon_out:
            - name: axout
              attributes:
                energy_message_out: 1.0
                latency_message_out: 1.0
)";

    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    const ryml::ConstNodeRef root = tree.rootref();

    ASSERT_TRUE(root.has_child("architecture"));
    sanafe::Architecture arch = sanafe::description_parse_arch_section_yaml(
            parser, root["architecture"]);

    // 1 tile with core[0..3] should expand to 4 cores
    EXPECT_EQ(arch.tiles.size(), 1);
    EXPECT_EQ(arch.core_count, 4);

    const auto &cores = arch.cores();
    EXPECT_EQ(cores.at(0).get().name, "core[0]");
    EXPECT_EQ(cores.at(1).get().name, "core[1]");
    EXPECT_EQ(cores.at(2).get().name, "core[2]");
    EXPECT_EQ(cores.at(3).get().name, "core[3]");
}

TEST(YamlArchTest, MissingTileSectionThrows)
{
    const std::string yaml = R"(
architecture:
  name: missing_tile_arch
  attributes:
    link_buffer_size: 1
    width: 1
    height: 1
)";

    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    const ryml::ConstNodeRef root = tree.rootref();

    EXPECT_THROW(sanafe::description_parse_arch_section_yaml(
                         parser, root["architecture"]),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlArchTest, MissingCoreSectionThrows)
{
    const std::string yaml = R"(
architecture:
  name: missing_core_arch
  attributes:
    link_buffer_size: 1
    width: 1
    height: 1
  tile:
    - name: tile0
      attributes:
        energy_north_hop: 1.0
        latency_north_hop: 1.0
        energy_east_hop: 1.0
        latency_east_hop: 1.0
        energy_south_hop: 1.0
        latency_south_hop: 1.0
        energy_west_hop: 1.0
        latency_west_hop: 1.0
)";

    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    const ryml::ConstNodeRef root = tree.rootref();

    EXPECT_THROW(sanafe::description_parse_arch_section_yaml(
                         parser, root["architecture"]),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlArchTest, MissingSomaSectionThrows)
{
    const std::string yaml = R"(
architecture:
  name: missing_soma_arch
  attributes:
    link_buffer_size: 1
    width: 1
    height: 1
  tile:
    - name: tile0
      attributes:
        energy_north_hop: 1.0
        latency_north_hop: 1.0
        energy_east_hop: 1.0
        latency_east_hop: 1.0
        energy_south_hop: 1.0
        latency_south_hop: 1.0
        energy_west_hop: 1.0
        latency_west_hop: 1.0
      core:
        - name: core0
          attributes:
            buffer_position: soma
            max_neurons_supported: 10
          axon_in:
            - name: axin
              attributes:
                energy_message_in: 0.0
                latency_message_in: 0.0
          synapse:
            - name: syn
              attributes:
                model: current_based
                energy_process_spike: 1.0
                latency_process_spike: 1.0
          dendrite:
            - name: dend
              attributes:
                model: accumulator
                energy_update: 0.0
                latency_update: 0.0
          axon_out:
            - name: axout
              attributes:
                energy_message_out: 1.0
                latency_message_out: 1.0
)";

    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    const ryml::ConstNodeRef root = tree.rootref();

    EXPECT_THROW(sanafe::description_parse_arch_section_yaml(
                         parser, root["architecture"]),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlArchTest, LoadArchFromFile_FileNotOpen)
{
    std::ifstream bad_stream;
    EXPECT_THROW(sanafe::description_parse_arch_file_yaml(bad_stream),
            std::runtime_error);
}

TEST(YamlArchTest, LoadArchFromFile_ValidFile)
{
    std::filesystem::path path(SANAFE_ROOT_PATH);
    sanafe::Architecture arch =
            sanafe::load_arch(path.string() + "/arch/example_chip.yaml");

    EXPECT_EQ(arch.name, "demo");
    EXPECT_EQ(arch.tiles.size(), 2);
    EXPECT_EQ(arch.noc_width_in_tiles, 2);
    EXPECT_EQ(arch.noc_height_in_tiles, 1);
}

TEST(YamlArchTest, LoadArchFromFile_VerifiesNestedStructure)
{
    std::filesystem::path path(SANAFE_ROOT_PATH);
    sanafe::Architecture arch =
            sanafe::load_arch(path.string() + "/arch/example_chip.yaml");

    // Tiles and cores expanded correctly + attributes
    EXPECT_EQ(arch.tiles.size(), 2);
    EXPECT_EQ(arch.core_count, 8);

    EXPECT_EQ(arch.tiles[0].name, "demo_tile[0]");
    EXPECT_EQ(arch.tiles[0].cores[0].name, "demo_core[0]");

    EXPECT_DOUBLE_EQ(arch.tiles[0].power_metrics.energy_north_hop, 2.0e-12);
    EXPECT_DOUBLE_EQ(arch.tiles[0].cores[0].pipeline_hw[0].model_info.model_attributes.at("energy_process_spike"), 20.0e-12); // Hardware pipeline units via pipeline_hw vector
    EXPECT_DOUBLE_EQ(arch.tiles[0].cores[0].axon_out[0].metrics.latency_message_out, 5.0e-9);
}