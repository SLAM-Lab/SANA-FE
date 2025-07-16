#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "arg_parsing.hpp"
#include "arch.hpp"
#include "yaml_common.hpp"
#include "yaml_arch.hpp"

// any helper functions go here
namespace {
    ryml::Tree parse_yaml_snippet(const std::string &yaml_text, ryml::Parser& parser) {
        ryml::Tree tree = ryml::parse_in_place(&parser, const_cast<char*>(yaml_text.c_str()));
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
    sanafe::AxonInPowerMetrics axon_in_metrics = sanafe::yaml_parse_axon_in_attributes(parser, node);
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
    EXPECT_THROW(sanafe::yaml_parse_axon_in_attributes(parser, node), sanafe::YamlDescriptionParsingError);
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
    sanafe::AxonOutPowerMetrics axon_out_metrics = sanafe::yaml_parse_axon_out_attributes(parser, node);
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
    EXPECT_THROW(sanafe::yaml_parse_axon_out_attributes(parser, node), sanafe::YamlDescriptionParsingError);
}

TEST(YamlArchTest, ParseProcessingUnitAttributesWithPlugin) {
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

TEST(YamlArchTest, DescriptionParseTileMetricsYaml_Valid) {
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
log_latency: false
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
    EXPECT_FALSE(result.log_latency);
}