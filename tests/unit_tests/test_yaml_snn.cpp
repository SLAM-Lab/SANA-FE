#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "arch.hpp"
#include "yaml_common.hpp"
#include "yaml_snn.hpp"
#include "network.hpp"

// any helper functions go here
namespace {
    ryml::Tree parse_yaml_snippet(const std::string &yaml_text, ryml::Parser& parser) {
        ryml::Tree tree = ryml::parse_in_place(&parser, const_cast<char*>(yaml_text.c_str()));
        return tree;
    }
} 

TEST(YamlSnnTest, ParseEdgeDescription_Valid) {
    auto [src, tgt] = sanafe::description_parse_edge_description("A.1 -> B.2");

    EXPECT_EQ(src.group_name, "A");
    EXPECT_EQ(src.neuron_offset, 1);
    EXPECT_EQ(tgt.group_name, "B");
    EXPECT_EQ(tgt.neuron_offset, 2);
}

TEST(YamlSnnTest, ParseEdgeDescription_MissingDotThrows) {
    EXPECT_THROW(sanafe::description_parse_edge_description("A -> B.2"), std::runtime_error);
    EXPECT_THROW(sanafe::description_parse_edge_description("A.1 -> B"), std::runtime_error);
}

TEST(YamlSnnTest, CountNeurons_WithRangesAndSingles) {
    const std::string yaml = R"(
- 0..2
- 5
- 10..12
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    ryml::Tree tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    size_t count = sanafe::description_count_neurons(parser, node);
    EXPECT_EQ(count, 3 + 1 + 3);  // 0,1,2 + 5 + 10,11,12
}

TEST(YamlSnnTest, CountNeurons_InvalidFormatThrows) {
    const std::string yaml = R"(
invalid: stuff
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();
    EXPECT_THROW(sanafe::description_count_neurons(parser, node), sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseFullNetworkSection) {
    const std::string yaml = R"(
  name: example
  groups:
    - name: Input
      neurons:
        - 0..1
    - name: Output
      neurons:
        - 0..1
  edges:
    - Input.0 -> Output.0: [weight: -1.0]
    - Input.1 -> Output.1: [weight: -2.0]
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    sanafe::SpikingNetwork net = sanafe::yaml_parse_network_section(parser, node);
    ASSERT_EQ(net.groups.size(), 2);
    ASSERT_TRUE(net.groups.find("Input") != net.groups.end());
    ASSERT_TRUE(net.groups.find("Output") != net.groups.end());

    EXPECT_EQ(net.groups.at("Input").neurons.size(), 2);
    EXPECT_EQ(net.groups.at("Output").neurons.size(), 2);

    const auto& input0 = net.groups.at("Input").neurons[0];
    const auto& input1 = net.groups.at("Input").neurons[1];

    ASSERT_EQ(input0.edges_out.size(), 1);
    ASSERT_EQ(input1.edges_out.size(), 1);

    EXPECT_EQ(input0.edges_out[0].post_neuron.group_name, "Output");
    EXPECT_EQ(input1.edges_out[0].post_neuron.group_name, "Output");

    EXPECT_EQ(input0.edges_out[0].post_neuron.neuron_offset, 0);
    EXPECT_EQ(input1.edges_out[0].post_neuron.neuron_offset, 1);
}