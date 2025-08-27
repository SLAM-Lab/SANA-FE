#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>

#include "arch.hpp"
#include "network.hpp"
#include "yaml_common.hpp"
#include "yaml_snn.hpp"

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

TEST(YamlSnnTest, ParseEdgeDescription_Valid)
{
    auto [src, tgt] = sanafe::description_parse_edge_description("A.1 -> B.2");

    EXPECT_EQ(src.group_name, "A");
    EXPECT_EQ(src.neuron_offset, 1);
    EXPECT_EQ(tgt.group_name, "B");
    EXPECT_EQ(tgt.neuron_offset, 2);
}

TEST(YamlSnnTest, ParseEdgeDescription_MissingDotThrows)
{
    EXPECT_THROW(sanafe::description_parse_edge_description("A -> B.2"),
            std::runtime_error);
    EXPECT_THROW(sanafe::description_parse_edge_description("A.1 -> B"),
            std::runtime_error);
}

TEST(YamlSnnTest, CountNeurons_WithRangesAndSingles)
{
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
    EXPECT_EQ(count, 3 + 1 + 3); // 0,1,2 + 5 + 10,11,12
}

TEST(YamlSnnTest, CountNeurons_InvalidFormatThrows)
{
    const std::string yaml = R"(
invalid: stuff
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();
    EXPECT_THROW(sanafe::description_count_neurons(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseNeuronSimAttributesListOfMapsFlow)
{
    const std::string yaml = R"(
- log_spikes: True
- log_potential: True
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();
    auto config = sanafe::yaml_parse_neuron_attributes(parser, node);
    EXPECT_EQ(config.log_spikes, true);
    EXPECT_EQ(config.log_potential, true);
    // Other attributes should be unset
    EXPECT_EQ(config.default_synapse_hw_name, std::nullopt);
}

TEST(YamlSnnTest, ParseNeuronSimAttributesMapFlow)
{
    const std::string yaml = R"(
log_spikes: True
log_potential: False
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();
    auto config = sanafe::yaml_parse_neuron_attributes(parser, node);
    EXPECT_EQ(config.log_spikes, true);
    EXPECT_EQ(config.log_potential, false);
    // Other attributes should be unset
    EXPECT_EQ(config.default_synapse_hw_name, std::nullopt);
}

TEST(YamlSnnTest, ParseNeuronSimAttributesListOfMapsInline)
{
    const std::string yaml = R"(
[log_spikes: True, log_potential: True]
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();
    auto config = sanafe::yaml_parse_neuron_attributes(parser, node);
    EXPECT_EQ(config.log_spikes, true);
    EXPECT_EQ(config.log_potential, true);
    // Other attributes should be unset
    EXPECT_EQ(config.default_synapse_hw_name, std::nullopt);
}

TEST(YamlSnnTest, ParseNeuronSimAttributesMapInline)
{
    const std::string yaml = R"(
{log_spikes: True, log_potential: False}
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();
    auto config = sanafe::yaml_parse_neuron_attributes(parser, node);
    EXPECT_EQ(config.log_spikes, true);
    EXPECT_EQ(config.log_potential, false);
    // Other attributes should be unset
    EXPECT_EQ(config.default_synapse_hw_name, std::nullopt);
}

TEST(YamlSnnTest, ParseFullNetworkSection)
{
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

    sanafe::SpikingNetwork net =
            sanafe::yaml_parse_network_section(parser, node);
    ASSERT_EQ(net.groups.size(), 2);
    ASSERT_TRUE(net.groups.find("Input") != net.groups.end());
    ASSERT_TRUE(net.groups.find("Output") != net.groups.end());

    EXPECT_EQ(net.groups.at("Input").neurons.size(), 2);
    EXPECT_EQ(net.groups.at("Output").neurons.size(), 2);

    const auto &input0 = net.groups.at("Input").neurons[0];
    const auto &input1 = net.groups.at("Input").neurons[1];

    ASSERT_EQ(input0.edges_out.size(), 1);
    ASSERT_EQ(input1.edges_out.size(), 1);

    EXPECT_EQ(input0.edges_out[0].post_neuron.group_name, "Output");
    EXPECT_EQ(input1.edges_out[0].post_neuron.group_name, "Output");

    EXPECT_EQ(input0.edges_out[0].post_neuron.neuron_offset, 0);
    EXPECT_EQ(input1.edges_out[0].post_neuron.neuron_offset, 1);
}

TEST(YamlSnnTest, ParseNetworkSection_InvalidFormatThrows)
{
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
    - oops: [weight: -3.0]
)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            std::runtime_error);
}

TEST(YamlSnnTest, ParseMultipleNetworks)
{
    const std::string yaml = R"(
  name: example[0..2]
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

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, WriteEdgeFormat)
{
    sanafe::NeuronAddress src{"A", 1};
    sanafe::NeuronAddress tgt{"B", 2};
    sanafe::Connection conn(0);
    conn.pre_neuron = src;
    conn.post_neuron = tgt;

    EXPECT_EQ(sanafe::write_edge_format(conn), "A.1 -> B.2");
}

TEST(YamlSnnTest, SerializeNetworkToYaml)
{
    std::filesystem::path path(SANAFE_ROOT_PATH);
    // FAIL() << "Current path: " << path.string() + "/arch/example.yaml";
    sanafe::Architecture arch =
            sanafe::load_arch(path.string() + "/arch/example.yaml");
    sanafe::SpikingNetwork net =
            sanafe::load_net(path.string() + "/snn/example.yaml", arch);
    // FAIL() << "loaded network\n";
    std::filesystem::path output_path = path / "tests/output.yaml";
    // FAIL() << "opening output path: " << output_path.string();
    net.save(output_path);

    // Now read back the file and check if it matches the original
    sanafe::SpikingNetwork loaded_net = sanafe::load_net(output_path, arch);
    ASSERT_EQ(loaded_net.groups.size(), 2);
    ASSERT_TRUE(loaded_net.groups.find("in") != loaded_net.groups.end());
    ASSERT_TRUE(loaded_net.groups.find("out") != loaded_net.groups.end());

    EXPECT_EQ(loaded_net.groups.at("in").neurons.size(), 3);
    EXPECT_EQ(loaded_net.groups.at("out").neurons.size(), 3);

    const auto &input0 = loaded_net.groups.at("in").neurons[0];
    const auto &input1 = loaded_net.groups.at("in").neurons[1];
    const auto &input2 = loaded_net.groups.at("in").neurons[2];

    ASSERT_EQ(input0.edges_out.size(), 1);
    ASSERT_EQ(input1.edges_out.size(), 1);
    ASSERT_EQ(input2.edges_out.size(), 1);

    EXPECT_EQ(input0.edges_out[0].post_neuron.group_name, "out");
    EXPECT_EQ(input1.edges_out[0].post_neuron.group_name, "out");
    EXPECT_EQ(input2.edges_out[0].post_neuron.group_name, "out");

    EXPECT_EQ(input0.edges_out[0].post_neuron.neuron_offset, 0);
    EXPECT_EQ(input1.edges_out[0].post_neuron.neuron_offset, 2);
    EXPECT_EQ(input2.edges_out[0].post_neuron.neuron_offset, 2);

    EXPECT_DOUBLE_EQ(input0.edges_out[0].synapse_attributes.at("weight"), -1.0);
    EXPECT_DOUBLE_EQ(input1.edges_out[0].synapse_attributes.at("weight"), -2.0);
    EXPECT_DOUBLE_EQ(input2.edges_out[0].synapse_attributes.at("weight"), 3.0);

    std::filesystem::remove(output_path); // Clean up the output file
}
