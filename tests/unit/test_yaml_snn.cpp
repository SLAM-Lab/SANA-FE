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
    const std::string edge_description = R"(A.1 -> B.2)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(edge_description, parser);
    auto node = tree.rootref();

    auto [src, tgt] = sanafe::description_parse_edge_description(
            edge_description, parser, node);

    EXPECT_EQ(src.group_name, "A");
    EXPECT_EQ(src.neuron_offset, 1);
    EXPECT_EQ(tgt.group_name, "B");
    EXPECT_EQ(tgt.neuron_offset, 2);
}

TEST(YamlSnnTest, ParseEdgeDescription_MissingDotThrows)
{
    const std::string placeholder_yaml = R"(p)";
    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(placeholder_yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::description_parse_edge_description(
                         "A -> B.2", parser, node),
            sanafe::YamlDescriptionParsingError);
    EXPECT_THROW(sanafe::description_parse_edge_description(
                         "A.1 -> B", parser, node),
            sanafe::YamlDescriptionParsingError);
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
            sanafe::YamlDescriptionParsingError);
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

    EXPECT_EQ(loaded_net.groups.at("in").neurons.size(), 2);
    EXPECT_EQ(loaded_net.groups.at("out").neurons.size(), 2);

    const auto &input0 = loaded_net.groups.at("in").neurons[0];
    const auto &input1 = loaded_net.groups.at("in").neurons[1];

    ASSERT_EQ(input0.edges_out.size(), 2);
    ASSERT_EQ(input1.edges_out.size(), 2);

    EXPECT_EQ(input0.edges_out[0].post_neuron.group_name, "out");
    EXPECT_EQ(input0.edges_out[1].post_neuron.group_name, "out");
    EXPECT_EQ(input1.edges_out[0].post_neuron.group_name, "out");
    EXPECT_EQ(input1.edges_out[1].post_neuron.group_name, "out");

    EXPECT_EQ(input0.edges_out[0].post_neuron.neuron_offset, 0);
    EXPECT_EQ(input0.edges_out[1].post_neuron.neuron_offset, 1);
    EXPECT_EQ(input1.edges_out[0].post_neuron.neuron_offset, 0);
    EXPECT_EQ(input1.edges_out[1].post_neuron.neuron_offset, 1);

    EXPECT_DOUBLE_EQ(input0.edges_out[0].synapse_attributes.at("weight"), -1.0);
    EXPECT_DOUBLE_EQ(input0.edges_out[1].synapse_attributes.at("weight"), 2.0);
    EXPECT_DOUBLE_EQ(input1.edges_out[0].synapse_attributes.at("weight"), 1.0);
    EXPECT_DOUBLE_EQ(input1.edges_out[1].synapse_attributes.at("weight"), 3.0);

    std::filesystem::remove(output_path); // Clean up the output file
}

TEST(YamlSnnTest, ParseEdgeDescription_NoArrowThrows)
{
    const std::string edge_description = R"(A.1 B.2)";

    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(edge_description, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::description_parse_edge_description(
                         edge_description, parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseEdgeDescription_HyperedgeNoNeuronOffset)
{
    const std::string edge_description = R"(A -> B)";

    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(edge_description, parser);
    auto node = tree.rootref();
    auto [src, tgt] = sanafe::description_parse_edge_description(
            edge_description, parser, node);

    EXPECT_EQ(src.group_name, "A");
    EXPECT_FALSE(src.neuron_offset.has_value());
    EXPECT_EQ(tgt.group_name, "B");
    EXPECT_FALSE(tgt.neuron_offset.has_value());
}

TEST(YamlSnnTest, ParseEdgeDescription_WithWhitespace)
{
    const std::string edge_description = R"(A.1  ->  B.2)";

    // NOLINTBEGIN(misc-include-cleaner)
    ryml::EventHandlerTree event_handler = {};
    // Enable location tracking for helpful error prints
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(edge_description, parser);
    auto node = tree.rootref();
    auto [src, tgt] = sanafe::description_parse_edge_description(
            edge_description, parser, node);

    EXPECT_EQ(src.group_name, "A");
    EXPECT_EQ(src.neuron_offset, 1);
    EXPECT_EQ(tgt.group_name, "B");
    EXPECT_EQ(tgt.neuron_offset, 2);
}

TEST(YamlSnnTest, CountNeurons_MapFormatThrows)
{
    const std::string yaml = R"(
0:
  1:
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::description_count_neurons(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, CountNeurons_NestedMapInList)
{
    const std::string yaml = R"(
- 0: {attr: value}
- 1: {attr: value}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    size_t count = sanafe::description_count_neurons(parser, node);
    EXPECT_EQ(count, 2);
}

TEST(YamlSnnTest, ParseNeuronAttributes_HardwareUnits)
{
    const std::string yaml = R"(
synapse_hw_name: syn_unit_1
dendrite_hw_name: dend_unit_1
soma_hw_name: soma_unit_1
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    auto config = sanafe::yaml_parse_neuron_attributes(parser, node);
    EXPECT_EQ(config.default_synapse_hw_name, "syn_unit_1");
    EXPECT_EQ(config.dendrite_hw_name, "dend_unit_1");
    EXPECT_EQ(config.soma_hw_name, "soma_unit_1");
}

TEST(YamlSnnTest, ParseNeuronAttributes_UnitSpecificModelAttributes)
{
    const std::string yaml = R"(
shared_attr: 1.0
dendrite:
  dend_specific: 2.0
soma:
  soma_specific: 3.0
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    auto config = sanafe::yaml_parse_neuron_attributes(parser, node);

    ASSERT_TRUE(config.model_attributes.find("shared_attr") !=
            config.model_attributes.end());
    EXPECT_TRUE(config.model_attributes.at("shared_attr").forward_to_dendrite);
    EXPECT_TRUE(config.model_attributes.at("shared_attr").forward_to_soma);

    ASSERT_TRUE(config.model_attributes.find("dend_specific") !=
            config.model_attributes.end());
    EXPECT_FALSE(
            config.model_attributes.at("dend_specific").forward_to_synapse);
    EXPECT_FALSE(config.model_attributes.at("dend_specific").forward_to_soma);

    ASSERT_TRUE(config.model_attributes.find("soma_specific") !=
            config.model_attributes.end());
    EXPECT_FALSE(
            config.model_attributes.at("soma_specific").forward_to_synapse);
    EXPECT_FALSE(
            config.model_attributes.at("soma_specific").forward_to_dendrite);
}

TEST(YamlSnnTest, ParseNeuronSection_InvalidNeuronId)
{
    const std::string yaml = R"(
  name: test
  groups:
    - name: Input
      neurons:
        - 0..1
        - 5: {weight: 1.0}
  edges: []
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            std::out_of_range);
}

TEST(YamlSnnTest, ParseNetworkSection_MissingGroupsThrows)
{
    const std::string yaml = R"(
  name: example
  edges: []
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseNetworkSection_MissingEdgesThrows)
{
    const std::string yaml = R"(
  name: example
  groups:
    - name: Input
      neurons:
        - 0
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseNeuronConnection_InvalidSourceGroup)
{
    const std::string yaml = R"(
  name: test
  groups:
    - name: Output
      neurons:
        - 0
  edges:
    - Invalid.0 -> Output.0: {}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseNeuronConnection_InvalidTargetGroup)
{
    const std::string yaml = R"(
  name: test
  groups:
    - name: Input
      neurons:
        - 0
  edges:
    - Input.0 -> Invalid.0: {}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseNeuronConnection_OutOfBoundsNeuronId)
{
    const std::string yaml = R"(
  name: test
  groups:
    - name: Input
      neurons:
        - 0
    - name: Output
      neurons:
        - 0
  edges:
    - Input.5 -> Output.0: {}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseHyperedge_NoTypeThrows)
{
    const std::string yaml = R"(
  name: test
  groups:
    - name: Input
      neurons:
        - 0..1
    - name: Output
      neurons:
        - 0..1
  edges:
    - Input -> Output: {weight: 1.0}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseHyperedge_InvalidTypeThrows)
{
    const std::string yaml = R"(
  name: test
  groups:
    - name: Input
      neurons:
        - 0..1
    - name: Output
      neurons:
        - 0..1
  edges:
    - Input -> Output: {type: invalid_type}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    EXPECT_THROW(sanafe::yaml_parse_network_section(parser, node),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseEdgeAttributes_UnitSpecific)
{
    const std::string yaml = R"(
  name: test
  groups:
    - name: Input
      neurons:
        - 0
    - name: Output
      neurons:
        - 0
  edges:
    - Input.0 -> Output.0:
        synapse:
          weight: 1.5
        dendrite:
          delay: 2
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);
    auto node = tree.rootref();

    sanafe::SpikingNetwork net =
            sanafe::yaml_parse_network_section(parser, node);

    const auto &conn = net.groups.at("Input").neurons[0].edges_out[0];
    EXPECT_DOUBLE_EQ(conn.synapse_attributes.at("weight"), 1.5);
    EXPECT_EQ(static_cast<int>(conn.dendrite_attributes.at("delay")), 2);
}

TEST(YamlSnnTest, ParseMappingSection_InvalidNeuronGroup)
{
    std::filesystem::path path(SANAFE_ROOT_PATH);
    sanafe::Architecture arch =
            sanafe::load_arch(path.string() + "/arch/example.yaml");

    const std::string yaml = R"(
network:
  name: test
  groups:
    - name: Input
      neurons:
        - 0
  edges: []
mappings:
  - InvalidGroup.0: {core: 0.0}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);

    auto snn = sanafe::yaml_parse_network_section(parser, tree["network"]);
    EXPECT_THROW(sanafe::description_parse_mapping_section_yaml(
                         parser, tree["mappings"], arch, snn),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseMappingSection_OutOfBoundsTile)
{
    std::filesystem::path path(SANAFE_ROOT_PATH);
    sanafe::Architecture arch =
            sanafe::load_arch(path.string() + "/arch/example.yaml");

    const std::string yaml = R"(
network:
  name: test
  groups:
    - name: Input
      neurons:
        - 0
  edges: []
mappings:
  - Input.0: {core: 999.0}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);

    auto snn = sanafe::yaml_parse_network_section(parser, tree["network"]);
    EXPECT_THROW(sanafe::description_parse_mapping_section_yaml(
                         parser, tree["mappings"], arch, snn),
            sanafe::YamlDescriptionParsingError);
}

TEST(YamlSnnTest, ParseMappingSection_NeuronRange)
{
    std::filesystem::path path(SANAFE_ROOT_PATH);
    sanafe::Architecture arch =
            sanafe::load_arch(path.string() + "/arch/example.yaml");

    const std::string yaml = R"(
network:
  name: test
  groups:
    - name: Input
      neurons:
        - 0..2
  edges: []
mappings:
  - Input.0..2: {core: 0.0}
)";
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    auto tree = parse_yaml_snippet(yaml, parser);

    sanafe::SpikingNetwork net =
            sanafe::yaml_parse_network_section(parser, tree["network"]);

    EXPECT_NO_THROW(sanafe::description_parse_mapping_section_yaml(
            parser, tree["mappings"], arch, net));

    EXPECT_TRUE(net.groups.at("Input").neurons[0].core_address.has_value());
    EXPECT_TRUE(net.groups.at("Input").neurons[1].core_address.has_value());
    EXPECT_TRUE(net.groups.at("Input").neurons[2].core_address.has_value());
}
