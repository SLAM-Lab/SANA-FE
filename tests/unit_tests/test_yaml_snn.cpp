#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "arg_parsing.hpp"
#include "arch.hpp"
#include "yaml_common.hpp"
#include "yaml_arch.hpp"
#include "yaml_snn.hpp"

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

TEST(YamlSnnTest, CountNeurons_Empty)
{

}