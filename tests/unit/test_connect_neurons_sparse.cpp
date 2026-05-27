// test_network_sparse.cpp
// Regression tests for sanafe::NeuronGroup::connect_neurons_sparse
//
// Bug being regressed:
//   connect_neurons_sparse indexed the per-edge attribute list by source_id
//   (a neuron offset) instead of by edge index. This caused edges to receive
//   the wrong attribute values whenever the source_id of an edge didn't match
//   its position in source_dest_id_pairs, and out-of-bounds reads when
//   source_id >= source_dest_id_pairs.size().

#include "attribute.hpp"
#include "network.hpp"
#include <gtest/gtest.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

sanafe::ModelAttribute make_attr_double(double val)
{
    sanafe::ModelAttribute attr;
    attr.value = val;
    // Bug-1 fix may honor these flags; the original sparse path ignored them.
    // Set both true so the test passes under either behavior.
    attr.forward_to_synapse = true;
    attr.forward_to_dendrite = true;
    return attr;
}

// Helper: pull a double weight out of a connection's synapse_attributes.
double weight_of(const sanafe::Connection &con)
{
    const auto it = con.synapse_attributes.find("weight");
    EXPECT_NE(it, con.synapse_attributes.end())
            << "connection missing 'weight' attribute";
    return std::get<double>(it->second.value);
}

} // namespace

// ---------------------------------------------------------------------------
// Regression test 1: test per-edge attributes
// ---------------------------------------------------------------------------
TEST(ConnectNeuronsSparseTest, AttributesIndexedByEdgePositionNotSourceId)
{
    sanafe::SpikingNetwork net;
    sanafe::NeuronConfiguration cfg{};

    auto &src = net.create_neuron_group("src", 3, cfg);
    auto &dst = net.create_neuron_group("dst", 3, cfg);

    // Edges defined in an order where source_id never equals edge index:
    //   edge 0: src[2] -> dst[0]   weight should be 10.0
    //   edge 1: src[0] -> dst[1]   weight should be 20.0
    //   edge 2: src[1] -> dst[2]   weight should be 30.0
    const std::vector<std::pair<size_t, size_t>> pairs = {
            {2, 0}, {0, 1}, {1, 2}};

    std::map<std::string, std::vector<sanafe::ModelAttribute>> attrs;
    attrs["weight"] = {make_attr_double(10.0), make_attr_double(20.0),
            make_attr_double(30.0)};

    src.connect_neurons_sparse(dst, attrs, pairs);

    // Verify each source neuron has exactly the expected outgoing edge with
    // the correct weight.
    ASSERT_EQ(src.neurons[2].edges_out.size(), 1u);
    ASSERT_EQ(src.neurons[0].edges_out.size(), 1u);
    ASSERT_EQ(src.neurons[1].edges_out.size(), 1u);

    EXPECT_DOUBLE_EQ(weight_of(src.neurons[2].edges_out[0]), 10.0)
            << "edge 0 (src[2]->dst[0]) got wrong weight; "
               "likely indexed by source_id=2 instead of edge=0";
    EXPECT_DOUBLE_EQ(weight_of(src.neurons[0].edges_out[0]), 20.0)
            << "edge 1 (src[0]->dst[1]) got wrong weight; "
               "likely indexed by source_id=0 instead of edge=1";
    EXPECT_DOUBLE_EQ(weight_of(src.neurons[1].edges_out[0]), 30.0)
            << "edge 2 (src[1]->dst[2]) got wrong weight; "
               "likely indexed by source_id=1 instead of edge=2";
}

// ---------------------------------------------------------------------------
// Regression test 2: multiple edges from the same source neuron must each
// receive their own distinct attribute, not the same one.
// ---------------------------------------------------------------------------
TEST(ConnectNeuronsSparseTest, MultipleEdgesFromSameSourceGetDistinctAttributes)
{
    sanafe::SpikingNetwork net;
    sanafe::NeuronConfiguration cfg{};

    auto &src = net.create_neuron_group("src", 2, cfg);
    auto &dst = net.create_neuron_group("dst", 3, cfg);

    // Two edges from src[0], one from src[1]; three distinct weights.
    const std::vector<std::pair<size_t, size_t>> pairs = {
            {0, 0}, {0, 1}, {1, 2}};

    std::map<std::string, std::vector<sanafe::ModelAttribute>> attrs;
    attrs["weight"] = {make_attr_double(1.0), make_attr_double(2.0),
            make_attr_double(3.0)};

    src.connect_neurons_sparse(dst, attrs, pairs);

    ASSERT_EQ(src.neurons[0].edges_out.size(), 2u);
    ASSERT_EQ(src.neurons[1].edges_out.size(), 1u);

    // edges_out is appended in iteration order of pairs, so:
    //   src[0].edges_out[0] is edge 0 -> weight 1.0
    //   src[0].edges_out[1] is edge 1 -> weight 2.0
    //   src[1].edges_out[0] is edge 2 -> weight 3.0
    EXPECT_DOUBLE_EQ(weight_of(src.neurons[0].edges_out[0]), 1.0);
    EXPECT_DOUBLE_EQ(weight_of(src.neurons[0].edges_out[1]), 2.0)
            << "second edge from src[0] got the same weight as the first; "
               "bug indexes by source_id, so both reuse value_list[0]";
    EXPECT_DOUBLE_EQ(weight_of(src.neurons[1].edges_out[0]), 3.0);
}

// ---------------------------------------------------------------------------
// Regression test 3: out-of-bounds protection.
// ---------------------------------------------------------------------------
TEST(ConnectNeuronsSparseTest, LargeSourceIdSmallEdgeCountDoesNotOverrun)
{
    sanafe::SpikingNetwork net;
    sanafe::NeuronConfiguration cfg{};

    // 10 source neurons but only 2 edges, with source_ids 5 and 7.
    // Buggy code would access value_list[5] and value_list[7] on a
    // length-2 vector.
    auto &src = net.create_neuron_group("src", 10, cfg);
    auto &dst = net.create_neuron_group("dst", 10, cfg);

    const std::vector<std::pair<size_t, size_t>> pairs = {{5, 0}, {7, 1}};

    std::map<std::string, std::vector<sanafe::ModelAttribute>> attrs;
    attrs["weight"] = {make_attr_double(100.0), make_attr_double(200.0)};

    ASSERT_NO_THROW(src.connect_neurons_sparse(dst, attrs, pairs));

    ASSERT_EQ(src.neurons[5].edges_out.size(), 1u);
    ASSERT_EQ(src.neurons[7].edges_out.size(), 1u);

    EXPECT_DOUBLE_EQ(weight_of(src.neurons[5].edges_out[0]), 100.0);
    EXPECT_DOUBLE_EQ(weight_of(src.neurons[7].edges_out[0]), 200.0);
}
