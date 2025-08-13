#include <gtest/gtest.h>
#define private public
#define protected public
#include "models.hpp"
#undef private
#undef protected
#include "attribute.hpp"
#include "network.hpp"
#include "mapped.hpp"
#include "core.hpp"
#include "arch.hpp"

using namespace sanafe;

namespace {
class TestLoihiSynapseModel : public ::testing::Test {
protected:
    LoihiSynapseModel model;

    ModelAttribute make_attr_double(double val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    ModelAttribute make_attr_int(int val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    ModelAttribute make_attr_bool(bool val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    void set_weight(size_t idx, double val) {
        model.set_attribute_edge(idx, "weight", make_attr_double(val));
    }
};
}

TEST_F(TestLoihiSynapseModel, WeightIsReturnedWhenReadIsTrue) {
    set_weight(0, 2.5);
    auto result = model.update(0, true);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 2.5);
}

TEST_F(TestLoihiSynapseModel, CurrentIsZeroWhenReadIsFalse) {
    set_weight(0, 3.3);
    auto result = model.update(0, false);
    EXPECT_FALSE(result.current.has_value());
}

TEST_F(TestLoihiSynapseModel, LatencyAttributeIsReturned) {
    model.set_attribute_edge(0, "latency", make_attr_double(4.0));
    set_weight(0, 1.0);

    auto result = model.update(0, true);
    ASSERT_TRUE(result.latency.has_value());
    EXPECT_DOUBLE_EQ(result.latency.value(), 4.0);
}

TEST_F(TestLoihiSynapseModel, InvalidIndexDoesNotCrash) {
    EXPECT_NO_THROW({
        auto result = model.update(999, false);
        EXPECT_FALSE(result.current.has_value());
    });
}

TEST_F(TestLoihiSynapseModel, MixedSignModeCanBeSet) {
    EXPECT_NO_THROW({
        model.set_attribute_edge(0, "mixed", make_attr_bool(true));
        model.set_attribute_edge(1, "mixed", make_attr_bool(false));
    });
}

TEST_F(TestLoihiSynapseModel, ResetIsSafeAndSilent) {
    set_weight(0, 2.0);
    model.reset();
    auto result = model.update(0, true);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 2.0);
}

TEST_F(TestLoihiSynapseModel, GroupAttributeIsSetCorrectly) {
    EXPECT_NO_THROW({
        model.set_attribute_edge(0, "g", make_attr_int(5));
    });
}

TEST_F(TestLoihiSynapseModel, ReadSynapseWithLatencyFirstAccessThrows) {
    set_weight(0, 1.0);
    model.set_attribute_hw("latency_concurrent_access", make_attr_double(2.0));

    EXPECT_THROW(model.update(0, true), std::out_of_range);
}

TEST_F(TestLoihiSynapseModel, ReadSynapseSubsequentAccessThrows) {
    set_weight(0, 1.0);
    model.set_attribute_hw("latency_concurrent_access", make_attr_double(1.0));

    EXPECT_THROW(model.update(0, true), std::out_of_range); // first access
    EXPECT_THROW(model.update(0, true), std::out_of_range); // second access
}

TEST_F(TestLoihiSynapseModel, ReadSynapseClearsConcurrentAccessesThrows) {
    model.set_attribute_hw("latency_concurrent_access", make_attr_double(1.0));

    set_weight(0, 1.0);
    EXPECT_THROW(model.update(0, true), std::out_of_range);

    set_weight(1, 2.0);
    EXPECT_THROW(model.update(1, true), std::out_of_range);
}

TEST_F(TestLoihiSynapseModel, ConcurrentLatency_FirstThenSubsequent) {
    LoihiSynapseModel m;
    m.set_attribute_hw("latency_concurrent_access", make_attr_double(2.0));
    m.set_attribute_edge(0, "weight", make_attr_double(1.0));
    m.synapse_to_pre_neuron[0] = 42;   // sender A

    auto r1 = m.update(0, /*read=*/true);  // first -> 3 * 2.0 = 6
    ASSERT_TRUE(r1.latency.has_value());
    EXPECT_DOUBLE_EQ(r1.latency.value(), 6.0);

    auto r2 = m.update(0, /*read=*/true);
    ASSERT_TRUE(r2.latency.has_value());
    EXPECT_DOUBLE_EQ(r2.latency.value(), 0.0);
}

TEST_F(TestLoihiSynapseModel, ConcurrentLatency_ClearsOnMaxParallel) {
    LoihiSynapseModel m;
    m.set_attribute_hw("latency_concurrent_access", make_attr_double(1.0));
    m.set_attribute_edge(0, "weight", make_attr_double(1.0));
    m.synapse_to_pre_neuron[0] = 7;    // same sender for all calls

    auto a1 = m.update(0, true);  // first -> 3
    (void)m.update(0, true); 
    (void)m.update(0, true);    

    // size() is now 3; next call triggers clear -> treated as first again
    auto a4 = m.update(0, true);
    ASSERT_TRUE(a1.latency.has_value());
    ASSERT_TRUE(a4.latency.has_value());
    EXPECT_DOUBLE_EQ(a1.latency.value(), 3.0);
    EXPECT_DOUBLE_EQ(a4.latency.value(), 3.0);
}

TEST_F(TestLoihiSynapseModel, ConcurrentLatency_ClearsOnDifferentSender) {
    LoihiSynapseModel m;
    m.set_attribute_hw("latency_concurrent_access", make_attr_double(2.0));
    m.set_attribute_edge(0, "weight", make_attr_double(1.0));
    m.set_attribute_edge(1, "weight", make_attr_double(2.0));

    m.synapse_to_pre_neuron[0] = 100;  // sender A
    m.synapse_to_pre_neuron[1] = 200;  // sender B

    auto a1 = m.update(0, true);
    auto b1 = m.update(1, true);  

    ASSERT_TRUE(a1.latency.has_value());
    ASSERT_TRUE(b1.latency.has_value());
    EXPECT_DOUBLE_EQ(a1.latency.value(), 6.0);
    EXPECT_DOUBLE_EQ(b1.latency.value(), 6.0);
}

TEST_F(TestLoihiSynapseModel, UsesPerSynapseCostWhenNoConcurrentLatency) {
    LoihiSynapseModel m;
    m.set_attribute_edge(5, "weight", make_attr_double(1.0));
    m.set_attribute_edge(5, "latency", make_attr_double(4.5)); // costs[5] = 4.5

    auto r = m.update(5, /*read=*/true);
    ASSERT_TRUE(r.latency.has_value());
    EXPECT_DOUBLE_EQ(r.latency.value(), 4.5);
}

TEST_F(TestLoihiSynapseModel, CoversMapConnection) {
    using namespace sanafe;

    CoreAddress core_addr{0,0,0};
    CorePipelineConfiguration pipeline_cfg{};
    CoreConfiguration core_cfg("dummy_core", core_addr, pipeline_cfg);
    Core core(core_cfg);

    AxonOutPowerMetrics metrics{};
    AxonOutConfiguration axon_out_cfg(metrics, "axon");
    AxonOutUnit axon_out(axon_out_cfg);

    SpikingNetwork net("net");
    NeuronConfiguration neuron_cfg{};
    Neuron pre_neuron_obj(0, net, "g", neuron_cfg);
    MappedNeuron mpre(pre_neuron_obj, 123, &core, nullptr, 0, &axon_out, nullptr);
    MappedNeuron mpost(pre_neuron_obj, 456, &core, nullptr, 1, &axon_out, nullptr);

    // Put them in a container so addresses are stable
    std::vector<MappedNeuron> neurons;
    neurons.push_back(mpre);
    neurons.push_back(mpost);

    // Create connection and map
    MappedConnection con(neurons[0], neurons[1]);
    con.synapse_address = 0;

    LoihiSynapseModel model;
    model.set_attribute_edge(0, "weight", make_attr_double(1.0));
    model.set_attribute_hw("latency_concurrent_access", make_attr_double(1.0));

    EXPECT_NO_THROW(model.map_connection(con));

    // Triggers read_synapse
    auto r1 = model.update(0, true);
    EXPECT_TRUE(r1.latency.has_value());
}
