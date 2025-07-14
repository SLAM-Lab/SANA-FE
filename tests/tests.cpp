#include <gtest/gtest.h>
#include "models.hpp"
#include "attribute.hpp"

using namespace sanafe;

TEST(InitialTest, CheckTestFunctionality) {
    // making sure that the test framework is working correctly
    EXPECT_EQ(0, 0);
}

TEST(InitialTest, CheckChipInitialization) {
    char *a = new char[100];
    (void)a;
}

class TestCurrentBasedSynapseModel : public ::testing::Test {
protected:
    CurrentBasedSynapseModel model;

    ModelAttribute make_attr_double(double val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }
};

TEST_F(TestCurrentBasedSynapseModel, ReadReturnsCorrectWeight) {
    model.set_attribute_edge(0, "weight", make_attr_double(1.23));
    PipelineResult result = model.update(0, true);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 1.23);

}

TEST_F(TestCurrentBasedSynapseModel, WriteReturnsZero) {
    model.set_attribute_edge(0, "w", make_attr_double(2.5));
    PipelineResult result = model.update(0, false);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 0.0);
}

TEST_F(TestCurrentBasedSynapseModel, ResizesCorrectlyOnLargeIndex) {
    model.set_attribute_edge(100, "weight", make_attr_double(3.14));
    PipelineResult result = model.update(100, true);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 3.14);
}

TEST_F(TestCurrentBasedSynapseModel, MultipleWeightsMaintainValues) {
    model.set_attribute_edge(0, "w", make_attr_double(1.0));
    model.set_attribute_edge(1, "w", make_attr_double(2.0));
    model.set_attribute_edge(2, "w", make_attr_double(3.0));

    auto result0 = model.update(0, true);
    ASSERT_TRUE(result0.current.has_value());
    EXPECT_DOUBLE_EQ(result0.current.value(), 1.0);

    auto result1 = model.update(1, true);
    ASSERT_TRUE(result1.current.has_value());
    EXPECT_DOUBLE_EQ(result1.current.value(), 2.0);

    auto result2 = model.update(2, true);
    ASSERT_TRUE(result2.current.has_value());
}

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

TEST_F(TestLoihiSynapseModel, WeightIsReturnedWhenReadIsTrue) {
    set_weight(0, 2.5);
    auto result = model.update(0, true);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 2.5);
}

TEST_F(TestLoihiSynapseModel, CurrentIsZeroWhenReadIsFalse) {
    set_weight(0, 3.3);
    auto result = model.update(0, false);  // simulate write mode
    EXPECT_FALSE(result.current.has_value());  // not reading current
}

TEST_F(TestLoihiSynapseModel, LatencyAttributeIsReturned) {
    model.set_attribute_edge(0, "latency", make_attr_double(4.0));
    set_weight(0, 1.0);

    auto result = model.update(0, true);
    ASSERT_TRUE(result.latency.has_value());
    EXPECT_DOUBLE_EQ(result.latency.value(), 4.0);
}

TEST_F(TestLoihiSynapseModel, InvalidIndexDoesNotCrash) {
    // No weights set — should handle gracefully
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
    // Behavior can't be verified here, but flag can be toggled
}

TEST_F(TestLoihiSynapseModel, ResetIsSafeAndSilent) {
    set_weight(0, 2.0);
    model.reset();  // should not crash or throw
    auto result = model.update(0, true);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 2.0);  // still there after reset
}

class TestLoihiLifModel : public sanafe::LoihiLifModel {
public:
    void set_simulation_time(long int t) {
        simulation_time = t;
    }
};

sanafe::ModelAttribute make_attr_double(double val) {
    sanafe::ModelAttribute attr;
    attr.value = val;
    return attr;
}

sanafe::ModelAttribute make_attr_string(const std::string &val) {
    sanafe::ModelAttribute attr;
    attr.value = val;
    return attr;
}

sanafe::ModelAttribute make_attr_bool(bool val) {
    sanafe::ModelAttribute attr;
    attr.value = val;
    return attr;
}

TEST(LoihiLifModelTest, FiresWhenAboveThreshold) {
    TestLoihiLifModel neuron;

    neuron.set_attribute_neuron(0, "threshold", make_attr_double(64.0));
    neuron.set_attribute_neuron(0, "reset", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "reset_mode", make_attr_string("hard"));
    neuron.set_attribute_neuron(0, "leak_decay", make_attr_double(1.0));
    neuron.set_attribute_neuron(0, "input_decay", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "bias", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "force_update", make_attr_bool(false));

    neuron.reset();
    neuron.set_simulation_time(1);

    auto result = neuron.update(0, 80.0);
    EXPECT_EQ(result.status, NeuronStatus::fired);
    EXPECT_NEAR(neuron.get_potential(0), 0.0, 1e-6);
}

TEST(LoihiLifModelTest, DoesNotFireBelowThreshold) {
    TestLoihiLifModel neuron;

    neuron.set_attribute_neuron(0, "threshold", make_attr_double(64.0));
    neuron.set_attribute_neuron(0, "reset", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "reset_mode", make_attr_string("hard"));
    neuron.set_attribute_neuron(0, "leak_decay", make_attr_double(1.0));
    neuron.set_attribute_neuron(0, "input_decay", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "bias", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "force_update", make_attr_bool(false));

    neuron.reset();
    neuron.set_simulation_time(1);

    auto result = neuron.update(0, 50.0);
    EXPECT_EQ(result.status, NeuronStatus::updated);
    EXPECT_NEAR(neuron.get_potential(0), 50.0, 1e-6);
}

TEST(LoihiLifModelTest, StableWithoutInput) {
    TestLoihiLifModel neuron;

    neuron.set_attribute_neuron(0, "threshold", make_attr_double(64.0));
    neuron.set_attribute_neuron(0, "reset", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "reset_mode", make_attr_string("hard"));
    neuron.set_attribute_neuron(0, "leak_decay", make_attr_double(1.0));
    neuron.set_attribute_neuron(0, "input_decay", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "bias", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "force_update", make_attr_bool(false));

    neuron.reset();
    neuron.set_simulation_time(1);
    neuron.update(0, 50.0);

    neuron.set_simulation_time(2);
    auto result = neuron.update(0, std::nullopt);

    EXPECT_EQ(result.status, NeuronStatus::updated);
    EXPECT_NEAR(neuron.get_potential(0), 50.0, 1e-6);
}