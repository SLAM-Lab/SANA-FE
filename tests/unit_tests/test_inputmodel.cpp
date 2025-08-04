#include <gtest/gtest.h>
#include "models.hpp"
#include "attribute.hpp"

using namespace sanafe;

namespace {
class TestInputModel : public ::testing::Test {
protected:
    InputModel model;

    ModelAttribute make_attr_bool(bool val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    ModelAttribute make_attr_spike_vector(std::vector<bool> vals) {
        ModelAttribute attr;
        std::vector<ModelAttribute> elements;
        for (bool v : vals) {
            ModelAttribute e;
            e.value = v;
            elements.push_back(e);
        }
        attr.value = elements;
        return attr;
    }
};
}

TEST_F(TestInputModel, GeneratesSpikeWhenSpikeValueSet) {
    model.set_attribute_neuron(0, "spikes", make_attr_spike_vector({true}));

    auto result = model.update(0, std::nullopt);
    EXPECT_EQ(result.status, sanafe::fired);  // Check status instead of current
}

TEST_F(TestInputModel, NoSpikeWhenSpikeValueZero) {
    model.set_attribute_neuron(0, "spikes", make_attr_spike_vector({false}));

    auto result = model.update(0, std::nullopt);
    EXPECT_EQ(result.status, sanafe::idle);
}

TEST_F(TestInputModel, ResetClearsState) {
    model.set_attribute_neuron(0, "spikes", make_attr_spike_vector({true}));
    model.update(0, std::nullopt);  // consume the spike

    model.reset();
    auto result = model.update(0, std::nullopt);
    EXPECT_EQ(result.status, sanafe::idle);
}

TEST_F(TestInputModel, ExternalCurrentThrows) {
    EXPECT_THROW(model.update(0, 3.5), std::runtime_error);
}

TEST_F(TestInputModel, SetsPoissonProbability) {
    ModelAttribute attr;
    attr.value = 0.8;  // high probability
    model.set_attribute_neuron(0, "poisson", attr);
}

TEST_F(TestInputModel, SetsRate) {
    ModelAttribute attr;
    attr.value = 1.0;  // 1 Hz
    model.set_attribute_neuron(0, "rate", attr);
}

TEST_F(TestInputModel, GeneratesSpikeWithPoisson) {
    ModelAttribute attr;
    attr.value = 1.0;  // guaranteed spike (since uniform_distribution ∈ [0,1))
    model.set_attribute_neuron(0, "poisson", attr);
    auto result = model.update(0, std::nullopt);
    EXPECT_EQ(result.status, sanafe::fired);
}

TEST_F(TestInputModel, GeneratesSpikeWithRate) {
    ModelAttribute attr;
    attr.value = 1.0; // 1 spike per timestep
    model.set_attribute_neuron(0, "rate", attr);
    auto result = model.update(0, std::nullopt);
    EXPECT_EQ(result.status, sanafe::fired);
}

TEST(ModelParseResetMode, ReturnsCorrectModes) {
    using namespace sanafe;
    EXPECT_EQ(model_parse_reset_mode("none"), neuron_no_reset);
    EXPECT_EQ(model_parse_reset_mode("soft"), neuron_reset_soft);
    EXPECT_EQ(model_parse_reset_mode("hard"), neuron_reset_hard);
    EXPECT_EQ(model_parse_reset_mode("saturate"), neuron_reset_saturate);
    EXPECT_THROW(model_parse_reset_mode("invalid"), std::invalid_argument);
}

TEST(ModelGetPipelineUnit, ReturnsCorrectModels) {
    using namespace sanafe;
    EXPECT_TRUE(std::dynamic_pointer_cast<CurrentBasedSynapseModel>(model_get_pipeline_unit("current_based")));
    EXPECT_TRUE(std::dynamic_pointer_cast<LoihiSynapseModel>(model_get_pipeline_unit("loihi")));
    EXPECT_TRUE(std::dynamic_pointer_cast<AccumulatorModel>(model_get_pipeline_unit("accumulator")));
    EXPECT_TRUE(std::dynamic_pointer_cast<MultiTapModel1D>(model_get_pipeline_unit("taps")));
    EXPECT_TRUE(std::dynamic_pointer_cast<InputModel>(model_get_pipeline_unit("input")));
    EXPECT_TRUE(std::dynamic_pointer_cast<LoihiLifModel>(model_get_pipeline_unit("leaky_integrate_fire")));
    EXPECT_TRUE(std::dynamic_pointer_cast<TrueNorthModel>(model_get_pipeline_unit("truenorth")));
    EXPECT_THROW(model_get_pipeline_unit("invalid_model"), std::invalid_argument);
}