#include <gtest/gtest.h>

#include "attribute.hpp"
#include "models.hpp"

namespace
{
class TestInputModel : public ::testing::Test
{
protected:
    sanafe::InputModel model;

    sanafe::ModelAttribute make_attr_bool(bool val)
    {
        sanafe::ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    sanafe::ModelAttribute make_attr_spike_vector(std::vector<bool> vals)
    {
        sanafe::ModelAttribute attr;
        std::vector<sanafe::ModelAttribute> elements;
        for (bool v : vals)
        {
            sanafe::ModelAttribute e;
            e.value = v;
            elements.push_back(e);
        }
        attr.value = elements;
        return attr;
    }
};
}

TEST_F(TestInputModel, GeneratesSpikeWhenSpikeValueSet)
{
    model.set_attribute_neuron(0, "spikes", make_attr_spike_vector({true}));

    auto result = model.update(0UL, std::nullopt, 1L);
    EXPECT_EQ(result.status, sanafe::fired); // Check status instead of current
}

TEST_F(TestInputModel, NoSpikeWhenSpikeValueZero)
{
    model.set_attribute_neuron(0, "spikes", make_attr_spike_vector({false}));

    auto result = model.update(0UL, std::nullopt, 1L);
    EXPECT_EQ(result.status, sanafe::idle);
}

TEST_F(TestInputModel, ResetClearsState)
{
    model.set_attribute_neuron(0, "spikes", make_attr_spike_vector({true}));
    model.update(0UL, std::nullopt, 1L); // consume the spike

    model.reset();
    auto result = model.update(0UL, std::nullopt, 1L);
    EXPECT_EQ(result.status, sanafe::idle);
}

TEST_F(TestInputModel, ExternalCurrentThrows)
{
    EXPECT_THROW(model.update(0UL, 3.5, 1L), std::runtime_error);
}

TEST_F(TestInputModel, SetsPoissonProbability)
{
    sanafe::ModelAttribute attr;
    attr.value = 0.8; // high probability
    model.set_attribute_neuron(0UL, "poisson", attr);
}

TEST_F(TestInputModel, SetsRate)
{
    sanafe::ModelAttribute attr;
    attr.value = 1.0; // 1 Hz
    model.set_attribute_neuron(0UL, "rate", attr);
}

TEST_F(TestInputModel, GeneratesSpikeWithPoisson)
{
    sanafe::ModelAttribute attr;
    attr.value = 1.0; // guaranteed spike (since uniform_distribution ∈ [0,1))
    model.set_attribute_neuron(0UL, "poisson", attr);
    auto result = model.update(0UL, std::nullopt, 1L);
    EXPECT_EQ(result.status, sanafe::fired);
}

TEST_F(TestInputModel, GeneratesSpikeWithRate)
{
    sanafe::ModelAttribute attr;
    attr.value = 1.0; // 1 spike per timestep
    model.set_attribute_neuron(0UL, "rate", attr);
    auto result = model.update(0UL, std::nullopt, 1L);
    EXPECT_EQ(result.status, sanafe::fired);
}

TEST(ModelParseResetMode, ReturnsCorrectModes)
{
    EXPECT_EQ(sanafe::model_parse_reset_mode("none"), sanafe::neuron_no_reset);
    EXPECT_EQ(
            sanafe::model_parse_reset_mode("soft"), sanafe::neuron_reset_soft);
    EXPECT_EQ(
            sanafe::model_parse_reset_mode("hard"), sanafe::neuron_reset_hard);
    EXPECT_EQ(sanafe::model_parse_reset_mode("saturate"),
            sanafe::neuron_reset_saturate);
    EXPECT_THROW(
            sanafe::model_parse_reset_mode("invalid"), std::invalid_argument);
}

TEST(ModelGetPipelineUnit, ReturnsCorrectModels)
{
    EXPECT_TRUE(std::dynamic_pointer_cast<sanafe::CurrentBasedSynapseModel>(
            sanafe::model_get_pipeline_unit("current_based")));
    EXPECT_TRUE(std::dynamic_pointer_cast<sanafe::AccumulatorModel>(
            sanafe::model_get_pipeline_unit("accumulator")));
    EXPECT_TRUE(std::dynamic_pointer_cast<sanafe::MultiTapModel1D>(
            sanafe::model_get_pipeline_unit("taps")));
    EXPECT_TRUE(std::dynamic_pointer_cast<sanafe::InputModel>(
            sanafe::model_get_pipeline_unit("input")));
    EXPECT_TRUE(std::dynamic_pointer_cast<sanafe::LoihiLifModel>(
            sanafe::model_get_pipeline_unit("leaky_integrate_fire")));
    EXPECT_TRUE(std::dynamic_pointer_cast<sanafe::TrueNorthModel>(
            sanafe::model_get_pipeline_unit("truenorth")));
    EXPECT_THROW(sanafe::model_get_pipeline_unit("invalid_model"),
            std::invalid_argument);
}
