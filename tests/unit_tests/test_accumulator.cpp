#include <gtest/gtest.h>

#include "attribute.hpp"
#include "models.hpp"

namespace
{
class TestAccumulatorModel : public ::testing::Test
{
protected:
    class AccumulatorModelExposed : public sanafe::AccumulatorModel
    {
    public:
        void set_simulation_time(long int t)
        {
            simulation_time = t;
        }
    };

    AccumulatorModelExposed model;

    sanafe::ModelAttribute make_attr_double(double val)
    {
        sanafe::ModelAttribute attr;
        attr.value = val;
        return attr;
    }
};
}

TEST_F(TestAccumulatorModel, IntegratesCurrent)
{
    model.set_attribute_neuron(0, "dendrite_leak_decay", make_attr_double(1.0));
    model.update(0, 5.0, std::nullopt);
    auto result = model.update(0, std::nullopt, std::nullopt);

    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 5.0);
}

TEST_F(TestAccumulatorModel, AppliesLeakDecay)
{
    model.set_attribute_neuron(0, "dendrite_leak_decay", make_attr_double(0.5));
    model.update(0, 4.0, std::nullopt);

    model.set_simulation_time(2); // advances two timesteps
    auto result = model.update(0, std::nullopt, std::nullopt);

    double expected = 4.0 * std::pow(0.5, 2);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_NEAR(result.current.value(), expected, 1e-9);
}

TEST_F(TestAccumulatorModel, AccumulatesChargeOverTime)
{
    model.set_attribute_neuron(0, "dendrite_leak_decay", make_attr_double(1.0));
    model.update(0, 2.0, std::nullopt);
    model.update(0, 3.0, std::nullopt);

    auto result = model.update(0, std::nullopt, std::nullopt);

    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 5.0);
}

TEST_F(TestAccumulatorModel, UnknownAttributeDoesNotThrow)
{
    sanafe::ModelAttribute attr;
    attr.value = 42.0;
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "unknown_attribute", attr));
}

TEST_F(TestAccumulatorModel, LeakDecayChangeAffectsNextUpdate)
{
    model.set_attribute_neuron(0, "dendrite_leak_decay", make_attr_double(1.0));
    model.update(0, 8.0, std::nullopt);

    model.set_attribute_neuron(
            0, "dendrite_leak_decay", make_attr_double(0.25));
    model.set_simulation_time(2);

    auto result = model.update(0, std::nullopt, std::nullopt);
    double expected = 8.0 * std::pow(0.25, 2);

    ASSERT_TRUE(result.current.has_value());
    EXPECT_NEAR(result.current.value(), expected, 1e-9);
}
