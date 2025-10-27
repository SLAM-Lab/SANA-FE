#include <gtest/gtest.h>

#include "attribute.hpp"
#include "models.hpp"

namespace
{
class TestAccumulatorModel : public ::testing::Test
{
protected:
    sanafe::AccumulatorModel model;

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
    model.update(0UL, 5.0, std::nullopt, 1L);
    auto result = model.update(0UL, std::nullopt, std::nullopt, 1L);

    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 5.0);
}

TEST_F(TestAccumulatorModel, AccumulatesChargeOverTime)
{
    model.update(0UL, 2.0, std::nullopt, 1L);
    model.update(0UL, 3.0, std::nullopt, 1L);
    auto result = model.update(0UL, std::nullopt, std::nullopt, 1L);

    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 5.0);
}

TEST_F(TestAccumulatorModel, UnknownAttributeDoesNotThrow)
{
    sanafe::ModelAttribute attr;
    attr.value = 42.0;
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "unknown_attribute", attr));
}
