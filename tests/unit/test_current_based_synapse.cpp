#include <gtest/gtest.h>

#include "attribute.hpp"
#include "models.hpp"

namespace
{
class TestCurrentBasedSynapseModel : public ::testing::Test
{
protected:
    sanafe::CurrentBasedSynapseModel model;

    sanafe::ModelAttribute make_attr_double(double val)
    {
        sanafe::ModelAttribute attr;
        attr.value = val;
        return attr;
    }
};
}

TEST_F(TestCurrentBasedSynapseModel, ReadReturnsCorrectWeight)
{
    model.set_attribute_edge(0UL, "weight", make_attr_double(1.23));
    sanafe::PipelineResult result = model.update(0UL, true, 1L);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_NEAR(result.current.value(), 1.23, 1e-6);
}

TEST_F(TestCurrentBasedSynapseModel, WriteReturnsZero)
{
    model.set_attribute_edge(0, "w", make_attr_double(2.5));
    sanafe::PipelineResult result = model.update(0UL, false, 1L);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_NEAR(result.current.value(), 0.0, 1e-6);
}

TEST_F(TestCurrentBasedSynapseModel, ResizesCorrectlyOnLargeIndex)
{
    model.set_attribute_edge(100UL, "weight", make_attr_double(3.14));
    sanafe::PipelineResult result = model.update(100UL, true, 1L);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_NEAR(result.current.value(), 3.14, 1e-6);
}

TEST_F(TestCurrentBasedSynapseModel, MultipleWeightsMaintainValues)
{
    model.set_attribute_edge(0, "w", make_attr_double(1.0));
    model.set_attribute_edge(1, "w", make_attr_double(2.0));
    model.set_attribute_edge(2, "w", make_attr_double(3.0));

    auto result0 = model.update(0UL, true, 1L);
    ASSERT_TRUE(result0.current.has_value());
    EXPECT_NEAR(result0.current.value(), 1.0, 1e-6);

    auto result1 = model.update(1UL, true, 1L);
    ASSERT_TRUE(result1.current.has_value());
    EXPECT_NEAR(result1.current.value(), 2.0, 1e-6);

    auto result2 = model.update(2UL, true, 1L);
    ASSERT_TRUE(result2.current.has_value());
}

TEST_F(TestCurrentBasedSynapseModel, TestReset)
{
    model.set_attribute_edge(0UL, "weight", make_attr_double(1.23));
    model.set_attribute_edge(1UL, "w", make_attr_double(4.56));

    EXPECT_NO_THROW(model.reset());

    auto result = model.update(0UL, true, 1L);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_NEAR(result.current.value(), 1.23, 1e-6);
}