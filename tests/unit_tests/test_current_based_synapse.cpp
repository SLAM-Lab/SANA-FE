#include <gtest/gtest.h>
#include "models.hpp"
#include "attribute.hpp"

using namespace sanafe;

namespace {
    class TestCurrentBasedSynapseModel : public ::testing::Test {
    protected:
        CurrentBasedSynapseModel model;
    
        ModelAttribute make_attr_double(double val) {
            ModelAttribute attr;
            attr.value = val;
            return attr;
        }
    };
}

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

TEST_F(TestCurrentBasedSynapseModel, TestReset) {
    model.set_attribute_edge(0, "weight", make_attr_double(1.23));
    model.set_attribute_edge(1, "w", make_attr_double(4.56));

    EXPECT_NO_THROW(model.reset());

    auto result = model.update(0, true);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 1.23);
}