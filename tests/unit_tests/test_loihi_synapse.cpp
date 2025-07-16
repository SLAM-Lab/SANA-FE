#include <gtest/gtest.h>
#include "models.hpp"
#include "attribute.hpp"

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