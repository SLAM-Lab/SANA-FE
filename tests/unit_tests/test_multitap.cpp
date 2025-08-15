#include <gtest/gtest.h>
#include "models.hpp"
#include "attribute.hpp"

using namespace sanafe;

namespace {
    class TestMultiTapModel1D : public ::testing::Test {
    protected:
        class MultiTapModelExposed : public MultiTapModel1D {
        public:
            void set_simulation_time(long int t) { simulation_time = t; }
        };

        MultiTapModelExposed model;

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

        ModelAttribute make_attr_vec(const std::vector<double> &vals) {
            ModelAttribute attr;
            std::vector<ModelAttribute> inner;
            inner.reserve(vals.size());
            for (double v : vals) {
                ModelAttribute a;
                a.value = v;
                inner.push_back(a);
            }
            attr.value = inner;
            return attr;
        }
    };
}

TEST_F(TestMultiTapModel1D, TapsZeroThrows) {
    ModelAttribute attr;
    attr.value = static_cast<int>(0);
    EXPECT_THROW(model.set_attribute_neuron(0, "taps", attr), std::invalid_argument);
}

TEST_F(TestMultiTapModel1D, TapsResizeValid) {
    ModelAttribute attr;
    attr.value = static_cast<int>(3);
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "taps", attr));
}

TEST_F(TestMultiTapModel1D, TimeConstantsResizing) {
    ModelAttribute attr;
    attr.value = static_cast<int>(2);
    model.set_attribute_neuron(0, "taps", attr);

    EXPECT_NO_THROW(model.set_attribute_neuron(0, "time_constants", make_attr_vec({0.9, 0.8})));
}

TEST_F(TestMultiTapModel1D, TimeConstantsResizeLargerVector) {
    model.set_attribute_neuron(0, "taps", make_attr_int(2));
    auto attr = make_attr_vec({0.5, 0.5, 0.5});
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "time_constants", attr));
}

TEST_F(TestMultiTapModel1D, SpaceConstantsResizing) {
    ModelAttribute attr;
    attr.value = static_cast<int>(3);
    model.set_attribute_neuron(0, "taps", attr);

    EXPECT_NO_THROW(model.set_attribute_neuron(0, "space_constants", make_attr_vec({0.5, 0.5})));
}

TEST_F(TestMultiTapModel1D, SpaceConstantsResizeLargerVector) {
    model.set_attribute_neuron(0, "taps", make_attr_int(2));
    auto attr = make_attr_vec({0.4, 0.4, 0.4});
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "space_constants", attr));
}

TEST_F(TestMultiTapModel1D, UnknownAttributeDoesNotThrow) {
    ModelAttribute attr;
    attr.value = 1.0;
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "unknown_attribute", attr));
}

TEST_F(TestMultiTapModel1D, InputCurrentAdds) {
    ModelAttribute taps_attr;
    taps_attr.value = static_cast<int>(2);
    model.set_attribute_neuron(0, "taps", taps_attr);

    model.set_attribute_neuron(0, "time_constants", make_attr_vec({1.0, 1.0}));
    model.set_attribute_neuron(0, "space_constants", make_attr_vec({0.0}));

    auto result_before = model.update(0, 1.5, std::nullopt);
    ASSERT_TRUE(result_before.current.has_value());
    EXPECT_DOUBLE_EQ(result_before.current.value(), 1.5);
}

TEST_F(TestMultiTapModel1D, InputCurrentToMappedTap) {
    ModelAttribute taps_attr;
    taps_attr.value = static_cast<int>(2);
    model.set_attribute_neuron(0, "taps", taps_attr);

    ModelAttribute tap_attr;
    tap_attr.value = static_cast<int>(1);
    model.set_attribute_edge(0, "tap", tap_attr);

    model.update(0, std::nullopt, std::nullopt);

    EXPECT_NO_THROW(model.update(0, 2.0, 0));
}

TEST_F(TestMultiTapModel1D, InvalidTapThrows) {
    ModelAttribute taps_attr;
    taps_attr.value = static_cast<int>(1);
    model.set_attribute_neuron(0, "taps", taps_attr);

    ModelAttribute tap_attr;
    tap_attr.value = static_cast<int>(5);
    model.set_attribute_edge(0, "tap", tap_attr);

    EXPECT_THROW(model.update(0, 1.0, 0), std::logic_error);
}

TEST_F(TestMultiTapModel1D, ResetClearsVoltages) {
    ModelAttribute taps_attr;
    taps_attr.value = static_cast<int>(1);
    model.set_attribute_neuron(0, "taps", taps_attr);

    model.update(0, 3.0, std::nullopt);
    model.reset();

    auto result = model.update(0, std::nullopt, std::nullopt);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_DOUBLE_EQ(result.current.value(), 0.0);
}

TEST_F(TestMultiTapModel1D, CalculateNextStateChangesVoltages) {
    ModelAttribute taps_attr;
    taps_attr.value = static_cast<int>(2);
    model.set_attribute_neuron(0, "taps", taps_attr);

    model.set_attribute_neuron(0, "time_constants", make_attr_vec({0.5, 0.5}));
    model.set_attribute_neuron(0, "space_constants", make_attr_vec({0.0}));

    model.update(0, 2.0, std::nullopt);
    model.set_simulation_time(1);

    auto result = model.update(0, std::nullopt, std::nullopt);
    ASSERT_TRUE(result.current.has_value());
    EXPECT_LT(result.current.value(), 2.0);
}

TEST_F(TestMultiTapModel1D, ReduceNumberOfTapsTriggersWarningPath) {
    ModelAttribute attr_large;
    attr_large.value = static_cast<int>(4);
    model.set_attribute_neuron(0, "taps", attr_large);

    ModelAttribute attr_small;
    attr_small.value = static_cast<int>(2);
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "taps", attr_small));
}

TEST_F(TestMultiTapModel1D, TimeConstantsTooFewThrows) {
    model.set_attribute_neuron(0, "taps", make_attr_int(3));
    EXPECT_THROW(
        model.set_attribute_neuron(0, "time_constants", make_attr_vec({0.9, 0.8})),
        std::invalid_argument
    );
}

TEST_F(TestMultiTapModel1D, SpaceConstantsTooFewThrows) {
    model.set_attribute_neuron(0, "taps", make_attr_int(3));
    EXPECT_THROW(
        model.set_attribute_neuron(0, "space_constants", make_attr_vec({0.5})),
        std::invalid_argument
    );
}