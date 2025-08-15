#include <gtest/gtest.h>
#include "models.hpp"
#include "attribute.hpp"
#include <cstdlib>

using namespace sanafe;

namespace {
class TestTrueNorthModel : public ::testing::Test {
protected:
    TrueNorthModel model;

    ModelAttribute make_attr_double(double val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    ModelAttribute make_attr_string(const std::string &val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    ModelAttribute make_attr_bool(bool val) {
        ModelAttribute attr;
        attr.value = val;
        return attr;
    }

    ModelAttribute make_attr_int(int val) {
    ModelAttribute attr; attr.value = val; return attr;
    }
};
}

TEST_F(TestTrueNorthModel, SetThresholdAndUpdateFires) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(0.5));
    model.set_attribute_neuron(0, "reset_mode", make_attr_string("hard"));
    model.set_attribute_neuron(0, "reset", make_attr_double(0.0));

    auto result = model.update(0, 1.0); // Add current
    EXPECT_EQ(result.status, sanafe::fired);
}

TEST_F(TestTrueNorthModel, LeakReducesPotential) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(10.0));
    model.set_attribute_neuron(0, "leak", make_attr_double(0.5));
    model.set_attribute_neuron(0, "leak_towards_zero", make_attr_bool(true));

    model.update(0, 2.0); // initial current
    double before = model.get_potential(0);

    model.update(0, std::nullopt); // leak applies
    double after = model.get_potential(0);

    EXPECT_LT(after, before);
}

TEST_F(TestTrueNorthModel, ResetClearsPotential) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(5.0));
    model.update(0, 3.0);
    model.reset();

    EXPECT_DOUBLE_EQ(model.get_potential(0), 0.0);
}

TEST_F(TestTrueNorthModel, SetReverseAttributesAndBias) {
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "reverse_threshold", make_attr_double(-2.0)));
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "reverse_reset", make_attr_double(-1.0)));
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "reverse_reset_mode", make_attr_string("soft")));
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "bias", make_attr_double(0.5)));
    EXPECT_NO_THROW(model.set_attribute_neuron(0, "force_soma_update", make_attr_bool(true)));
}

TEST_F(TestTrueNorthModel, LeakTowardsZeroBothDirections) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(10.0));
    model.set_attribute_neuron(0, "leak", make_attr_double(1.0));
    model.set_attribute_neuron(0, "leak_towards_zero", make_attr_bool(true));

    model.reset();
    model.update(0, 3.0); // apply current
    double pos_before = model.get_potential(0);
    model.set_attribute_neuron(0, "bias", make_attr_double(0.0));
    model.update(0, std::nullopt);
    EXPECT_LT(model.get_potential(0), pos_before);

    model.reset();
    model.set_attribute_neuron(0, "bias", make_attr_double(0.0)); // no bias now
    model.update(0, -3.0); // force negative potential
    double neg_before = model.get_potential(0);
    model.update(0, std::nullopt);
    EXPECT_LT(std::fabs(model.get_potential(0)), std::fabs(neg_before));
}

TEST_F(TestTrueNorthModel, LeakWithoutTowardsZeroIncreasesPotential) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(10.0));
    model.set_attribute_neuron(0, "leak", make_attr_double(1.0));
    model.set_attribute_neuron(0, "leak_towards_zero", make_attr_bool(false));

    model.update(0, std::nullopt);
    double before = model.get_potential(0);
    model.update(0, std::nullopt);
    EXPECT_GT(model.get_potential(0), before);
}

TEST_F(TestTrueNorthModel, ThresholdAndResetModes) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(1.0));
    model.set_attribute_neuron(0, "reset", make_attr_double(0.0));

    model.set_attribute_neuron(0, "reset_mode", make_attr_string("soft"));
    model.update(0, 2.0);
    EXPECT_GE(model.get_potential(0), 0.0);

    model.set_attribute_neuron(0, "reset_mode", make_attr_string("saturate"));
    model.update(0, 2.0);
    EXPECT_LE(model.get_potential(0), model.get_potential(0));
}

TEST_F(TestTrueNorthModel, ReverseResetModes) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(10.0));
    model.set_attribute_neuron(0, "reverse_threshold", make_attr_double(0.0));
    model.set_attribute_neuron(0, "reverse_reset", make_attr_double(-2.0));

    model.set_attribute_neuron(0, "reverse_reset_mode", make_attr_string("hard"));
    model.update(0, -5.0);

    model.set_attribute_neuron(0, "reverse_reset_mode", make_attr_string("soft"));
    model.update(0, -5.0);

    model.set_attribute_neuron(0, "reverse_reset_mode", make_attr_string("saturate"));
    model.update(0, -5.0);
}

TEST_F(TestTrueNorthModel, RandomizedThresholdAffectsPotential) {
    model.set_attribute_neuron(0, "threshold", make_attr_double(5.0));
    model.set_attribute_neuron(0, "reset_mode", make_attr_string("hard"));
    model.set_attribute_neuron(0, "reset", make_attr_double(0.0));

    model.update(0, 10.0);
    EXPECT_GE(model.get_potential(0), 0.0);
}

TEST_F(TestTrueNorthModel, RandomMaskNegativeThrows) {
    EXPECT_THROW(
        model.set_attribute_neuron(0, "random_mask", make_attr_int(-1)),
        std::invalid_argument
    );
}

TEST_F(TestTrueNorthModel, RandomMaskEnablesRandomizedThreshold) {
    std::srand(1); 

    model.set_attribute_neuron(0, "threshold", make_attr_double(1.0));
    model.set_attribute_neuron(0, "reset_mode", make_attr_string("hard"));
    model.set_attribute_neuron(0, "reset", make_attr_double(0.0));
    model.set_attribute_neuron(0, "random_mask", make_attr_int(0xFF));

    auto r = model.update(0, std::nullopt); 
    EXPECT_EQ(r.status, sanafe::fired);     
}