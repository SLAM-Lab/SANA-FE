#include <gtest/gtest.h>
#include "models.hpp"
#include "attribute.hpp"

using namespace sanafe;

namespace {
    class TestLoihiLifModel : public sanafe::LoihiLifModel {
    public:
        void set_simulation_time(long int t) {
            simulation_time = t;
        }
    };
}

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