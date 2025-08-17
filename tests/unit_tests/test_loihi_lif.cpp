#include "attribute.hpp"
#include "models.hpp"
#include <fstream>
#include <gtest/gtest.h>

namespace
{
class TestLoihiLifModel : public sanafe::LoihiLifModel
{
public:
    void set_simulation_time(long int t)
    {
        simulation_time = t;
    }
};
}

sanafe::ModelAttribute make_attr_double(double val)
{
    sanafe::ModelAttribute attr;
    attr.value = val;
    return attr;
}

sanafe::ModelAttribute make_attr_string(const std::string &val)
{
    sanafe::ModelAttribute attr;
    attr.value = val;
    return attr;
}

sanafe::ModelAttribute make_attr_bool(bool val)
{
    sanafe::ModelAttribute attr;
    attr.value = val;
    return attr;
}

TEST(LoihiLifModelTest, FiresWhenAboveThreshold)
{
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
    EXPECT_EQ(result.status, sanafe::NeuronStatus::fired);
    EXPECT_NEAR(neuron.get_potential(0), 0.0, 1e-6);
}

TEST(LoihiLifModelTest, DoesNotFireBelowThreshold)
{
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
    EXPECT_EQ(result.status, sanafe::NeuronStatus::updated);
    EXPECT_NEAR(neuron.get_potential(0), 50.0, 1e-6);
}

TEST(LoihiLifModelTest, StableWithoutInput)
{
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

    EXPECT_EQ(result.status, sanafe::NeuronStatus::updated);
    EXPECT_NEAR(neuron.get_potential(0), 50.0, 1e-6);
}

TEST(LoihiLifModelTest, NoiseFileFailsToOpen)
{
    TestLoihiLifModel neuron;
    EXPECT_THROW(neuron.set_attribute_hw(
                         "noise", make_attr_string("nonexistent.txt")),
            std::runtime_error);
}

TEST(LoihiLifModelTest, SetReverseAttributesAndBias)
{
    TestLoihiLifModel neuron;

    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0, "reverse_threshold", make_attr_double(-10.0)));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0, "reverse_reset", make_attr_double(-5.0)));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0, "reverse_reset_mode", make_attr_string("hard")));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0, "input_decay", make_attr_double(0.5)));
    EXPECT_NO_THROW(
            neuron.set_attribute_neuron(0, "bias", make_attr_double(1.5)));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0, "force_update", make_attr_bool(true)));
}

TEST(LoihiLifModelTest, LeakAndQuantizeReducesPotential)
{
    TestLoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "leak_decay", make_attr_double(0.5));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));

    neuron.reset();
    neuron.set_simulation_time(1);
    neuron.update(0, 80.0);

    neuron.set_simulation_time(2);
    double before = neuron.get_potential(0);
    neuron.update(0, std::nullopt);
    double after = neuron.get_potential(0);
    EXPECT_LT(after, before);
}

TEST(LoihiLifModelTest, FiresWithSoftReset)
{
    TestLoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(20.0));
    neuron.set_attribute_neuron(0, "reset_mode", make_attr_string("soft"));
    neuron.set_attribute_neuron(0, "reset", make_attr_double(5.0));

    neuron.reset();
    neuron.set_simulation_time(1);
    auto result = neuron.update(0, 25.0);
    EXPECT_EQ(result.status, sanafe::NeuronStatus::fired);
    EXPECT_GT(neuron.get_potential(0), 0.0); // soft reset subtracts threshold
}

TEST(LoihiLifModelTest, ReverseThresholdBranches)
{
    TestLoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));
    neuron.set_attribute_neuron(0, "reverse_threshold", make_attr_double(0.0));
    neuron.set_attribute_neuron(0, "reset", make_attr_double(0.0));

    neuron.reset();

    neuron.set_attribute_neuron(
            0, "reverse_reset_mode", make_attr_string("soft"));
    neuron.set_simulation_time(1);
    neuron.update(0, -10.0);

    neuron.set_attribute_neuron(
            0, "reverse_reset_mode", make_attr_string("hard"));
    neuron.set_simulation_time(2); // increment time
    neuron.update(0, -10.0);

    neuron.set_attribute_neuron(
            0, "reverse_reset_mode", make_attr_string("saturate"));
    neuron.set_simulation_time(3); // increment time
    neuron.update(0, -10.0);
}

TEST(LoihiLifModelTest, GenerateNoiseFromFile)
{
    std::ofstream file("noise_test.txt");
    file << "10\ninvalid\n20\n";
    file.close();

    TestLoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_test.txt"));

    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));
    neuron.reset();
    neuron.set_simulation_time(1);

    double before = neuron.get_potential(0);
    neuron.update(0, 10.0); // update will internally add noise
    double after = neuron.get_potential(0);

    EXPECT_NE(before, after); // potential should change due to noise or input
}

TEST(LoihiLifModelTest, NoiseFileNotOpenThrows)
{
    TestLoihiLifModel neuron;
    EXPECT_THROW(neuron.set_attribute_hw(
                         "noise", make_attr_string("nonexistent.txt")),
            std::runtime_error);
}

TEST(LoihiLifModelTest, ThrowsWhenUpdatingTwiceSameTimeStep)
{
    TestLoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.set_simulation_time(1);
    neuron.update(0, 5.0);

    EXPECT_THROW(
            neuron.update(0, 5.0), std::runtime_error); // same simulation_time
}

TEST(LoihiLifModelTest, ThrowsWhenSkippingTimestep)
{
    TestLoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.set_simulation_time(1);
    neuron.update(0, 5.0);

    neuron.set_simulation_time(3); // skip one timestep
    EXPECT_THROW(neuron.update(0, 5.0), std::runtime_error);
}

TEST(LoihiLifModelTest, AddsInputCurrentWhenProvided)
{
    TestLoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));

    neuron.reset();
    neuron.set_simulation_time(1);
    neuron.update(0, 2.0); // provide input current
    EXPECT_GT(neuron.get_potential(0), 0.0);
}

TEST(LoihiLifModelTest, ResetClearsState)
{
    TestLoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.set_simulation_time(1);
    neuron.update(0, 5.0);

    neuron.reset(); // clear potential
    EXPECT_DOUBLE_EQ(neuron.get_potential(0), 0.0);
}

TEST(LoihiLifModelTest, NoiseStreamEOFTriggersResetAndInvalidEntry)
{
    std::ofstream file("noise_test_eof.txt");
    file << "12\nbad_value\n"; // valid then invalid
    file.close();

    TestLoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_test_eof.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));
    neuron.reset();

    for (int i = 1; i <= 3; i++)
    {
        neuron.set_simulation_time(i);
        neuron.update(0, 5.0); // triggers file reads, EOF, reset
    }
}

TEST(LoihiLifModelTest, NoiseFileEmptyThrows)
{
    std::ofstream file("noise_empty.txt");
    file.close(); // empty file

    TestLoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_empty.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.set_simulation_time(1);

    EXPECT_THROW(neuron.update(0, 5.0), std::runtime_error);
}

TEST(LoihiLifModelTest, NoiseGeneratesSignBit)
{
    std::ofstream file("noise_signbit.txt");
    file << "256\n"; // triggers sign_bit != 0
    file.close();

    TestLoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_signbit.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.set_simulation_time(1);
    neuron.update(0, 1.0); // should read 256 and sign-extend
}

TEST(LoihiLifModelTest, NoiseEOFTriggersReset)
{
    std::ofstream file("noise_eof_trigger.txt");
    file << "5\n"; // Only one valid line
    file.close();

    TestLoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_eof_trigger.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    for (int i = 1; i <= 3; i++)
    {
        neuron.set_simulation_time(i);
        neuron.update(
                0, 1.0); // by 2nd or 3rd iteration, EOF will trigger reset
    }
}

TEST(LoihiLifModelTest, NoiseStreamNotOpenOnUpdateTriggersInfo)
{
    TestLoihiLifModel neuron;

    try
    {
        neuron.set_attribute_hw(
                "noise", make_attr_string("definitely_missing_noise_file.txt"));
        FAIL() << "set_attribute_hw should have thrown";
    }
    catch (const std::runtime_error &)
    {
    }

    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));
    neuron.reset();
    neuron.set_simulation_time(1);

    EXPECT_THROW(neuron.update(0, 1.0), std::runtime_error);
}

TEST(LoihiLifModelTest, SetForceSomaUpdate)
{
    TestLoihiLifModel neuron;
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0, "force_soma_update", make_attr_bool(true)));

    neuron.reset();
    neuron.set_simulation_time(1);
    auto result = neuron.update(0, std::nullopt);
    EXPECT_EQ(result.status, sanafe::NeuronStatus::updated);
}
