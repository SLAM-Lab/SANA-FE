#include "attribute.hpp"
#include "models.hpp"
#include <fstream>
#include <gtest/gtest.h>


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
    sanafe::LoihiLifModel neuron;

    neuron.set_attribute_neuron(0UL, "threshold", make_attr_double(64.0));
    neuron.set_attribute_neuron(0UL, "reset", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "reset_mode", make_attr_string("hard"));
    neuron.set_attribute_neuron(0UL, "leak_decay", make_attr_double(1.0));
    neuron.set_attribute_neuron(0UL, "input_decay", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "bias", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "force_update", make_attr_bool(false));

    neuron.reset();

    auto result = neuron.update(0UL, 80.0, 1L);
    EXPECT_EQ(result.status, sanafe::NeuronStatus::fired);
    EXPECT_NEAR(neuron.get_potential(0UL), 0.0, 1e-6);
}

TEST(LoihiLifModelTest, DoesNotFireBelowThreshold)
{
    sanafe::LoihiLifModel neuron;

    neuron.set_attribute_neuron(0UL, "threshold", make_attr_double(64.0));
    neuron.set_attribute_neuron(0UL, "reset", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "reset_mode", make_attr_string("hard"));
    neuron.set_attribute_neuron(0UL, "leak_decay", make_attr_double(1.0));
    neuron.set_attribute_neuron(0UL, "input_decay", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "bias", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "force_update", make_attr_bool(false));

    neuron.reset();

    auto result = neuron.update(0UL, 50.0, 1L);
    EXPECT_EQ(result.status, sanafe::NeuronStatus::updated);
    EXPECT_NEAR(neuron.get_potential(0UL), 50.0, 1e-6);
}

TEST(LoihiLifModelTest, StableWithoutInput)
{
    sanafe::LoihiLifModel neuron;

    neuron.set_attribute_neuron(0UL, "threshold", make_attr_double(64.0));
    neuron.set_attribute_neuron(0UL, "reset", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "reset_mode", make_attr_string("hard"));
    neuron.set_attribute_neuron(0UL, "leak_decay", make_attr_double(1.0));
    neuron.set_attribute_neuron(0UL, "input_decay", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "bias", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "force_update", make_attr_bool(false));

    neuron.reset();
    neuron.update(0UL, 50.0, 1L);

    auto result = neuron.update(0UL, std::nullopt, 2L);

    EXPECT_EQ(result.status, sanafe::NeuronStatus::updated);
    EXPECT_NEAR(neuron.get_potential(0UL), 50.0, 1e-6);
}

TEST(LoihiLifModelTest, NoiseFileFailsToOpen)
{
    sanafe::LoihiLifModel neuron;
    EXPECT_THROW(neuron.set_attribute_hw(
                         "noise", make_attr_string("nonexistent.txt")),
            std::runtime_error);
}

TEST(LoihiLifModelTest, SetReverseAttributesAndBias)
{
    sanafe::LoihiLifModel neuron;

    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0UL, "reverse_threshold", make_attr_double(-10.0)));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0UL, "reverse_reset", make_attr_double(-5.0)));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0UL, "reverse_reset_mode", make_attr_string("hard")));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0UL, "input_decay", make_attr_double(0.5)));
    EXPECT_NO_THROW(
            neuron.set_attribute_neuron(0UL, "bias", make_attr_double(1.5)));
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0UL, "force_update", make_attr_bool(true)));
}

TEST(LoihiLifModelTest, LeakAndQuantizeReducesPotential)
{
    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_neuron(0UL, "leak_decay", make_attr_double(0.5));
    neuron.set_attribute_neuron(0UL, "threshold", make_attr_double(100.0));

    neuron.reset();
    neuron.update(0, 80.0, 1L);

    double before = neuron.get_potential(0UL);
    neuron.update(0UL, std::nullopt, 2L);
    double after = neuron.get_potential(0UL);
    EXPECT_LT(after, before);
}

TEST(LoihiLifModelTest, FiresWithSoftReset)
{
    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(20.0));
    neuron.set_attribute_neuron(0, "reset_mode", make_attr_string("soft"));
    neuron.set_attribute_neuron(0, "reset", make_attr_double(5.0));

    neuron.reset();
    auto result = neuron.update(0UL, 25.0, 1L);
    EXPECT_EQ(result.status, sanafe::NeuronStatus::fired);
    EXPECT_GT(neuron.get_potential(0UL), 0.0); // soft reset subtracts threshold
}

TEST(LoihiLifModelTest, ReverseThresholdBranches)
{
    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_neuron(0UL, "threshold", make_attr_double(100.0));
    neuron.set_attribute_neuron(0UL, "reverse_threshold", make_attr_double(0.0));
    neuron.set_attribute_neuron(0UL, "reset", make_attr_double(0.0));

    neuron.reset();

    neuron.set_attribute_neuron(
            0UL, "reverse_reset_mode", make_attr_string("soft"));
    neuron.update(0UL, -10.0, 1L);

    neuron.set_attribute_neuron(
            0UL, "reverse_reset_mode", make_attr_string("hard"));
    neuron.update(0UL, -10.0, 2L);

    neuron.set_attribute_neuron(
            0UL, "reverse_reset_mode", make_attr_string("saturate"));
    neuron.update(0UL, -10.0, 3L);
}

TEST(LoihiLifModelTest, GenerateNoiseFromFile)
{
    std::ofstream file("noise_test.txt");
    file << "10\ninvalid\n20\n";
    file.close();

    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_test.txt"));

    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));
    neuron.reset();

    double before = neuron.get_potential(0UL);
    neuron.update(0UL, 10.0, 1L); // update will internally add noise
    double after = neuron.get_potential(0UL);

    EXPECT_NE(before, after); // potential should change due to noise or input
}

TEST(LoihiLifModelTest, NoiseFileNotOpenThrows)
{
    sanafe::LoihiLifModel neuron;
    EXPECT_THROW(neuron.set_attribute_hw(
                         "noise", make_attr_string("nonexistent.txt")),
            std::runtime_error);
}

TEST(LoihiLifModelTest, ThrowsWhenUpdatingTwiceSameTimeStep)
{
    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_neuron(0UL, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.update(0UL, 5.0, 1L);

    EXPECT_THROW(
            neuron.update(0UL, 5.0, 1L), std::runtime_error);  // Same sim time
}

TEST(LoihiLifModelTest, ThrowsWhenSkippingTimestep)
{
    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.update(0UL, 5.0, 1L);

    EXPECT_THROW(neuron.update(0UL, 5.0, 3L), std::runtime_error);
}

TEST(LoihiLifModelTest, AddsInputCurrentWhenProvided)
{
    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));

    neuron.reset();
    neuron.update(0UL, 2.0, 1L); // provide input current
    EXPECT_GT(neuron.get_potential(0UL), 0.0);
}

TEST(LoihiLifModelTest, ResetClearsState)
{
    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.update(0UL, 5.0, 1L);

    neuron.reset(); // clear potential
    EXPECT_DOUBLE_EQ(neuron.get_potential(0UL), 0.0);
}

TEST(LoihiLifModelTest, NoiseStreamEOFTriggersResetAndInvalidEntry)
{
    std::ofstream file("noise_test_eof.txt");
    file << "12\nbad_value\n"; // valid then invalid
    file.close();

    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_test_eof.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(100.0));
    neuron.reset();

    for (int i = 1; i <= 3; i++)
    {
        neuron.update(0UL, 5.0, i); // triggers file reads, EOF, reset
    }
}

TEST(LoihiLifModelTest, NoiseFileEmptyThrows)
{
    std::ofstream file("noise_empty.txt");
    file.close(); // empty file

    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_empty.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();

    EXPECT_THROW(neuron.update(0UL, 5.0, 1L), std::runtime_error);
}

TEST(LoihiLifModelTest, NoiseGeneratesSignBit)
{
    std::ofstream file("noise_signbit.txt");
    file << "256\n"; // triggers sign_bit != 0
    file.close();

    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_signbit.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    neuron.update(0UL, 1.0, 1L); // should read 256 and sign-extend
}

TEST(LoihiLifModelTest, NoiseEOFTriggersReset)
{
    std::ofstream file("noise_eof_trigger.txt");
    file << "5\n"; // Only one valid line
    file.close();

    sanafe::LoihiLifModel neuron;
    neuron.set_attribute_hw("noise", make_attr_string("noise_eof_trigger.txt"));
    neuron.set_attribute_neuron(0, "threshold", make_attr_double(10.0));

    neuron.reset();
    for (int i = 1; i <= 3; i++)
    {
        // by 2nd or 3rd iteration, EOF will trigger reset
        neuron.update(0UL, 1.0, i);
    }
}

TEST(LoihiLifModelTest, NoiseStreamNotOpenOnUpdateTriggersInfo)
{
    sanafe::LoihiLifModel neuron;

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

    EXPECT_THROW(neuron.update(0UL, 1.0, 1L), std::runtime_error);
}

TEST(LoihiLifModelTest, SetForceSomaUpdate)
{
    sanafe::LoihiLifModel neuron;
    EXPECT_NO_THROW(neuron.set_attribute_neuron(
            0, "force_update", make_attr_bool(true)));

    neuron.reset();
    auto result = neuron.update(0UL, std::nullopt, 1L);
    EXPECT_EQ(result.status, sanafe::NeuronStatus::updated);
}
