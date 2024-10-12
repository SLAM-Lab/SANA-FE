#include <cassert>
#include <cmath>
#include <optional>
#include <sstream>
#include <vector>

#include "description.hpp"
#include "models.hpp"
#include "network.hpp"
#include "plugins.hpp"
#include "print.hpp"

// *** Synapse models ***
sanafe::SynapseUnit::SynapseResult sanafe::CurrentBasedSynapseModel::update(
        const size_t synapse_address, const bool read)
{
    if (read)
    {
        TRACE1("w:%lf\n", weights[synapse_address]);
        return {weights[synapse_address], std::nullopt, std::nullopt};
    }
    return {0.0, std::nullopt, std::nullopt};
}

void sanafe::CurrentBasedSynapseModel::set_attribute(
        const size_t synapse_address, const std::string &param_name,
        const ModelParam &param)
{
    if (weights.size() <= synapse_address)
    {
        TRACE1("Resizing weights to: %zu\n", synapse_address + 1);
        weights.resize(std::max(weights.size() * 2, synapse_address + 1));
    }

    if ((param_name == "w") || (param_name == "weight"))
    {
        TRACE1("Setting weight at address:%zu = %lf\n", synapse_address,
                static_cast<double>(param));
        weights[synapse_address] = static_cast<double>(param);
    }

    min_synaptic_resolution = (1.0 / weight_bits);
}

// *** Dendrite models ***
sanafe::DendriteUnit::DendriteResult sanafe::AccumulatorModel::update(
        const size_t neuron_address, const std::optional<Synapse> synapse_in)
{
    while (timesteps_simulated[neuron_address] < simulation_time)
    {
        // Apply leak for 1 or more timesteps
        ++(timesteps_simulated[neuron_address]);
        accumulated_charges[neuron_address] *= leak_decay;
    }
    if (synapse_in.has_value())
    {
        // Integrate input charges
        accumulated_charges[neuron_address] += synapse_in.value().current;
    }

    return {accumulated_charges[neuron_address], std::nullopt, std::nullopt};
}

void sanafe::AccumulatorModel::set_attribute(const size_t neuron_address,
        const std::string &param_name, const ModelParam &param)
{
    if (param_name == "dendrite_leak_decay")
    {
        leak_decay = static_cast<double>(param);
    }
}

sanafe::DendriteUnit::DendriteResult sanafe::MultiTapModel1D::update(
        const size_t neuron_address, const std::optional<Synapse> synapse_in)
{
    while (timesteps_simulated < simulation_time)
    {
        ++timesteps_simulated;
        const size_t taps = tap_voltages.size();
        for (size_t t = 0; t < taps; t++)
        {
            next_voltages[t] = tap_voltages[t] * time_constants[t];
        }
        for (size_t src_tap = 0; src_tap < taps; src_tap++)
        {
            if (src_tap > 0)
            {
                const double proximal_current =
                        tap_voltages[src_tap] * space_constants[src_tap - 1];
                next_voltages[src_tap - 1] += proximal_current;
                next_voltages[src_tap] -= proximal_current;
                //INFO("t%ld proximal current:%lf ", src_tap, proximal_current);
            }
            if (src_tap < (taps - 1))
            {
                const double distal_current =
                        tap_voltages[src_tap] * space_constants[src_tap];
                next_voltages[src_tap + 1] += distal_current;
                next_voltages[src_tap] -= distal_current;
            }
        }
        for (size_t t = 0; t < taps; t++)
        {
            tap_voltages[t] = next_voltages[t];
        }

        INFO("***\n");
        for (size_t i = 0; i < tap_voltages.size(); i++)
        {
            printf("\tv[%zu]: %lf\n", i, tap_voltages[i]);
        }
        INFO("***\n");
    }

    if (synapse_in.has_value())
    {
        const Synapse &syn = synapse_in.value();
        int compartment = 0;
        const auto &tap = syn.con.dendrite_params.find("tap");
        if (tap != syn.con.dendrite_params.end())
        {
            compartment = static_cast<int>(tap->second);
        }
        assert(compartment >= 0);
        assert((size_t) compartment < tap_voltages.size());
        tap_voltages[compartment] += synapse_in.value().current;
    }

    // Return current for most proximal tap (which is always the first tap)
    return {tap_voltages[0], std::nullopt, std::nullopt};
}

void sanafe::MultiTapModel1D::set_attribute(const size_t neuron_address,
        const std::string &param_name, const ModelParam &param)
{
    if (param_name == "taps")
    {
        if (tap_voltages.size() > 1)
        {
            INFO("Warning: Redefining number of taps, constants might be "
                 "wrong.\n");
        }
        const size_t n_taps = static_cast<int>(param);
        if (n_taps == 0)
        {
            throw std::invalid_argument("Error: Number of taps must be > 0\n");
        }
        tap_voltages.resize(n_taps);
        next_voltages.resize(n_taps);
        time_constants.resize(n_taps);
        space_constants.resize(n_taps - 1);
    }
    else if (param_name == "time_constants")
    {
        time_constants = static_cast<std::vector<double>>(param);
        const size_t n_taps = tap_voltages.size();
        if (time_constants.size() != n_taps)
        {
            std::string error = "Error: Expected " + std::to_string(n_taps) +
                    " but received " + std::to_string(time_constants.size()) +
                    "time constants.";
            throw std::invalid_argument(error);
        }
    }
    else if (param_name == "space_constants")
    {
        space_constants = static_cast<std::vector<double>>(param);
        const size_t n_taps = tap_voltages.size();
        if (space_constants.size() != (n_taps - 1))
        {
            std::string error = "Error: Expected " +
                    std::to_string(n_taps - 1) + " but received " +
                    std::to_string(time_constants.size()) + "time constants.";
            throw std::invalid_argument(error);
        }
    }
    else
    {
        INFO("Warning: attribute '%s' not recognized.\n", param_name.c_str());
    }
}

// **** Soma models ****
void sanafe::LoihiLifModel::set_attribute(const size_t neuron_address,
        const std::string &param_name, const ModelParam &param)
{
    LoihiCompartment &cx = compartments[neuron_address];

    if (param_name == "threshold")
    {
        cx.threshold = static_cast<double>(param);
    }
    else if (param_name == "reverse_threshold")
    {
        cx.reverse_threshold = static_cast<double>(param);
    }
    else if (param_name == "reset")
    {
        cx.reset = static_cast<double>(param);
    }
    else if (param_name == "reverse_reset")
    {
        cx.reverse_reset = static_cast<double>(param);
    }
    else if (param_name == "reset_mode")
    {
        const std::string reset_mode_str = static_cast<std::string>(param);
        cx.reset_mode = model_parse_reset_mode(reset_mode_str);
    }
    else if (param_name == "reverse_reset_mode")
    {
        const std::string reverse_reset_mode_str =
                static_cast<std::string>(param);
        cx.reverse_reset_mode = model_parse_reset_mode(reverse_reset_mode_str);
    }
    else if (param_name == "leak_decay")
    {
        cx.leak_decay = static_cast<double>(param);
    }
    else if (param_name == "bias")
    {
        cx.bias = static_cast<double>(param);
        INFO("setting bias of %zu=%lf\n", neuron_address, cx.bias);
    }
    else if ((param_name == "force_update") ||
            (param_name == "force_soma_update"))
    {
        cx.force_update = static_cast<bool>(param);
    }

    INFO("Set param: %s\n", param_name.c_str());
}

sanafe::SomaUnit::SomaResult sanafe::LoihiLifModel::update(
        const size_t neuron_address, const std::optional<double> current_in)
{
    LoihiCompartment &cx = compartments[neuron_address];
    if (cx.timesteps_simulated == simulation_time)
    {
        throw std::runtime_error("Error: This Loihi model does not support "
                                 "multiple updates to the same compartment "
                                 "in one time-step.");
    }
    else if (cx.timesteps_simulated < (simulation_time - 1))
    {
        throw std::runtime_error(
                "Error: This Loihi model must update every time-step.\n");
    }

    // Calculate the change in potential since the last update e.g.
    //  integate inputs and apply any potential leak
    TRACE1("Updating potential (nid:%d.%d), before:%lf\n", group_id, neuron_id,
            potential);
    sanafe::NeuronStatus state = sanafe::IDLE;
    // Update soma, if there are any received spikes, there is a non-zero
    //  bias or we force the neuron to update every time-step
    if ((std::fabs(cx.potential) > 0.0) || current_in.has_value() ||
            (std::fabs(cx.bias) > 0.0) || cx.force_update)
    {
        // Neuron is turned on and potential write
        state = sanafe::UPDATED;
    }

    cx.potential *= cx.leak_decay;
    // Add randomized noise to potential if enabled
    /*
    if (noise_type == NOISE_FILE_STREAM)
    {
        // TODO: fix noise generation. This depends on which core
        //  is simulating the neuron. So somehow the neuron can still
        //  need information about which core it is executing on
        double random_potential = sim_generate_noise(n);
        n->potential += random_potential;
    }
    */
    // Add the synaptic / dendrite current to the potential
    TRACE1("bias:%lf potential before:%lf\n", bias, potential);
    cx.potential += cx.bias;

    if (current_in.has_value())
    {
        cx.potential += current_in.value();
    }
    TRACE1("Updating potential (nid:%d.%d), after:%lf\n", group_id, neuron_id,
            potential);

    // Check against threshold potential (for spiking)
    if (cx.potential > cx.threshold)
    {
        if (cx.reset_mode == sanafe::NEURON_RESET_HARD)
        {
            cx.potential = cx.reset;
        }
        else if (cx.reset_mode == sanafe::NEURON_RESET_SOFT)
        {
            cx.potential -= cx.threshold;
        }
        state = sanafe::FIRED;
        TRACE1("Neuron fired.\n");
    }
    // Check against reverse threshold
    if (cx.potential < cx.reverse_threshold)
    {
        if (cx.reverse_reset_mode == sanafe::NEURON_RESET_SOFT)
        {
            cx.potential -= cx.reverse_threshold;
        }
        else if (cx.reverse_reset_mode == sanafe::NEURON_RESET_HARD)
        {
            cx.potential = cx.reverse_reset;
        }
        else if (cx.reverse_reset_mode == sanafe::NEURON_RESET_SATURATE)
        {
            cx.potential = cx.reverse_threshold;
        }
    }

    ++(cx.timesteps_simulated);
    return {state, std::nullopt, std::nullopt};
}

void sanafe::TrueNorthModel::set_attribute(const size_t neuron_address,
        const std::string &param_name, const ModelParam &param)
{
    TrueNorthNeuron &n = neurons[neuron_address];
    if (param_name == "threshold")
    {
        n.threshold = static_cast<double>(param);
    }
    else if (param_name == "reverse_threshold")
    {
        n.reverse_threshold = static_cast<double>(param);
    }
    else if (param_name == "reset")
    {
        n.reset = static_cast<double>(param);
    }
    else if (param_name == "reverse_reset")
    {
        n.reverse_reset = static_cast<double>(param);
    }
    else if (param_name == "reset_mode")
    {
        const std::string reset_mode_str = static_cast<std::string>(param);
        n.reset_mode = model_parse_reset_mode(reset_mode_str);
    }
    else if (param_name == "reverse_reset_mode")
    {
        const std::string reverse_reset_mode_str =
                static_cast<std::string>(param);
        n.reverse_reset_mode = model_parse_reset_mode(reverse_reset_mode_str);
    }
    else if (param_name == "leak")
    {
        n.leak = static_cast<double>(param);
    }
    else if (param_name == "bias")
    {
        n.bias = static_cast<double>(param);
    }
    else if (param_name == "force_soma_update")
    {
        n.force_update = static_cast<bool>(param);
    }
    else if (param_name == "leak_towards_zero")
    {
        n.leak_towards_zero = static_cast<bool>(param);
    }
}

sanafe::SomaUnit::SomaResult sanafe::TrueNorthModel::update(
        const size_t neuron_address, const std::optional<double> current_in)
{
    bool randomize_threshold;
    sanafe::NeuronStatus state = sanafe::IDLE;
    TrueNorthNeuron &n = neurons[neuron_address];

    if ((std::fabs(n.potential) > 0.0) || current_in.has_value() ||
            (std::fabs(n.bias) > 0.0) || n.force_update)
    {
        // Neuron is turned on and potential write
        state = sanafe::UPDATED;
    }

    // Apply leak
    if (n.leak_towards_zero)
    {
        // TODO: what happens if we're above zero but by less
        //  than the leak amount (for convergent), will we
        //  oscillate between the two? Does it matter
        if (n.potential > 0.0)
        {
            // TrueNorth uses additive leak
            n.potential -= n.leak;
        }
        else if (n.potential < 0.0)
        {
            n.potential += n.leak;
        }
        // else equals zero, so no leak is applied
    }
    else
    {
        n.potential += n.leak;
    }

    n.potential += n.bias;

    if (current_in.has_value())
    {
        n.potential += current_in.value();
    }
    // Apply thresholding and reset
    double v = n.potential;
    randomize_threshold = (n.random_range_mask != 0);
    if (randomize_threshold)
    {
        // rand() generates values using a Linear Feedback Shift Register (LFSR)
        //  which are often implemented in H/W to generate psuedorandom numbers.
        //  Checkers will complain this function is not very 'random'
        //  but here we care about emulating hardware behavior over randomness.
        // codechecker_suppress [cert-msc30-c, cert-msc50-cpp]
        const unsigned int r = std::rand() & n.random_range_mask;
        v += static_cast<double>(r);
    }
    TRACE2("v:%lf +vth:%lf mode:%d -vth:%lf mode:%d\n", v, threshold,
            reset_mode, reverse_threshold, reverse_reset_mode);
    if (v >= n.threshold)
    {
        if (n.reset_mode == NEURON_RESET_HARD)
        {
            n.potential = n.reset;
        }
        else if (n.reset_mode == NEURON_RESET_SOFT)
        {
            n.potential -= n.threshold;
        }
        else if (n.reset_mode == NEURON_RESET_SATURATE)
        {
            n.potential = n.threshold;
        }
        state = sanafe::FIRED;
    }
    else if (v <= n.reverse_threshold)
    {
        if (n.reverse_reset_mode == NEURON_RESET_HARD)
        {
            n.potential = n.reverse_reset;
        }
        else if (n.reverse_reset_mode == NEURON_RESET_SOFT)
        {
            n.potential += n.reverse_threshold;
        }
        else if (n.reverse_reset_mode == NEURON_RESET_SATURATE)
        {
            n.potential = n.reverse_threshold;
        }
        // No spike is generated
    }
    TRACE2("potential:%lf threshold %lf\n", potential, threshold);
    return {state, std::nullopt, std::nullopt};
}

void sanafe::InputModel::set_attribute(const size_t neuron_address,
        const std::string &param_name, const ModelParam &param)
{
    if (param_name == "spikes")
    {
        spikes = static_cast<std::vector<bool>>(param);
        curr_spike = spikes.begin();
    }
    else if (param_name == "poisson")
    {
        poisson_probability = static_cast<double>(param);
    }
}

sanafe::SomaUnit::SomaResult sanafe::InputModel::update(
        const size_t neuron_address, std::optional<double> current_in)
{
    // This models a dummy input node
    if (current_in.has_value())
    {
        throw std::runtime_error("Error: Sending current to input node.\n");
    }

    bool send_spike = false;
    if (curr_spike != spikes.end())
    {
        send_spike = *curr_spike;
        curr_spike = std::next(curr_spike);
    }

    if (poisson_probability > uniform_distribution(gen))
    {
        send_spike = true;
    }

    const NeuronStatus status = send_spike ? FIRED : IDLE;
    return {status, std::nullopt, std::nullopt};
}

sanafe::NeuronResetModes sanafe::model_parse_reset_mode(const std::string &str)
{
    sanafe::NeuronResetModes reset_mode;

    if (str == "none")
    {
        reset_mode = sanafe::NEURON_NO_RESET;
    }
    else if (str == "soft")
    {
        reset_mode = sanafe::NEURON_RESET_SOFT;
    }
    else if (str == "hard")
    {
        reset_mode = sanafe::NEURON_RESET_HARD;
    }
    else if (str == "saturate")
    {
        reset_mode = sanafe::NEURON_RESET_SATURATE;
    }
    else
    {
        throw std::invalid_argument("Error: reset mode not recognized");
    }

    return reset_mode;
}

std::shared_ptr<sanafe::SynapseUnit> sanafe::model_get_synapse(
        const std::string &model_name)
{
    if (model_name == "current_based")
    {
        return std::shared_ptr<SynapseUnit>(new CurrentBasedSynapseModel());
    }
    const std::string error =
            "Synapse model not supported (" + model_name + ")\n";
    throw std::invalid_argument(error);
}

std::shared_ptr<sanafe::DendriteUnit> sanafe::model_get_dendrite(
        const std::string &model_name)
{
    if (model_name == "accumulator")
    {
        return std::shared_ptr<DendriteUnit>(new AccumulatorModel());
    }
    else if (model_name == "taps")
    {
        return std::shared_ptr<DendriteUnit>(new MultiTapModel1D());
    }

    throw std::invalid_argument("Model not supported.");
}

std::shared_ptr<sanafe::SomaUnit> sanafe::model_get_soma(
        const std::string &model_name)
{
    if (model_name == "input")
    {
        return std::shared_ptr<SomaUnit>(new InputModel());
    }
    if (model_name == "leaky_integrate_fire")
    {
        return std::shared_ptr<SomaUnit>(new LoihiLifModel());
    }
    else if (model_name == "truenorth")
    {
        return std::shared_ptr<SomaUnit>(new TrueNorthModel());
    }

    throw std::invalid_argument("Model not supported.");
}