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
sanafe::SynapseStatus sanafe::CurrentBasedSynapseModel::update(
        const bool read, const bool step)
{
    if (step)
    {
        current *= synaptic_current_decay;
        if (fabs(current) < min_synaptic_resolution)
        {
            current = 0.0;
        }
    }
    if (read)
    {
        current += weight;
    }

    return {current, 0.0, 0.0};
}

void sanafe::CurrentBasedSynapseModel::set_attributes(
        const std::map<std::string, ModelParam> &attr)
{
    for (const auto &a : attr)
    {
        const std::string &key = a.first;
        const ModelParam &value = a.second;
        if ((key == "w") || (key == "weight"))
        {
            weight = static_cast<double>(value);
        }
    }

    min_synaptic_resolution = (1.0 / weight_bits);
}

// *** Dendrite models ***
sanafe::SingleCompartmentModel::SingleCompartmentModel()
{
}

sanafe::DendriteStatus sanafe::SingleCompartmentModel::update(
        const std::optional<Synapse> synapse_in, const bool step)
{
    if (step)
    {
        accumulated_charge *= leak_decay;
    }
    if (synapse_in.has_value())
    {
        // Integrate any input current
        accumulated_charge += synapse_in.value().current;
    }

    return {accumulated_charge, 0.0, 0.0};
}

void sanafe::SingleCompartmentModel::set_attributes(
        const std::map<std::string, ModelParam> &attr)
{
    for (const auto &a : attr)
    {
        const std::string &key = a.first;
        const ModelParam &value = a.second;
        if (key == "dendrite_leak_decay")
        {
            leak_decay = static_cast<double>(value);
        }
    }
}

sanafe::MultiTapModel1D::MultiTapModel1D()
        : tap_voltages(std::vector<double>(1, 0.0))
        , next_voltages(std::vector<double>(1, 0.0))
        , space_constants(std::vector<double>(0))
        , time_constants(std::vector<double>(1, 0.0))
{
}

sanafe::DendriteStatus sanafe::MultiTapModel1D::update(
        const std::optional<Synapse> synapse_in, const bool step)
{
    if (step)
    {
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

    INFO("***\n");
    for (size_t i = 0; i < tap_voltages.size(); i++)
    {
        printf("\tv[%zu]: %lf\n", i, tap_voltages[i]);
    }
    INFO("***\n");
    // Return current for most proximal tap (which is always the first tap)
    return {tap_voltages[0], 0.0, 0.0};
}

void sanafe::MultiTapModel1D::set_attributes(
        const std::map<std::string, ModelParam> &attr)
{
    if (attr.find("taps") != attr.end())
    {
        if (tap_voltages.size() > 1)
        {
            INFO("Warning: Redefining number of taps, constants might be "
                 "wrong.\n");
        }
        const size_t n_taps = static_cast<int>(attr.at("taps"));
        if (n_taps == 0)
        {
            throw std::invalid_argument("Error: Number of taps must be > 0\n");
        }
        tap_voltages.resize(n_taps);
        next_voltages.resize(n_taps);
        time_constants.resize(n_taps);
        space_constants.resize(n_taps - 1);
    }

    for (const auto &a : attr)
    {
        const std::string &key = a.first;
        const ModelParam &params = a.second;

        if (key == "time_constants")
        {
            time_constants = static_cast<std::vector<double>>(params);
            const size_t n_taps = tap_voltages.size();
            if (time_constants.size() != n_taps)
            {
                std::string error = "Error: Expected " +
                        std::to_string(n_taps) + " but received " +
                        std::to_string(time_constants.size()) +
                        "time constants.";
                throw std::invalid_argument(error);
            }
        }
        else if (key == "space_constants")
        {
            space_constants = static_cast<std::vector<double>>(params);
            const size_t n_taps = tap_voltages.size();
            if (space_constants.size() != (n_taps - 1))
            {
                std::string error = "Error: Expected " +
                        std::to_string(n_taps - 1) + " but received " +
                        std::to_string(time_constants.size()) +
                        "time constants.";
                throw std::invalid_argument(error);
            }
        }
        else if (key != "taps")
        {
            INFO("Warning: attribute '%s' not recognized.\n", key.c_str());
        }
    }
}

// **** Soma models ****
sanafe::LoihiLifModel::LoihiLifModel(
        std::string gid, size_t nid)
        : sanafe::SomaModel(std::move(gid), nid)
{
}

void sanafe::LoihiLifModel::set_attributes(
        const std::map<std::string, ModelParam> &attr)
{
    for (const auto &a : attr)
    {
        const std::string &key = a.first;
        const ModelParam &value = a.second;
        if (key == "threshold")
        {
            threshold = static_cast<double>(value);
        }
        else if (key == "reverse_threshold")
        {
            reverse_threshold = static_cast<double>(value);
        }
        else if (key == "reset")
        {
            reset = static_cast<double>(value);
        }
        else if (key == "reverse_reset")
        {
            reverse_reset = static_cast<double>(value);
        }
        else if (key == "reset_mode")
        {
            const std::string reset_mode_str = static_cast<std::string>(value);
            reset_mode = model_parse_reset_mode(reset_mode_str);
        }
        else if (key == "reverse_reset_mode")
        {
            const std::string reverse_reset_mode_str =
                    static_cast<std::string>(value);
            reverse_reset_mode = model_parse_reset_mode(reverse_reset_mode_str);
        }
        else if (key == "leak_decay")
        {
            leak_decay = static_cast<double>(value);
        }
        else if (key == "bias")
        {
            bias = static_cast<double>(value);
        }
        else if (key == "force_soma_update")
        {
            force_update = static_cast<bool>(value);
        }
    }
}

sanafe::SomaStatus sanafe::LoihiLifModel::update(
        const std::optional<double> current_in, const bool step)
{
    // Calculate the change in potential since the last update e.g.
    //  integate inputs and apply any potential leak
    TRACE1("Updating potential (nid:%d.%d), before:%lf\n", group_id, neuron_id,
            potential);
    sanafe::NeuronStatus state = sanafe::IDLE;
    // Update soma, if there are any received spikes, there is a non-zero
    //  bias or we force the neuron to update every time-step
    if ((std::fabs(potential) > 0.0) || current_in.has_value() ||
            (std::fabs(bias) > 0.0) || force_update)
    {
        // Neuron is turned on and potential write
        state = sanafe::UPDATED;
    }

    if (step)
    {
        potential *= leak_decay;
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
        potential += bias;
    }
    if (current_in.has_value())
    {
        potential += current_in.value();
    }
    TRACE1("Updating potential (nid:%d.%d), after:%lf\n", group_id, neuron_id,
            potential);

    // Check against threshold potential (for spiking)
    if (potential > threshold)
    {
        if (reset_mode == sanafe::NEURON_RESET_HARD)
        {
            potential = reset;
        }
        else if (reset_mode == sanafe::NEURON_RESET_SOFT)
        {
            potential -= threshold;
        }
        state = sanafe::FIRED;
        TRACE1("Neuron fired.\n");
    }
    // Check against reverse threshold
    if (potential < reverse_threshold)
    {
        if (reverse_reset_mode == sanafe::NEURON_RESET_SOFT)
        {
            potential -= reverse_threshold;
        }
        else if (reverse_reset_mode == sanafe::NEURON_RESET_HARD)
        {
            potential = reverse_reset;
        }
        else if (reverse_reset_mode == sanafe::NEURON_RESET_SATURATE)
        {
            potential = reverse_threshold;
        }
    }
    return {state, 0.0, 0.0};
}

sanafe::TrueNorthModel::TrueNorthModel(
        const std::string &gid, const size_t nid)
        : sanafe::SomaModel(gid, nid)
{
}

void sanafe::TrueNorthModel::set_attributes(
        const std::map<std::string, ModelParam> &attr)
{
    for (const auto &a : attr)
    {
        const std::string &key = a.first;
        const ModelParam &value = a.second;

        if (key == "threshold")
        {
            threshold = static_cast<double>(value);
        }
        else if (key == "reverse_threshold")
        {
            reverse_threshold = static_cast<double>(value);
        }
        else if (key == "reset")
        {
            reset = static_cast<double>(value);
        }
        else if (key == "reverse_reset")
        {
            reverse_reset = static_cast<double>(value);
        }
        else if (key == "reset_mode")
        {
            const std::string reset_mode_str = static_cast<std::string>(value);
            reset_mode = model_parse_reset_mode(reset_mode_str);
        }
        else if (key == "reverse_reset_mode")
        {
            const std::string reverse_reset_mode_str =
                    static_cast<std::string>(value);
            reverse_reset_mode = model_parse_reset_mode(reverse_reset_mode_str);
        }
        else if (key == "leak")
        {
            leak = static_cast<double>(value);
        }
        else if (key == "bias")
        {
            bias = static_cast<double>(value);
        }
        else if (key == "force_soma_update")
        {
            force_update = static_cast<bool>(value);
        }
        else if (key == "leak_towards_zero")
        {
            leak_towards_zero = static_cast<bool>(value);
        }
    }
}

sanafe::SomaStatus sanafe::TrueNorthModel::update(
        const std::optional<double> current_in, const bool step)
{
    bool randomize_threshold;
    sanafe::NeuronStatus state = sanafe::IDLE;

    if ((std::fabs(potential) > 0.0) || current_in.has_value() ||
            (std::fabs(bias) > 0.0) || force_update)
    {
        // Neuron is turned on and potential write
        state = sanafe::UPDATED;
    }
    if (step)
    {
        // Apply leak
        if (leak_towards_zero)
        {
            // TODO: what happens if we're above zero but by less
            //  than the leak amount (for convergent), will we
            //  oscillate between the two? Does it matter
            if (potential > 0.0)
            {
                // TrueNorth uses additive leak
                potential -= leak;
            }
            else if (potential < 0.0)
            {
                potential += leak;
            }
            // else equals zero, so no leak is applied
        }
        else
        {
            potential += leak;
        }

        potential += bias;
    }
    if (current_in.has_value())
    {
        potential += current_in.value();
    }
    // Apply thresholding and reset
    double v = potential;
    randomize_threshold = (random_range_mask != 0);
    if (randomize_threshold)
    {
        // rand() generates values using a Linear Feedback Shift Register (LFSR)
        //  which are often implemented in H/W to generate psuedorandom numbers.
        //  Checkers will complain this function is not very 'random'
        //  but here we care about emulating hardware behavior over randomness.
        // codechecker_suppress [cert-msc30-c, cert-msc50-cpp]
        const unsigned int r = std::rand() & random_range_mask;
        v += static_cast<double>(r);
    }
    TRACE2("v:%lf +vth:%lf mode:%d -vth:%lf mode:%d\n", v, threshold,
            reset_mode, reverse_threshold, reverse_reset_mode);
    if (v >= threshold)
    {
        if (reset_mode == NEURON_RESET_HARD)
        {
            potential = reset;
        }
        else if (reset_mode == NEURON_RESET_SOFT)
        {
            potential -= threshold;
        }
        else if (reset_mode == NEURON_RESET_SATURATE)
        {
            potential = threshold;
        }
        state = sanafe::FIRED;
    }
    else if (v <= reverse_threshold)
    {
        if (reverse_reset_mode == NEURON_RESET_HARD)
        {
            potential = reverse_reset;
        }
        else if (reverse_reset_mode == NEURON_RESET_SOFT)
        {
            potential += reverse_threshold;
        }
        else if (reverse_reset_mode == NEURON_RESET_SATURATE)
        {
            potential = reverse_threshold;
        }
        // No spike is generated
    }
    TRACE2("potential:%lf threshold %lf\n", potential, threshold);
    return {state, 0.0, 0.0};
}

void sanafe::InputModel::set_attributes(const std::map<std::string, ModelParam> &attr)
{
    for (const auto &a : attr)
    {
        const std::string &key = a.first;
        const ModelParam &value = a.second;
        if (key == "spikes")
        {
            spikes = static_cast<std::vector<bool>>(value);
            curr_spike = spikes.begin();
        }
    }
}

sanafe::SomaStatus sanafe::InputModel::update(
        std::optional<double> current_in, bool step)
{
    // This models a dummy input node; all input currents are ignored
    if (step)
    {
        if (curr_spike != spikes.end())
        {
            send_spike = *curr_spike;
            curr_spike = std::next(curr_spike);
        }
        else
        {
            send_spike = false;
        }
    }

    NeuronStatus status = IDLE;
    if (send_spike)
    {
        status = FIRED;
    }

    return {status, 0.0, 0.0};
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

std::shared_ptr<sanafe::SynapseModel> sanafe::model_get_synapse(
        const std::string &model_name, const size_t synapse_address)
{
    if (model_name == "current_based")
    {
        return std::shared_ptr<SynapseModel>(
                new CurrentBasedSynapseModel(synapse_address));
    }
    const std::string error =
            "Synapse model not supported (" + model_name + ")\n";
    throw std::invalid_argument(error);
}

std::shared_ptr<sanafe::DendriteModel> sanafe::model_get_dendrite(
        const std::string &model_name)
{
    if (model_name == "single_compartment")
    {
        return std::shared_ptr<DendriteModel>(new SingleCompartmentModel());
    }
    else if (model_name == "taps")
    {
        return std::shared_ptr<DendriteModel>(new MultiTapModel1D());
    }

    throw std::invalid_argument("Model not supported.");
}

std::shared_ptr<sanafe::SomaModel> sanafe::model_get_soma(
        const std::string &model_name, const std::string &group_id,
        const size_t id)
{
    if (model_name == "input")
    {
        return std::shared_ptr<SomaModel>(new InputModel(group_id, id));
    }
    if (model_name == "leaky_integrate_fire")
    {
        return std::shared_ptr<SomaModel>(new LoihiLifModel(group_id, id));
    }
    else if (model_name == "truenorth")
    {
        return std::shared_ptr<SomaModel>(new TrueNorthModel(group_id, id));
    }

    throw std::invalid_argument("Model not supported.");
}
