// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ios>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "arch.hpp"
#include "attribute.hpp"
#include "mapped.hpp"
#include "models.hpp"
#include "pipeline.hpp"
#include "print.hpp"

// *** Synapse hardware unit models ***
sanafe::PipelineResult sanafe::CurrentBasedSynapseModel::update(
        const size_t synapse_address, const bool read)
{
    PipelineResult output;
    if (read)
    {
        TRACE1(MODELS, "w:%lf\n", weights[synapse_address]);
        output.current = weights[synapse_address];
    }
    else
    {
        output.current = 0.0;
    }
    return output;
}

void sanafe::CurrentBasedSynapseModel::set_attribute_edge(
        const size_t synapse_address, const std::string &attribute_name,
        const ModelAttribute &param)
{
    if (weights.size() <= synapse_address)
    {
        TRACE1(MODELS, "Resizing weights to: %zu\n", synapse_address + 1);
        weights.resize(std::max(weights.size() * 2, synapse_address + 1));
    }

    if ((attribute_name == "w") || (attribute_name == "weight"))
    {
        TRACE1(MODELS, "Setting weight at address:%zu = %lf\n", synapse_address,
                static_cast<double>(param));
        weights[synapse_address] = static_cast<double>(param);
    }

    min_synaptic_resolution = (1.0 / weight_bits);
}

void sanafe::CurrentBasedSynapseModel::reset()
{
}

sanafe::PipelineResult sanafe::LoihiSynapseModel::read_synapse(
        const size_t synapse_address)
{
    constexpr size_t max_parallel_accesses = 3;
    // Model the concurrent synaptic weights being access
    double latency = std::numeric_limits<double>::quiet_NaN();
    if (concurrent_access_latency.has_value())
    {
        const size_t current_sending_neuron =
                synapse_to_pre_neuron.at(synapse_address);
        if (!concurrent_accesses.empty())
        {
            const auto [first_access, first_sending_neuron] =
                    concurrent_accesses.front();
            TRACE2(MODELS, "first access address:%zu\n", first_access);
            if (concurrent_accesses.size() >= max_parallel_accesses ||
                    first_sending_neuron != current_sending_neuron)
            {
                concurrent_accesses.clear();
            }
        }

        const bool is_first_access = concurrent_accesses.empty();
        if (is_first_access)
        {
            latency = max_parallel_accesses * concurrent_access_latency.value();
        }
        else
        {
            latency = 0.0;
        }
        concurrent_accesses.emplace_back(
                synapse_address, current_sending_neuron);
    }
    else // Used latency tagged cost in file
    {
        latency = costs[synapse_address];
    }

    TRACE1(MODELS, "w:%lf\n", weights[synapse_address]);
    PipelineResult result{};
    result.current = weights[synapse_address];
    result.latency = latency;

    return result;
}

sanafe::PipelineResult sanafe::LoihiSynapseModel::update(
        const size_t synapse_address, const bool read)
{
    // TODO: add lookup table with different synapse read costs?
    // More detailed Loihi synaptic model
    // Either use a latency cost per synapse from earlier conversion/mapping stages
    //  (most detailed) e.g., using additional info from the Loihi synapse
    //  compiler to tag each edge. Or, if this isn't available, use a slightly
    //  more detailed h/w model that considers the parallel 4-way processing
    //  in the synapse unit, noting that we only parallelize reads from a single
    //  spike message
    if (read)
    {
        return read_synapse(synapse_address);
    }

    // No read, don't generate current but return 0 latency
    PipelineResult result{};
    result.latency = 0.0;
    return result;
}

void sanafe::LoihiSynapseModel::reset()
{
    concurrent_accesses.clear();
}

void sanafe::LoihiSynapseModel::set_attribute_hw(
        const std::string &attribute_name, const ModelAttribute &param)
{
    if (attribute_name == "latency_concurrent_access")
    {
        concurrent_access_latency = static_cast<double>(param);
    }
}

void sanafe::LoihiSynapseModel::set_attribute_edge(const size_t synapse_address,
        const std::string &attribute_name, const ModelAttribute &param)
{
    if (weights.size() <= synapse_address)
    {
        TRACE1(MODELS, "Resizing weights to: %zu\n", synapse_address + 1);
        weights.resize(std::max(weights.size() * 2, synapse_address + 1));
    }
    if (costs.size() <= synapse_address)
    {
        costs.resize(synapse_address + 1, 0.0);
    }
    if (groups.size() <= synapse_address)
    {
        groups.resize(synapse_address + 1, -1);
    }

    if ((attribute_name == "w") || (attribute_name == "weight"))
    {
        TRACE1(MODELS, "Setting weight at address:%zu = %lf\n", synapse_address,
                static_cast<double>(param));
        weights[synapse_address] = static_cast<double>(param);
    }
    else if (attribute_name == "mixed")
    {
        mixed_sign_mode = static_cast<bool>(param);
    }
    else if (attribute_name == "latency")
    {
        costs[synapse_address] = static_cast<double>(param);
    }
    else if (attribute_name == "g")
    {
        groups[synapse_address] = static_cast<int>(param);
    }

    min_synaptic_resolution = (1.0 / weight_bits);
}

void sanafe::LoihiSynapseModel::map_connection(MappedConnection &con)
{
    const size_t synapse_address = con.synapse_address;
    // Get unique identifier for spiking neuron
    const MappedNeuron &pre_neuron = con.pre_neuron_ref;
    synapse_to_pre_neuron[synapse_address] = pre_neuron.id;
}

// *** Dendrite models ***
sanafe::PipelineResult sanafe::AccumulatorModel::update(size_t neuron_address,
        std::optional<double> current,
        std::optional<size_t> /*synapse_address*/)
{
    PipelineResult output;

    while (timesteps_simulated[neuron_address] < simulation_time)
    {
        // Apply leak for 1 or more timesteps
        ++(timesteps_simulated[neuron_address]);
        accumulated_charges[neuron_address] *= leak_decay;
    }
    if (current.has_value())
    {
        // Integrate input charges
        accumulated_charges[neuron_address] += current.value();
    }
    output.current = accumulated_charges[neuron_address];

    return output;
}

void sanafe::AccumulatorModel::set_attribute_neuron(
        const size_t /*neuron_address*/, const std::string &attribute_name,
        const ModelAttribute &param)
{
    if (attribute_name == "dendrite_leak_decay")
    {
        leak_decay = static_cast<double>(param);
    }
}

void sanafe::MultiTapModel1D::calculate_next_state()
{
    const size_t taps = tap_voltages.size();

    for (size_t i = 0; i < tap_voltages.size(); i++)
    {
        TRACE1(MODELS, "\tv[%zu]: %lf\n", i, tap_voltages[i]);
    }
    TRACE1(MODELS, "^^^\n");

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

void sanafe::MultiTapModel1D::print_taps()
{
    for (size_t i = 0; i < tap_voltages.size(); i++)
    {
        TRACE1(MODELS, "\tv[%zu]: %lf\n", i, tap_voltages[i]);
    }
}

void sanafe::MultiTapModel1D::input_current(const double current,
        const std::optional<size_t> synapse_address)
{
    constexpr int default_tap = 0;

    int tap = default_tap;
    if (synapse_address.has_value() &&
            synapse_address.value() < synapse_to_tap.size())
    {
        tap = synapse_to_tap[synapse_address.value()];
    }
    if ((tap < 0) || (static_cast<size_t>(tap) >= tap_voltages.size()))
    {
        const std::string error(
                "Error: tap should be >= 0 and less than taps.\n");
        throw std::logic_error(error);
    }

    tap_voltages[tap] += current;
    TRACE2(MODELS, "Adding current:%lf to tap %d\n", current, tap);
}

sanafe::PipelineResult sanafe::MultiTapModel1D::update(size_t neuron_address,
        std::optional<double> current, std::optional<size_t> synapse_address)
{
    while (timesteps_simulated < simulation_time)
    {
        TRACE1(MODELS, "***\n");
        ++timesteps_simulated;
        calculate_next_state();
    }

    if (current.has_value())
    {
        input_current(current.value(), synapse_address);
    }

    TRACE1(MODELS, "Tap voltages after update for address:%zu\n",
            neuron_address);
    print_taps();
    TRACE1(MODELS, "***\n");

    // Return current for most proximal tap (which is always the first tap)
    PipelineResult output{};
    output.current = tap_voltages[0];
    return output;
}

void sanafe::MultiTapModel1D::set_attribute_neuron(const size_t /*address*/,
        const std::string &attribute_name, const ModelAttribute &param)
{
    if (attribute_name == "taps")
    {
        const size_t n_taps = static_cast<int>(param);
        if (n_taps == 0)
        {
            throw std::invalid_argument("Error: Number of taps must be > 0\n");
        }
        if (tap_voltages.size() > n_taps)
        {
            INFO("Warning: Reducing number of taps, constants might be "
                 "wrong.\n");
        }
        tap_voltages.resize(n_taps);
        next_voltages.resize(n_taps);
        time_constants.resize(n_taps);
        space_constants.resize(n_taps - 1);
    }
    else if (attribute_name == "time_constants")
    {
        time_constants = static_cast<std::vector<double>>(param);
        const size_t n_taps = time_constants.size();
        if (time_constants.size() < n_taps)
        {
            const std::string error = "Error: Expected " +
                    std::to_string(n_taps) + " but received " +
                    std::to_string(time_constants.size()) + "time constants.";
            throw std::invalid_argument(error);
        }
        if (time_constants.size() > n_taps)
        {
            // Extend the other defined constants
            tap_voltages.resize(n_taps);
            next_voltages.resize(n_taps);
            space_constants.resize(n_taps - 1);
        }
    }
    else if (attribute_name == "space_constants")
    {
        space_constants = static_cast<std::vector<double>>(param);
        const size_t n_taps = space_constants.size() + 1;
        if (space_constants.size() < (n_taps - 1))
        {
            const std::string error = "Error: Expected " +
                    std::to_string(n_taps - 1) + " but received " +
                    std::to_string(time_constants.size()) + "time constants.";
            throw std::invalid_argument(error);
        }
        if (space_constants.size() > (n_taps - 1))
        {
            // Extend the other defined constants
            tap_voltages.resize(n_taps);
            next_voltages.resize(n_taps);
            time_constants.resize(n_taps);
        }
    }
    else
    {
        INFO("Warning: attribute '%s' not recognized.\n",
                attribute_name.c_str());
    }
}

void sanafe::MultiTapModel1D::set_attribute_edge(const size_t address,
        const std::string &attribute_name, const ModelAttribute &param)
{
    if (attribute_name == "tap")
    {
        // Each connection/synapse will specify a destination tap
        if (synapse_to_tap.size() <= address)
        {
            synapse_to_tap.resize(address + 1, 0);
        }
        synapse_to_tap[address] = static_cast<int>(param);
    }
}

void sanafe::MultiTapModel1D::reset()
{
    const size_t n_taps = tap_voltages.size();
    for (size_t tap = 0; tap < n_taps; ++tap)
    {
        tap_voltages[tap] = 0.0;
        next_voltages[tap] = 0.0;
    }
}

// **** Soma hardware unit models ****
void sanafe::LoihiLifModel::set_attribute_hw(
        const std::string &attribute_name, const ModelAttribute &param)
{
    if (attribute_name == "noise")
    {
        const std::string noise_filename = static_cast<std::string>(param);
        noise_type = noise_file_stream;
        noise_stream.open(noise_filename);
        TRACE1(MODELS, "Opening noise str: %s\n", noise_filename.c_str());
        if (!noise_stream.is_open())
        {
            INFO("Error: Failed to open noise stream: %s.\n",
                    noise_filename.c_str());
            throw std::runtime_error("Error: Failed to open noise stream");
        }
    }
}

void sanafe::LoihiLifModel::set_attribute_neuron(const size_t neuron_address,
        const std::string &attribute_name, const ModelAttribute &param)
{
    LoihiCompartment &cx = compartments[neuron_address];

    if (attribute_name == "threshold")
    {
        cx.threshold = static_cast<double>(param);
    }
    else if (attribute_name == "reverse_threshold")
    {
        cx.reverse_threshold = static_cast<double>(param);
    }
    else if (attribute_name == "reset")
    {
        cx.reset = static_cast<double>(param);
    }
    else if (attribute_name == "reverse_reset")
    {
        cx.reverse_reset = static_cast<double>(param);
    }
    else if (attribute_name == "reset_mode")
    {
        const std::string reset_mode_str = static_cast<std::string>(param);
        cx.reset_mode = model_parse_reset_mode(reset_mode_str);
    }
    else if (attribute_name == "reverse_reset_mode")
    {
        const std::string reverse_reset_mode_str =
                static_cast<std::string>(param);
        cx.reverse_reset_mode = model_parse_reset_mode(reverse_reset_mode_str);
    }
    else if (attribute_name == "leak_decay")
    {
        cx.leak_decay = static_cast<double>(param);
    }
    else if (attribute_name == "input_decay")
    {
        cx.input_decay = static_cast<double>(param);
    }
    else if (attribute_name == "bias")
    {
        cx.bias = static_cast<double>(param);
        TRACE2(MODELS, "Setting bias of %zu=%lf\n", neuron_address, cx.bias);
    }
    else if ((attribute_name == "force_update") ||
            (attribute_name == "force_soma_update"))
    {
        cx.force_update = static_cast<bool>(param);
    }

    TRACE1(MODELS, "Set attribute: %s\n", attribute_name.c_str());
}

void sanafe::LoihiLifModel::loihi_leak_and_quantize(LoihiCompartment &cx)
{
    cx.input_current *= cx.input_decay;
    cx.potential *= cx.leak_decay;
    // TODO: formalize quantization for Loihi and remove hack
    //  Make sure we multiple all biases and thresholds by 64 in the snn
    //  description for loihi benchmarks, it shouldn't be managed here
    //  Then we can simply apply quantization without multiplication first..
    //  Again, this has a pretty tiny effect on simulator predictions
    cx.potential = static_cast<int>(cx.potential * 64.0) / 64.0;
}

bool sanafe::LoihiLifModel::loihi_threshold_and_reset(LoihiCompartment &cx)
{
    bool fired = false;

    if (cx.potential > cx.threshold)
    {
        if (cx.reset_mode == sanafe::neuron_reset_hard)
        {
            cx.potential = cx.reset;
        }
        else if (cx.reset_mode == sanafe::neuron_reset_soft)
        {
            cx.potential -= cx.threshold;
        }
        TRACE1(MODELS, "Neuron fired.\n");
        fired = true;
    }
    // Check against reverse threshold
    if (cx.potential < cx.reverse_threshold)
    {
        if (cx.reverse_reset_mode == sanafe::neuron_reset_soft)
        {
            cx.potential -= cx.reverse_threshold;
        }
        else if (cx.reverse_reset_mode == sanafe::neuron_reset_hard)
        {
            cx.potential = cx.reverse_reset;
        }
        else if (cx.reverse_reset_mode == sanafe::neuron_reset_saturate)
        {
            cx.potential = cx.reverse_threshold;
        }
    }

    return fired;
}

sanafe::PipelineResult sanafe::LoihiLifModel::update(
        const size_t neuron_address, const std::optional<double> current_in)
{
    LoihiCompartment &cx = compartments[neuron_address];
    if (cx.timesteps_simulated == simulation_time)
    {
        throw std::runtime_error("Error: This Loihi model does not support "
                                 "multiple updates to the same compartment "
                                 "in one time-step.");
    }
    if (cx.timesteps_simulated < (simulation_time - 1))
    {
        throw std::runtime_error(
                "Error: This Loihi model must update every time-step.\n");
    }

    // Calculate the change in potential since the last update e.g.
    //  integate inputs and apply any potential leak
    TRACE1(MODELS, "Updating potential (cx:%zu), before:%lf\n", neuron_address,
            cx.potential);
    sanafe::NeuronStatus state = sanafe::idle;
    // Update soma, if there are any received spikes, there is a non-zero
    //  bias or we force the neuron to update every time-step
    if ((std::fabs(cx.potential) > 0.0) || current_in.has_value() ||
            (std::fabs(cx.bias) > 0.0) || cx.force_update)
    {
        // Neuron is turned on and potential write
        state = sanafe::updated;
    }

    loihi_leak_and_quantize(cx);
    // Add randomized noise to potential if enabled
    if (noise_type == noise_file_stream)
    {
        cx.potential += loihi_generate_noise();
    }
    // Add the synaptic / dendrite current to the potential
    TRACE1(MODELS, "bias:%lf potential before:%lf\n", cx.bias, cx.potential);
    cx.potential += cx.bias;

    if (current_in.has_value())
    {
        cx.input_current += current_in.value();
    }
    cx.potential += cx.input_current;
    TRACE1(MODELS, "Updating potential (nid:%zu), after:%lf\n", neuron_address,
            cx.potential);

    // Check against threshold potential (for spiking)
    if (loihi_threshold_and_reset(cx))
    {
        state = sanafe::fired;
    }

    ++(cx.timesteps_simulated);
    PipelineResult output{};
    output.status = state;
    return output;
}

void sanafe::LoihiLifModel::reset()
{
    for (LoihiCompartment &cx : compartments)
    {
        cx.input_current = 0.0;
        cx.potential = 0.0;
    }
}

void sanafe::LoihiLifModel::loihi_reset_noise_stream()
{
    // If we get to the end of the stream, by default reset it.
    //  However, it is unlikely the stream will be correct at this
    //  point
    INFO("Warning: At the end of the noise stream. Random values are unlikely "
         "to be correct.\n");
    noise_stream.clear();
    noise_stream.seekg(0, std::ios::beg);
}

int sanafe::LoihiLifModel::loihi_read_noise_stream()
{
    int random_val = 0;
    // With a noise stream, we have a file containing a series of
    //  random values. This is useful if we want to exactly
    //  replicate h/w without knowing how the stream is generated.
    //  We can record the random sequence and replicate it here
    if (!noise_stream.is_open())
    {
        INFO("Error: Noise stream is not open.\n");
        throw std::runtime_error("Noise stream is not open");
    }

    // Peek ahead to see if we're at the end of the file without consuming
    if (noise_stream.eof() ||
            (noise_stream.peek() == std::ifstream::traits_type::eof()))
    {
        loihi_reset_noise_stream();
    }

    std::string noise_str;
    if (std::getline(noise_stream, noise_str))
    {
        std::istringstream iss(noise_str);
        if (!(iss >> random_val))
        {
            INFO("Error: invalid noise stream entry: %s.\n", noise_str.c_str());
        }

        TRACE2(MODELS, "Generated random val: %d\n", random_val);
    }
    else
    {
        INFO("Error: Couldn't read noise entry from file\n");
        throw std::runtime_error("Couldn't read noise entry");
    }

    return random_val;
}

double sanafe::LoihiLifModel::loihi_generate_noise()
{
    int random_val = 0;

    if (noise_type == noise_file_stream)
    {
        random_val = loihi_read_noise_stream();
    }
    // else, don't generate any noise and return 0.0

    // Get the number of noise bits required TODO: generalize
    const int sign_bit = random_val & 0x100;
    random_val &= 0x7f; // TODO: hack, fixed for 8 bits
    if (sign_bit != 0)
    {
        // Sign extend
        random_val |= ~(0x7f);
    }

    return static_cast<double>(random_val);
}

void sanafe::TrueNorthModel::set_attribute_neuron(const size_t neuron_address,
        const std::string &attribute_name, const ModelAttribute &param)
{
    TrueNorthNeuron &n = neurons[neuron_address];
    if (attribute_name == "threshold")
    {
        n.threshold = static_cast<double>(param);
    }
    else if (attribute_name == "reverse_threshold")
    {
        n.reverse_threshold = static_cast<double>(param);
    }
    else if (attribute_name == "reset")
    {
        n.reset = static_cast<double>(param);
    }
    else if (attribute_name == "reverse_reset")
    {
        n.reverse_reset = static_cast<double>(param);
    }
    else if (attribute_name == "reset_mode")
    {
        const std::string reset_mode_str = static_cast<std::string>(param);
        n.reset_mode = model_parse_reset_mode(reset_mode_str);
    }
    else if (attribute_name == "reverse_reset_mode")
    {
        const std::string reverse_reset_mode_str =
                static_cast<std::string>(param);
        n.reverse_reset_mode = model_parse_reset_mode(reverse_reset_mode_str);
    }
    else if (attribute_name == "leak")
    {
        n.leak = static_cast<double>(param);
    }
    else if (attribute_name == "bias")
    {
        n.bias = static_cast<double>(param);
    }
    else if (attribute_name == "force_soma_update")
    {
        n.force_update = static_cast<bool>(param);
    }
    else if (attribute_name == "leak_towards_zero")
    {
        n.leak_towards_zero = static_cast<bool>(param);
    }
}

void sanafe::TrueNorthModel::truenorth_leak(TrueNorthNeuron &n)
{
    if (n.leak_towards_zero)
    {
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
}

bool sanafe::TrueNorthModel::truenorth_threshold_and_reset(TrueNorthNeuron &n)
{
    // Apply thresholding and reset
    double v = n.potential;
    const bool randomize_threshold = (n.random_range_mask != 0);
    if (randomize_threshold)
    {
        // rand() generates values using a Linear Feedback Shift Register (LFSR)
        //  which are often implemented in H/W to generate psuedorandom numbers.
        //  Checkers will complain this function is not very 'random'
        //  but here we care about emulating hardware behavior over randomness.
        // NOLINTNEXTLINE(cert-msc30-c, cert-msc50-cpp)
        const unsigned int r = std::rand() & n.random_range_mask;
        v += static_cast<double>(r);
    }
    TRACE2(MODELS, "v:%lf +vth:%lf mode:%d -vth:%lf mode:%d\n", v, n.threshold,
            n.reset_mode, n.reverse_threshold, n.reverse_reset_mode);
    bool fired = false;
    if (v >= n.threshold)
    {
        if (n.reset_mode == neuron_reset_hard)
        {
            n.potential = n.reset;
        }
        else if (n.reset_mode == neuron_reset_soft)
        {
            n.potential -= n.threshold;
        }
        else if (n.reset_mode == neuron_reset_saturate)
        {
            n.potential = n.threshold;
        }
        fired = true;
    }
    else if (v <= n.reverse_threshold)
    {
        if (n.reverse_reset_mode == neuron_reset_hard)
        {
            n.potential = n.reverse_reset;
        }
        else if (n.reverse_reset_mode == neuron_reset_soft)
        {
            n.potential += n.reverse_threshold;
        }
        else if (n.reverse_reset_mode == neuron_reset_saturate)
        {
            n.potential = n.reverse_threshold;
        }
        // No spike is generated
    }

    return fired;
}

sanafe::PipelineResult sanafe::TrueNorthModel::update(
        const size_t neuron_address, const std::optional<double> current_in)
{
    sanafe::NeuronStatus state = sanafe::idle;
    TrueNorthNeuron &n = neurons[neuron_address];

    if ((std::fabs(n.potential) > 0.0) || current_in.has_value() ||
            (std::fabs(n.bias) > 0.0) || n.force_update)
    {
        // Neuron is turned on and potential write
        state = sanafe::updated;
    }

    // Apply leak
    truenorth_leak(n);

    n.potential += n.bias;

    if (current_in.has_value())
    {
        n.potential += current_in.value();
    }
    if (truenorth_threshold_and_reset(n))
    {
        state = sanafe::fired;
    }
    TRACE2(MODELS, "potential:%lf threshold %lf\n", n.potential, n.threshold);
    PipelineResult output{};
    output.status = state;
    return output;
}

void sanafe::InputModel::set_attribute_neuron(const size_t /*neuron_address*/,
        const std::string &attribute_name, const ModelAttribute &param)
{
    TRACE1(MODELS, "Setting attribute:%s\n", attribute_name.c_str());
    if (attribute_name == "spikes")
    {
        spikes = static_cast<std::vector<bool>>(param);
        TRACE1(MODELS, "Setting input spike train (len:%zu)\n", spikes.size());
        curr_spike = spikes.begin();
    }
    else if (attribute_name == "poisson")
    {
        poisson_probability = static_cast<double>(param);
        TRACE2(MODELS, "Setting poisson probability:%lf\n",
                poisson_probability);
    }
    else if (attribute_name == "rate")
    {
        rate = static_cast<double>(param);
        TRACE2(MODELS, "Setting rate probability:%lf\n", rate);
    }
}

void sanafe::TrueNorthModel::reset()
{
    for (TrueNorthNeuron &neuron : neurons)
    {
        neuron.potential = 0.0;
    }
}

sanafe::PipelineResult sanafe::InputModel::update(
        const size_t /*neuron_address*/, std::optional<double> current_in)
{
    // This models a dummy input node
    if (current_in.has_value() && (current_in.value() != 0.0))
    {
        const std::string error =
                "Error: Current sent to input neuron which cannot "
                "be processed (" +
                std::to_string(current_in.value()) + ")";
        INFO("%s\n", error.c_str());
        throw std::runtime_error(error);
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
        TRACE2(MODELS, "Randomly generating spike (Poisson).\n");
    }

    TRACE1(MODELS, "Simulation time:%ld\n", simulation_time);
    if ((rate > 0.0) &&
            ((simulation_time % static_cast<long int>(1.0 / rate)) == 0))
    {
        send_spike = true;
        TRACE2(MODELS, "Randomly generating spikes (rate).\n");
    }

    const NeuronStatus status = send_spike ? fired : idle;
    PipelineResult output{};
    output.status = status;
    return output;
}

sanafe::NeuronResetModes sanafe::model_parse_reset_mode(const std::string &str)
{
    sanafe::NeuronResetModes reset_mode{sanafe::neuron_no_reset};

    if (str == "none")
    {
        reset_mode = sanafe::neuron_no_reset;
    }
    else if (str == "soft")
    {
        reset_mode = sanafe::neuron_reset_soft;
    }
    else if (str == "hard")
    {
        reset_mode = sanafe::neuron_reset_hard;
    }
    else if (str == "saturate")
    {
        reset_mode = sanafe::neuron_reset_saturate;
    }
    else
    {
        throw std::invalid_argument("Error: reset mode not recognized");
    }

    return reset_mode;
}

std::shared_ptr<sanafe::PipelineUnit> sanafe::model_get_pipeline_unit(
        const std::string &model_name)
{
    if (model_name == "current_based")
    {
        return std::shared_ptr<PipelineUnit>(new CurrentBasedSynapseModel());
    }
    if (model_name == "loihi")
    {
        return std::shared_ptr<PipelineUnit>(new LoihiSynapseModel());
    }
    if (model_name == "accumulator")
    {
        return std::shared_ptr<PipelineUnit>(new AccumulatorModel());
    }
    if (model_name == "taps")
    {
        return std::shared_ptr<PipelineUnit>(new MultiTapModel1D());
    }
    if (model_name == "input")
    {
        return std::shared_ptr<PipelineUnit>(new InputModel());
    }
    if (model_name == "leaky_integrate_fire")
    {
        return std::shared_ptr<PipelineUnit>(new LoihiLifModel());
    }
    if (model_name == "truenorth")
    {
        return std::shared_ptr<PipelineUnit>(new TrueNorthModel());
    }
    const std::string error =
            "Pipeline model not supported (" + model_name + ")\n";
    throw std::invalid_argument(error);
}
