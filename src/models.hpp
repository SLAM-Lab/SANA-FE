// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// models.hpp
#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <memory> // For shared_ptr<T>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "arch.hpp"
#include "fwd.hpp"
#include "pipeline.hpp"

namespace sanafe
{
constexpr int default_weight_bits = 8; // Based on real-world H/W e.g., Loihi
constexpr int loihi_max_compartments{1024};

class CurrentBasedSynapseModel : public SynapseUnit
{
public:
    CurrentBasedSynapseModel() { register_attributes(current_based_synapse_attributes); }
    CurrentBasedSynapseModel(const CurrentBasedSynapseModel &copy) = default;
    CurrentBasedSynapseModel(CurrentBasedSynapseModel &&other) = default;
    ~CurrentBasedSynapseModel() override = default;
    CurrentBasedSynapseModel &operator=(const CurrentBasedSynapseModel &other) = delete;
    CurrentBasedSynapseModel &operator=(CurrentBasedSynapseModel &&other) = delete;

    PipelineResult update(size_t synapse_address, bool read) override;
    void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) override {};
    void set_attribute_edge(size_t synapse_address, const std::string &attribute_name, const ModelAttribute &param) override;
    void reset() override;

private:
    static inline const std::set<std::string> current_based_synapse_attributes{"weight", "w"};
    std::vector<double> weights;
    double min_synaptic_resolution{0.0};
    int weight_bits{default_weight_bits};
};

// More detailed Loihi-specific synapse model
class LoihiSynapseModel : public SynapseUnit
{
public:
    LoihiSynapseModel() { register_attributes(loihi_synapse_attributes); }
    LoihiSynapseModel(const LoihiSynapseModel &copy) = default;
    LoihiSynapseModel(LoihiSynapseModel &&other) = default;
    ~LoihiSynapseModel() override = default;
    LoihiSynapseModel &operator=(const LoihiSynapseModel &other) = delete;
    LoihiSynapseModel &operator=(LoihiSynapseModel &&other) = delete;

    PipelineResult update(size_t synapse_address, bool read) override;
    void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) override;
    void set_attribute_edge(size_t synapse_address, const std::string &attribute_name, const ModelAttribute &param) override;
    void reset() override;
    void map_connection(MappedConnection &con) override;

private:
    static inline const std::set<std::string> loihi_synapse_attributes{
            "weight", "w", "latency", "latency_concurrent_access"};
    std::vector<double> weights;
    std::vector<double> groups;
    std::vector<double> costs;
    std::vector<std::pair<size_t, size_t>> concurrent_accesses;
    std::map<size_t, size_t> synapse_to_pre_neuron;
    std::optional<double> concurrent_access_latency{std::nullopt};
    double min_synaptic_resolution{0.0};
    int weight_bits{default_weight_bits};
    bool mixed_sign_mode{true};

    PipelineResult read_synapse(size_t synapse_address);
};

class AccumulatorModel : public DendriteUnit
{
public:
    AccumulatorModel() { register_attributes(accumulator_attributes); }
    AccumulatorModel(const AccumulatorModel &copy) = default;
    AccumulatorModel(AccumulatorModel &&other) = default;
    ~AccumulatorModel() override = default;
    AccumulatorModel &operator=(const AccumulatorModel &other) = delete;
    AccumulatorModel &operator=(AccumulatorModel &&other) = delete;

    PipelineResult update(size_t neuron_address, std::optional<double> current, std::optional<size_t> synapse_address) override;
    void reset() override {
    }
    void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) override {};
    void set_attribute_neuron(size_t neuron_address, const std::string &attribute_name, const ModelAttribute &param) override;
    void set_attribute_edge(size_t neuron_address, const std::string &attribute_name, const ModelAttribute &param) override {};

private:
    // Technically soma attributes, but due to common usage, suppres warnings
    //  for these in the dendrite H/W too
    static inline const std::set<std::string> accumulator_attributes{
            "reset_mode",
            "reverse_reset_mode",
            "reset",
            "reverse_reset",
            "bias",
            "threshold",
            "reverse_threshold",
            "leak_decay",
            "noise",
            "weight",
            "w",
            "latency"
    };
    std::vector<double> accumulated_charges{std::vector<double>(loihi_max_compartments, 0.0)};
    std::vector<long int> timesteps_simulated{std::vector<long int>(loihi_max_compartments, 0)};
    double leak_decay{0.0};
};

class MultiTapModel1D : public DendriteUnit
{
public:
    MultiTapModel1D() = default;
    MultiTapModel1D(const MultiTapModel1D &copy) = default;
    MultiTapModel1D(MultiTapModel1D &&other) = default;
    ~MultiTapModel1D() override = default;
    MultiTapModel1D &operator=(const MultiTapModel1D &other) = delete;
    MultiTapModel1D &operator=(MultiTapModel1D &&other) = delete;

    PipelineResult update(size_t neuron_address, std::optional<double> current, std::optional<size_t> synapse_address) override;
    void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) override {};
    void set_attribute_neuron(size_t neuron_address, const std::string &attribute_name, const ModelAttribute &param) override;
    void set_attribute_edge(size_t synapse_address, const std::string &attribute_name, const ModelAttribute &param) override;
    void reset() override;

private:
    static inline const std::set<std::string> multitap_attributes{
            "taps", "time_constants", "space_constants"};
    // Modeling a 1D dendrite with taps
    std::vector<double> tap_voltages{std::vector<double>(1, 0.0)};
    std::vector<double> next_voltages{std::vector<double>(1, 0.0)};
    std::vector<double> space_constants{std::vector<double>(0)};
    std::vector<double> time_constants{std::vector<double>(1, 0.0)};
    std::vector<int> synapse_to_tap;
    long int timesteps_simulated{0L};

    void calculate_next_state();
    void input_current(double current, std::optional<size_t> synapse_address);
    void print_taps();
};

class LoihiLifModel : public SomaUnit
{
public:
    LoihiLifModel() { register_attributes(loihi_lif_attributes); }
    LoihiLifModel(const LoihiLifModel &copy) = delete;
    LoihiLifModel(LoihiLifModel &&other) = delete;
    ~LoihiLifModel() override = default;
    LoihiLifModel &operator=(const LoihiLifModel &other) = delete;
    LoihiLifModel &operator=(LoihiLifModel &&other) = delete;

    void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) override;
    void set_attribute_neuron(size_t neuron_address, const std::string &attribute_name, const ModelAttribute &param) override;
    PipelineResult update(size_t neuron_address, std::optional<double> current_in) override;
    void reset() override;
    double get_potential(const size_t neuron_address) override
    {
        return compartments[neuron_address].potential;
    }

    struct LoihiCompartment
    {
        std::vector<double> biases;
        bool force_update{false};
        long int timesteps_simulated{0L};
        int reset_mode{neuron_reset_hard};
        int reverse_reset_mode{neuron_no_reset};
        double input_current{0.0};
        double potential{0.0};
        double leak_decay{1.0};
        double input_decay{0.0};
        double bias{0.0};
        double threshold{0.0};
        double reverse_threshold{0.0};
        double reset{0.0};
        double reverse_reset{0.0};
    };

    enum NoiseType : uint8_t
    {
        noise_none = 0U,
        noise_file_stream = 1U,
    };

private:
    static inline const std::set<std::string> loihi_lif_attributes{"reset_mode",
            "reverse_reset_mode", "reset", "reverse_reset", "bias", "threshold",
            "reverse_threshold", "leak_decay", "noise"};
    std::vector<LoihiCompartment> compartments{loihi_max_compartments};
    NoiseType noise_type{noise_none};
    std::ifstream noise_stream;

    static void loihi_leak_and_quantize(LoihiCompartment &cx);
    static bool loihi_threshold_and_reset(LoihiCompartment &cx);

    double loihi_generate_noise();
    int loihi_read_noise_stream();
    void loihi_reset_noise_stream();
};

constexpr int truenorth_max_neurons{4096};
class TrueNorthModel : public SomaUnit
{
public:
    TrueNorthModel() { register_attributes(truenorth_attributes); }
    TrueNorthModel(const TrueNorthModel &copy) = default;
    TrueNorthModel(TrueNorthModel &&other) = default;
    ~TrueNorthModel() override = default;
    TrueNorthModel &operator=(const TrueNorthModel &other) = delete;
    TrueNorthModel &operator=(TrueNorthModel &&other) = delete;

    void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) override {};
    void set_attribute_neuron(size_t neuron_address,
            const std::string &attribute_name,
            const ModelAttribute &param) override;
    PipelineResult update(size_t neuron_address,
            std::optional<double> current_in = std::nullopt) override;
    void reset() override;
    double get_potential(const size_t neuron_address) override
    {
        return neurons[neuron_address].potential;
    }

    struct TrueNorthNeuron
    {
        bool force_update{false};
        unsigned int random_range_mask{0U};
        int reset_mode{neuron_reset_hard};
        int reverse_reset_mode{neuron_no_reset};
        bool leak_towards_zero{true};
        // int noise_type;
        double potential{0.0};
        double leak{0.0}; // Default is no leak (potential decay)
        double bias{0.0};
        double threshold{0.0};
        double reverse_threshold{0.0};
        double reset{0.0};
        double reverse_reset{0.0};
    };

private:
    static inline const std::set<std::string> truenorth_attributes{
            "reset_mode",
            "reverse_reset_mode",
            "reset",
            "reverse_reset",
            "bias",
            "threshold",
            "reverse_threshold",
            "leak",
    };
    std::vector<TrueNorthNeuron> neurons{truenorth_max_neurons};

    static void truenorth_leak(TrueNorthNeuron &n);
    static bool truenorth_threshold_and_reset(TrueNorthNeuron &n);
};

class InputModel : public SomaUnit
{
public:
    InputModel() { register_attributes(input_attributes); }
    InputModel(const InputModel &copy) = delete; // Because of random_device
    InputModel(InputModel &&other) = delete; // Because of random_device
    ~InputModel() override = default;
    InputModel &operator=(const InputModel &other) = delete;
    InputModel &operator=(InputModel &&other) = delete;

    void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) override {};
    void set_attribute_neuron(size_t neuron_address, const std::string &attribute_name, const ModelAttribute &param) override;
    PipelineResult update(size_t neuron_address, std::optional<double> current_in = std::nullopt) override;
    void reset() override { send_spike = false;
    }

private:
    static inline const std::set<std::string> input_attributes{"rate", "poisson", "spikes"};
    std::vector<bool> spikes;
    std::vector<bool>::const_iterator curr_spike{spikes.begin()};
    std::uniform_real_distribution<double> uniform_distribution{0.0, 1.0};
    // Intentionally use a fixed seed to get deterministic results across runs
    // To truly randomize, use:
    //std::random_device rd;
    //std::mt19937 gen{rd()};
    std::mt19937 gen{1};
    double poisson_probability{0.0};
    double rate{0.0};
    bool send_spike{false};
};

std::shared_ptr<PipelineUnit> model_get_pipeline_unit(const std::string &model_name);
NeuronResetModes model_parse_reset_mode(const std::string &str);
}

#endif
