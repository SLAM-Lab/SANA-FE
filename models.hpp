// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// models.hpp
#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include <map>
#include <memory> // For shared_ptr<T>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "network.hpp"
#include "print.hpp"

namespace sanafe
{

// Forward declarations
struct Synapse;
enum NeuronStatus : int;

enum NeuronResetModes
{
    NEURON_NO_RESET,
    NEURON_RESET_SOFT,
    NEURON_RESET_HARD,
    NEURON_RESET_SATURATE,
    NEURON_RESET_MODE_COUNT,
};

// An attribute can contain a scalar value, or either a list or named set of
//  attributes i.e., attributes are recursively defined attributes. However,
//  in C++, variants cannot be defined recursively.
struct ModelParam
{
    operator bool() const
    {
        if (std::holds_alternative<bool>(value))
        {
            return std::get<bool>(value);
        }
        else if (std::holds_alternative<int>(value))
        {
            TRACE1("Warning: Casting integer value to bool type.\n");
            return (std::get<int>(value) != 0);
        }

        std::string error = "Error: Attribute ";
        if (name.has_value())
        {
            error += name.value();
        }
        error += " cannot be cast to a bool";
        throw std::runtime_error(error);
    }
    operator int() const
    {
        return std::get<int>(value);
    }
    operator double() const
    {
        if (std::holds_alternative<double>(value))
        {
            return std::get<double>(value);
        }
        else if (std::holds_alternative<int>(value))
        {
            // Assume it is safe to convert from any integer to double
            TRACE1("Warning: Casting integer value to double type.\n");
            return static_cast<double>(std::get<int>(value));
        }

        std::string error = "Error: Attribute ";
        if (name.has_value())
        {
            error += name.value();
        }
        error += " cannot be cast to a double";
        throw std::runtime_error(error);
    }

    operator std::string() const
    {
        return std::get<std::string>(value);
    }
    template <typename T> operator std::vector<T>() const
    {
        std::vector<T> cast_vector;
        const auto &value_vector = std::get<std::vector<ModelParam>>(value);
        cast_vector.reserve(value_vector.size());

        for (const auto &element : value_vector)
        {
            cast_vector.push_back(static_cast<T>(element));
        }
        return cast_vector;
    }
    template <typename T> operator std::map<std::string, ModelParam>() const
    {
        std::map<std::string, ModelParam> cast_map;
        const auto &value_vector = std::get<std::vector<ModelParam>>(value);
        for (const auto &element : value_vector)
        {
            cast_map[element.name.value()] = static_cast<T>(element);
        }
        return cast_map;
    }
    bool operator==(const ModelParam &rhs) const
    {
        return (value == rhs.value);
    }

    std::variant<bool, int, double, std::string, std::vector<ModelParam>> value;
    std::optional<std::string> name;
    // In C++17, we cannot use std::map (which would be the natural choice) with
    //  incomplete types i.e., cannot use std::map in such a recursive
    //  structure. Considering this, and the fact that performance is not as
    //  important for this struct, label every attribute with a name and if the
    //  user wants to use "map" style lookups e.g., foo = attribute["key"]
    //  then support casting the struct to a std::map.
    //  There have been other discussions on this topic e.g., for implementing
    //  JSON and YAML parsers, but they end up either requiring Boost or other
    //  dependencies, and / or rely on undefined C++ behavior and generally
    //  require complex solutions.
};

struct ModelInfo
{
    std::map<std::string, ModelParam> model_parameters{};
    std::optional<std::filesystem::path> plugin_library_path{};
    std::string name;
};

class SynapseModel
{
public:
    struct SynapseModelResult
    {
        double current;
        std::optional<double> energy{std::nullopt};
        std::optional<double> latency{std::nullopt};
    };

    SynapseModel() = default;
    SynapseModel(const SynapseModel &copy) = default;
    SynapseModel(SynapseModel &&other) = default;
    virtual ~SynapseModel() = default;
    SynapseModel &operator=(const SynapseModel &other) = default;
    SynapseModel &operator=(SynapseModel &&other) = default;

    virtual SynapseModelResult update(bool read = false) = 0;
    virtual void set_attributes(
            const std::map<std::string, ModelParam> &attr) = 0;

    // Additional helper functions
    void set_time(const long int timestep)
    {
        simulation_time = timestep;
    }

protected:
    long int simulation_time{0L};
};

class DendriteModel
{
public:
    struct DendriteModelResult
    {
        double current;
        std::optional<double> energy{std::nullopt};
        std::optional<double> latency{std::nullopt};
    };

    DendriteModel(const DendriteModel &copy) = default;
    DendriteModel(DendriteModel &&other) = default;
    virtual ~DendriteModel() = default;
    DendriteModel &operator=(const DendriteModel &other) = default;
    DendriteModel &operator=(DendriteModel &&other) = default;

    virtual void set_attributes(
            const std::map<std::string, ModelParam> &attributes) = 0;
    virtual DendriteModelResult update(std::optional<Synapse> synapse_in) = 0;

    // Additional helper functions
    void set_time(const long int timestep)
    {
        simulation_time = timestep;
    }

protected:
    DendriteModel() = default; // Abstract base class; do not instantiate

    long int simulation_time{0L};
};

class SomaModel
{
public:
    struct SomaModelResult
    {
        NeuronStatus status;
        std::optional<double> energy{std::nullopt};
        std::optional<double> latency{std::nullopt};
    };

    SomaModel() = default;
    SomaModel(const SomaModel &copy) = default;
    SomaModel(SomaModel &&other) = default;
    virtual ~SomaModel() = default;
    SomaModel &operator=(const SomaModel &other) = delete;
    SomaModel &operator=(SomaModel &&other) = delete;

    void set_time(const long int timestep)
    {
        simulation_time = timestep;
    }

    virtual SomaModelResult update(std::optional<double> current_in) = 0;
    virtual void set_attributes(
            const std::map<std::string, ModelParam> &attr) = 0;
    virtual double get_potential()
    {
        return 0.0;
    }

protected:
    long int simulation_time{0L};
};

constexpr int default_weight_bits = 8; // Based on real-world H/W e.g., Loihi
class CurrentBasedSynapseModel : public SynapseModel
{
public:
    CurrentBasedSynapseModel() = default;
    CurrentBasedSynapseModel(const CurrentBasedSynapseModel &copy) = default;
    CurrentBasedSynapseModel(CurrentBasedSynapseModel &&other) = default;
    ~CurrentBasedSynapseModel() override = default;
    CurrentBasedSynapseModel &operator=(
            const CurrentBasedSynapseModel &other) = default;
    CurrentBasedSynapseModel &operator=(
            CurrentBasedSynapseModel &&other) = default;

    SynapseModelResult update(bool read) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;

private:
    double weight{0.0};
    double min_synaptic_resolution{0.0};
    int weight_bits{default_weight_bits};
};

class SingleCompartmentModel : public DendriteModel
{
public:
    SingleCompartmentModel() = default;
    SingleCompartmentModel(const SingleCompartmentModel &copy) = default;
    SingleCompartmentModel(SingleCompartmentModel &&other) = default;
    ~SingleCompartmentModel() override = default;
    SingleCompartmentModel &operator=(
            const SingleCompartmentModel &other) = default;
    SingleCompartmentModel &operator=(SingleCompartmentModel &&other) = default;

    DendriteModelResult update(std::optional<Synapse> synapse_in) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;

private:
    double accumulated_charge{0.0};
    double leak_decay{0.0};
    long int timesteps_simulated{0L};
};

class MultiTapModel1D : public DendriteModel
{
public:
    MultiTapModel1D() = default;
    MultiTapModel1D(const MultiTapModel1D &copy) = default;
    MultiTapModel1D(MultiTapModel1D &&other) = default;
    ~MultiTapModel1D() override = default;
    MultiTapModel1D &operator=(const MultiTapModel1D &other) = default;
    MultiTapModel1D &operator=(MultiTapModel1D &&other) = default;

    DendriteModelResult update(std::optional<Synapse> synapse_in) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;

private:
    // Modeling a 1D dendrite with taps
    std::vector<double> tap_voltages{std::vector<double>(1, 0.0)};
    std::vector<double> next_voltages{std::vector<double>(1, 0.0)};
    std::vector<double> space_constants{std::vector<double>(0)};
    std::vector<double> time_constants{std::vector<double>(1, 0.0)};
    long int timesteps_simulated{0L};
};

class LoihiLifModel : public SomaModel
{
public:
    LoihiLifModel() = default;
    LoihiLifModel(const LoihiLifModel &copy) = default;
    LoihiLifModel(LoihiLifModel &&other) = default;
    ~LoihiLifModel() override = default;
    LoihiLifModel &operator=(const LoihiLifModel &other) = delete;
    LoihiLifModel &operator=(LoihiLifModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    SomaModelResult update(std::optional<double> current_in) override;
    double get_potential() override
    {
        return potential;
    }

private:
    bool force_update{false};
    long int timesteps_simulated{0L};
    int reset_mode{NEURON_RESET_HARD};
    int reverse_reset_mode{NEURON_NO_RESET};
    //int noise_type;
    double potential{0.0};
    double leak_decay{1.0};
    double bias{0.0};
    double threshold{0.0};
    double reverse_threshold{0.0};
    double reset{0.0};
    double reverse_reset{0.0};
};

class TrueNorthModel : public SomaModel
{
public:
    TrueNorthModel() = default;
    TrueNorthModel(const TrueNorthModel &copy) = default;
    TrueNorthModel(TrueNorthModel &&other) = default;
    ~TrueNorthModel() override = default;
    TrueNorthModel &operator=(const TrueNorthModel &other) = delete;
    TrueNorthModel &operator=(TrueNorthModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    SomaModelResult update(
            std::optional<double> current_in = std::nullopt) override;
    double get_potential() override
    {
        return potential;
    }

private:
    bool force_update{false};
    unsigned int random_range_mask{0U};
    int reset_mode{NEURON_RESET_HARD};
    int reverse_reset_mode{NEURON_NO_RESET};
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

class InputModel : public SomaModel
{
public:
    InputModel() = default;
    InputModel(const InputModel &copy) = default;
    InputModel(InputModel &&other) = default;
    ~InputModel() override = default;
    InputModel &operator=(const InputModel &other) = delete;
    InputModel &operator=(InputModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    SomaModelResult update(
            std::optional<double> current_in = std::nullopt) override;

private:
    std::vector<bool> spikes{};
    std::vector<bool>::const_iterator curr_spike{spikes.begin()};
    std::uniform_real_distribution<double> uniform_distribution{0.0, 1.0};
    std::random_device rd;
    std::mt19937 gen{rd()};
    double poisson_probability{0.0};
    bool send_spike{false};
};

NeuronResetModes model_parse_reset_mode(const std::string &str);
std::shared_ptr<SynapseModel> model_get_synapse(const std::string &model_name);
std::shared_ptr<DendriteModel> model_get_dendrite(
        const std::string &model_name);
std::shared_ptr<SomaModel> model_get_soma(const std::string &model_name);

}

#endif
