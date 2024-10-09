// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.hpp: Create a neuromorphic design based on a set of commands
//  In this simulator an architecture is a represented as a set of different
//  hardware blocks. The design (chip) is a set of tiles, connected by NoC
//  interconnect. Within each tile is one or more cores. Each core contains
//  neuromorphic computation. The neuromorphic pipeline (which seems sufficient
//  for any design) is a series of elements:

#ifndef ARCH_HEADER_INCLUDED_
#define ARCH_HEADER_INCLUDED_

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem> // For std::filesystem::path
#include <functional> // For std::reference_wrapper
#include <list>
#include <map>
#include <optional>
#include <vector>

#include "network.hpp"

namespace sanafe
{

class Neuron;
struct MappedConnection;
struct Message;

class SpikingHardware;
class Tile;
class Core;
struct AxonInUnit;
struct SynapseUnit;
struct DendriteUnit;
struct SomaUnit;
struct AxonOutUnit;

struct TileConfiguration;
struct CoreConfiguration;
struct AxonInConfiguration;
struct SynapseConfiguration;
struct DendriteConfiguration;
struct SomaConfiguration;
struct AxonOutConfiguration;
struct TilePowerMetrics;
struct AxonInPowerMetrics;
struct AxonOutPowerMetrics;
struct CorePipelineConfiguration;
struct NetworkOnChipConfiguration;

struct AxonInModel;
struct AxonOutModel;

enum BufferPosition : int
{
    BUFFER_BEFORE_DENDRITE_UNIT,
    BUFFER_BEFORE_SOMA_UNIT,
    BUFFER_BEFORE_AXON_OUT_UNIT,
    BUFFER_POSITIONS,
};

enum NoiseType
{
    NOISE_NONE = -1,
    NOISE_FILE_STREAM,
};

struct ModelInfo
{
    std::map<std::string, ModelParam> model_parameters{};
    std::optional<std::filesystem::path> plugin_library_path{};
    std::string name;
};

enum NeuronResetModes
{
    NEURON_NO_RESET,
    NEURON_RESET_SOFT,
    NEURON_RESET_HARD,
    NEURON_RESET_SATURATE,
    NEURON_RESET_MODE_COUNT,
};

// An attribute can contain a scalar value, or either a list or named set of
//  attributes i.e., attributes can be recursively defined. However,
//  in C++, variants cannot be defined recursively, so create this new class.
struct ModelParam
{
    operator bool() const
    {
        if (std::holds_alternative<bool>(value))
        {
            return std::get<bool>(value);
        }
        if (std::holds_alternative<int>(value))
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
        if (std::holds_alternative<int>(value))
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
    std::variant<bool, int, double, std::string, std::vector<ModelParam>> value;
    std::optional<std::string> name;

    // Filters control which hardware units can receive this parameter
    bool forward_to_synapse{true};
    bool forward_to_dendrite{true};
    bool forward_to_soma{true};
};

Architecture load_arch(const std::filesystem::path &path);
constexpr size_t default_max_neurons = 1024;  // The same as Loihi 1


class Architecture
{
public:
    std::vector<TileConfiguration> tiles{};
    std::string name;
    size_t core_count{0UL};
    int noc_width;
    int noc_height;
    int noc_buffer_size;
    int max_cores_per_tile{0};

    Architecture(std::string name, const NetworkOnChipConfiguration &noc);
    std::vector<std::reference_wrapper<CoreConfiguration>> cores();
    TileConfiguration &create_tile(std::string name, const TilePowerMetrics &power_metrics);
    CoreConfiguration &create_core(std::string name, size_t parent_tile_id, const CorePipelineConfiguration &pipeline_config);
    [[nodiscard]] std::string info() const;
};

struct NetworkOnChipConfiguration
{
    int width_in_tiles{1};
    int height_in_tiles{1};
    int link_buffer_size{0};
};

struct TilePowerMetrics
{
    double energy_north_hop{0.0};
    double latency_north_hop{0.0};
    double energy_east_hop{0.0};
    double latency_east_hop{0.0};
    double energy_south_hop{0.0};
    double latency_south_hop{0.0};
    double energy_west_hop{0.0};
    double latency_west_hop{0.0};
};

struct TileConfiguration
{
    std::vector<CoreConfiguration> cores{};
    TilePowerMetrics power_metrics{};
    std::string name{};
    size_t id;
    size_t x{};
    size_t y{};

    TileConfiguration(std::string name, const size_t id, const TilePowerMetrics &metrics);
};

struct CorePipelineConfiguration
{
    BufferPosition buffer_position{BUFFER_BEFORE_SOMA_UNIT};
    size_t max_neurons_supported{default_max_neurons};
};

struct CoreAddress
{
    size_t parent_tile_id;
    size_t offset_within_tile;
    size_t id;
};

struct CoreConfiguration
{
    CorePipelineConfiguration pipeline{};
    std::string name{};
    CoreAddress address{};

    std::vector<AxonInConfiguration> axon_in;
    std::vector<SynapseConfiguration> synapses;
    std::vector<DendriteConfiguration> dendrites;
    std::vector<SomaConfiguration> somas;
    std::vector<AxonOutConfiguration> axon_out;

    AxonInConfiguration &create_axon_in(const std::string &name, const AxonInPowerMetrics &power_metrics);
    SynapseConfiguration &create_synapse(const std::string &name, const ModelInfo &model_details);
    DendriteConfiguration &create_dendrite(const std::string &name, const ModelInfo &model_details);
    SomaConfiguration &create_soma(std::string name, const ModelInfo &model_details);
    AxonOutConfiguration &create_axon_out(const std::string &name, const AxonOutPowerMetrics &power_metrics);

    CoreConfiguration(std::string name, const CoreAddress &address, const CorePipelineConfiguration &pipeline);
};

struct AxonInPowerMetrics
{
    double energy_message_in{0.0};
    double latency_message_in{0.0};
};


struct AxonInConfiguration
{
    AxonInPowerMetrics metrics;
    std::string name;
    AxonInConfiguration(const AxonInPowerMetrics &metrics, std::string name) : metrics(metrics), name(name) {}
};

struct SynapseConfiguration
{
    ModelInfo model_info;
    std::string name;
    SynapseConfiguration(const ModelInfo &model_info, std::string name) : model_info(model_info), name(name) {}
};

struct DendriteConfiguration
{
    ModelInfo model_info;
    std::string name;
    DendriteConfiguration(const ModelInfo &model_info, std::string name) : model_info(model_info), name(name) {}

};

struct SomaConfiguration
{
    ModelInfo model_info;
    std::string name;
    SomaConfiguration(const ModelInfo &model_info, std::string name) : model_info(model_info), name(name) {}
};

struct AxonOutPowerMetrics
{
    double energy_message_out{0.0};
    double latency_message_out{0.0};
};

struct AxonOutConfiguration
{
    AxonOutPowerMetrics metrics;
    std::string name;
    AxonOutConfiguration(const AxonOutPowerMetrics &metrics, std::string name) : metrics(metrics), name(name) {}
};
}

#endif
