// Copyright (c) 2025 - The University of Texas at Austin
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

#include "param.hpp"

namespace sanafe
{

class Neuron;
struct MappedConnection;
struct Message;

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

Architecture load_arch(const std::filesystem::path &path);

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
