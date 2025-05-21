// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.hpp
//
//  Classes to specify different neuromorphic (spiking) architectures.
//  In SANA-FE, an architecture is a represented as a hierarchy of different
//  hardware tiles, cores, and spike pipeline hardware units. Within a
//  neuromorphic chip, cores share some network resources (tiles) and
//  communicate over a Network-on-Chip (NoC). An architecture contains one or
//  more network tiles, where each tile contains one or more cores. Each core
//  has a neuromorphic pipeline composed of a sequence of neural-inspired
//  hardware units: implementing synaptic, dendritic and somatic models. SANA-FE
//  uses these classes to represent an abstract architecture. The Architecture
//  TileConfiguration, and CoreConfiguration classes are later used by SANA-FE
//  to construct a SpikingChip, which is a simulation of the realized hardware.
//  One Architecture can later be used to define multiple SpikingChip
//  simulations.

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

#include "attribute.hpp"

namespace sanafe
{

class Neuron;
class MappedConnection;
struct Message;

struct TileConfiguration;
struct CoreConfiguration;
struct AxonInConfiguration;
struct PipelineUnitConfiguration;
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
    BUFFER_BEFORE_DENDRITE_UNIT = 0,
    BUFFER_INSIDE_DENDRITE_UNIT = 1,
    BUFFER_BEFORE_SOMA_UNIT = 2,
    BUFFER_INSIDE_SOMA_UNIT = 3,
    BUFFER_BEFORE_AXON_OUT_UNIT = 4,
    BUFFER_POSITIONS = 5,
};

struct ModelInfo
{
    std::map<std::string, ModelAttribute> model_attributes{};
    std::optional<std::filesystem::path> plugin_library_path{};
    std::string name{};
    bool log_energy{false};
    bool log_latency{false};
};

enum NeuronResetModes
{
    NEURON_NO_RESET = 0,
    NEURON_RESET_SOFT = 1,
    NEURON_RESET_HARD = 2,
    NEURON_RESET_SATURATE = 3,
    NEURON_RESET_MODE_COUNT = 4,
};

class Architecture
{
public:
    std::vector<TileConfiguration> tiles{};
    std::string name{};
    size_t core_count{0UL};
    int noc_width_in_tiles{1};
    int noc_height_in_tiles{1};
    int noc_buffer_size{0};
    int max_cores_per_tile{0};

    Architecture(std::string name, const NetworkOnChipConfiguration &noc);
    [[nodiscard]] std::vector<std::reference_wrapper<CoreConfiguration>> cores();
    TileConfiguration &create_tile(std::string name, const TilePowerMetrics &power_metrics);
    CoreConfiguration &create_core(std::string name, size_t parent_tile_id, const CorePipelineConfiguration &pipeline_config);
    [[nodiscard]] std::string info() const noexcept;

private:
    std::pair<int, int> calculate_tile_coordinates(const size_t tile_id);
};

Architecture load_arch(const std::filesystem::path path);

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
    bool log_energy{false};
    bool log_latency{false};
};

struct TileConfiguration
{
    std::vector<CoreConfiguration> cores{};
    TilePowerMetrics power_metrics{};
    std::string name{};
    size_t id{};
    size_t x{};
    size_t y{};

    TileConfiguration(std::string name, const size_t id,
            const TilePowerMetrics &metrics);
};

constexpr size_t default_max_neurons = 1024; // The same as Loihi 1
struct CorePipelineConfiguration
{
    BufferPosition buffer_position{BUFFER_BEFORE_SOMA_UNIT};
    size_t max_neurons_supported{default_max_neurons};
    bool log_energy{false};
    bool log_latency{false};
};

struct CoreAddress
{
    size_t parent_tile_id{};
    size_t offset_within_tile{};
    size_t id{};
};

struct CoreConfiguration
{
    CorePipelineConfiguration pipeline{};
    std::string name{};
    CoreAddress address{};

    std::vector<AxonInConfiguration> axon_in{};
    std::vector<PipelineUnitConfiguration> pipeline_hw{};
    std::vector<AxonOutConfiguration> axon_out{};

    AxonInConfiguration &create_axon_in(std::string name, const AxonInPowerMetrics &power_metrics);
    PipelineUnitConfiguration &create_hardware_unit(std::string name, const ModelInfo &model_details);
    AxonOutConfiguration &create_axon_out(std::string name, const AxonOutPowerMetrics &power_metrics);

    CoreConfiguration(std::string name, const CoreAddress &address, const CorePipelineConfiguration &pipeline);
};

struct AxonInPowerMetrics
{
    double energy_message_in{0.0};
    double latency_message_in{0.0};
};

struct AxonInConfiguration
{
    AxonInPowerMetrics metrics{};
    std::string name{};
    AxonInConfiguration(const AxonInPowerMetrics &metrics, std::string name) : metrics(metrics), name(std::move(name)) {}
};

struct PipelineUnitConfiguration
{
    ModelInfo model_info{};
    std::string name{};
    size_t tile_id{};
    size_t core_offset{};
    size_t core_id{};
    bool implements_synapse{false};
    bool implements_dendrite{false};
    bool implements_soma{false};

    PipelineUnitConfiguration(const ModelInfo &model_info, const std::string name) : model_info(model_info), name(std::move(name)) {}
};

struct AxonOutPowerMetrics
{
    double energy_message_out{0.0};
    double latency_message_out{0.0};
};

struct AxonOutConfiguration
{
    AxonOutPowerMetrics metrics{};
    std::string name{};
    AxonOutConfiguration(const AxonOutPowerMetrics &metrics, const std::string name) : metrics(metrics), name(std::move(name)) {}
};
}

#endif
