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
struct Connection;
struct Message;

class Tile;
class Core;
struct AxonInUnit;
struct SynapseUnit;
struct DendriteUnit;
struct SomaUnit;
struct AxonOutUnit;

struct TilePowerMetrics;
struct AxonOutPowerMetrics;
struct CorePipelineConfiguration;
struct NetworkOnChipConfiguration;

struct AxonInModel;
struct AxonOutModel;

enum BufferPosition
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

class Architecture
{
public:
    std::vector<Tile> tiles{};
    std::string name;
    size_t core_count{0UL};
    int noc_width;
    int noc_height;
    int noc_buffer_size;
    int max_cores_per_tile{0};

    Architecture(std::string name, const NetworkOnChipConfiguration &noc);
    std::vector<std::reference_wrapper<Core>> cores();
    Tile &create_tile(std::string name, const TilePowerMetrics &power_metrics);
    Core &create_core(std::string name, size_t parent_tile_id, const CorePipelineConfiguration &pipeline_config);
    [[nodiscard]] std::string info() const;
};

Architecture load_arch(const std::filesystem::path &path);

struct Message
{
    double generation_delay{0.0};
    double network_delay{0.0};
    double receive_delay{0.0};
    double blocked_delay{0.0};
    double sent_timestamp{-std::numeric_limits<double>::infinity()};
    double received_timestamp{-std::numeric_limits<double>::infinity()};
    double processed_timestamp{-std::numeric_limits<double>::infinity()};
    long int timestep;
    int spikes{0};
    size_t hops{0UL};
    size_t src_neuron_id;
    std::string src_neuron_group_id;
    int src_x;
    int dest_x{0};
    int src_y;
    int dest_y{0};
    int src_tile_id;
    int src_core_id;
    int src_core_offset;
    int dest_tile_id{0};
    int dest_core_id{0};
    int dest_core_offset{0};
    int dest_axon_hw{0};
    int dest_axon_id{0};
    bool placeholder{true};
    bool in_noc{false};

    explicit Message(const Architecture &arch, const Neuron &n, long int timestep);
    explicit Message(const Architecture &arch, const Neuron &n, long int timestep, int axon_address);
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

struct AxonInPowerMetrics
{
    double energy_message_in{0.0};
    double latency_message_in{0.0};
};

struct AxonOutPowerMetrics
{
    double energy_message_out{0.0};
    double latency_message_out{0.0};
};

constexpr size_t default_max_neurons = 1024;  // The same as Loihi 1
struct CorePipelineConfiguration
{
    BufferPosition buffer_position{BUFFER_BEFORE_SOMA_UNIT};
    size_t max_neurons_supported{default_max_neurons};
};

struct NetworkOnChipConfiguration
{
    int width_in_tiles{1};
    int height_in_tiles{1};
    int link_buffer_size{0};
};

class Tile
{
public:
    std::vector<Core> cores{};
    std::string name;
    double energy{0.0};
    double energy_north_hop;
    double latency_north_hop;
    double energy_east_hop;
    double latency_east_hop;
    double energy_south_hop;
    double latency_south_hop;
    double energy_west_hop;
    double latency_west_hop;
    long int hops{0L};
    long int messages_received{0L};
    long int total_neurons_fired{0L};
    long int north_hops{0L};
    long int east_hops{0L};
    long int south_hops{0L};
    long int west_hops{0L};
    size_t id;
    size_t x{0};
    size_t y{0};

    explicit Tile(std::string name, size_t tile_id, const TilePowerMetrics &power_metrics);
    [[nodiscard]] int get_id() const { return id; }
    [[nodiscard]] std::string info() const;
};

struct CoreAddress
{
    size_t parent_tile_id;
    size_t offset_within_tile;
    size_t id;
};

class Core
{
public:
    std::vector<AxonInUnit> axon_in_hw;
    std::vector<std::shared_ptr<SynapseUnit>> synapse;
    std::vector<std::shared_ptr<DendriteUnit>> dendrite;
    std::vector<std::shared_ptr<SomaUnit>> soma;
    std::vector<AxonOutUnit> axon_out_hw;

    std::vector<Message *> messages_in;
    std::vector<AxonInModel> axons_in;
    std::vector<Neuron *> neurons;
    std::vector<Connection *> connections_in;
    std::vector<AxonOutModel> axons_out;

    std::list<BufferPosition> neuron_processing_units{};
    std::list<BufferPosition> message_processing_units{};
    CorePipelineConfiguration pipeline_config{};
    std::string name;
    double energy{0.0};
    double next_message_generation_delay{0.0};
    size_t id;
    size_t offset;
    size_t parent_tile_id;
    int message_count{0};

    explicit Core(std::string name, const CoreAddress &address, const CorePipelineConfiguration &pipeline);
    AxonInUnit &create_axon_in(const std::string &name, const AxonInPowerMetrics &power_metrics);
    SynapseUnit &create_synapse(const std::string &name, const ModelInfo &model_details);
    DendriteUnit &create_dendrite(const std::string &name, const ModelInfo &model_details);
    SomaUnit &create_soma(std::string name, const ModelInfo &model_details);
    AxonOutUnit &create_axon_out(const std::string &name, const AxonOutPowerMetrics &power_metrics);
    void map_neuron(Neuron &n);
    [[nodiscard]] int get_id() const { return id; }
    [[nodiscard]] int get_offset() const { return offset; }
    [[nodiscard]] std::string info() const;
};

struct AxonInUnit
{
    std::string name;
    long int spike_messages_in{0L};
    double energy{0.0};
    double time{0.0};
    double energy_spike_message;
    double latency_spike_message;

    explicit AxonInUnit(std::string axon_in_name, const AxonInPowerMetrics &power_metrics);
};

class SynapseUnit
{
public:
    struct SynapseResult
    {
        double current;
        std::optional<double> energy{std::nullopt};
        std::optional<double> latency{std::nullopt};
    };

    std::map<std::string, ModelParam> model_parameters{};
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;
    std::optional<double> default_energy_process_spike{std::nullopt};
    std::optional<double> default_latency_process_spike{std::nullopt};

    long int spikes_processed{0L};
    double energy{0.0};
    double time{0.0};
    size_t mapped_connections{0UL};

    //explicit SynapseUnit(std::string synapse_name, const CoreAddress &parent_core, const ModelInfo &model_details);
    SynapseUnit() = default;
    SynapseUnit(const SynapseUnit &copy) = default;
    SynapseUnit(SynapseUnit &&other) = default;
    virtual ~SynapseUnit() = default;
    SynapseUnit &operator=(const SynapseUnit &other) = default;
    SynapseUnit &operator=(SynapseUnit &&other) = default;

    virtual SynapseResult update(size_t synapse_address, bool read = false) = 0;
    virtual void set_attribute(size_t synapse_address, const std::string &param_name, const ModelParam &param) = 0;

    // Additional helper functions
    void set_time(const long int timestep)
    {
        simulation_time = timestep;
    }
    void configure(std::string synapse_name, const ModelInfo &model);

protected:
    long int simulation_time{0L};
};

class DendriteUnit
{
public:
    struct DendriteResult
    {
        double current;
        std::optional<double> energy{std::nullopt};
        std::optional<double> latency{std::nullopt};
    };

    std::map<std::string, ModelParam> model_parameters{};
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;
    std::optional<double> default_energy_update{std::nullopt};
    std::optional<double> default_latency_update{std::nullopt};
    double energy{0.0};
    double time{0.0};

    DendriteUnit(const DendriteUnit &copy) = default;
    DendriteUnit(DendriteUnit &&other) = default;
    virtual ~DendriteUnit() = default;
    DendriteUnit &operator=(const DendriteUnit &other) = default;
    DendriteUnit &operator=(DendriteUnit &&other) = default;

    virtual void set_attribute(size_t neuron_address, const std::string &param_name, const ModelParam &param) = 0;
    virtual DendriteResult update(size_t neuron_address, std::optional<Synapse> synapse_in) = 0;

    // Additional helper functions
    void set_time(const long int timestep)
    {
        simulation_time = timestep;
    }
    void configure(std::string dendrite_name, const ModelInfo &model_details);

protected:
    // Abstract base class; do not instantiate
    DendriteUnit() = default;

    long int simulation_time{0L};
};

class SomaUnit
{
public:
    struct SomaEnergyMetrics
    {
        double energy_update_neuron{0.0};
        double energy_access_neuron{0.0};
        double energy_spike_out{0.0};
    };

    struct SomaLatencyMetrics
    {
        double latency_update_neuron{0.0};
        double latency_access_neuron{0.0};
        double latency_spike_out{0.0};
    };

    struct SomaResult
    {
        NeuronStatus status;
        std::optional<double> energy{std::nullopt};
        std::optional<double> latency{std::nullopt};
    };

    //explicit SomaUnit(std::string soma_name, const CoreAddress &parent_core, const ModelInfo &model_details);
    SomaUnit() = default;
    SomaUnit(const SomaUnit &copy) = default;
    SomaUnit(SomaUnit &&other) = default;
    virtual ~SomaUnit() = default;
    SomaUnit &operator=(const SomaUnit &other) = delete;
    SomaUnit &operator=(SomaUnit &&other) = delete;

    virtual SomaResult update(size_t neuron_address, std::optional<double> current_in) = 0;
    virtual void set_attribute(size_t neuron_address, const std::string &param_name, const ModelParam &param) = 0;
    virtual double get_potential(size_t neuron_address)
    {
        return 0.0;
    }

    void set_time(const long int timestep)
    {
        simulation_time = timestep;
    }
    void configure(const std::string &soma_name, const ModelInfo &model_details);

    std::map<std::string, ModelParam> model_parameters{};
    FILE *noise_stream{nullptr};
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;
    long int neuron_updates{0L};
    long int neurons_fired{0L};
    long int neuron_count{0L};
    double energy{0.0};
    double time{0.0};
    std::optional<SomaEnergyMetrics> default_energy_metrics;
    std::optional<SomaLatencyMetrics> default_latency_metrics;
    int noise_type{NOISE_NONE};

protected:
    long int simulation_time{0L};
};

struct AxonOutUnit
{
    // The axon output points to a number of axons, stored at the
    //  post-synaptic core. A neuron can point to a number of these
    std::string name;
    long int packets_out{0L};
    double energy{0.0};
    double time{0.0};
    double energy_access;
    double latency_access;

    explicit AxonOutUnit(std::string axon_out_name, const CoreAddress &parent_core, const AxonOutPowerMetrics &power_metrics);
    //std::string description() const;
};

struct AxonInModel
{
    // List of all neuron connections to send spikes to
    std::vector<int> synapse_addresses{};
    Message *message{nullptr};
    int spikes_received{0};
    int active_synapses{0};
};

struct AxonOutModel
{
    int dest_axon_id{-1};
    int dest_tile_id{-1};
    int dest_core_offset{-1};
    size_t src_neuron_id{};
};

void arch_create_axons(Architecture &arch);
void arch_print_axon_summary(Architecture &arch);
void arch_map_neuron_connections(Neuron &n);
void arch_allocate_axon(Neuron &pre_neuron, Core &post_core);
void arch_add_connection_to_axon(Connection &con, Core &post_core);
}

#endif
