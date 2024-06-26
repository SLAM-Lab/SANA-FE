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
#include <vector>

#define ARCH_INVALID_ID (-1)

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

struct AxonInModel;
class SynapseModel;
class DendriteModel;
class SomaModel;
struct AxonOutModel;

struct TilePowerMetrics;
struct SynapsePowerMetrics;
struct DendritePowerMetrics;
struct SomaPowerMetrics;
struct AxonOutPowerMetrics;
struct CorePipelineConfiguration;
struct NetworkOnChipConfiguration;

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
    // TODO: implement different random noise generators
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
    //std::string description() const;
    //void save_arch_description(const std::filesystem::path &path);
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
    int hops{0};
    std::string src_neuron_id;
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

struct SynapsePowerMetrics
{
    double energy_process_spike{0.0};
    double latency_process_spike{0.0};
};

struct DendritePowerMetrics
{
    double energy_access{0.0};
    double latency_access{0.0};
};

struct SomaPowerMetrics
{
    double energy_update_neuron{0.0};
    double latency_update_neuron{0.0};
    double energy_access_neuron{0.0};
    double latency_access_neuron{0.0};
    double energy_spike_out{0.0};
    double latency_spike_out{0.0};
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
    //std::string description() const;
};

struct CoreAddress
{
    size_t parent_tile_id;
    size_t offset_within_tile;
    size_t id;
};

struct ModelInfo
{
    std::optional<std::filesystem::path> plugin_library_path{};
    std::string name;
};

class Core
{
public:
    std::vector<AxonInUnit> axon_in_hw;
    std::vector<SynapseUnit> synapse;
    std::vector<DendriteUnit> dendrite;
    std::vector <SomaUnit> soma;
    std::vector<AxonOutUnit> axon_out_hw;

    std::vector<Message *> messages_in;
    std::vector<AxonInModel> axons_in;
    std::vector<Neuron *> neurons;
    std::vector<Connection *> connections_in;
    std::vector<AxonOutModel> axons_out;

    std::list<BufferPosition> neuron_processing_units;
    std::list<BufferPosition> message_processing_units;
    CorePipelineConfiguration pipeline_config;
    std::string name;
    double energy{0.0};
    double next_message_generation_delay{0.0};
    size_t id;
    size_t offset;
    size_t parent_tile_id;
    int message_count{0};

    explicit Core(std::string name, const CoreAddress &address, const CorePipelineConfiguration &pipeline);
    AxonInUnit &create_axon_in(const std::string &name, const AxonInPowerMetrics &power_metrics);
    SynapseUnit &create_synapse(const std::string &name, const SynapsePowerMetrics &power_metrics, const ModelInfo &model);
    DendriteUnit &create_dendrite(const std::string &name, const DendritePowerMetrics &power_metrics, const ModelInfo &model);
    SomaUnit &create_soma(std::string name, const SomaPowerMetrics &power_metrics, const ModelInfo &model);
    AxonOutUnit &create_axon_out(const std::string &name, const AxonOutPowerMetrics &power_metrics);
    void map_neuron(Neuron &n);
    [[nodiscard]] int get_id() const { return id; }
    [[nodiscard]] int get_offset() const { return offset; }
    [[nodiscard]] std::string info() const;
    //std::string description() const;
};

struct AxonInUnit
{
    std::string name;
    CoreAddress parent_core_address;
    long int spike_messages_in{0L};
    double energy{0.0};
    double time{0.0};
    double energy_spike_message;
    double latency_spike_message;

    explicit AxonInUnit(std::string axon_in_name, const CoreAddress &parent_core_address, const AxonInPowerMetrics &power_metrics);
    //std::string description() const;
};

struct SynapseUnit
{
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;
    CoreAddress parent_core_address;
    long int spikes_processed{0L};
    double energy{0.0};
    double time{0.0};
    //double energy_memory_access, latency_memory_access;
    double energy_spike_op;
    double latency_spike_op;

    explicit SynapseUnit(std::string synapse_name, const CoreAddress &parent_core, const SynapsePowerMetrics &power_metrics,  const ModelInfo &model);
    //std::string description() const;
};

struct DendriteUnit
{
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;
    CoreAddress parent_core_address;
    double energy{0.0};
    double time{0.0};
    double energy_access;
    double latency_access;
    explicit DendriteUnit(std::string dendrite_name, const CoreAddress &parent_core, const DendritePowerMetrics &power_metrics, const ModelInfo &model);
    //std::string description() const;
};

struct SomaUnit
{
    FILE *noise_stream{nullptr};
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;
    CoreAddress parent_core_address;
    long int neuron_updates{0L};
    long int neurons_fired{0L};
    long int neuron_count{0L};
    double energy{0.0};
    double time{0.0};
    double energy_update_neuron;
    double latency_update_neuron;
    double energy_access_neuron;
    double latency_access_neuron;
    double energy_spiking;
    double latency_spiking;
    int noise_type{NOISE_NONE};

    explicit SomaUnit(std::string soma_name, const CoreAddress &parent_core, const SomaPowerMetrics &power_metrics, const ModelInfo &model);
    //std::string description() const;
};

struct AxonOutUnit
{
    // The axon output points to a number of axons, stored at the
    //  post-synaptic core. A neuron can point to a number of these
    std::string name;
    CoreAddress parent_core_address;
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
    long int last_updated{0L};
    int spikes_received{0};
    int active_synapses{0};
};

struct AxonOutModel
{
    int dest_axon_id{-1};
    int dest_tile_id{-1};
    int dest_core_offset{-1};
    std::string src_neuron_id{};
};

void arch_create_axons(Architecture &arch);
void arch_print_axon_summary(Architecture &arch);
void arch_map_neuron_connections(Neuron &n);
void arch_allocate_axon(Neuron &pre_neuron, Core &post_core);
void arch_add_connection_to_axon(Connection &con, Core &post_core);
}

#endif
