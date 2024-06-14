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
struct SynapseModel;
struct DendriteModel;
struct SomaModel;
struct AxonOutModel;

struct TilePowerMetrics;
struct SynapsePowerMetrics;
struct SomaPowerMetrics;
struct CorePipelineConfiguration;
struct NetworkOnChipConfiguration;

enum BufferPosition
{
    BUFFER_BEFORE_DENDRITE_UNIT,
    BUFFER_BEFORE_SOMA_UNIT,
    BUFFER_BEFORE_AXON_OUT_UNIT,
    BUFFER_POSITIONS,
};

class Architecture
{
public:
    std::vector<Tile> tiles;
    std::string name;
    size_t core_count;
    int noc_width, noc_height, noc_buffer_size, max_cores_per_tile;

    Architecture(const std::string &name, const NetworkOnChipConfiguration &noc);
    std::vector<std::reference_wrapper<Core>> cores();
    Tile &create_tile(const std::string &name, const TilePowerMetrics &power_metrics);
    Core &create_core(const std::string &name, const size_t parent_tile_id, const CorePipelineConfiguration &pipeline_config);
    std::string info();
    std::string description() const;
    void save_arch_description(const std::filesystem::path &path);
};

Architecture load_arch_description(const std::filesystem::path &filename);

struct Message
{
    double generation_delay, network_delay, receive_delay, blocked_delay;
    double sent_timestamp, received_timestamp, processed_timestamp;
    long int timestep;
    int spikes, hops;
    int src_neuron_id, src_neuron_group_id;
    int src_x, dest_x, src_y, dest_y;
    int src_tile_id, src_core_id, src_core_offset;
    int dest_tile_id, dest_core_id, dest_core_offset;
    int dest_axon_hw, dest_axon_id;
    bool placeholder, in_noc;

    explicit Message(const Architecture &arch, const Neuron &n, const int timestep);
    explicit Message(const Architecture &arch, const Neuron &n, const int timestep, const int axon_address);
};

struct TilePowerMetrics
{
    double energy_north_hop, latency_north_hop;
    double energy_east_hop, latency_east_hop;
    double energy_south_hop, latency_south_hop;
    double energy_west_hop, latency_west_hop;

    TilePowerMetrics(const double energy_north = 0.0, const double latency_north = 0.0, const double energy_east = 0.0, const double latency_east = 0.0, const double energy_south = 0.0, const double latency_south = 0.0, const double energy_west = 0.0, const double latency_west = 0.0);
};

struct SynapsePowerMetrics
{
    double energy_memory_access, latency_memory_access;
    double energy_spike_op, latency_spike_op;

    SynapsePowerMetrics(const double energy_memory = 0.0, const double latency_memory = 0.0, const double energy_spike = 0.0, const double latency_spike = 0.0);
};

struct SomaPowerMetrics
{
    double energy_update_neuron, latency_update_neuron;
    double energy_access_neuron, latency_access_neuron;
    double energy_spiking, latency_spiking;

    SomaPowerMetrics(const double energy_update = 0.0, const double latency_update = 0.0, const double energy_access = 0.0, const double latency_access = 0.0, const double energy_spiking = 0.0, const double latency_spiking = 0.0);
};

struct CorePipelineConfiguration
{
    BufferPosition timestep_buffer_pos;
    size_t max_neurons_supported;

    CorePipelineConfiguration(const std::string &buffer_pos="soma", const size_t max_neurons_supported=1024);
};

struct NetworkOnChipConfiguration
{
    int width_in_tiles, height_in_tiles, buffer_size;

    NetworkOnChipConfiguration(const int width = 1, const int height = 1, const int buffer_size = 0);
};

class Tile
{
public:
    // TODO: this may need to be a vector of core pointers, to avoid reallocation messing up the python
    //  version of the objects
    std::vector<Core> cores;
    std::string name;
    double energy;
    double energy_north_hop, latency_north_hop;
    double energy_east_hop, latency_east_hop;
    double energy_south_hop, latency_south_hop;
    double energy_west_hop, latency_west_hop;
    long int hops, messages_received, total_neurons_fired;
    long int north_hops, east_hops, south_hops, west_hops;
    size_t id, x, y;

    explicit Tile(const std::string &name, const size_t tile_id, const TilePowerMetrics &power_metrics);
    int get_id() { return id; }
    std::string info() const;
    std::string description() const;
};

struct CoreAddress
{
    size_t id, parent_tile_id, offset_within_tile;
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
    double energy, next_message_generation_delay;
    size_t id, offset, parent_tile_id;
    int message_count;

    explicit Core(const std::string &name, const CoreAddress &address, const CorePipelineConfiguration &pipeline);
    AxonInUnit &create_axon_in(const std::string &name, const double energy_message, const double latency_message);
    SynapseUnit &create_synapse(const std::string &name, const std::string &model_str, const SynapsePowerMetrics &power_metrics);
    DendriteUnit &create_dendrite(const std::string &name, const std::string &model_str, const double energy_access, const double latency_access);
    SomaUnit &create_soma(const std::string &name, const std::string &model_str, const SomaPowerMetrics &power_metrics);
    AxonOutUnit &create_axon_out(const std::string &name, const double energy_access, const double latency_access);
    void map_neuron(Neuron &n);
    int get_id() { return id; }
    int get_offset() { return offset; }
    std::string info() const;
    std::string description() const;
};

struct AxonInUnit
{
    std::string name;
    CoreAddress parent_core_address;
    long int spike_messages_in;
    double energy, time;
    double energy_spike_message, latency_spike_message;

    explicit AxonInUnit(const std::string &axon_in_name, const CoreAddress &parent_core_address, const double energy_message=0.0, const double latency_message=0.0);
    std::string description() const;
};

struct SynapseUnit
{
    std::filesystem::path plugin_lib;
    std::string name, model;
    CoreAddress parent_core_address;
    long int spikes_processed;
    double energy, time;
    double energy_memory_access, latency_memory_access;
    double energy_spike_op, latency_spike_op;

    explicit SynapseUnit(const std::string &synapse_name, const std::string &model_str, const CoreAddress &parent_core, const SynapsePowerMetrics &power_metrics);
    std::string description() const;
};

struct DendriteUnit
{
    std::filesystem::path plugin_lib;
    std::string name, model;
    CoreAddress parent_core_address;
    double energy, time;
    double energy_access, latency_access;
    explicit DendriteUnit(const std::string &dendrite_name, const std::string &model_str, const CoreAddress &parent_core, const double energy_cost=0.0, const double latency_cost=0.0);
    std::string description() const;
};

struct SomaUnit
{
    FILE *noise_stream;
    std::filesystem::path plugin_lib;
    std::string name, model;
    CoreAddress parent_core_address;
    long int neuron_updates, neurons_fired, neuron_count;
    double energy, time;
    double energy_update_neuron, latency_update_neuron;
    double energy_access_neuron, latency_access_neuron;
    double energy_spiking, latency_spiking;
    int noise_type;

    explicit SomaUnit(const std::string &soma_name, const std::string &model_str, const CoreAddress &parent_core, const SomaPowerMetrics &power_metrics);
    std::string description() const;
};

struct AxonOutUnit
{
    // The axon output points to a number of axons, stored at the
    //  post-synaptic core. A neuron can point to a number of these
    std::string name;
    CoreAddress parent_core_address;
    long int packets_out;
    double energy, time;
    double energy_access, latency_access;

    explicit AxonOutUnit(const std::string &axon_out_name, const CoreAddress &parent_core, const double energy_access, const double latency_access);
    std::string description() const;
};

struct AxonInModel
{
    // List of all neuron connections to send spikes to
    std::vector<int> synapse_addresses;

    Message *message;
    long int last_updated;
    int spikes_received, active_synapses;
};

struct AxonOutModel
{
    // List of all neuron connections to send spike to
    int dest_axon_id, dest_tile_id, dest_core_offset;
    int src_neuron_id;
};

enum NoiseType
{
    NOISE_NONE = -1,
    NOISE_FILE_STREAM,
    // TODO: implement different random noise generators
};

void arch_create_axons(Architecture &arch);
void arch_print_axon_summary(Architecture &arch);
void arch_map_neuron_connections(Neuron &n);
void arch_allocate_axon(Neuron &pre_neuron, Core &post_core);
void arch_add_connection_to_axon(Connection &con, Core &post_core);
int arch_parse_synapse_model(const std::string &model_str);
}

#endif
