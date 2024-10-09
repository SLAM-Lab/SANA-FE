#ifndef HARDWARE_HEADER_INCLUDED_
#define HARDWARE_HEADER_INCLUDED_

#include <cassert>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include "arch.hpp"

namespace sanafe
{
// Forward declarations
struct AxonInModel;
struct AxonOutModel;
struct MappedNeuron;
struct TilePowerMetrics;
struct AxonInPowerMetrics;

struct MappedConnection
{
    std::map<std::string, ModelParam> dendrite_params{};
    MappedNeuron *post_neuron{nullptr};
    MappedNeuron *pre_neuron{nullptr};
    SynapseUnit *synapse_hw{nullptr};
    size_t synapse_address{0UL};
    int id;

    explicit MappedConnection(int connection_id);
};

struct MappedNeuron
{
    std::vector<MappedConnection> connections_out;
    std::vector<int> axon_out_addresses;
    std::string parent_group_name;
    size_t id;

    // Internal pointers to mapped hardware
    Core *core{nullptr};
    Core *post_synaptic_cores{nullptr};
    DendriteUnit *dendrite_hw{nullptr};
    SomaUnit *soma_hw{nullptr};
    AxonOutUnit *axon_out_hw{nullptr};

    size_t mapped_address{};
    size_t mapping_order;
    int spike_count{0};
    int maps_in_count{0};
    int maps_out_count{0};

    // Flags and traces
    bool force_synapse_update{false};
    bool force_dendrite_update{false};
    bool force_soma_update{false};
    bool log_spikes{false};
    bool log_potential{false};

    // Inputs to H/W units
    NeuronStatus status{sanafe::IDLE};
    std::vector<Synapse> dendrite_input_synapses{};
    double soma_input_charge{0.0};
    bool axon_out_input_spike{false};

    void set_attributes(const NeuronTemplate &attributes);
    MappedNeuron(size_t id, size_t mapping_order) : id(id), mapping_order{mapping_order} {}
    // TODO: allow user to manually set mapping order
};

struct Synapse
{
    double current;
    MappedConnection &con;
};

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

    explicit Message(const SpikingHardware &hw, const MappedNeuron &n, long int timestep);
    explicit Message(const SpikingHardware &hw, const MappedNeuron &n, long int timestep, int axon_address);
};

struct AxonInUnit
{
    std::string name;
    long int spike_messages_in{0L};
    double energy{0.0};
    double time{0.0};
    double energy_spike_message;
    double latency_spike_message;

    explicit AxonInUnit(const AxonInConfiguration &config);
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

    explicit AxonOutUnit(const AxonOutConfiguration &config);
    //std::string description() const;
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
    std::vector<MappedNeuron> neurons;
    std::vector<MappedConnection *> connections_in;
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

    explicit Core(const CoreConfiguration &config);
    MappedNeuron &map_neuron(const Neuron &n);
    AxonInUnit &create_axon_in(const AxonInConfiguration &config);
    SynapseUnit &create_synapse(const SynapseConfiguration &config);
    DendriteUnit &create_dendrite(const DendriteConfiguration &config);
    SomaUnit &create_soma(const SomaConfiguration &config);
    AxonOutUnit &create_axon_out(const AxonOutConfiguration &config);
    [[nodiscard]] int get_id() const { return id; }
    [[nodiscard]] int get_offset() const { return offset; }
    [[nodiscard]] std::string info() const;
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

    explicit Tile(const TileConfiguration &config);
    [[nodiscard]] int get_id() const { return id; }
    [[nodiscard]] std::string info() const;
};

enum ProgramArgs
{
    ARCH_FILENAME = 0,
    NETWORK_FILENAME,
    TIMESTEPS,
    PROGRAM_NARGS,
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
}

#endif
