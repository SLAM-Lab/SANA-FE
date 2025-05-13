// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// chip.hpp - Neuromorphic simulator kernel
//
// Time-step based simulation, based on loop:
// 1) seed any input spikes
// 2) route spikes
// 3) update neurons and check firing
#ifndef SIM_HEADER_INCLUDED_
#define SIM_HEADER_INCLUDED_

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#include <atomic>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <list>
#include <memory>
#include <queue>
#include <variant>

struct BookSimConfig;

#include "arch.hpp"
#include "print.hpp"

namespace sanafe
{
struct Timestep;
struct RunData;
struct Scheduler;
struct NocInfo;

class Architecture;
struct TileConfiguration;
struct CoreConfiguration;
struct CorePipelineConfiguration;
struct AxonInConfiguration;
struct PipelineUnitConfiguration;
struct AxonOutConfiguration;

class SpikingNetwork;
struct Message;
class Neuron;
struct Connection;
struct MappedNeuron;
class MappedConnection;
struct Synapse;
struct NeuronConfiguration;
struct ModelInfo;
struct PipelineResult;

class Tile;
class Core;
class AxonInUnit;
class AxonOutUnit;
class PipelineUnit;

struct AxonInModel;
struct AxonOutModel;
enum BufferPosition : int;
enum NeuronStatus : int
{
    INVALID_NEURON_STATE,
    IDLE,
    UPDATED,
    FIRED
};

enum TimingModel : int
{
    TIMING_MODEL_SIMPLE, // analytical model
    TIMING_MODEL_DETAILED, // semi-analytical model
    TIMING_MODEL_CYCLE_ACCURATE, // Booksim2 simulator
};

constexpr long int default_heartbeat_timesteps = 100L;
class SpikingChip
{
public:
    std::vector<Tile> tiles{};
    // Keep a reference to the different neuron groups mapped to the H/W
    std::map<std::string, std::vector<MappedNeuron *>> mapped_neuron_groups{};

    SpikingChip(const Architecture &arch, const std::filesystem::path &output_dir = ".", bool record_spikes = false, bool record_potentials = false, bool record_perf = false, bool record_messages = false);
    ~SpikingChip();
    // Do not allow copying
    SpikingChip(const SpikingChip &copy) = delete;
    SpikingChip(SpikingChip &&other) = delete;
    SpikingChip &operator=(const SpikingChip &copy) = delete;
    SpikingChip &operator=(SpikingChip &&other) = delete;
    RunData sim(long int timesteps = 1, long int heartbeat = default_heartbeat_timesteps, const TimingModel timing_model = TIMING_MODEL_DETAILED);
    void load(const SpikingNetwork &net);
    double get_power() const;
    RunData get_run_summary() const;
    void sim_output_run_summary(const std::filesystem::path &output_dir) const;
    void reset();

    std::vector<std::reference_wrapper<Core>> cores();

    size_t core_count{0UL};
    int noc_width;
    int noc_height;
    int noc_buffer_size;
    int max_cores_per_tile{0};

private:
    std::unique_ptr<BookSimConfig> booksim_config{};
    std::string out_dir;
    size_t total_neurons_mapped{0UL};
    long int total_neurons_fired{0L};
    long int total_timesteps{0L};
    long int total_spikes{0L};
    long int total_messages_sent{0L};
    double wall_time{0.0};
    double total_sim_time{0.0};
    double total_energy{0.0};
    double synapse_energy{0.0};
    double dendrite_energy{0.0};
    double soma_energy{0.0};
    double network_energy{0.0};

    // Performance and other misc tracking
    static std::atomic<int> chip_count;
    double neuron_processing_wall{0.0};
    double message_processing_wall{0.0};
    double scheduler_wall{0.0};
    double other_stats_wall{0.0};

    // Flags and filestreams
    bool spike_trace_enabled{false};
    bool potential_trace_enabled{false};
    bool perf_trace_enabled{false};
    bool message_trace_enabled{false};
    std::ofstream spike_trace{};
    std::ofstream potential_trace{};
    std::ofstream message_trace{};
    std::ofstream perf_trace{};

    Timestep step(TimingModel timing_model=TIMING_MODEL_DETAILED);
    void map_neurons(const SpikingNetwork &net);
    void map_connections(const SpikingNetwork &net);
    MappedConnection &map_connection(const Connection &con);
    void map_axons();
    void sim_reset_measurements();
    void sim_timestep(Timestep &ts, const TimingModel timing_model = TIMING_MODEL_DETAILED);
    double sim_estimate_network_costs(const Tile &src, Tile &dest);
    void sim_calculate_energy(Timestep &ts);

    void sim_format_run_summary(std::ostream &out, const RunData &run_data) const;

    void sim_print_axon_summary();
    void sim_create_neuron_axons(MappedNeuron &pre_neuron);
    void sim_allocate_axon(MappedNeuron &pre_neuron, Core &post_core);
    void sim_add_connection_to_axon(MappedConnection &con, Core &post_core);

    void process_neurons(Timestep &ts);
    void process_messages(Timestep &ts);
    void forced_updates(const Timestep &ts);


    void process_neuron(Timestep &ts, MappedNeuron &n);
    void receive_message(Message &m);
    double process_message(Timestep &ts, Core &c, Message &m);
    PipelineResult execute_pipeline(const std::vector<PipelineUnit *> &pipeline, Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);

    double pipeline_process_axon_in(Core &core, const Message &m);
    PipelineResult pipeline_process_axon_out(Timestep &ts, MappedNeuron &n);

    std::ofstream sim_trace_open_perf_trace(const std::filesystem::path &out_dir);
    std::ofstream sim_trace_open_spike_trace(const std::filesystem::path &out_dir);
    std::ofstream sim_trace_open_potential_trace(const std::filesystem::path &out_dir);
    std::ofstream sim_trace_open_message_trace(const std::filesystem::path &out_dir);
    void sim_trace_write_spike_header(std::ofstream &spike_trace_file);
    void sim_trace_write_potential_header(std::ofstream &potential_trace_file);
    void sim_trace_write_perf_header(std::ofstream &perf_trace_file);
    void sim_trace_write_message_header(std::ofstream &message_trace_file);
    void sim_trace_record_spikes(std::ofstream &spike_trace_file, long int timesteps);
    void sim_trace_record_potentials(std::ofstream &potential_trace_file, long int timestep);
    void sim_trace_record_message(std::ofstream &message_trace_file, const Message &m);
    void sim_trace_record_perf(std::ofstream &out, const Timestep &ts);
    std::map<std::string, double> sim_trace_get_optional_traces();
};

struct RunData
{
    double wall_time{0.0};
    long int spikes{0L};
    long int packets_sent{0L};
    long int neurons_fired{0L};
    long int timestep_start;
    long int timesteps_executed;

    double total_energy{0.0};
    double synapse_energy{0.0};
    double dendrite_energy{0.0};
    double soma_energy{0.0};
    double network_energy{0.0};
    double sim_time{0.0};

    RunData(long int start,  long int steps);
};


struct Timestep
{
    std::vector<std::list<Message>> messages;
    long int timestep;
    long int spike_count{0L};
    long int total_hops{0L};
    long int packets_sent{0L};
    long int neurons_fired{0L};

    double total_energy{0.0};
    double synapse_energy{0.0};
    double dendrite_energy{0.0};
    double soma_energy{0.0};
    double network_energy{0.0};
    double sim_time{0.0};

    Timestep(long int ts, int core_count);
};

class MappedConnection
{
public:
    MappedNeuron *post_neuron{nullptr};
    MappedNeuron *pre_neuron{nullptr};
    PipelineUnit *synapse_hw{nullptr};
    std::vector<PipelineUnit *> message_processing_pipeline{};
    size_t synapse_address{0UL};

    explicit MappedConnection();
    void build_message_processing_pipeline();
};

class MappedNeuron
{
public:
    std::vector<MappedConnection> connections_out;
    std::vector<int> axon_out_addresses;
    std::string parent_group_name;
    size_t offset;
    size_t id;

    // Internal pointers to mapped hardware
    Core *core{nullptr};
    Core *post_synaptic_cores{nullptr};
    PipelineUnit *dendrite_hw{nullptr};
    PipelineUnit *soma_hw{nullptr};
    AxonOutUnit *axon_out_hw{nullptr};
    std::vector<PipelineUnit *> neuron_processing_pipeline{};

    size_t mapped_address{-1ULL};
    size_t mapping_order;
    int spike_count{0};
    int maps_in_count{0};
    int maps_out_count{0};
    NeuronStatus status{INVALID_NEURON_STATE};

    // Flags and traces
    bool force_synapse_update{false};
    bool force_dendrite_update{false};
    bool force_soma_update{false};
    bool log_spikes{false};
    bool log_potential{false};

    // Track spikes
    bool axon_out_input_spike{false};

    MappedNeuron(const Neuron &neuron_to_map, Core *mapped_core, const size_t nid, const size_t mapped_address, PipelineUnit *mapped_dendrite, PipelineUnit *mapped_soma, AxonOutUnit *mapped_axon_out);
    MappedNeuron(const MappedNeuron &copy) = default;
    MappedNeuron& operator=(const MappedNeuron& other) = default;
    MappedNeuron(MappedNeuron&& other) = default;
    MappedNeuron& operator=(MappedNeuron&& other) = default;
    void set_model_attributes(const std::map<std::string, ModelParam> &model_parameters);

private:
    void build_neuron_processing_pipeline();
};

struct Synapse
{
    double current;
    MappedConnection *con;
};

constexpr long int placeholder_mid = -1UL; // An invalid message id for placeholders
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
    const long int mid;
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

    explicit Message(const long int id, const SpikingChip &hw, const MappedNeuron &n, long int timestep);
    explicit Message(const long int id, const SpikingChip &hw, const MappedNeuron &n, long int timestep, int axon_address);
};

struct PipelineResult
{
    // Hardware outputs
    std::optional<double> current{std::nullopt};
    NeuronStatus status{INVALID_NEURON_STATE};
    // Optionally simulate energy and/or latency
    std::optional<double> energy{std::nullopt};
    std::optional<double> latency{std::nullopt};
};

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

class AxonInUnit
{
public:
    std::string name;
    long int spike_messages_in{0L};
    double energy{0.0};
    double time{0.0};
    double energy_spike_message;
    double latency_spike_message;

    explicit AxonInUnit(const AxonInConfiguration &config);
};

class PipelineUnit
{
public:
    PipelineUnit(const PipelineUnit &copy) = default;
    PipelineUnit(PipelineUnit &&other) = default;
    virtual ~PipelineUnit() = default;
    PipelineUnit &operator=(const PipelineUnit &other) = default;
    PipelineUnit &operator=(PipelineUnit &&other) = default;

    // Virtual member functions
    virtual void set_attribute_hw(const std::string &param_name, const ModelParam &param) = 0;
    virtual void set_attribute_neuron(size_t neuron_address, const std::string &param_name, const ModelParam &param) = 0;
    virtual void set_attribute_edge(size_t synapse_address, const std::string &param_name, const ModelParam &param) = 0;
    virtual void reset() = 0;

    // Optional virtual functions
    // The user of this class must implement the interfaces they wish to support
    //  Depending on whether you want to support Synapse, Dendrite, Soma
    //  operations, or a combination of the three.
    // If using synaptic inputs (address and read vs. synapse update without read)
    virtual PipelineResult update(size_t synapse_address, bool read = false) { throw std::logic_error("Error: Synapse input not implemented"); }
    // If using dendritic inputs (neuron address, synaptic current and synaptic address for additional info)
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) { throw std::logic_error("Error: Dendrite input not implemented"); }
    // If using somatic inputs (address and current in)
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in) { throw std::logic_error("Error: Soma input not implemented"); }
    virtual double get_potential(size_t neuron_address) { return 0.0; }
    virtual void map_connection(MappedConnection &con) {}

    // Normal member functions and function pointers
    void set_time(const long int timestep) { simulation_time = timestep; }
    void configure(std::string unit_name, const ModelInfo &model);
    PipelineResult process(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);

    using InputInterfaceFunc = PipelineResult (PipelineUnit:: *)(Timestep &, MappedNeuron &, std::optional<MappedConnection*>, const PipelineResult &);
    using OutputInterfaceFunc = PipelineResult (PipelineUnit:: *)(MappedNeuron &, std::optional<MappedConnection *>, const PipelineResult &);
    InputInterfaceFunc process_input_fn{nullptr};
    OutputInterfaceFunc process_output_fn{nullptr};

    // Model information
    std::map<std::string, ModelParam> model_parameters{};
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;

    // Performance metrics
    std::optional<double> default_energy_process_spike{std::nullopt};
    std::optional<double> default_latency_process_spike{std::nullopt};
    std::optional<double> default_energy_update{std::nullopt};
    std::optional<double> default_latency_update{std::nullopt};
    std::optional<SomaEnergyMetrics> default_soma_energy_metrics;
    std::optional<SomaLatencyMetrics> default_soma_latency_metrics;
    double energy{0.0};
    double time{0.0};

    // Performance counters
    long int spikes_processed{0L};
    long int neuron_updates{0L};
    long int neurons_fired{0L};
    long int neuron_count{0L};

    // Implementation flags, set whichever operations your derived unit supports
    //  to 'true'. Note that a hardware unit must support one or more of these
    const bool implements_synapse;
    const bool implements_dendrite;
    const bool implements_soma;

    // Performance monitoring flags
    bool log_energy{false};
    bool log_latency{false};

protected:
    long int simulation_time{0L};
    PipelineUnit(const bool implements_synapse, const bool implements_dendrite, const bool implements_soma);

    PipelineResult process_synapse_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_dendrite_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_soma_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);

    PipelineResult process_synapse_output(MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &output);
    PipelineResult process_dendrite_output(MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &output);
    PipelineResult process_soma_output(MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &output);

private:
    PipelineResult calculate_synapse_default_energy_latency(MappedConnection &con, const PipelineResult &simulation_result);
    PipelineResult calculate_dendrite_default_energy_latency(MappedNeuron &n, const PipelineResult &simulation_result);
    PipelineResult calculate_soma_default_energy_latency(MappedNeuron &n, const PipelineResult &simulation_result);
    void update_soma_activity(MappedNeuron &n, const PipelineResult &simulation_result);
    void check_outputs(const MappedNeuron &n, const PipelineResult &result);
};

// Specific unit base classes, for the normal use case where the model
//  implements only synaptic, dendritic or somatic functionality (and not a)
//  combination
class SynapseUnit : public PipelineUnit
{
public:
    SynapseUnit() : PipelineUnit(true, false, false) {};
    virtual PipelineResult update(size_t synapse_address, bool read = false) = 0;
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) override final { throw std::logic_error("Error: Synapse H/W called with dendrite inputs"); }
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in) override final { throw std::logic_error("Error: Synapse H/W called with soma inputs"); }
    virtual void set_attribute_neuron(size_t neuron_address,  const std::string &param_name, const ModelParam &param) override final {};
};

class DendriteUnit : public PipelineUnit
{
public:
    DendriteUnit() : PipelineUnit(false, true, false) {};
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) override = 0;
    virtual PipelineResult update(size_t synapse_address, bool read = false) override final { throw std::logic_error("Error: Dendrite H/W called with synapse inputs"); }
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in) override final { throw std::logic_error("Error: Dendrite H/W called with soma inputs"); }
};

class SomaUnit : public PipelineUnit
{
public:
    SomaUnit() : PipelineUnit(false, false, true) {};
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in) = 0;
    virtual PipelineResult update(size_t synapse_address, bool read = false) override final { throw std::logic_error("Error: Soma H/W called with synapse inputs"); }
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) override final { throw std::logic_error("Error: Soma H/W called with dendrite inputs"); }
    virtual void set_attribute_edge(size_t synapse_address, const std::string &param_name, const ModelParam &param) override final {};
    virtual void map_connection(MappedConnection &con) override final {};
};

class AxonOutUnit
{
public:
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
    std::vector<std::shared_ptr<PipelineUnit>> pipeline_hw;
    std::vector<AxonOutUnit> axon_out_hw;

    std::vector<Message *> messages_in;
    std::vector<AxonInModel> axons_in;
    std::vector<MappedNeuron> neurons;
    std::vector<MappedConnection *> connections_in;
    std::vector<AxonOutModel> axons_out;
    std::vector<PipelineResult> timestep_buffer{};

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
    bool log_energy{false};
    bool log_latency{false};

    explicit Core(const CoreConfiguration &config);
    void map_neuron(const Neuron &n, const size_t neuron_id);
    AxonInUnit &create_axon_in(const AxonInConfiguration &config);
    PipelineUnit &create_pipeline_unit(const PipelineUnitConfiguration &config);
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
    bool log_energy{false};
    bool log_latency{false};

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

double calculate_elapsed_time(const std::chrono::time_point<std::chrono::high_resolution_clock> &ts_start, const std::chrono::time_point<std::chrono::high_resolution_clock> &ts_end);
size_t abs_diff(size_t a, size_t b);
BufferPosition pipeline_parse_buffer_pos_str(const std::string &buffer_pos_str, const bool buffer_inside_unit);
}

#endif
