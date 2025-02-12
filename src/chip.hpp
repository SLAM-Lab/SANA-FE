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

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <list>
#include <memory>
#include <queue>
#include <variant>

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
struct MappedConnection;
struct Synapse;
struct NeuronConfiguration;
struct ModelInfo;

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
    TIMING_MODEL_SIMPLE,
    TIMING_MODEL_DETAILED
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
    void reset();

    std::vector<std::reference_wrapper<Core>> cores();

    size_t core_count{0UL};
    int noc_width;
    int noc_height;
    int noc_buffer_size;
    int max_cores_per_tile{0};

private:
    std::string out_dir;
    long int total_neurons_fired{0L};
    long int total_timesteps{0L};
    long int total_spikes{0L};
    long int total_messages_sent{0L};
    double total_energy{0.0};
    double total_sim_time{0.0};
    double wall_time{0.0};
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
};

struct RunData
{
    double energy{0.0};
    double sim_time{0.0};
    double wall_time{0.0};
    long int spikes{0L};
    long int packets_sent{0L};
    long int neurons_fired{0L};
    long int timestep_start;
    long int timesteps_executed;

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
    double energy{0.0};
    double sim_time{0.0};

    Timestep(long int ts, int core_count);
};

struct MappedConnection
{
    std::map<std::string, ModelParam> dendrite_params{};
    MappedNeuron *post_neuron{nullptr};
    MappedNeuron *pre_neuron{nullptr};
    PipelineUnit *synapse_hw{nullptr};
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
    PipelineUnit *dendrite_hw{nullptr};
    PipelineUnit *soma_hw{nullptr};
    AxonOutUnit *axon_out_hw{nullptr};

    size_t mapped_address{-1ULL};
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

    void configure_models(const std::map<std::string, ModelParam> &model_parameters);
    MappedNeuron(const Neuron &neuron_to_map, Core *mapped_core, const size_t address, PipelineUnit *mapped_dendrite, PipelineUnit *mapped_soma, AxonOutUnit *mapped_axon_out);
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

    explicit Message(const SpikingChip &hw, const MappedNeuron &n, long int timestep);
    explicit Message(const SpikingChip &hw, const MappedNeuron &n, long int timestep, int axon_address);
};

struct PipelineResult
{
    std::variant<double, NeuronStatus> pipeline_output;
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
    virtual void set_attribute(size_t address, const std::string &param_name, const ModelParam &param) = 0;
    virtual void reset() = 0;

    // The user of this class must implement the interface they wish to support
    //  Depending on whether you want to support Synapse, Dendrite, Soma or a
    //  combination of the three types in a PipelineUnit implementation
    //
    // If using synaptic inputs (address and read/update without read)
    virtual PipelineResult update(size_t synapse_address, bool read = false) { throw std::logic_error("Error: Synapse input not implemented"); }
    // If using dendritic inputs (address and synaptic information)
    virtual PipelineResult update(size_t neuron_address, std::optional<Synapse> synapse_in) { throw std::logic_error("Error: Dendrite input not implemented"); }
    // If using somatic inputs (address and current in)
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in) { throw std::logic_error("Error: Soma input not implemented"); }

    // Optional virtual functions that may be useful
    virtual double get_potential(size_t neuron_address) { return 0.0; }

    // Normal member functions
    void set_time(const long int timestep) { simulation_time = timestep; }
    void add_connection(MappedConnection &con);
    void configure(std::string unit_name, const ModelInfo &model);
    double process_connection(const Timestep &ts, MappedConnection &con);
    double process_neuron(const Timestep &ts, MappedNeuron &n);

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

    // Implementation flags, set whichever to your derived unit supports 'true'
    //  Note that a unit can support one or more of these
    bool implements_synapse{false};
    bool implements_dendrite{false};
    bool implements_soma{false};

protected:
    long int simulation_time{0L};
    std::vector<MappedConnection *> mapped_connections_in{};
    PipelineUnit() = default;
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
    void map_neuron(const Neuron &n);
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

void sim_timestep(Timestep &ts, SpikingChip &hw, const TimingModel timing_model = TIMING_MODEL_DETAILED);
double sim_estimate_network_costs(const Tile &src, Tile &dest);
void sim_reset_measurements(SpikingChip &hw);
double sim_calculate_energy(const SpikingChip &hw);

std::ofstream sim_trace_open_perf_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_spike_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_potential_trace(const std::filesystem::path &out_dir, const SpikingChip &hw);
std::ofstream sim_trace_open_message_trace(const std::filesystem::path &out_dir);
void sim_trace_write_spike_header(std::ofstream &spike_trace_file);
void sim_trace_write_potential_header(std::ofstream &potential_trace_file, const SpikingChip &hw);
void sim_trace_write_perf_header(std::ofstream &perf_trace_file);
void sim_trace_write_message_header(std::ofstream &message_trace_file);
void sim_trace_record_spikes(std::ofstream &spike_trace_file, long int timesteps, const SpikingChip &hw);
void sim_trace_record_potentials(std::ofstream &potential_trace_file, long int timestep, const SpikingChip &hw);
void sim_trace_record_message(std::ofstream &message_trace_file, const Message &m);
void sim_trace_perf_log_timestep(std::ofstream &out, const Timestep &ts);

void sim_output_run_summary(const std::filesystem::path &output_dir, const RunData &run_data);
void sim_format_run_summary(std::ostream &out, const RunData &run_data);

void sim_print_axon_summary(SpikingChip &hw);
void sim_create_neuron_axons(MappedNeuron &pre_neuron);
void sim_allocate_axon(MappedNeuron &pre_neuron, Core &post_core);
void sim_add_connection_to_axon(MappedConnection &con, Core &post_core);

void pipeline_process_neurons(Timestep &ts, SpikingChip &hw);
void pipeline_process_messages(Timestep &ts, SpikingChip &hw);

void pipeline_process_neuron(Timestep &ts, const SpikingChip &arch, MappedNeuron &n);
void pipeline_receive_message(SpikingChip &arch, Message &m);
double pipeline_process_message(const Timestep &ts, Core &c, Message &m);

double pipeline_process_axon_in(Core &core, const Message &m);
double pipeline_process_synapse(const Timestep &ts, MappedConnection &con);
double pipeline_process_dendrite(const Timestep &ts, MappedNeuron &n);
double pipeline_process_soma(const Timestep &ts, MappedNeuron &n);
double pipeline_process_axon_out(Timestep &ts, const SpikingChip &arch, MappedNeuron &n);
BufferPosition pipeline_parse_buffer_pos_str(const std::string &buffer_pos_str);

std::pair<double, double> pipeline_apply_default_dendrite_power_model(MappedNeuron &n, std::optional<double> energy, std::optional<double> latency);

timespec calculate_elapsed_time(const timespec &ts_start, const timespec &ts_end);
size_t abs_diff(size_t a, size_t b);
}

#endif
