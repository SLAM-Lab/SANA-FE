// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// sim.hpp - Neuromorphic simulator kernel
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

class SpikingNetwork;
struct Message;
class Neuron;
struct MappedNeuron;
struct MappedConnection;

class Tile;
class Core;
struct AxonInUnit;
struct AxonOutUnit;

struct AxonInModel;
struct AxonOutModel;
enum BufferPosition : int;

const long int default_heartbeat_timesteps = 100L;
class SpikingHardware
{
public:
    std::vector<Tile> tiles{};
    // Keep a reference to the different neuron groups mapped to the H/W
    std::map<std::string, std::vector<MappedNeuron *>> mapped_neuron_groups{};

    SpikingHardware(const Architecture &arch, const std::filesystem::path &output_dir = ".", bool record_spikes = false, bool record_potentials = false, bool record_perf = false, bool record_messages = false);
    ~SpikingHardware();
    // Do not allow copying
    SpikingHardware(const SpikingHardware &copy) = delete;
    SpikingHardware(SpikingHardware &&other) = delete;
    SpikingHardware &operator=(const SpikingHardware &copy) = delete;
    SpikingHardware &operator=(SpikingHardware &&other) = delete;
    RunData sim(long int timesteps = 1, long int heartbeat = default_heartbeat_timesteps);
    void load(const SpikingNetwork &net);
    double get_power() const;
    RunData get_run_summary() const;

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

    Timestep step();
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

void sim_timestep(Timestep &ts, SpikingHardware &hw);
double sim_estimate_network_costs(const Tile &src, Tile &dest);
void sim_reset_measurements(SpikingHardware &hw);
double sim_calculate_energy(const SpikingHardware &hw);

std::ofstream sim_trace_open_perf_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_spike_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_potential_trace(const std::filesystem::path &out_dir, const SpikingHardware &hw);
std::ofstream sim_trace_open_message_trace(const std::filesystem::path &out_dir);
void sim_trace_write_spike_header(std::ofstream &spike_trace_file);
void sim_trace_write_potential_header(std::ofstream &potential_trace_file, const SpikingHardware &hw);
void sim_trace_write_perf_header(std::ofstream &perf_trace_file);
void sim_trace_write_message_header(std::ofstream &message_trace_file);
void sim_trace_record_spikes(std::ofstream &spike_trace_file, long int timesteps, const SpikingHardware &hw);
void sim_trace_record_potentials(std::ofstream &potential_trace_file, long int timestep, const SpikingHardware &hw);
void sim_trace_record_message(std::ofstream &message_trace_file, const Message &m);
void sim_trace_perf_log_timestep(std::ofstream &out, const Timestep &ts);

void sim_output_run_summary(const std::filesystem::path &output_dir, const RunData &run_data);
void sim_format_run_summary(std::ostream &out, const RunData &run_data);

void sim_print_axon_summary(SpikingHardware &hw);
void sim_create_neuron_axons(MappedNeuron &pre_neuron);
void sim_allocate_axon(MappedNeuron &pre_neuron, Core &post_core);
void sim_add_connection_to_axon(MappedConnection &con, Core &post_core);

timespec calculate_elapsed_time(const timespec &ts_start, const timespec &ts_end);
}

#endif
