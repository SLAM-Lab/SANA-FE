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

#define RAND_SEED 0xbeef // For srand()
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
class Network;
struct Message;
class Neuron;

class Tile;
class Core;

//const long int default_heartbeat_timesteps = 100L;
// TODO: revert
const long int default_heartbeat_timesteps = 1L;
class Simulation
{
public:
    Simulation(Architecture &arch, Network &net, const std::filesystem::path &output_dir = ".", bool record_spikes = false, bool record_potentials = false, bool record_perf = false, bool record_messages = false);
    ~Simulation();
    // Do not allow copying of Simulation object
    Simulation(const Simulation &copy) = delete;
    Simulation(Simulation &&other) = delete;
    Simulation &operator=(const Simulation &copy) = delete;
    Simulation &operator=(Simulation &&other) = delete;
    RunData run(long int timesteps = 1, long int heartbeat = default_heartbeat_timesteps);
    //int update_neuron(std::vector<NeuronGroup>::size_type group_id, std::vector<Neuron>::size_type n_id, std::vector<std::string> kwargs, int count);
    double get_power() const;
    RunData get_run_summary() const;


private:
    Architecture &arch;
    Network &net;
    std::string out_dir;
    long int total_neurons_fired, total_timesteps, total_spikes;
    long int total_messages_sent;
    double total_energy, total_sim_time, wall_time;
    bool spike_trace_enabled, potential_trace_enabled;
    bool perf_trace_enabled, message_trace_enabled;
    std::ofstream spike_trace, potential_trace, message_trace, perf_trace;

    Timestep step();
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

    Timestep(long int ts,  int core_count);
};

enum ProgramArgs
{
    ARCH_FILENAME = 0,
    NETWORK_FILENAME,
    TIMESTEPS,
    PROGRAM_NARGS,
};

void sim_timestep(Timestep &ts, Architecture &arch, Network &net);
double sim_estimate_network_costs(const Tile &src, Tile &dest);
void sim_reset_measurements(Network &net, Architecture &arch);
double sim_calculate_energy(const Architecture &arch);

std::ofstream sim_trace_open_perf_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_spike_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_potential_trace(const std::filesystem::path &out_dir, const Network &net);
std::ofstream sim_trace_open_message_trace(const std::filesystem::path &out_dir);
void sim_trace_write_spike_header(std::ofstream &spike_trace_file);
void sim_trace_write_potential_header(std::ofstream &potential_trace_file, const Network &net);
void sim_trace_write_perf_header(std::ofstream &perf_trace_file);
void sim_trace_write_message_header(std::ofstream &message_trace_file);
void sim_trace_record_spikes(std::ofstream &spike_trace_file, long int timesteps, const Network &net);
void sim_trace_record_potentials(std::ofstream &potential_trace_file, long int timestep, const Network &net);
void sim_trace_record_message(std::ofstream &message_trace_file, const Message &m);
void sim_trace_perf_log_timestep(std::ofstream &out, const Timestep &ts);

void sim_output_run_summary(const std::filesystem::path &output_dir, const RunData &run_data);
void sim_format_run_summary(std::ostream &out, const RunData &run_data);

//double sim_generate_noise(Neuron *n);
timespec calculate_elapsed_time(const timespec &ts_start, const timespec &ts_end);

//int sim_poisson_input( double firing_probability);
//int sim_rate_input( double firing_rate, double *spike_val);
}

#endif
