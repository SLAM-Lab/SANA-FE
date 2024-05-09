// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// sim.h - Neuromorphic simulator kernel
//
// Time-step based simulation, based on loop:
// 1) seed any input spikes
// 2) route spikes
// 3) update neurons and check firing
#ifndef SIM_HEADER_INCLUDED_
#define SIM_HEADER_INCLUDED_

#define RAND_SEED 0xbeef // For srand()
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#include <memory>
#include <list>
#include <filesystem>
#include "arch.hpp"
#include "network.hpp"
#include "stdio.h"

namespace sanafe
{
struct Timestep;
struct RunData;

class Simulation
{
public:
	Simulation(Architecture &arch, Network &net, const std::string &output_dir, const bool record_spikes, const bool record_potentials, const bool record_perf, const bool record_messages);
	~Simulation();
	RunData run(const long int timesteps=1, const long int heartbeat=100);
	//int update_neuron(std::vector<NeuronGroup>::size_type group_id, std::vector<Neuron>::size_type n_id, std::vector<std::string> kwargs, int count);
	double get_power();
	RunData get_run_summary();

private:
	Architecture &arch;
	Network &net;
	std::string out_dir;
	long int total_neurons_fired, total_timesteps, total_spikes;
	long int total_messages_sent;
	double total_energy, total_sim_time, wall_time;
	bool spike_trace_enabled, potential_trace_enabled;
	bool perf_trace_enabled, message_trace_enabled;
	FILE *stats_fp;
	std::ofstream spike_trace, potential_trace, message_trace, perf_trace;

	Timestep step();
	Simulation(const Simulation &copy);
};

struct RunData
{
	long int timestep_start, timesteps_executed;
	double energy, sim_time, wall_time;
	long int spikes, packets_sent, neurons_fired;

	RunData(const long int start, const long int steps);
};

struct Timestep
{
	std::vector<std::list<Message> > messages;
	std::vector<MessageFifo> message_queues;
	long int timestep, spike_count, total_hops, packets_sent;
	long int neurons_fired, spikes;
	double energy, sim_time;

	Timestep(const long int ts, const int core_count);
};

enum ProgramArgs
{
	ARCH_FILENAME = 0,
	NETWORK_FILENAME,
	TIMESTEPS,
	PROGRAM_NARGS,
};

struct NocInfo
{
	MessageFifo messages_received[ARCH_MAX_CORES];
	int noc_width, noc_height;
	double noc_messages_in[ARCH_MAX_X][ARCH_MAX_Y][4+ARCH_MAX_CORES_PER_TILE];
	double core_finished_receiving[ARCH_MAX_CORES];
	double mean_in_flight_receive_delay;
	long int messages_in_noc;
};

struct Scheduler
{
	int noc_width, noc_height, buffer_size;
};

enum Direction
{
	NORTH = 0,
	EAST,
	SOUTH,
	WEST
};

void print_run_data(FILE *fp, RunData &data);

std::unique_ptr<Simulation> sim_init_sim(void);
void sim_init_timestep(Timestep &ts, Architecture &arch);
void sim_timestep(Timestep &ts, Architecture &arch, Network &net);

void sim_process_neurons(Timestep &ts, Network &net, Architecture &arch);
void sim_receive_messages(Timestep &sim, Architecture &arch);
double sim_schedule_messages(std::vector<MessageFifo> &messages_sent, const Scheduler &scheduler);
void sim_update_noc_message_counts(const Message &m, NocInfo &noc, const int message_in);
double sim_calculate_messages_along_route(Message &m, NocInfo &noc);
void sim_update_noc(const double t, NocInfo &noc);

// TODO: reimplement

void sim_process_neuron(Timestep &ts, Architecture &arch, Neuron &n);
double sim_pipeline_receive(Timestep &ts, Architecture &arch, Core &c, Message &m);
double sim_update_synapse(Timestep &ts, Architecture &arch, Core &c, const int synapse_address, const bool synaptic_lookup);
double sim_update_dendrite(Timestep &ts, Architecture &arch, Neuron &n, const double charge);
double sim_update_soma(Timestep &ts, Architecture &arch, Neuron &n, const double current_in);
double sim_estimate_network_costs(Tile &src, Tile &dest);
void sim_neuron_send_spike_message(Timestep &ts, Architecture &arch, Neuron &n);

double sim_update_soma_lif(Timestep &ts, Neuron &n, const double current_in);
double sim_update_soma_truenorth(Timestep &ts, Neuron &n, const double current_in);

void sim_reset_measurements(Network &net, Architecture &arch);
double sim_calculate_energy(const Architecture &arch);
double sim_calculate_time(const Architecture &arch);
long int sim_calculate_packets(const Architecture &arch);

std::ofstream sim_trace_open_perf_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_spike_trace(const std::filesystem::path &out_dir);
std::ofstream sim_trace_open_potential_trace(const std::filesystem::path &out_dir, const Network &net);
std::ofstream sim_trace_open_message_trace(const std::filesystem::path &out_dir);
void sim_trace_write_spike_header(std::ofstream &spike_trace_file);
void sim_trace_write_potential_header(std::ofstream &potential_trace_file, const Network &net);
void sim_trace_write_perf_header(std::ofstream &perf_trace_file);
void sim_trace_write_message_header(std::ofstream &message_trace_file);
void sim_trace_record_spikes(std::ofstream &spike_trace_file, const long int timesteps, const Network &net);
void sim_trace_record_potentials(std::ofstream &potential_trace_file, const int timestep, const Network &net);
void sim_trace_record_message(std::ofstream &perf_trace_file, const Message &m);
void sim_trace_perf_log_timestep(std::ofstream &out, const Timestep &ts);

void sim_output_run_summary(const std::filesystem::path &output_file, const RunData &run_data);
void sim_format_run_summary(std::ostream &out, const RunData &run_data);

void sim_write_summary(FILE *fp, const Simulation &stats);

int sim_poisson_input(const double firing_probability);
int sim_rate_input(const double firing_rate, double *spike_val);

MessageFifo *sim_init_timing_priority(std::vector<MessageFifo> &message_queues);
void sim_insert_priority_queue(MessageFifo **priority_queue, MessageFifo &c);
MessageFifo *sim_pop_priority_queue(MessageFifo **priority_queue);

void sim_init_fifo(MessageFifo &f);
void sim_message_fifo_push(MessageFifo &queue, Message &m);
Message *sim_message_fifo_pop(MessageFifo *queue);

double sim_generate_noise(Neuron *n);
timespec calculate_elapsed_time(const timespec &ts_start, const timespec &ts_end);

}

#endif