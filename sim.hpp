// Copyright (c) 2023 - The University of Texas at Austin
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
#include "arch.hpp"
#include "network.hpp"
#include "stdio.h"

struct Timestep
{
	std::vector<std::list<Message> > messages;
	std::vector<MessageFifo> message_queues;
	long int timestep, spike_count, total_hops, packets_sent;
	long int total_neurons_fired, spikes;
	double energy, sim_time;

	Timestep(const long int ts, const int core_count);
};

struct NocInfo
{
	struct MessageFifo messages_received[ARCH_MAX_CORES];
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

struct Simulation
{
	long int timesteps, total_spikes, total_messages_sent;
	long int total_neurons_fired;
	double total_energy, total_sim_time, wall_time;
	int log_perf, log_spikes, log_potential, log_messages;
	int gui_on;
	FILE *spike_trace_fp, *potential_trace_fp, *message_trace_fp, *perf_fp;
	FILE *stats_fp;

	Simulation();
};

//#include "pybind11/pybind11.h"
//#include "pybind11/stl.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

enum program_args
{
	ARCH_FILENAME,
	NETWORK_FILENAME,
	TIMESTEPS,
	PROGRAM_NARGS,
};

struct run_ts_data
{
	double energy;
	double time;
	double wall_time;
	long int spikes, packets, neurons;
	long int timestep_start, timesteps;
};

class SanaFe
{
public:
	std::unique_ptr<Simulation> sim;
	Network net;
	Architecture *arch;
	int timesteps;
	FILE *input_fp;
	struct run_ts_data run_data;
	std::string out_dir;

        SanaFe();
	void init();
	int update_neuron(std::vector<NeuronGroup>::size_type group_id, std::vector<Neuron>::size_type n_id, std::vector<string> kwargs, int count);
	void run_timesteps(int timesteps = 1);
	void set_spike_trace(const bool enable=true);
	void set_potential_trace(const bool enable=true);
	void set_perf_trace(const bool enable=true);
	void set_message_trace(const bool enable=true);
	void set_out_dir(const std::string dir);
	void open_perf_trace(void);
	void open_spike_trace(void);
	void open_potential_trace(void);
	void open_message_trace(void);
	void set_gui_flag(bool flag = true);
	void set_arch(const char *filename);
	void set_net(const char *filename);
        double get_power();
	std::vector<int> get_status(const std::vector<NeuronGroup>::size_type gid);
        void sim_summary();
	std::vector<std::vector<int>> run_summary();
	void clean_up(int ret = RET_OK);
	~SanaFe() { clean_up(); };
};

class Vector_Cleanup_Class
{
public:
	Architecture* arch;
	Vector_Cleanup_Class(){}
	Vector_Cleanup_Class(Architecture* arch)
	{
		this->arch = arch;
	}
	~Vector_Cleanup_Class()
	{
		//arch->spike_vector.clear();
	}
};

void run(Simulation &sim, Network &net, Architecture &arch);
struct timespec calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);
void store_data_init(run_ts_data *data, Simulation &sim, int timesteps);
void store_data(run_ts_data *data, Simulation &sim);
void print_run_data(FILE *fp, run_ts_data* data);

void sim_timestep(Timestep &ts, Network &net, Architecture &arch);
std::unique_ptr<Simulation> sim_init_sim(void);
void sim_init_timestep(Timestep &ts, Architecture &arch);

void sim_process_neurons(Timestep &ts, Network &net, Architecture &arch);
void sim_receive_messages(Timestep &sim, Architecture &arch);
double sim_schedule_messages(std::vector<MessageFifo> &messages_sent, const Scheduler &scheduler);
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

void sim_write_summary(FILE *fp, const Simulation &stats);
void sim_spike_trace_write_header(const Simulation &sim);
void sim_potential_trace_write_header(const Simulation &sim, const Network &net);
void sim_message_trace_write_header(const Simulation &sim);
void sim_trace_record_spikes(const Simulation &sim, const Network &net);
void sim_trace_record_potentials(const Simulation &sim, const Network &net);
void sim_trace_record_message(const Simulation &sim, const Message &m);
void sim_perf_write_header(FILE *perf_fp);
void sim_perf_log_timestep(const Timestep &ts, FILE *fp);
int sim_poisson_input(const double firing_probability);
int sim_rate_input(const double firing_rate, double *spike_val);

MessageFifo *sim_init_timing_priority(std::vector<MessageFifo> &message_queues);
void sim_insert_priority_queue(MessageFifo **priority_queue, MessageFifo &c);
MessageFifo *sim_pop_priority_queue(MessageFifo **priority_queue);

void sim_message_fifo_push(struct MessageFifo &queue, Message &m);
Message *sim_message_fifo_pop(struct MessageFifo &queue);

#endif