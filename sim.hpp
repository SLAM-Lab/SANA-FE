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

#include "arch.hpp"
#include "network.hpp"
#include "stdio.h"

struct timestep
{
	struct message messages[ARCH_MAX_CORES][ARCH_MAX_CONNECTION_MAP+1];
	struct message_fifo message_queues[ARCH_MAX_CORES];
	long int timestep, spike_count, total_hops, packets_sent;
	long int total_neurons_fired, spikes;
	int message_counts[ARCH_MAX_CORES];
	double energy, sim_time;
};

struct noc_info
{
	struct message_fifo messages_received[ARCH_MAX_CORES];
	int noc_width, noc_height;
	double noc_messages_in[ARCH_MAX_X][ARCH_MAX_Y][4+ARCH_MAX_CORES_PER_TILE];
	double core_finished_receiving[ARCH_MAX_CORES];
	double mean_in_flight_receive_delay;
	long int messages_in_noc;
};

struct scheduler
{
	int noc_width, noc_height, buffer_size;
};

enum direction
{
	NORTH = 0,
	EAST,
	SOUTH,
	WEST
};

struct simulation
{
	struct timestep ts;
	long int timesteps, total_spikes, total_messages_sent;
	long int total_neurons_fired;
	double total_energy, total_sim_time, wall_time;
	int log_perf, log_spikes, log_potential, log_messages;
	int gui_on;
	FILE *spike_trace_fp, *potential_trace_fp, *message_trace_fp, *perf_fp;
	FILE *stats_fp;
};

//#include "pybind11/pybind11.h"
//#include "pybind11/stl.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

using namespace std;

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


class sana_fe
{
public:
	struct simulation *sim;
	struct network net;
	struct architecture *arch;
	int timesteps;
	FILE *input_fp;
	struct run_ts_data run_data;

        sana_fe();
	void init();
	int update_neuron(int group_id, int n_id, vector<string> kwargs, int count);
	void run_timesteps(int timesteps = 1);
	void set_input(char *filename);
	void open_perf_trace(void);
	void open_spike_trace(void);
	void open_potential_trace(void);
	void open_message_trace(void);
	void set_gui_flag(bool flag = true);
	void set_arch(const char *filename);
	void set_net(const char *filename);
        double get_power();
	vector<int> get_status(int gid);
        void sim_summary();
	vector<vector<int>> run_summary();
	void clean_up(int ret = RET_OK);
	~sana_fe() { clean_up(); };
};

class Vector_Cleanup_Class
{
public:
	struct architecture* arch;
	Vector_Cleanup_Class(){}
	Vector_Cleanup_Class(struct architecture* arch)
	{
		this->arch = arch;
	}
	~Vector_Cleanup_Class()
	{
		arch->spike_vector.clear();
	}
};

void run(struct simulation *sim, struct network *net, struct architecture *arch);
struct timespec calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);
void store_data_init(run_ts_data* data, simulation* sim, int timesteps);
void store_data(run_ts_data* data, simulation* sim);
void print_run_data(FILE *fp, run_ts_data* data);


void sim_timestep(struct timestep *const ts, struct network *const net, struct architecture *const arch);
struct simulation *sim_init_sim(void);
void sim_init_timestep(struct timestep *const ts);

void sim_process_neurons(struct timestep *const ts, struct network *net, struct architecture *arch);
void sim_receive_messages(struct timestep *const sim, struct architecture *arch);
double sim_schedule_messages(struct message_fifo *const messages_sent, const struct scheduler *const scheduler);
// TODO: reimplement
int sim_input_spikes(struct network *net);

void sim_process_neuron(struct timestep *const ts, struct neuron *n);
double sim_pipeline_receive(struct timestep *const ts, struct core *c, struct connection_map *axon);
double sim_update_synapse(struct timestep *const ts, struct connection_map *axon, const int synaptic_lookup);
double sim_update_dendrite(struct timestep *const ts, struct neuron *n, const double charge);
double sim_update_soma(struct timestep *const ts, struct neuron *n, const double current_in);
double sim_update_axon(struct neuron *n);
double sim_estimate_network_costs(struct tile *const src, struct tile *const dest);
void sim_neuron_send_spike_message(struct timestep *const ts, struct neuron *n);

double sim_update_soma_latency(struct timestep *const ts, struct neuron *n);
double sim_update_soma_lif(struct timestep *const ts, struct neuron *n, const double current_in);
double sim_update_soma_truenorth(struct timestep *const ts, struct neuron *n, const double current_in);

void sim_reset_measurements(struct network *net, struct architecture *arch);
double sim_calculate_energy(const struct architecture *const arch);
double sim_calculate_time(const struct architecture *const arch);
long int sim_calculate_packets(const struct architecture *arch);

void sim_write_summary(FILE *fp, const struct simulation *stats);
void sim_spike_trace_write_header(const struct simulation *const sim);
void sim_potential_trace_write_header(const struct simulation *const sim, const struct network *net);
void sim_message_trace_write_header(const struct simulation *const sim);
void sim_trace_record_spikes(const struct simulation *const sim, const struct network *net);
void sim_trace_record_potentials(const struct simulation *const sim, const struct network *net);
void sim_trace_record_message(const struct simulation *const sim, const struct message *const m);
void sim_perf_write_header(FILE *perf_fp);
void sim_perf_log_timestep(const struct timestep *const ts, FILE *fp);
int sim_poisson_input(const double firing_probability);
int sim_rate_input(const double firing_rate, double *spike_val);

struct message_fifo *sim_init_timing_priority(struct message_fifo *const send_queues);
void sim_insert_priority_queue(struct message_fifo **priority_queue, struct message_fifo *c);
struct message_fifo *sim_pop_priority_queue(struct message_fifo **priority_queue);

void sim_message_fifo_push(struct message_fifo *queue, struct message *m);
struct message *sim_message_fifo_pop(struct message_fifo *queue);


#endif