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
#define MAX_NOISE_FILE_ENTRY 128

#include "arch.h"
#include "network.h"
#include "stdio.h"

struct timestep
{
	struct message messages[ARCH_MAX_CORES][ARCH_MAX_CONNECTION_MAP+1];
	struct message_fifo message_queues[ARCH_MAX_CORES];
	long int timestep, spike_count, total_hops, packets_sent;
	long int total_neurons_fired, spikes;
	double energy, sim_time;
};

struct noc_timings
{
	struct message_fifo messages_received[ARCH_MAX_CORES];
	size_t noc_width, noc_height;
	double messages_in_flight[ARCH_MAX_X][ARCH_MAX_Y];
	double noc_messages_in[ARCH_MAX_X][ARCH_MAX_Y][4+ARCH_MAX_CORES_PER_TILE];
	double noc_messages_out[ARCH_MAX_X][ARCH_MAX_Y][4+ARCH_MAX_CORES_PER_TILE];
	long int messages_in_noc, noc_capacity;
	double mean_in_flight_receive_delay, mean_hops, utilization;
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
	FILE *spike_trace_fp, *potential_trace_fp, *message_trace_fp, *perf_fp;
	FILE *stats_fp;
};

void sim_timestep(struct timestep *const ts, struct network *const net, struct architecture *const arch);
struct simulation *sim_init_sim(void);
void sim_init_timestep(struct timestep *const ts);

void sim_process_neurons(struct timestep *const ts, struct network *net, struct architecture *arch);
void sim_receive_messages(struct timestep *const sim, struct architecture *arch);
double sim_schedule_messages(struct message_fifo *const messages_sent, const int noc_width, const int noc_height);
double sim_schedule_messages_blocking(struct message_fifo *const messages_sent);
double sim_schedule_messages_simplest(struct message_fifo *const messages_sent);
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
