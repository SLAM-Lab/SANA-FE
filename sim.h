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

struct simulation
{
	int log_perf, log_spikes, log_potential;
	long int timesteps, total_spikes, total_packets_sent;
	double total_energy, total_sim_time, wall_time;
	FILE *spike_trace_fp, *potential_trace_fp, *perf_fp, *stats_fp;
};

struct timestep
{
	long int timestep, spike_count, packets_sent, total_hops;
	long int total_neurons_fired, spikes;
	double energy, sim_time;
};

struct timing
{
	struct core *c;
	int next;
};

#include "arch.h"
#include "network.h"
#include "stdio.h"

struct timestep sim_timestep(struct simulation *const sim, struct network *const net, struct architecture *const arch);
void sim_init_sim(struct simulation *sim);
void sim_init_timestep(struct timestep *const ts);

void sim_send_messages(struct timestep *const ts, struct network *net, struct architecture *arch);
void sim_receive_messages(struct timestep *const ts, struct network *net, struct architecture *arch);
int sim_input_spikes(struct network *net);
int sim_schedule_messages(struct timestep *const ts, struct network *net, struct architecture *arch);

void sim_process_neuron(struct timestep *const ts, struct network *net, struct neuron *n);
double sim_pipeline_receive(struct timestep *const ts, struct core *c, struct axon_map *axon);
double sim_neuron_send_spike(struct neuron *n);
double sim_update_synapse(struct timestep *const ts, struct axon_map *axon);
double sim_update_dendrite(struct timestep *const ts, struct neuron *n, const double charge);
double sim_update_soma(struct timestep *const ts, struct neuron *n, const double current_in);
double sim_update_axon(struct neuron *n);
double sim_estimate_network_costs(struct timestep *const ts, struct tile *const src, struct tile *const dest);

double sim_update_soma_lif(struct timestep *const ts, struct neuron *n, const double current_in);
double sim_update_soma_truenorth(struct timestep *const ts, struct neuron *n, const double current_in);

void sim_reset_measurements(struct network *net, struct architecture *arch);
double sim_calculate_energy(const struct architecture *const arch);
double sim_calculate_time(const struct architecture *const arch);
double sim_calculate_time_old(const struct architecture *const arch);
long int sim_calculate_packets(const struct architecture *arch);

void sim_write_summary(FILE *fp, const struct architecture *arch, const struct simulation *stats);
void sim_probe_write_header(const struct simulation *sim, const struct network *net);
void sim_probe_log_timestep(const struct simulation *sim, const struct network *net);
void sim_perf_write_header(FILE *perf_fp);
void sim_perf_log_timestep(const struct timestep *const ts, FILE *fp);
int sim_poisson_input(const double firing_probability);
int sim_rate_input(const double firing_rate, double *spike_val);

struct core *sim_init_timing_priority(struct architecture *arch);
struct core *sim_update_timing_queue(struct architecture *arch,
						struct core *top_priority);

#endif
