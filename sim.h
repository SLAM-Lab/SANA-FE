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

struct sim_stats
{
	int time_steps;
	long int total_spikes, total_packets_sent;
	double total_energy, total_sim_time, wall_time;
	double network_time;
};

#include "arch.h"
#include "network.h"
#include "stdio.h"
struct sim_stats sim_timestep(struct network *const net, struct architecture *const arch, FILE *probe_spike_fp, FILE *probe_potential_fp, FILE *perf_fp);

int sim_route_spikes(struct network *net);
int sim_input_spikes(struct network *net);

void sim_update(struct network *net);
void sim_update_neuron(struct neuron *n);
void sim_update_synapse_cuba(struct neuron *n);
void sim_update_dendrite(struct neuron *n);
void sim_update_potential(struct neuron *n);
void sim_update_axon(struct neuron *n);

void sim_reset_measurements(struct network *net, struct architecture *arch);
double sim_calculate_energy(const struct architecture *const arch, const double time);
double sim_calculate_time(const struct architecture *const arch, double *network_time);
long int sim_calculate_packets(const struct architecture *arch);

void sim_write_summary(FILE *fp, const struct sim_stats *stats);
void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp, const struct network *net);
void sim_probe_log_timestep(FILE *spike_fp, FILE *potential_fp, const struct network *net);
void sim_perf_write_header(FILE *perf_fp, const struct architecture *arch);
void sim_perf_log_timestep(FILE *fp, const struct architecture *arch);
int sim_poisson_input(const double firing_probability);
int sim_rate_input(const double firing_rate, double *spike_val);
#endif
