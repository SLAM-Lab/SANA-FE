#ifndef SIM_HEADER_INCLUDED
#define SIM_HEADER_INCLUDED

// DEBUG and TRACE print with source annotations, TRACE is only enabled for
//  verbose debug printing
#define INFO(fmt, args...) fprintf(stdout, "[%s:%d:%s()] " fmt, __FILE__, __LINE__, __func__, ##args)
#ifdef DEBUG
#define TRACE(fmt, args...) fprintf(stdout, "[%s:%d:%s()] " fmt, __FILE__, __LINE__, __func__, ##args)
#else
//#define NDEBUG 1 // Turn off assert()
#define TRACE(fmt, args...) do {} while (0)
#endif

#define RAND_SEED 0xbeef // For srand()
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct sim_stats
{
	int time_steps;
	long int total_spikes, total_packets_sent;
	double total_energy, total_sim_time, wall_time;
};

#include "arch.h"
#include "network.h"
#include "stdio.h"
struct sim_stats sim_timestep(struct network *net, struct architecture *arch, FILE *probe_spike_fp, FILE *probe_potential_fp, FILE *perf_fp);

int sim_route_spikes(struct network *net);
int sim_input_spikes(struct network *net);

void sim_update(struct network *net);
void sim_update_neuron(struct neuron *n);
void sim_update_synapse_cuba(struct neuron *n);
void sim_update_dendrite(struct neuron *n);
void sim_update_potential(struct neuron *n);
void sim_update_axon(struct neuron *n);

void sim_reset_measurements(struct network *net, struct architecture *arch);
double sim_calculate_energy(const struct architecture *arch);
double sim_calculate_time(const struct architecture *arch);
long int sim_calculate_packets(const struct architecture *arch);

void sim_write_summary(FILE *fp, const struct sim_stats *stats);
void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp, const struct network *net);
void sim_probe_log_timestep(FILE *spike_fp, FILE *potential_fp, const struct network *net);
void sim_perf_write_header(FILE *perf_fp, const struct architecture *arch);
void sim_perf_log_timestep(FILE *fp, const struct architecture *arch);
int sim_input(const double firing_probability);
#endif
