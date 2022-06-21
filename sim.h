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

struct neuron
{
	int fired, post_connection_count, core_id, id, compartment;
	int log_spikes, log_voltage;
	double potential, current, bias, threshold, reset;
	double potential_decay, current_decay, potential_time_const;
	double current_time_const;
};

struct synapse
{
	struct neuron *post_neuron, *pre_neuron;
	float weight;
};

//struct spike
//{
//	struct synapse *synapse;
//};

struct core
{
	int id, x, y, spike_count, compartments;
	unsigned int *packets_sent;
	double energy, time;
	struct neuron *neurons;
	struct synapse **synapses;
};

struct input
{
	int send_spike, post_connection_count;
	struct synapse *synapses;
};

struct sim_results
{
	double total_energy, total_sim_time, wall_time;
	int time_steps;
	long int total_spikes;
};

#include "tech.h"
struct sim_results sim_timestep(const struct technology *tech, struct core *cores, struct input *inputs, FILE *probe_spike_fp, FILE *probe_potential_fp);

int sim_route_spikes(const struct technology *tech, struct core *cores);
int sim_input_spikes(const struct technology *tech, struct core *cores, struct input *inputs);

void sim_update_neurons(const struct technology *tech, struct core *cores);
void sim_update_potential(const struct technology *tech, struct neuron *n, struct core *c);
void sim_update_synapse_cuba(const struct technology *tech, struct core *c, struct neuron *n);
void sim_update_dendrite(const struct technology *tech, struct core *c, struct neuron *n);
void sim_update_lif(const struct technology *tech, struct core *c, struct neuron *n);
void sim_update_axon(const struct technology *tech, struct core *c, struct neuron *n);

void sim_reset_measurements(struct core *cores, const int max_cores);
double sim_calculate_energy(struct core *cores, const int max_cores);
double sim_calculate_time(const struct technology *tech, struct core *cores);

void sim_write_results(FILE *fp, const struct sim_results *results);
void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp, const struct core *cores, const int max_cores);
void sim_probe_log_timestep(FILE *spike_fp, FILE *potential_fp, const struct core *cores, const int max_cores);
int sim_input(const double firing_probability);
#endif
