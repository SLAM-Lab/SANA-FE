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

struct compartment
{
	int id, compartment_used, fired, post_connection_count, spike_count;
	int log_spikes, log_voltage, update_needed, force_update;
	double *time;
	double potential, current,  bias, reset, threshold;
	double potential_decay, potential_time_const;
	double current_decay, current_time_const;
	double energy;

	struct synapse *synapses;
	struct axon_output *axon_out;
	struct axon_input *axon_in;
};

struct synapse
{
	int id;
	struct compartment *post_neuron, *pre_neuron;
	double energy, weight;
	struct mem *memory;
};

struct mem
{
	int id, memory_size_bytes, weight_bits;
};

struct axon_output
{
	int id, fan_out;
	double energy;
	struct router *r;
	unsigned int *packets_sent;
};

struct axon_input
{
	int id, fan_in;
	struct router *r;
};

struct router
{
	int x, y, id;
	int max_dimensions, width; // For now just support 2 dimensions
	double energy;
	struct router *east, *west, *north, *south;
};

struct input
{
	int send_spike, post_connection_count;
	struct synapse *synapses;
};

struct sim_results
{
	int time_steps;
	long int total_spikes;
	double total_energy, total_sim_time, wall_time;
};

#include "tech.h"
#include "arch.h"
struct sim_results sim_timestep(const struct technology *tech, struct architecture *arch, FILE *probe_spike_fp, FILE *probe_potential_fp);

int sim_route_spikes(const struct technology *tech, struct architecture *arch);
int sim_input_spikes(const struct technology *tech, struct architecture *arch);

void sim_update(const struct technology *tech, struct architecture *arch);
void sim_update_compartment(const struct technology *tech, struct compartment *c);
void sim_update_synapse_cuba(const struct technology *tech, struct compartment *c);
void sim_update_dendrite(const struct technology *tech, struct compartment *c);
void sim_update_potential(const struct technology *tech, struct compartment *c);
void sim_update_axon(const struct technology *tech, struct compartment *c);

void sim_reset_measurements(struct architecture *arch);
double sim_calculate_energy(struct architecture *arch);
double sim_calculate_time(const struct technology *tech, struct architecture *arch);

void sim_write_summary(FILE *fp, const struct sim_results *results);
void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp, const struct architecture *arch);
void sim_probe_log_timestep(FILE *spike_fp, FILE *potential_fp, const struct architecture *arch);
int sim_input(const double firing_probability);
#endif
