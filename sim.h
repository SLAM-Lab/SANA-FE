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

// Some max parameters specific to Loihi
#define MAX_CORES_LOIHI 128
#define CORES_PER_TILE 4
#define MAX_COMPARTMENTS 1024
#define FAN_OUT 4096
#define MAX_NEURONS (MAX_CORES_LOIHI * MAX_COMPARTMENTS)
#define MAX_SPIKES (MAX_COMPARTMENTS * FAN_OUT)
#define RAND_SEED 0xbeef // For srand()

// Energy and time estimates of different events, generated from SPICE
//  simulations of Loihi.  All numbers were taken from:
//  "Loihi: A Neuromorphic Manycore Processor with On-Chip Learning" (2018)
//  M. Davies et al

// Energy estimates
#define ACTIVE_NEURON_UPDATE_ENERGY 81.0e-12 // J
#define INACTIVE_NEURON_UPDATE_ENERGY 52.0e-12 // J
#define SPIKE_OP_ENERGY 23.6e-12 // J
#define SPIKE_WITHIN_TILE_ENERGY 1.7e-12 // J
#define EAST_WEST_HOP_ENERGY 3.0e-12 // J
#define NORTH_SOUTH_HOP_ENERGY 4.0e-12 // J

// Time estimates
#define ACTIVE_NEURON_UPDATE_TIME 8.4e-9 // s
#define INACTIVE_NEURON_UPDATE_TIME 5.3e-9 // s
#define SPIKE_OP_TIME 3.5e-9 // s
#define SPIKE_WITHIN_TILE_TIME 2.1e-9 // s
//#define SPIKE_OP_TIME 0.7e-9 // TODO: just experimenting to see trends remove this
#define EAST_WEST_HOP_TIME 4.1e-9// s
#define NORTH_SOUTH_HOP_TIME 6.5e-9 // s
// TODO: replcace with a function that calculates based on the number of tiles
//  used, the range is 113-465 ns (1-32 tiles).  At the moment most of my
//  experiments are just using a single tile
#define MESH_BARRIER_TIME 113.0e-9 // s

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct neuron
{
	int fired, post_connection_count, core_id, id, compartment;
	int log_spikes, log_voltage;
	double potential, current, bias, threshold, reset, energy, time;
	double potential_decay, current_decay, potential_time_const;
	double current_time_const, input_rate;
};

struct synapse
{
	struct neuron *post_neuron, *pre_neuron;
	float weight;
};

struct spike
{
	struct synapse *synapse;
};

struct core
{
	int id, x, y, spike_count, compartments;
	struct neuron neurons[MAX_COMPARTMENTS];
	struct synapse synapses[MAX_COMPARTMENTS][FAN_OUT];
};

struct sim_results
{
	double total_energy, total_sim_time, wall_time;
	int time_steps;
	long int total_spikes;
};

void sim_run(const int timesteps, struct core *cores, const int max_cores, struct sim_results *results);
void sim_update_neurons(struct core *cores, const int max_cores);
int sim_route_spikes(struct core *cores, const int max_cores);
void sim_update_potential(struct neuron *n);
void sim_send_spike(struct synapse *s);
void sim_seed_input_spikes(struct core *cores, const int max_cores);
double sim_calculate_time(struct core *cores, const int max_cores);
void sim_reset_measurements(struct core *cores, const int max_cores);
double sim_calculate_energy(struct core *cores, const int max_cores);
struct timespec sim_calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);
void sim_write_results(FILE *fp, struct sim_results *results);
int sim_input(const double firing_probability);
void sim_timestep(struct sim_results *results, struct core *cores, const int max_cores);

#endif
