// main.c
// Performance simulation for neuromorphic architectures
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sim.h"
#include "tech.h"
#include "network.h"

void init_stats(struct sim_stats *stats);
void run(struct technology *tech, struct architecture *arch, struct sim_stats *stats, FILE *probe_spikes_fp, FILE *probe_potential_fp, FILE *perf_fp);
//void next_inputs(char *buffer, struct core *cores, const int max_cores, struct neuron **neuron_ptrs);
struct timespec calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);
int parse_dvs(FILE *fp, struct architecture *arch);

enum program_args
{
	TECHNOLOGY_FILENAME,
	ARCH_FILENAME,
	NETWORK_FILENAME,
	TIMESTEPS,
	PROGRAM_NARGS,
};

int main(int argc, char *argv[])
{
	FILE *input_fp, *network_fp, *stats_fp, *tech_fp, *arch_fp;
	FILE *probe_spikes_fp, *probe_potential_fp, *perf_fp;
	struct architecture arch;
	struct technology tech;
	struct neuron **neuron_ptrs;
	struct sim_stats stats;
	char *filename, *input_buffer;
	int timesteps, max_input_line;

	filename = NULL;
	input_fp = NULL;
	network_fp = NULL;
	stats_fp = NULL;
	input_buffer = NULL;
	probe_spikes_fp = NULL;
	probe_potential_fp = NULL;

	if (argc < 1)
	{
		INFO("Error: No program arguments.\n");
		exit(1);
	}
	// First arg is always program name, skip
	argc--;
	argv++;

	// Parse optional args
	if (argc > 2)
	{
		if ((argv[0][0] == '-') && (argv[0][1] == 'i'))
		{
			// Optional input vector argument
			filename = argv[1];
			argc--;
			argv += 2;
		}
	}

	if (filename)
	{
		input_fp = fopen(filename, "r");
		if (input_fp == NULL)
		{
			INFO("Error: Couldn't open inputs %s.\n", filename);
			exit(1);
		}
	}

	if (argc < PROGRAM_NARGS)
	{
		INFO("Usage: ./sim [-i <input vectors>] <tech file> <arch description>"
				" <neuron config> <timesteps>\n");
		exit(1);
	}

	// Initialize the technology parameters - the chip parameters and key
	//  metrics
	tech_init(&tech);

	// Read in program args, sanity check and parse inputs
	filename = argv[TECHNOLOGY_FILENAME];
	tech_fp = fopen(filename, "r");
	if (tech_fp == NULL)
	{
		INFO("Error: Technology file failed to open.\n");
		exit(1);
	}
	tech_read_file(&tech, tech_fp);

	filename = argv[ARCH_FILENAME];
	arch_fp = fopen(filename, "r");
	if (arch_fp == NULL)
	{
		INFO("Error: Architecture file failed to open.\n");
		exit(1);
	}
	arch = arch_read_file(arch_fp);

	sscanf(argv[TIMESTEPS], "%d", &timesteps);
	if (timesteps <= 0)
	{
		INFO("Time-steps must be > 0 (%d)\n", timesteps);
		exit(1);
	}

	//max_neurons = max_cores * tech.max_neurons;
	// Input line must be long enough to encode inputs for all neurons
	//  simultaneously
	max_input_line = 32 + (arch.max_neurons*32);

	filename = argv[NETWORK_FILENAME];
	// Create the network
	network_fp = fopen(filename, "r");
	if (network_fp == NULL)
	{
		INFO("Neuron data (%s) failed to open.\n", filename);
		exit(1);
	}

	// Open the probe output files for writing, for now hard code filenames
	//  TODO: add option for command line argument
	probe_potential_fp = fopen("probe_potential.csv", "w");
	probe_spikes_fp = fopen("probe_spikes.csv", "w");
	if ((probe_potential_fp == NULL) || (probe_spikes_fp == NULL))
	{
		INFO("Error: Couldn't open probe output files for writing.\n");
		exit(1);
	}

	// TODO: add option for command line argument
	// TODO: also want the option of enabling and disabling this
	//  for quick runs we probably don't want it?
	perf_fp = fopen("perf.csv", "w");
	if (perf_fp == NULL)
	{
		INFO("Error: Couldn't open perf output files for writing.\n");
		exit(1);
	}

	INFO("Allocating memory to track %d neurons.\n",
							arch.max_neurons);
	neuron_ptrs = (struct neuron **)
		malloc(arch.max_neurons * sizeof(struct neuron *));
	if (neuron_ptrs == NULL)
	{
		INFO("Error: failed to allocate memory.\n");
		exit(1);
	}
	for (int i = 0; i < arch.max_neurons; i++)
	{
		neuron_ptrs[i] = NULL;
	}

	INFO("Reading network from file.\n");
	network_read_csv(network_fp, &tech, &arch, neuron_ptrs);
	fclose(network_fp);

	init_stats(&stats);
	INFO("Creating probe and perf data files.\n");
	sim_probe_write_header(probe_spikes_fp, probe_potential_fp, &arch);
	sim_perf_write_header(perf_fp, &arch);
	if (input_fp != NULL)
	{

		// Allocate a buffer for the inputs
		input_buffer = (char *) malloc(sizeof(char) * max_input_line);
		if (input_buffer == NULL)
		{
			INFO("Error: Couldn't allocate memory for inputs.\n");
			exit(1);
		}
		// Parse DVS128 gesture input
		// TODO: make this parameterised so we can parse different input
		//  formats
		//while (fgets(input_buffer, max_input_line, input_fp))

		while (parse_dvs(input_fp, &arch))
		{
			run(&tech, &arch, &stats, probe_spikes_fp,
						probe_potential_fp, perf_fp);
		}
	}
	// TODO: another option to generate rate based / poisson spike trains
	else
	{
		// Single step simulation, based on initial network state and
		//  no inputs
		for (int i = 0; i < timesteps; i++)
		{
			INFO("*** Time-step %d ***\n", i+1);
			run(&tech, &arch, &stats, probe_spikes_fp,
						probe_potential_fp, perf_fp);
		}
	}

	INFO("***** Run Summary *****\n");
	sim_write_summary(stdout, &stats);
	INFO("Average power consumption: %f W.\n",
				stats.total_energy / stats.total_sim_time);
	stats_fp = fopen("stats.yaml", "w");
	if (stats_fp != NULL)
	{
		sim_write_summary(stats_fp, &stats);
	}
	fclose(stats_fp);
	INFO("Run finished.\n");

	// Cleanup
	arch_free(&arch);
	free(neuron_ptrs);

	// Close any open files
	fclose(probe_potential_fp);
	fclose(probe_spikes_fp);
	fclose(perf_fp);

	return 0;
}

void run(struct technology *tech, struct architecture *arch,
			struct sim_stats *stats, FILE *probe_spikes_fp,
					FILE *probe_potential_fp, FILE *perf_fp)
{
	// Run neuromorphic hardware simulation for one timestep
	//  Measure the CPU time it takes and accumulate the stats
	struct timespec ts_start, ts_end, ts_elapsed;
	struct sim_stats timestep_stats;

	// Measure the wall-clock time taken to run the simulation
	//  on the host machine
	clock_gettime(CLOCK_MONOTONIC, &ts_start);

	timestep_stats = sim_timestep(tech, arch, probe_spikes_fp,
					probe_potential_fp, perf_fp);
	// Accumulate totals for the entire simulation
	// TODO: make a function
	stats->total_energy += timestep_stats.total_energy;
	stats->total_sim_time += timestep_stats.total_sim_time;
	stats->total_spikes += timestep_stats.total_spikes;
	stats->total_packets_sent += timestep_stats.total_packets_sent;

	// Calculate elapsed time
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	ts_elapsed = calculate_elapsed_time(ts_start, ts_end);
	stats->wall_time +=
		(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9);
	INFO("Time-step took: %fs.\n",
		(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9));
}

/*
void next_inputs(char *buffer, struct core *cores, const int max_cores,
						struct neuron **neuron_ptrs)
{
	char *token;
	int neuron_count;

	neuron_count = 0;
	token = strtok(buffer, ",");
	while (token != NULL)
	{
		// This time read all the fields in the line, we're
		//  interested in the synapse data
		int ret = sscanf(token, "%lf", &firing_rate);
		if (ret <= 0)
		{
			INFO("Error: invalid input format (%s)", buffer);
			exit(1);
		}

		token = strtok(NULL, ",");
		neuron_count++;
	}
}
*/

void init_stats(struct sim_stats *stats)
{
	stats->total_energy = 0.0; // Joules
	stats->total_sim_time = 0.0; // Seconds
	stats->wall_time = 0.0; // Seconds
	stats->time_steps = 0;
	stats->total_spikes = 0;
	stats->total_packets_sent = 0;
}

struct timespec calculate_elapsed_time(struct timespec ts_start,
							struct timespec ts_end)
{
	// Calculate elapsed wall-clock time between ts_start and ts_end
	struct timespec ts_elapsed;

	ts_elapsed.tv_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	ts_elapsed.tv_sec = ts_end.tv_sec - ts_start.tv_sec;
	if (ts_end.tv_nsec < ts_start.tv_nsec)
	{
		ts_elapsed.tv_sec--;
		ts_elapsed.tv_nsec += 1000000000UL;
	}

	return ts_elapsed;
}

int parse_dvs(FILE *fp, struct architecture *arch)
{
	// TODO
	return 0;
}
