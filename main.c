// main.c
// Performance simulation for neuromorphic architectures
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sim.h"
#include "network.h"

void init_results(struct sim_results *results);
void next_inputs(char *buffer, struct core *cores, const int max_cores, struct neuron **neuron_ptrs);

enum program_args
{
	TIMESTEPS,
	N_CORES,
	NETWORK_FILENAME,
	PROGRAM_NARGS,
};

#define MAX_INPUT_LINE (32 + MAX_NEURONS*32)

int main(int argc, char *argv[])
{
	FILE *input_fp, *network_fp, *results_fp;
	struct core *cores;
	struct neuron **neuron_ptrs;
	struct sim_results results;
	char *filename, *input_buffer;
	int timesteps, max_cores, max_neurons;

	filename = NULL;
	input_fp = NULL;
	network_fp = NULL;
	results_fp = NULL;
	input_buffer = NULL;

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

		input_buffer = (char *) malloc(sizeof(char) * MAX_INPUT_LINE);
		if (input_buffer == NULL)
		{
			INFO("Error: Couldn't allocate memory for inputs.\n");
			exit(1);
		}
	}

	if (argc < PROGRAM_NARGS)
	{
		INFO("Usage: ./sim [-i <input vectors>] <timesteps> <cores> "
							"<neuron config>\n");
		exit(1);
	}

	// Read in program args, sanity check and parse inputs
	sscanf(argv[TIMESTEPS], "%d", &timesteps);
	if (timesteps <= 0)
	{
		INFO("Time-steps must be > 0 (%d)\n", timesteps);
		exit(1);
	}

	sscanf(argv[N_CORES], "%d", &max_cores);
	if ((max_cores <= 0) || (max_cores > MAX_CORES_LOIHI))
	{
		INFO("Cores must be > 0 and < %d (%d)\n", MAX_CORES_LOIHI,
								max_cores);
		exit(1);
	}
	max_neurons = max_cores * MAX_COMPARTMENTS;

	filename = argv[NETWORK_FILENAME];

	// Create the network
	network_fp = fopen(filename, "r");
	if (network_fp == NULL)
	{
		INFO("Neuron data (%s) failed to open.\n", filename);
		exit(1);
	}

	INFO("Allocating memory for %d cores.\n", max_cores);
	cores = (struct core *) malloc(max_cores * sizeof(struct core));

	INFO("Allocating memory to track %d neurons", max_neurons);
	neuron_ptrs = (struct neuron **)
				malloc(MAX_NEURONS * sizeof(struct neuron *));
	if ((neuron_ptrs == NULL) || (cores == NULL))
	{
		INFO("Error: failed to allocate memory.\n");
		exit(1);
	}
	INFO("Allocated %ld bytes\n", max_cores * sizeof(struct core));
        for (int i = 0; i < max_neurons; i++)
        {
            neuron_ptrs[i] = NULL;
        }
	network_read_csv(network_fp, neuron_ptrs, cores, max_cores);
	fclose(network_fp);

	// TODO: eventually we could have some simple commands like
	//  run <n timesteps>
	//  set rate[neuron #] <firing rate>
	//  set threshold[neuron #] <threshold>
	//  That can either be input from the command line or a file
	//  This could replace having a separate input vector file format
	//  It would also be more general / powerful
	init_results(&results);
	if (input_fp != NULL)
	{
		// Run set of input vectors, each one is presented for the
		//  same number of timesteps
		while (fgets(input_buffer, MAX_INPUT_LINE, input_fp))
		{
			next_inputs(input_buffer, cores, max_cores,
								neuron_ptrs);
			INFO("Next inputs set.\n");
			sim_run(timesteps, cores, max_cores, &results);
		}
	}
	else
	{
		// Single step simulation, based on initial state of network
		sim_run(timesteps, cores, max_cores, &results);
	}

	INFO("Total simulated time: %es.\n", results.total_sim_time);
	INFO("Total energy calculated: %eJ.\n", results.total_energy);
	INFO("Average power consumption: %fW.\n",
				results.total_energy / results.total_sim_time);
	INFO("Run finished.\n");

	results_fp = fopen("results.yaml", "w");
	if (results_fp != NULL)
	{
		sim_write_results(results_fp, &results);
	}
	fclose(results_fp);

	// Cleanup
        free(cores);
	free(neuron_ptrs);
	free(input_buffer);

	return 0;
}

void next_inputs(char *buffer, struct core *cores, const int max_cores,
						struct neuron **neuron_ptrs)
{
	char *token;
	double firing_rate;
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
		if ((firing_rate < 0) || (firing_rate > 1))
		{
			INFO("Warning: input rate not in range [0,1] (%f)",
								firing_rate);
		}

		neuron_ptrs[neuron_count]->input_rate = firing_rate;
		token = strtok(NULL, ",");
		neuron_count++;
	}
}

void init_results(struct sim_results *results)
{
	results->total_energy = 0.0; // Joules
	results->total_sim_time = 0.0; // Seconds
	results->wall_time = 0.0; // Seconds
	results->time_steps = 0;
	results->total_spikes = 0;
}
