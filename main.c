// main.c
// Performance simulation for neuromorphic architectures
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sim.h"
#include "tech.h"
#include "network.h"

void init_results(struct sim_results *results);
void run(struct technology *tech, struct core *cores, struct sim_results *results, struct input *inputs, FILE *probe_spikes_fp, FILE *probe_potential_fp);
//void next_inputs(char *buffer, struct core *cores, const int max_cores, struct neuron **neuron_ptrs);
struct timespec calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);
int parse_dvs(FILE *fp, struct input *inputs, const int max_inputs);

enum program_args
{
	TECHNOLOGY_FILENAME,
	TIMESTEPS,
	NETWORK_FILENAME,
	PROGRAM_NARGS,
};

int main(int argc, char *argv[])
{
	FILE *input_fp, *network_fp, *results_fp, *tech_fp;
	FILE *probe_spikes_fp, *probe_potential_fp;
	struct technology tech;
	struct core *cores;
	struct input *inputs;
	struct neuron **neuron_ptrs;
	struct sim_results results;
	char *filename, *input_buffer;
	int timesteps, max_neurons, max_input_line, max_cores;

	filename = NULL;
	input_fp = NULL;
	network_fp = NULL;
	results_fp = NULL;
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
		INFO("Usage: ./sim [-i <input vectors>] <tech file> <timesteps>"
						" <neuron config>\n");
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
		INFO("Error: Tech file failed to open.\n");
		exit(1);
	}
	tech_read_file(&tech, tech_fp);
	max_cores = tech.max_cores;

	sscanf(argv[TIMESTEPS], "%d", &timesteps);
	if (timesteps <= 0)
	{
		INFO("Time-steps must be > 0 (%d)\n", timesteps);
		exit(1);
	}

	max_neurons = max_cores * tech.max_compartments;
	// Input line must be long enough to encode inputs for all neurons
	//  simultaneously
	max_input_line = 32 + (max_neurons*32);

	filename = argv[NETWORK_FILENAME];
	// Create the network
	network_fp = fopen(filename, "r");
	if (network_fp == NULL)
	{
		INFO("Neuron data (%s) failed to open.\n", filename);
		exit(1);
	}

	// Create all the input nodes, which are the interface between
	//  externally generated spike trains and our spiking network
	INFO("Allocating memory for %d inputs.\n", tech.max_inputs);
	inputs = (struct input *)
				malloc(tech.max_inputs * sizeof(struct input));
	if (inputs == NULL)
	{
		INFO("Failed to allocate input memory.\n");
		exit(1);
	}
	// Zero initialize the input nodes
	for (int i = 0; i < tech.max_inputs; i++)
	{
		struct input *in = &(inputs[i]);

		in->synapses = (struct synapse *) malloc(tech.fan_out *
							sizeof(struct synapse));
		in->post_connection_count = 0;
		in->send_spike = 0;
		if (in->synapses == NULL)
		{
			INFO("Failed to allocate input synapse memory.\n");
			exit(1);
		}
	}

	// Open the probe output files for writing, for now hard code filenames
	//  TODO: add option to command line argument
	probe_potential_fp = fopen("probe_potential.csv", "w");
	probe_spikes_fp = fopen("probe_spikes.csv", "w");
	if ((probe_potential_fp == NULL) || (probe_spikes_fp == NULL))
	{
		INFO("Error: Couldn't open probe output files for writing.\n");
		exit(1);
	}

	INFO("Allocating memory for %d cores.\n", max_cores);
	cores = (struct core *) malloc(max_cores * sizeof(struct core));

	for (int i = 0; i < max_cores; i++)
	{
		struct core *c = &(cores[i]);
		// Allocate each core, creating memory for the compartments i.e.
		//  neurons, and all the synaptic data. Since the parameters are
		//  defined at runtime, this must be done dynamically
		c->neurons = (struct neuron *) malloc(tech.max_compartments *
							sizeof(struct neuron));
		c->packets_sent = (unsigned int *)
				malloc(max_cores * sizeof(unsigned int));

		if ((c->neurons == NULL) || (c->packets_sent == NULL))
		{
			INFO("Error: failed to allocate neuron.\n");
			exit(1);
		}

		c->synapses = (struct synapse **) malloc(tech.max_compartments *
						sizeof(struct synapse *));
		if (c->synapses == NULL)
		{
			INFO("Error: failed to allocate synapse ptr.\n");
			exit(1);
		}
		// For each compartment, allocate a certain amount of synaptic
		//  memory
		// TODO: this isn't really how it works in the wild. What would
		//  be better is to allocate a single block of synaptic memory
		//  per core. Then each neuron compartment has an index into
		//  its memory within the block. This defines way too much
		//  synaptic memory (GB) when it should be MB.
		for (int j = 0; j < tech.max_compartments; j++)
		{
			c->synapses[j] =
				(struct synapse *) malloc(tech.fan_out *
							sizeof(struct synapse));
			if (c->synapses[j] == NULL)
			{
				INFO("Error: failed to allocate synapse.\n");
			}
		}
	}

	INFO("Allocating memory to track %d neurons.\n", max_neurons);
	neuron_ptrs = (struct neuron **)
				malloc(max_neurons * sizeof(struct neuron *));
	if ((neuron_ptrs == NULL) || (cores == NULL))
	{
		INFO("Error: failed to allocate memory.\n");
		exit(1);
	}
	//INFO("Allocated %ld bytes\n", max_cores * sizeof(struct core));
	for (int i = 0; i < max_neurons; i++)
	{
		neuron_ptrs[i] = NULL;
	}
	network_read_csv(network_fp, neuron_ptrs, cores, &tech, inputs);
	fclose(network_fp);

	init_results(&results);
	sim_probe_write_header(probe_spikes_fp, probe_potential_fp, cores,
								tech.max_cores);
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
		while (parse_dvs(input_fp, inputs, tech.max_inputs))
		{
			run(&tech, cores, &results, inputs, probe_spikes_fp,
							probe_potential_fp);
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
			run(&tech, cores, &results, inputs, probe_spikes_fp,
							probe_potential_fp);
		}
	}

	INFO("Total simulated time: %es.\n", results.total_sim_time);
	INFO("Total energy calculated: %eJ.\n", results.total_energy);
	INFO("Total spikes processed: %ld.\n", results.total_spikes);
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
	for (int i = 0; i < max_cores; i++)
	{
		struct core *c = &(cores[i]);

		free(c->neurons);
		free(c->packets_sent);
		for (int j = 0; j < tech.max_compartments; j++)
		{
			free(c->synapses[j]);
		}
		free(c->synapses);
	}
	for (int i = 0; i < tech.max_inputs; i++)
	{
		struct input *in = &(inputs[i]);

		free(in->synapses);
	}
	free(cores);
	free(neuron_ptrs);
	free(input_buffer);
	free(inputs);

	// Close any open files
	fclose(probe_potential_fp);
	fclose(probe_spikes_fp);

	return 0;
}

void run(struct technology *tech, struct core *cores,
			struct sim_results *results, struct input *inputs,
			FILE *probe_spikes_fp, FILE *probe_potential_fp)
{
	// Run neuromorphic hardware simulation for one timestep
	//  Measure the CPU time it takes and accumulate the results
	struct timespec ts_start, ts_end, ts_elapsed;
	struct sim_results timestep_results;

	// Measure the wall-clock time taken to run the simulation
	//  on the host machine
	clock_gettime(CLOCK_MONOTONIC, &ts_start);

	timestep_results = sim_timestep(tech, cores, inputs,
					probe_spikes_fp, probe_potential_fp);
	// Accumulate totals for the entire simulation
	// TODO: make a function
	results->total_energy += timestep_results.total_energy;
	results->total_sim_time += timestep_results.total_sim_time;
	results->total_spikes += timestep_results.total_spikes;

	// Calculate elapsed time
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	ts_elapsed = calculate_elapsed_time(ts_start, ts_end);
	results->wall_time +=
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

int sim_input(const double firing_probability)
{
	// Simulate a single external input (as one neuron) for a timestep
	//  Return 1 if the input fires, 0 otherwise
	double rand_uniform;
	int input_fired;

	rand_uniform = (double) rand() / RAND_MAX;
	input_fired = (rand_uniform < firing_probability);

	return input_fired;
}

void init_results(struct sim_results *results)
{
	results->total_energy = 0.0; // Joules
	results->total_sim_time = 0.0; // Seconds
	results->wall_time = 0.0; // Seconds
	results->time_steps = 0;
	results->total_spikes = 0;
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

int parse_dvs(FILE *fp, struct input *inputs, const int max_inputs)
{
	// TODO
	return 0;
}