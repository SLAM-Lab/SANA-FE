// main.c
// Performance simulation for neuromorphic architectures
#include <stdio.h>
#include <stdlib.h>

#include "sim.h"
#include "network.h"

enum program_args
{
	PROGRAM_NAME = 0,
	TIMESTEPS,
    N_CORES,
	NEURON_CONFIG_FILENAME,
	PROGRAM_NARGS,
};

int main(int argc, char *argv[])
{
	FILE *fp;
	struct sim_results results;
	struct core *cores;
	char *filename;
	int timesteps, max_cores;

	if (argc < PROGRAM_NARGS)
	{
		INFO("Usage: ./sim <timesteps> <cores> <neuron list>\n");
		exit(1);
	}

	// Read in program args, sanity check and parse inputs
	sscanf(argv[TIMESTEPS], "%d", &timesteps);
	printf(argv[TIMESTEPS]);
	if (timesteps <= 0)
	{
		INFO("Time-steps must be > 0 (%d)\n", timesteps);
		exit(1);
	}

    sscanf(argv[N_CORES], "%d", &max_cores);
    if ((max_cores <= 0) || (max_cores > MAX_CORES_LOIHI))
    {
        INFO("Cores must be > 0 and < %d (%d)\n", MAX_CORES_LOIHI, max_cores);
        exit(1);
    }

	filename = argv[NEURON_CONFIG_FILENAME];

	// Read in the configuration of all neurons in the network
	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		INFO("Neuron data (%s) failed to open.\n", filename);
		exit(1);
	}

	INFO("Allocating memory for %d cores.\n", max_cores);
	cores = (struct core *) malloc(max_cores * sizeof(struct core));
	if (cores == NULL)
	{
		INFO("Error: failed to allocate memory.\n");
		exit(1);
	}
	INFO("Allocated %ld bytes\n", max_cores * sizeof(struct core));

	sim_init_cores(cores, max_cores);
	network_read_csv(fp, cores);
	fclose(fp);

	results = sim_run(timesteps, cores, max_cores);
	fp = fopen("results.yaml", "w");
	if (fp != NULL)
	{
		sim_write_results(fp, &results);
	}
	fclose(fp);
        free(cores);

	return 0;
}

