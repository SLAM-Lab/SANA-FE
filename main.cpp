// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// main.c - Command line interface
// Performance simulation of neuromorphic architectures
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "print.hpp"
#include "sim.hpp"

using namespace sanafe;

int main(int argc, char *argv[])
{
	long int timesteps;
	int ret;

	// Assume that if we don't get to the point where we write this with
	//  a valid value, something went wrong and we errored out
	ret = RET_FAIL;

	if (argc < 1)
	{
		throw std::invalid_argument("Error: No program arguments.");
	}
	// First arg is always program name, skip
	argc--;
	argv++;

	// Parse optional args
	std::string output_dir(".");
	bool record_spikes, record_potentials, record_perf, record_messages;
	record_spikes = false;
	record_potentials = false;
	record_perf = false;
	record_messages = false;
	while (argc > 2)
	{
		if (argv[0][0] == '-')
		{
			switch (argv[0][1])
			{
			case 'o':
			{
				if (argc <= 0)
				{
					throw std::invalid_argument(
					"Error: No output dir given.\n");
				}
				argc--;
				argv++;
				std::string out_dir(argv[0]);
				// TODO: fix this for C++
				INFO("Writing output to %s\n", out_dir.c_str());;
				break;
			}
			case 's':
				record_spikes = true;
				break;
			case 'v':
				record_potentials = true;
				break;
			case 'p':
				record_perf = true;
				break;
			case 'm':
				record_messages = true;
				break;
			default:
				INFO("Error: Flag %c not recognized.\n",
								argv[0][1]);
				break;
			}
			argc--;
			argv++;
		}
		else
		{
			break;
		}
	}

	if (argc < PROGRAM_NARGS)
	{
		INFO("Usage: ./sim [-p<log perf> -s<spike trace> "
				"-v<potential trace> -m <message trace>] "
				"-o <output directory>"
				"<arch description> <network description> "
							"<timesteps>\n");
		return 0;
	}
	Simulation sim(output_dir, record_spikes, record_potentials,
		record_perf, record_messages);

	// Read in program args, sanity check and parse inputs
	sim.read_arch_file(argv[ARCH_FILENAME]);
	sim.read_net_file(argv[NETWORK_FILENAME]);

	timesteps = 0;
	ret = sscanf(argv[TIMESTEPS], "%ld", &timesteps);
	if (ret < 1)
	{
		INFO("Error: Time-steps must be integer > 0 (%s).\n",
							argv[TIMESTEPS]);
		return 1;
	}
	else if (timesteps <= 0)
	{
		INFO("Error: Time-steps must be > 0 (%ld)\n", timesteps);
		return 1;
	}

	// Step simulation
	INFO("Running simulation.\n");
	sim.run(timesteps);

	INFO("***** Run Summary *****\n");
	sim.get_run_summary(output_dir);
	double average_power = sim.get_power();
	INFO("Average power consumption: %f W.\n", average_power);
	INFO("Run finished.\n");

	return 0;
}

/*
void run(struct simulation *sim, struct network *net, Architecture *arch)
{
	// TODO: remove the need to pass the network struct, only the arch
	//  should be needed (since it links back to the net anyway)
	// Run neuromorphic hardware simulation for one timestep
	//  Measure the CPU time it takes and accumulate the stats
	struct timespec ts_start, ts_end, ts_elapsed;
	struct timestep *ts = &(sim->ts);

	// Measure the wall-clock time taken to run the simulation
	//  on the host machine
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	sim_timestep(ts, net, arch);
	// Calculate elapsed time
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

	sim->total_energy += ts->energy;
	sim->total_sim_time += ts->sim_time;
	sim->total_spikes += ts->spike_count;
	sim->total_neurons_fired += ts->total_neurons_fired;
	sim->total_messages_sent += ts->packets_sent;
	if (sim->log_spikes)
	{
		sim_trace_record_spikes(sim, net);
	}
	if (sim->log_potential)
	{
		sim_trace_record_potentials(sim, net);
	}
	if (sim->log_perf)
	{
		sim_perf_log_timestep(ts, sim->perf_fp);
	}
	if (sim->log_messages)
	{
		for (int i = 0; i < ARCH_MAX_CORES; i++)
		{
			for (int j = 0; j < ts->message_counts[i]; j++)
			{
				Message *m = &(ts->messages[i][j]);
				sim_trace_record_message(sim, m);
			}
		}
	}
	sim->timesteps = ts->timestep;
	sim->wall_time += (double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9);
	TRACE1("Time-step took: %fs.\n",
		(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9));

}
*/

struct timespec sanafe::calculate_elapsed_time(struct timespec ts_start,
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
