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
#include "module.hpp"

int main(int argc, char *argv[])
{
	SANA_FE sana_fe;
	int timesteps, ret;

	// Assume that if we don't get to the point where we write this with
	//  a valid value, something went wrong and we errored out
	ret = RET_FAIL;

	if (argc < 1)
	{
		INFO("Error: No program arguments.\n");
		sana_fe.clean_up(RET_FAIL);
	}
	// First arg is always program name, skip
	argc--;
	argv++;

	// Parse optional args
	while (argc > 2)
	{
		if (argv[0][0] == '-')
		{
			switch (argv[0][1])
			{
			case 'i':
				sana_fe.set_input(argv[1]);
				argv++;
				argc--;
				break;
			case 'p':
				sana_fe.set_perf_flag();
				break;
			case 's':
				sana_fe.set_spike_flag();
				break;
			case 'v':
				sana_fe.set_pot_flag();
				break;
			case 'm':
				sana_fe.set_mess_flag();
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
				"<arch description> <network description> "
							"<timesteps>\n");
		sana_fe.clean_up(RET_FAIL);
		goto clean_up;
	}

	if (sim->log_potential)
	{
		sim->potential_trace_fp = fopen("potential.csv", "w");
		if (sim->potential_trace_fp == NULL)
		{
			INFO("Error: Couldn't open trace file for writing.\n");
			goto clean_up;
		}
	}
	if (sim->log_spikes)
	{
		sim->spike_trace_fp = fopen("spikes.csv", "w");
		if (sim->spike_trace_fp == NULL)
		{
			INFO("Error: Couldn't open trace file for writing.\n");
			goto clean_up;
		}
	}
	if (sim->log_messages)
	{
		sim->message_trace_fp = fopen("messages.csv", "w");
		if (sim->message_trace_fp == NULL)
		{
			INFO("Error: Couldn't open trace file for writing.\n");
			goto clean_up;
		}

	}
	if (sim->log_perf)
	{
		sim->perf_fp = fopen("perf.csv", "w");
		if (sim->perf_fp == NULL)
		{
			INFO("Error: Couldn't open perf file for writing.\n");
			sana_fe.clean_up(RET_FAIL);
			goto clean_up;
		}
	}

	// Read in program args, sanity check and parse inputs
	sana_fe.set_arch(argv[ARCH_FILENAME]);
	sana_fe.set_net(argv[NETWORK_FILENAME]);

	timesteps = 0;
	ret = sscanf(argv[TIMESTEPS], "%d", &timesteps);
	if (ret < 1)
	{
		INFO("Error: Time-steps must be integer > 0 (%s).\n",
							argv[TIMESTEPS]);
		sana_fe.clean_up(RET_FAIL);
	}
	else if (timesteps <= 0)
	{
		INFO("Error: Time-steps must be > 0 (%d)\n", timesteps);
		sana_fe.clean_up(RET_FAIL);
	}

	// Step simulation
	INFO("Running simulation.\n");
	for (long timestep = 1; timestep <= timesteps; timestep++)
	{
		if ((timestep % 100) == 0)
		{
			// Print heart-beat every hundred timesteps
			INFO("*** Time-step %ld ***\n", timestep);
		}
		sana_fe.run_timesteps();
	}

	INFO("***** Run Summary *****\n");
	sana_fe.sim_summary();
	double average_power = sana_fe.get_power();
	INFO("Average power consumption: %f W.\n", average_power);
	INFO("Run finished.\n");

<<<<<<< HEAD
clean_up:
	// Free any larger structures here
	network_free(&net);
	arch_free(arch);
	// Free any locally allocated memory here
	free(input_buffer);

	// Close any open files here
	if (sim->potential_trace_fp != NULL)
	{
		fclose(sim->potential_trace_fp);
	}
	if (sim->spike_trace_fp != NULL)
	{
		fclose(sim->spike_trace_fp);
	}
	if (sim->message_trace_fp != NULL)
	{
		fclose(sim->message_trace_fp);
	}
	if (sim->perf_fp != NULL)
	{
		fclose(sim->perf_fp);
	}
	if (sim->stats_fp != NULL)
	{
		fclose(sim->stats_fp);
	}

	// Free the simulation structure only after we close all files
	free(sim);

	if (ret == RET_FAIL)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void run(struct simulation *sim, struct network *net, struct architecture *arch)
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
				struct message *m = &(ts->messages[i][j]);
				sim_trace_record_message(sim, m);
			}
		}
	}
	sim->timesteps = ts->timestep;
	sim->wall_time += (double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9);
	TRACE1("Time-step took: %fs.\n",
		(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9));

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
=======
	sana_fe.clean_up(RET_OK);
>>>>>>> Switching to class-based execution
}
