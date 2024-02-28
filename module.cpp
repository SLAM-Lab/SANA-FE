// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// module.cpp - Pybind module interface
// Performance simulation of neuromorphic architectures
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "print.hpp"
#include "sim.hpp"
#include "network.hpp"
#include "arch.hpp"
#include "description.hpp"
#include "command.hpp"
#include "main.hpp"
#include "pybind11/pybind11.h"

class SANA_FE{
	public:
		struct simulation *sim;
		struct network net;
		struct architecture *arch;
		double average_power;
		int timesteps;
		FILE *input_fp;

		SANA_FE(){
			init();
		}
		void init(){
			arch = arch_init();
			network_init(&net);
			INFO("Initializing simulation.\n");
			sim = sim_init_sim();
		}
		int forceSpike(int group_id, int n_id, int num_spikes){
			if (num_spikes < 0)
				return -1;
			if (group_id >= net.neuron_group_count)
				return -1;
			if (n_id >= net.groups[group_id].neuron_count)
				return -1;
			net.groups[group_id].neurons[n_id].forced_spikes = num_spikes;
			return num_spikes;
		}
		void run_timesteps(int timesteps = 1){
			for (int i = 0; i < timesteps; ++i){
				run(sim, &net, arch);
			}
		}

		// Temporary functions
		void parseArgs(int argc, char *argv[]){
			char *filename;
			int ret;
			FILE *network_fp, *arch_fp;

			if (argc < 1)
			{
				INFO("Error: No program arguments.\n");
				clean_up(RET_FAIL);
			}

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
						filename = argv[1];
						argv++;
						argc--;
						break;
					case 'p':
						sim->log_perf = 1;
						break;
					case 's':
						sim->log_spikes = 1;
						break;
					case 'v':
						sim->log_potential = 1;
						break;
					case 'm':
						sim->log_messages = 1;
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

			if (filename)
			{
				input_fp = fopen(filename, "r");
				if (input_fp == NULL)
				{
					INFO("Error: Couldn't open inputs %s.\n", filename);
					clean_up(RET_FAIL);
				}
			}

			if (argc < PROGRAM_NARGS)
			{
				INFO("Usage: ./sim [-p<log perf> -s<spike trace> "
						"-v<potential trace> -i <input vectors>] "
						"<arch description> <network description> "
									"<timesteps>\n");
				clean_up(RET_FAIL);
			}

			if (sim->log_potential)
			{
				sim->potential_trace_fp = fopen("potential.trace", "w");
				if (sim->potential_trace_fp == NULL)
				{
					INFO("Error: Couldn't open trace file for writing.\n");
					clean_up(RET_FAIL);
				}
			}
			if (sim->log_spikes)
			{
				sim->spike_trace_fp = fopen("spikes.trace", "w");
				if (sim->spike_trace_fp == NULL)
				{
					INFO("Error: Couldn't open trace file for writing.\n");
					clean_up(RET_FAIL);
				}
			}
			if (sim->log_messages)
			{
				sim->message_trace_fp = fopen("messages.trace", "w");
				if (sim->message_trace_fp == NULL)
				{
					INFO("Error: Couldn't open trace file for writing.\n");
					clean_up(RET_FAIL);
				}

			}
			if (sim->log_perf)
			{
				sim->perf_fp = fopen("perf.csv", "w");
				if (sim->perf_fp == NULL)
				{
					INFO("Error: Couldn't open perf file for writing.\n");
					clean_up(RET_FAIL);
				}
			}

			// Read in program args, sanity check and parse inputs
			filename = argv[ARCH_FILENAME];
			arch_fp = fopen(filename, "r");
			if (arch_fp == NULL)
			{
				INFO("Error: Architecture file %s failed to open.\n", filename);
				clean_up(RET_FAIL);
			}
			ret = description_parse_file(arch_fp, NULL, arch);
			//arch_print_description(&description, 0);
			fclose(arch_fp);
			if (ret == RET_FAIL)
			{
				clean_up(RET_FAIL);
			}
			//arch_init_message_scheduler(&scheduler, arch);

			timesteps = 0;
			ret = sscanf(argv[TIMESTEPS], "%d", &timesteps);
			if (ret < 1)
			{
				INFO("Error: Time-steps must be integer > 0 (%s).\n",
									argv[TIMESTEPS]);
				clean_up(RET_FAIL);
			}
			else if (timesteps <= 0)
			{
				INFO("Error: Time-steps must be > 0 (%d)\n", timesteps);
				clean_up(RET_FAIL);
			}

			filename = argv[NETWORK_FILENAME];
			// Create the network
			network_fp = fopen(filename, "r");
			if (network_fp == NULL)
			{
				INFO("Network data (%s) failed to open.\n", filename);
				clean_up(RET_FAIL);
			}
			INFO("Reading network from file.\n");
			ret = description_parse_file(network_fp, &net, arch);
			fclose(network_fp);
			if (ret == RET_FAIL)
			{
				clean_up(RET_FAIL);
			}
			network_check_mapped(&net);

			arch_create_connection_maps(arch);
			INFO("Creating probe and perf data files.\n");
			if (sim->spike_trace_fp != NULL)
			{
				sim_spike_trace_write_header(sim);
			}
			if (sim->potential_trace_fp != NULL)
			{
				sim_potential_trace_write_header(sim, &net);
			}
			if (sim->message_trace_fp != NULL)
			{
				sim_message_trace_write_header(sim);
			}
			if (sim->perf_fp != NULL)
			{
				sim_perf_write_header(sim->perf_fp);
			}
		}

		void clean_up(description_ret ret = RET_OK){
			// Free any larger structures here
			network_free(&net);
			arch_free(arch);

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
				exit(1);
			}
			else
			{
				exit(0);
			}
		}
};

void test_pybind(void){
	INFO("Printing through Pybind!\n");
}

PYBIND11_MODULE(simcpp, m) {
    m.doc() = R"pbdoc(
        SANA-FE Cpp Module with Pybind11 
        --------------------------------

        .. currentmodule:: simcpp

        .. autosummary::
           :toctree: _generate

           test_pybind
    )pbdoc";

    m.def("test_pybind", &test_pybind, R"pbdoc(
        test_pybind function from main.cpp

        Test pybind11 functionality.
    )pbdoc");
}
