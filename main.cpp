// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// main.c - Command line interface
// Performance simulation of neuromorphic architectures
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <filesystem>

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
	std::filesystem::path output_dir = std::filesystem::current_path();
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
				output_dir = std::filesystem::path(argv[0]);
				// TODO: fix this for C++
				INFO("Writing output to %s\n",
					output_dir.c_str());;
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

	// Read in program args, sanity check and parse inputs
	Architecture arch;
	arch.load_arch_file(argv[ARCH_FILENAME]);
	Network net;
	net.load_net_file(argv[NETWORK_FILENAME], arch);
	Simulation sim(
		arch, net, output_dir, record_spikes, record_potentials,
		record_perf, record_messages);

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
	const auto run_data = sim.get_run_summary();
	sim_output_run_summary(output_dir, run_data);
	double average_power = sim.get_power();
	INFO("Average power consumption: %f W.\n", average_power);
	INFO("Run finished.\n");

	return 0;
}
