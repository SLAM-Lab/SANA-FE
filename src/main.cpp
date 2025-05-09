// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// main.cpp - Command line interface
// Performance simulation of neuromorphic architectures
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem> // For std::filesystem::path
#include <optional>
#include <string>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "arch.hpp"
#include "description.hpp"
#include "models.hpp"
#include "network.hpp"
#include "print.hpp"
#include "chip.hpp"

int main(int argc, char *argv[])
{
    if (argc < 1)
    {
        INFO("Error: No program arguments.");
        return 1;
    }
    // First arg is always program name, skip
    argc--;
    argv++;

    // Parse optional args
    std::filesystem::path output_dir = std::filesystem::current_path();
    bool record_spikes{false};
    bool record_potentials{false};
    bool record_perf{false};
    bool record_messages{false};
    bool use_netlist_format{false};
    bool use_simple_timing_model{false}; // Default is the detailed model

    while (argc > 2)
    {
        if (argv[0][0] == '-')
        {
            switch (argv[0][1])
            {
            case 'o': {
                argc--;
                argv++;
                output_dir = std::filesystem::path(argv[0]);
                // TODO: fix this for C++
                INFO("Writing output to %s\n", output_dir.c_str());
                break;
            }

            case 'm':
                record_messages = true;
                break;
            case 'n':
                use_netlist_format = true;
            case 'p':
                record_perf = true;
                break;
            case 's':
                record_spikes = true;
                break;
            case 't':
                use_simple_timing_model = true;
            case 'v':
                record_potentials = true;
                break;

            default:
                INFO("Error: Flag %c not recognized.\n", argv[0][1]);
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

    if (argc < sanafe::PROGRAM_NARGS)
    {
        INFO("Usage: ./sim [-psvmo] <arch description> <network description> <timesteps>\n");
        return 0;
    }

    // Read in program args, sanity check and parse inputs
    try
    {
#ifndef GIT_COMMIT
#define GIT_COMMIT "git-hash-unknown"
#endif
        INFO("Running SANA-FE simulation (build:%s)\n", GIT_COMMIT);
#ifdef HAVE_OPENMP
        const int nthreads = omp_get_num_procs();
        INFO("OpenMP enabled, %d threads detected\n", nthreads);
#else
        INFO("No OpenMP multithreading enabled.\n");
#endif
        sanafe::Architecture arch =
                sanafe::load_arch(argv[sanafe::ARCH_FILENAME]);
        INFO("Architecture initialized.\n");
        sanafe::SpikingNetwork net = sanafe::load_net(
                argv[sanafe::NETWORK_FILENAME], arch, use_netlist_format);
        INFO("Network initialized.\n");

        sanafe::SpikingChip hw(arch, output_dir, record_spikes,
                record_potentials, record_perf, record_messages);
        hw.load(net);

        long int timesteps;
        try
        {
            timesteps = std::stol(argv[sanafe::TIMESTEPS]);
        }
        catch(const std::exception& e)
        {
            INFO("Error: Invalid time-step format: %s\n",
                    argv[sanafe::TIMESTEPS]);
            return 1;
        }

        if (timesteps <= 0)
        {
            INFO("Error: Time-steps must be > 0 (%ld)\n", timesteps);
            return 1;
        }

        const long int heartbeat = 100L;
        sanafe::TimingModel timing_model = sanafe::TIMING_MODEL_DETAILED;
        if (use_simple_timing_model)
        {
            timing_model = sanafe::TIMING_MODEL_SIMPLE;
        }

        INFO("Running simulation.\n");
        hw.sim(timesteps, heartbeat, timing_model);

        INFO("***** Run Summary *****\n");
        hw.sim_output_run_summary(output_dir);
        double average_power = hw.get_power();
        INFO("Average power consumption: %f W.\n", average_power);
        INFO("Run finished.\n");

        return 0;
    }
    catch (const sanafe::DescriptionParsingError &exc)
    {
        INFO("%s", exc.what());
        return 1;
    }
    catch (const std::runtime_error &exc)
    {
        INFO("Error: runtime exception thrown: %s\n", exc.what());
        return 1;
    }
    catch (const std::invalid_argument &exc)
    {
        INFO("Error: invalid argument thrown: %s\n", exc.what());
        return 1;
    }
}

// Project TODOs and wishlist roughly in priority order
//
// ** New simulator features **
// TODO: integrate Booksim2 as a 'cycle' accurate modeling option
// TODO: support saving SNNs to yaml and netlist
// TODO: extend the old netlist format to support embedded YAML, group mappings, and hyperedges
// ** Misc **
// TODO: update tutorials to use Jupyter notebook
// TODO: Improve plug-in documentation
