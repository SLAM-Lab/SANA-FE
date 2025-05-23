// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// main.cpp - Command line interface
// Performance simulation of neuromorphic architectures
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem> // For std::filesystem::path
#include <stdexcept>
#include <string>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <booksim_lib.hpp>

#include "arch.hpp"
#include "description.hpp"
#include "network.hpp"
#include "print.hpp"
#include "chip.hpp"

enum ProgramArgs : uint8_t
{
    ARCH_FILENAME = 0,
    NETWORK_FILENAME = 1,
    TIMESTEPS = 2,
    PROGRAM_NARGS = 3,
};

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
    int total_threads_available{1}; // NOLINT(misc-const-correctness)
#ifdef HAVE_OPENMP
    int processing_threads{1};
#endif
    int scheduler_threads{0};
    // Select the timing model on the command line. The default is the
    //  detailed build-in timing model (using a scheduler). Select only one of
    //  these two flags to enable either the simple analytical timing model or
    //  an external cycle-accurate timing model (Booksim2).
    sanafe::TimingModel timing_model = sanafe::TIMING_MODEL_DETAILED;
    std::string timing_model_str = "detailed";

    #ifdef HAVE_OPENMP
        total_threads_available = omp_get_num_procs();
        INFO("OpenMP enabled, %d threads detected\n", total_threads_available);
#else
        INFO("No OpenMP multithreading enabled.\n");
#endif

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
                INFO("Writing output to %s\n", output_dir.c_str());
                break;
            }
            case 'm':
                record_messages = true;
                break;
            case 'n':
                use_netlist_format = true;
                break;
            case 'p':
                record_perf = true;
                break;
            case 's':
                record_spikes = true;
                break;
            case 't':
                argc--;
                argv++;
                timing_model_str = std::string(argv[0]);
                break;
            case 'v':
                record_potentials = true;
                break;
            case 'N':
                argc--;
                argv++;
#ifndef HAVE_OPENMP
                INFO("Warning: multiple threads not supported; flag ignored");
#else
                processing_threads = std::min(total_threads_available, std::stoi(argv[0]));
                INFO("Setting processing threads to %d\n", processing_threads);
                omp_set_num_threads(processing_threads);
#endif
                break;
            case 'S':
                argc--;
                argv++;
                scheduler_threads =
                        std::min(total_threads_available, std::stoi(argv[0]));
                INFO("Setting scheduling threads to %d\n", scheduler_threads);

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

    if (argc < PROGRAM_NARGS)
    {
        INFO("Usage: ./sim [-psvmo] <arch description> <network description> <timesteps>\n");
        return 0;
    }

    if (timing_model_str == "simple")
    {
        timing_model = sanafe::TIMING_MODEL_SIMPLE;
    }
    else if (timing_model_str == "detailed")
    {
        timing_model = sanafe::TIMING_MODEL_DETAILED;

    }
    else if (timing_model_str == "cycle")
    {
        timing_model = sanafe::TIMING_MODEL_CYCLE_ACCURATE;
    }
    else
    {
        INFO("Error: Timing model %s not recognized.\n",
                timing_model_str.c_str());
    }

    // Read in program args, sanity check and parse inputs
    try
    {
#ifndef GIT_COMMIT
#define GIT_COMMIT "git-hash-unknown"
#endif
        INFO("Running SANA-FE simulation (build:%s)\n", GIT_COMMIT);
        INFO("Loading booksim2 library for cycle-accurate support\n");
        booksim_init();

        sanafe::Architecture arch =
                sanafe::load_arch(argv[ARCH_FILENAME]);
        INFO("Architecture initialized.\n");
        const sanafe::SpikingNetwork net = sanafe::load_net(
                argv[NETWORK_FILENAME], arch, use_netlist_format);
        INFO("Network initialized.\n");

        sanafe::SpikingChip hw(arch, output_dir, record_spikes,
                record_potentials, record_perf, record_messages);
        hw.load(net);

        long int timesteps{0L};
        try
        {
            timesteps = std::stol(argv[TIMESTEPS]);
        }
        catch(const std::exception& e)
        {
            INFO("Error: Invalid time-step format: %s\n",
                    argv[TIMESTEPS]);
            return 1;
        }

        if (timesteps <= 0)
        {
            INFO("Error: Time-steps must be > 0 (%ld)\n", timesteps);
            return 1;
        }

        INFO("Running simulation.\n");
        const sanafe::RunData run_summary =
                hw.sim(timesteps, timing_model, scheduler_threads);

        INFO("Closing Booksim2 library\n");
        booksim_close();

        INFO("***** Run Summary *****\n");
        hw.sim_output_run_summary(output_dir, run_summary);
        const double average_power = hw.get_power();
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
// TODO: Add support for Fugu and Lava frameworks as extra (optional) Python
//  dependencies
