// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// main.cpp - Command line interface
// Performance simulation of neuromorphic architectures
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem> // For std::filesystem::path
#include <stdexcept>
#include <string>
#include <vector>

#if HAVE_OPENMP
#include <omp.h>
#endif
#include <booksim_lib.hpp>

#include "arch.hpp"
#include "arg_parsing.hpp"
#include "chip.hpp"
#include "network.hpp"
#include "print.hpp"
#include "yaml_common.hpp"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,readability-function-size)
int main(int argc, const char *argv[])
{
    const std::vector<std::string> arg_vec = program_args_to_vector(argc, argv);
    const OptionalProgramFlags optional_flags =
            parse_command_line_flags(arg_vec);
    argc -= optional_flags.total_args_parsed;
    if (argc < RequiredProgramArgs::program_nargs)
    {
        INFO("Usage: ./sim [-psvmo] <arch description> <network description> "
             "<timesteps>\n");
        return 0;
    }

    // Read in program args, sanity check and parse inputs
    try
    {
#ifndef GIT_COMMIT
#define GIT_COMMIT "git-hash-unknown"
#endif
        INFO("Running SANA-FE simulation (build:%s)\n", GIT_COMMIT);
        booksim_init();

        const RequiredProgramArgs required_args =
                parse_required_args(arg_vec, optional_flags.total_args_parsed);
        sanafe::Architecture arch =
                sanafe::load_arch(required_args.arch_filename);
        INFO("Architecture initialized.\n");
        const sanafe::SpikingNetwork net =
                sanafe::load_net(required_args.network_filename, arch,
                        optional_flags.use_netlist_format);
        INFO("Network initialized.\n");

        sanafe::SpikingChip hw(arch);
        hw.load(net);

        INFO("Running simulation.\n");
        const sanafe::RunData run_summary = hw.sim(
                required_args.timesteps_to_execute, optional_flags.timing_model,
                optional_flags.scheduler_threads, optional_flags.record_spikes,
                optional_flags.record_potentials, optional_flags.record_perf,
                optional_flags.record_messages, optional_flags.output_dir);

        INFO("Closing Booksim2 library\n");
        booksim_close();

        INFO("***** Run Summary *****\n");
        hw.sim_output_run_summary(optional_flags.output_dir, run_summary);
        [[maybe_unused]] const double average_power = hw.get_power();
        INFO("Average power consumption: %f W.\n", average_power);
        INFO("Run finished.\n");

        return 0;
    }
    catch (const sanafe::YamlDescriptionParsingError &exc)
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

// SANA-FE tasklist roughly in priority order
//
// * Finish & integrate backends for Fugu, Lava, NIR, and SNNTorch
// * Support different multithreaded schedulers and clean-up interface
// * Implement multi-threaded BookSim 2
// * Create and push new docker images to DockerHub
