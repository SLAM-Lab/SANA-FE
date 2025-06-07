// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// main.cpp - Command line interface
// Performance simulation of neuromorphic architectures
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem> // For std::filesystem::path
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#if HAVE_OPENMP
#include <omp.h>
#endif
#include <booksim_lib.hpp>

#include "arch.hpp"
#include "chip.hpp"
#include "network.hpp"
#include "print.hpp"
#include "yaml_common.hpp"

namespace // anonymous to keep members private
{

struct OptionalProgramFlags
{
    std::filesystem::path output_dir{std::filesystem::current_path()};
    bool record_spikes{false};
    bool record_potentials{false};
    bool record_perf{false};
    bool record_messages{false};
    bool use_netlist_format{false};
    int total_args_parsed{0};
    int total_threads_available{1}; // NOLINT(misc-const-correctness)
    int processing_threads{1};
    int scheduler_threads{0};
    // Select the timing model on the command line. The default is the
    //  detailed build-in timing model (using a scheduler). Select only one of
    //  these two flags to enable either the simple analytical timing model or
    //  an external cycle-accurate timing model (Booksim2).
    sanafe::TimingModel timing_model = sanafe::timing_model_detailed;
};

struct RequiredProgramArgs
{
    enum ArgIdx : uint8_t
    {
        // Ignoring optional flags
        arch_filename_idx = 0,
        network_filename_idx = 1,
        timesteps_to_execute_idx = 2,
        program_nargs = 3,
    };
    std::string arch_filename;
    std::string network_filename;
    long int timesteps_to_execute{0L};
};

std::string_view get_next_arg(
        const std::vector<std::string> &args, const size_t current_idx)
{
    // Helper function to get next argument safely
    if (current_idx + 1 >= args.size())
    {
        throw std::runtime_error("Flag requires an argument but none provided");
    }
    return args[current_idx + 1];
}

// NOLINTNEXTLINE(readability-function-size)
int parse_flag(const std::vector<std::string> args, const size_t current_idx,
        OptionalProgramFlags &flags)
{
    const std::string_view arg = args[current_idx];
    if (arg.empty() || arg[0] != '-')
    {
        INFO("Invalid flag, should start with '-'");
        return 0;
    }
    if (arg.length() < 2)
    {
        INFO("Invalid flag, too short");
        return 0;
    }

    int args_consumed = 1;
    // Ignore the first character (arg[0]), which we know is a dash. Parse the
    //  character and fields immediately after
    switch (arg[1])
    {
    case 'o': {
        args_consumed = 2;
        flags.output_dir =
                std::filesystem::path(get_next_arg(args, current_idx));
        INFO("Writing output to %s\n", flags.output_dir.c_str());
        break;
    }
    case 'm':
        flags.record_messages = true;
        break;
    case 'n':
        flags.use_netlist_format = true;
        break;
    case 'p':
        flags.record_perf = true;
        break;
    case 's':
        flags.record_spikes = true;
        break;
    case 't':
        args_consumed = 2;
        flags.timing_model =
                sanafe::parse_timing_model(get_next_arg(args, current_idx));
        break;
    case 'v':
        flags.record_potentials = true;
        break;
    case 'N':
        args_consumed = 2;
#ifndef HAVE_OPENMP
        INFO("Warning: multiple threads not supported; flag ignored");
        INFO("Processing threads:%d (default)\n", flags.processing_threads);
        get_next_arg(
                args, current_idx); // consume the argument even if we ignore it
#else
        flags.processing_threads = std::min(flags.total_threads_available,
                std::stoi(std::string(get_next_arg(args, current_idx))));
        INFO("Setting processing threads to %d\n", flags.processing_threads);
        omp_set_num_threads(flags.processing_threads);
#endif
        break;
    case 'S':
        args_consumed = 2;
        flags.scheduler_threads = std::min(flags.total_threads_available,
                std::stoi(std::string(get_next_arg(args, current_idx))));
        INFO("Setting scheduling threads to %d\n", flags.scheduler_threads);

        break;

    default:
        INFO("Error: Flag %c not recognized.\n", arg[1]);
        break;
    }

    return args_consumed;
}

OptionalProgramFlags parse_command_line_flags(
        const std::vector<std::string> &args)
{
    OptionalProgramFlags flags;

    // Record the number of available threads to OpenMP, later we need this
    //  when configuring how many processing threads to use
#if HAVE_OPENMP
    flags.total_threads_available = omp_get_num_procs();
    INFO("OpenMP enabled, %d threads detected\n",
            flags.total_threads_available);
#else
    INFO("No OpenMP multithreading enabled.\n");
#endif

    // Should always be three arguments from the end
    size_t current_idx = 0;
    while (current_idx < args.size() &&
            (args.size() - current_idx) > (RequiredProgramArgs::program_nargs))
    {
        if (args[current_idx].empty() || args[current_idx][0] != '-')
        {
            break;
        }
        try
        {
            const int args_consumed = parse_flag(args, current_idx, flags);
            if (args_consumed == 0)
            {
                break; // Error while parsing flag
            }
            current_idx += args_consumed;
            flags.total_args_parsed += args_consumed;
        }
        catch (const std::exception &exc)
        {
            INFO("Error parsing flag: %s\n", exc.what());
            break;
        }
    }

    return flags;
}

RequiredProgramArgs parse_required_args(
        const std::vector<std::string> &args, const int optional_flags)
{
    RequiredProgramArgs required_args;

    const auto base_idx = static_cast<size_t>(optional_flags);
    if (base_idx + RequiredProgramArgs::program_nargs > args.size())
    {
        throw std::invalid_argument("Insufficient arguments provided");
    }

    required_args.arch_filename =
            args[optional_flags + RequiredProgramArgs::arch_filename_idx];
    required_args.network_filename =
            args[optional_flags + RequiredProgramArgs::network_filename_idx];

    try // Catch exceptions when converting from string to long
    {
        const auto timestep_arg = std::string(args[optional_flags +
                RequiredProgramArgs::timesteps_to_execute_idx]);
        required_args.timesteps_to_execute = std::stol(timestep_arg);
    }
    catch (const std::exception &e)
    {
        throw std::invalid_argument("Error: Invalid time-step arg:");
    }
    if (required_args.timesteps_to_execute <= 0)
    {
        throw std::invalid_argument("Error: Time-steps must be > 0");
    }

    return required_args;
}

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
std::vector<std::string> program_args_to_vector(
        const int argc, const char *argv[])
// NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
{
    std::vector<std::string> args;

    // Ignore the first argument, which should just be the program name 'sim'
    for (int i = 1; i < argc; i++)
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        args.emplace_back(argv[i]);
    }

    return args;
}

}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays,readability-function-size)
int main(int argc, const char *argv[])
{
    const std::vector<std::string> arg_vec =
            program_args_to_vector(argc, argv);
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

// Project TODOs and wishlist roughly in priority order
//
// ** New simulator features **
// TODO: Add support for Fugu and Lava frameworks as extra (optional) Python
//  dependencies
