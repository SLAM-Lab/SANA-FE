#include "arg_parsing.hpp"
#include "chip.hpp"
#include "print.hpp"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#if HAVE_OPENMP
#include <omp.h>
#endif

namespace
{

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

} // anonymous namespace

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
