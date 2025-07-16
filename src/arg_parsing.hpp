#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "chip.hpp"  // For sanafe::TimingModel and parse_timing_model

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

// Stores required arguments: architecture file, network file, and timestep count
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

// Parses optional command-line flags like -p, -s, -o, etc.
OptionalProgramFlags parse_command_line_flags(const std::vector<std::string>& args);

// Parses required positional arguments (architecture file, network file, timesteps)
RequiredProgramArgs parse_required_args(const std::vector<std::string> &args, const int optional_flags);

// Converts (argc, argv) to std::vector<std::string>
std::vector<std::string> program_args_to_vector(const int argc, const char *argv[]);