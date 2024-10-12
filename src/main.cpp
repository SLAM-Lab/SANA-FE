// Copyright (c) 2024 - The University of Texas at Austin
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

#include "arch.hpp"
#include "description.hpp"
#include "models.hpp"
#include "network.hpp"
#include "print.hpp"
#include "sim.hpp"

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
                ;
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
        sanafe::Architecture arch =
                sanafe::load_arch(argv[sanafe::ARCH_FILENAME]);
        INFO("Architecture initialized.\n");
        sanafe::SpikingNetwork net = sanafe::load_net(
                argv[sanafe::NETWORK_FILENAME], arch, use_netlist_format);
        INFO("Network initialized.\n");

        sanafe::SpikingHardware hw(arch, output_dir, record_spikes,
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

        // Step simulation
        INFO("Running simulation.\n");
        const auto run_data = hw.sim(timesteps);

        INFO("***** Run Summary *****\n");
        sim_output_run_summary(output_dir, run_data);
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
// ** New features **
// TODO: update tutorial to use Jupyter notebook
// TODO: add models.reset() so we can clear model state
// TODO: implement pipeline h/w units in pybind e.g., soma, dendrite
//
// ** Refactor and optimizing existing features **
// TODO: Parse model parameters after mapping, in two-pass approach
// TODO: Make the model parameter maps optional, revert to dendrite/soma/shared
//  maps