// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  sim.cpp
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <list>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <vector>

#include "arch.hpp"
#include "models.hpp"
#include "network.hpp"
#include "hardware.hpp"
#include "pipeline.hpp"
#include "plugins.hpp"
#include "print.hpp"
#include "schedule.hpp"
#include "sim.hpp"

sanafe::SpikingHardware::SpikingHardware(const Architecture &arch,
        const std::filesystem::path &output_dir, const bool record_spikes,
        const bool record_potentials, const bool record_perf,
        const bool record_messages)
        : core_count(arch.core_count)
        , noc_width(arch.noc_width)
        , noc_height(arch.noc_height)
        , noc_buffer_size(arch.noc_buffer_size)
        , max_cores_per_tile(arch.max_cores_per_tile)
        , out_dir(output_dir)
        , spike_trace_enabled(record_spikes)
        , potential_trace_enabled(record_potentials)
        , perf_trace_enabled(record_perf)
        , message_trace_enabled(record_messages)
{
    INFO("Initializing simulation.\n");
    for (const TileConfiguration &tile_config : arch.tiles)
    {
        tiles.emplace_back(tile_config);
        Tile &hardware_tile = tiles.back();
        for (const CoreConfiguration &core_config : tile_config.cores)
        {
            hardware_tile.cores.emplace_back(core_config);
            Core &hardware_core = hardware_tile.cores.back();

            for (const AxonInConfiguration &axon_config : core_config.axon_in)
            {
                hardware_core.create_axon_in(axon_config);
            }
            for (const SynapseConfiguration &synapse_config :
                    core_config.synapses)
            {
                hardware_core.create_synapse(synapse_config);
            }
            for (const DendriteConfiguration &dendrite_config :
                    core_config.dendrites)
            {
                hardware_core.create_dendrite(dendrite_config);
            }
            for (const SomaConfiguration &soma_config : core_config.somas)
            {
                hardware_core.create_soma(soma_config);
            }
            for (const AxonOutConfiguration &axon_config : core_config.axon_out)
            {
                hardware_core.create_axon_out(axon_config);
            }
        }
    }
}

sanafe::SpikingHardware::~SpikingHardware()
{
    // Close any open trace files
    spike_trace.close();
    potential_trace.close();
    perf_trace.close();
    message_trace.close();
}

void sanafe::SpikingHardware::load(const SpikingNetwork &net)
{
    map_neurons(net);
    map_connections(net);
}

void sanafe::SpikingHardware::map_neurons(const SpikingNetwork &net)
{
    auto list_of_cores = cores();
    // 1) map neurons to core
    // 2) set core execution order and attributes
    for (const auto &[name, group] : net.groups)
    {
        mapped_neuron_groups[name] =
                std::vector<MappedNeuron *>(group.neurons.size(), nullptr);
        for (const Neuron &neuron : group.neurons)
        {
            Core &mapped_core = list_of_cores[neuron.core_id];
            mapped_core.map_neuron(neuron);
        }
    }

    // Sort the neurons based on their mapped priority
    for (Tile &tile : tiles)
    {
        for (Core &core : tile.cores)
        {
            // TODO: it might be easier to temporarily reconstruct the mapping
            //  order in the spiking network, and actually just map them in
            //  order, rather than this complex 3-stage approach
            // Sort by mapping order, where a lower mapping order
            //  indicates the neuron should appear earlier in the vector
            std::sort(core.neurons.begin(), core.neurons.end(),
                    [](const MappedNeuron &a, const MappedNeuron &b) {
                        return a.mapping_order < b.mapping_order;
                    });
            // Store indexes inside all mapped neurons to allow us to
            //  interface with the hardware models
            for (size_t address = 0; address < core.neurons.size(); ++address)
            {
                MappedNeuron &mapped = core.neurons[address];
                mapped.mapped_address = address;
                const NeuronGroup &group =
                        net.groups.at(mapped.parent_group_name);
                const Neuron &neuron = group.neurons[mapped.id];

                for (auto &[name, param] : group.default_neuron_config.model_parameters)
                {
                    INFO("Setting group parameter:%s\n", name.c_str());
                }

                mapped.set_attributes(group.default_neuron_config);

                // TODO: Neuron object should keep this NeuronTemplate intact
                // Set attributes for this neuron
                NeuronTemplate neuron_specific_config;
                neuron_specific_config.default_synapse_hw_name =
                        neuron.default_synapse_hw_name;
                neuron_specific_config.dendrite_hw_name =
                        neuron.dendrite_hw_name;
                neuron_specific_config.force_dendrite_update =
                        neuron.force_dendrite_update;
                neuron_specific_config.force_synapse_update =
                        neuron.force_synapse_update;
                neuron_specific_config.force_soma_update =
                        neuron.force_soma_update;
                neuron_specific_config.log_potential = neuron.log_potential;
                neuron_specific_config.log_spikes = neuron.log_spikes;
                neuron_specific_config.soma_hw_name = neuron.soma_hw_name;
                neuron_specific_config.model_parameters =
                        neuron.model_parameters;
                mapped.set_attributes(neuron_specific_config);
                mapped_neuron_groups[group.name][neuron.id] = &mapped;
                for (auto &[name, param] : group.default_neuron_config.model_parameters)
                {
                    INFO("Setting neuron parameter:%s\n", name.c_str());
                }

                INFO("Set attributes of nid:%s.%zu at address cid:%zu[%zu]\n",
                        neuron.parent_group_id.c_str(), neuron.id, core.id,
                        address);
            }
        }
    }
}

void sanafe::SpikingHardware::map_connections(const SpikingNetwork &net)
{
    for (const auto &[name, group] : net.groups)
    {
        for (const Neuron &pre_neuron : group.neurons)
        {
            for (const Connection &curr_connection : pre_neuron.edges_out)
            {
                map_connection(curr_connection);
            }
        }
    }

    // TODO: it makes sense to build the axons as we go, rather than
    //  running in another loop, refactor and restructure this
    map_axons();

    // Set each connection attribute
    // TODO: refactor, at the high-level we want to
    // For all connections
    // 1) create a new mapped connection on the hardware
    // 2) create the axons and core-core connectivity
    // 3) set the attributes
    // Do this in one loop instead of three
    for (const auto &[name, group] : net.groups)
    {
        for (size_t nid = 0; nid < group.neurons.size(); ++nid)
        {
            const Neuron &pre_neuron = group.neurons[nid];
            MappedNeuron &mapped_neuron = *(mapped_neuron_groups[name][nid]);
            for (size_t idx = 0; idx < pre_neuron.edges_out.size(); ++idx)
            {
                const Connection &con = pre_neuron.edges_out[idx];
                MappedConnection &mapped_con =
                        mapped_neuron.connections_out[idx];

                for (auto &name_value_pair : con.synapse_params)
                {
                    if (name_value_pair.second.forward_to_synapse)
                    {
                        mapped_con.synapse_hw->set_attribute(
                                mapped_con.synapse_address,
                                name_value_pair.first, name_value_pair.second);
                    }
                }
            }
        }
    }

    return;
}

sanafe::MappedConnection &sanafe::SpikingHardware::map_connection(
        const Connection &con)
{
    auto list_of_cores = cores();

    auto &pre_group = mapped_neuron_groups.at(con.pre_neuron.group_name);
    MappedNeuron &pre_neuron =
            *(pre_group[con.pre_neuron.neuron_id.value()]);

    auto &post_group = mapped_neuron_groups.at(con.post_neuron.group_name);
    MappedNeuron &post_neuron =
            *(post_group[con.post_neuron.neuron_id.value()]);

    pre_neuron.connections_out.emplace_back(pre_neuron.connections_out.size());
    MappedConnection &mapped_con = pre_neuron.connections_out.back();
    mapped_con.pre_neuron = &pre_neuron;
    mapped_con.post_neuron = &post_neuron;

    // Map to synapse hardware unit
    Core &post_core = *(post_neuron.core);
    mapped_con.synapse_hw = post_core.synapse[0].get();

    if (con.synapse_hw_name.length() > 0)
    {
        bool synapse_found = false;
        for (auto &synapse_hw : post_core.synapse)
        {
            if (con.synapse_hw_name == synapse_hw->name)
            {
                mapped_con.synapse_hw = synapse_hw.get();
                synapse_found = true;
            }
        }
        if (!synapse_found)
        {
            INFO("Error: Could not map connection (hw:%s) "
                 "to any synapse h/w.\n",
                    con.synapse_hw_name.c_str());
            throw std::runtime_error(
                    "Error: Could not map connection to synapse h/w");
        }
    }

    return mapped_con;
}

void sanafe::SpikingHardware::map_axons()
{
    TRACE1("Creating all connection maps.\n");
    for (Tile &tile : tiles)
    {
        for (Core &core : tile.cores)
        {
            for (MappedNeuron &mapped_neuron : core.neurons)
            {
                sim_create_neuron_axons(mapped_neuron);
            }
        }
    }

    TRACE1("Finished creating connection maps.\n");
    sim_print_axon_summary(*this);
}

sanafe::RunData::RunData(const long int start, const long int steps)
        : timestep_start(start)
        , timesteps_executed(steps)
{
}

sanafe::RunData sanafe::SpikingHardware::sim(
        const long int timesteps, const long int heartbeat)
{
    RunData rd((total_timesteps + 1), timesteps);
    if (total_timesteps <= 0)
    {
        // If no timesteps have been simulated, open the trace files
        //  and simulate.
        if (spike_trace_enabled)
        {
            spike_trace = sim_trace_open_spike_trace(out_dir);
        }
        if (potential_trace_enabled)
        {
            potential_trace = sim_trace_open_potential_trace(out_dir, *this);
        }
        if (perf_trace_enabled)
        {
            perf_trace = sim_trace_open_perf_trace(out_dir);
        }
        if (message_trace_enabled)
        {
            message_trace = sim_trace_open_message_trace(out_dir);
        }
    }

    for (long int timestep = 1; timestep <= timesteps; timestep++)
    {
        if ((timestep % heartbeat) == 0)
        {
            // Print a heart-beat message to show that the simulation is running
            INFO("*** Time-step %ld ***\n", timestep);
        }
        const Timestep ts = step();
        INFO("neurons fired: %zu\n", ts.neurons_fired);
        rd.energy += ts.energy;
        rd.sim_time += ts.sim_time;
        rd.spikes += ts.spike_count;
        rd.packets_sent += ts.packets_sent;
        rd.neurons_fired += ts.neurons_fired;
        rd.wall_time = wall_time;
    }

    return rd;
}

sanafe::Timestep sanafe::SpikingHardware::step()
{
    // Run neuromorphic hardware simulation for one timestep
    //  Measure the CPU time it takes and accumulate the stats
    ++total_timesteps;
    Timestep ts = Timestep(total_timesteps, core_count);
    timespec ts_start;
    timespec ts_end;
    timespec ts_elapsed;

    // Run and measure the wall-clock time taken to run the simulation
    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    sim_timestep(ts, *this);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

    total_energy += ts.energy;
    total_sim_time += ts.sim_time;
    total_spikes += ts.spike_count;
    total_neurons_fired += ts.neurons_fired;
    total_messages_sent += ts.packets_sent;
    if (spike_trace_enabled)
    {
        sim_trace_record_spikes(spike_trace, total_timesteps, *this);
    }
    if (potential_trace_enabled)
    {
        sim_trace_record_potentials(potential_trace, total_timesteps, *this);
    }
    if (perf_trace_enabled)
    {
        sim_trace_perf_log_timestep(perf_trace, ts);
    }
    if (message_trace_enabled)
    {
        for (const auto &q : ts.messages)
        {
            for (const auto &m : q)
            {
                // Ignore dummy messages (without a destination). These account
                //  for processing that doesn't result in a spike being sent
                sim_trace_record_message(message_trace, m);
            }
        }
    }
    constexpr double ns_in_second = 1.0e9;
    wall_time +=
            (double) ts_elapsed.tv_sec + (ts_elapsed.tv_nsec / ns_in_second);
    SIM_TRACE1("Time-step took: %fs.\n",
            static_cast<double>(
                    ts_elapsed.tv_sec + (ts_elapsed.tv_nsec / ns_in_second)));

    return ts;
}

double sanafe::SpikingHardware::get_power() const
{
    double power; // Watts
    if (total_sim_time > 0.0)
    {
        power = total_energy / total_sim_time;
    }
    else
    {
        // Avoid divide by 0
        power = 0.0;
    }

    return power;
}

sanafe::RunData sanafe::SpikingHardware::get_run_summary() const
{
    // Store the summary data in a string to string mapping
    RunData run_data(0, total_timesteps);

    run_data.energy = total_energy;
    run_data.sim_time = total_sim_time;
    run_data.spikes = total_spikes;
    run_data.packets_sent = total_messages_sent;
    run_data.wall_time = wall_time;
    run_data.neurons_fired = total_neurons_fired;

    return run_data;
}

void sanafe::sim_output_run_summary(
        const std::filesystem::path &output_dir, const RunData &run_data)
{
    // Summarize and output the run data using a YAML format to the console
    sim_format_run_summary(std::cout, run_data);

    // Output the same YAML-formatted summary to the given output file
    const std::filesystem::path summary_filename("run_summary.yaml");
    const std::filesystem::path summary_path = output_dir / summary_filename;
    std::ofstream summary_file(summary_path);
    if (summary_file.is_open())
    {
        sim_format_run_summary(summary_file, run_data);
    }
    else
    {
        INFO("Summary file %s couldn't open.\n", summary_path.c_str());
    }
}

void sanafe::sim_format_run_summary(std::ostream &out, const RunData &run_data)
{
    out << "build_git_version: '" << GIT_COMMIT << "'" << std::endl;
    out << "energy: " << std::scientific << run_data.energy << std::endl;
    out << "sim_time: " << std::scientific << run_data.sim_time;
    out << std::endl;
    out << "total_spikes: " << run_data.spikes << std::endl;
    out << "total_messages_sent: " << run_data.packets_sent << std::endl;
    out << "wall_time: " << std::fixed << run_data.wall_time << std::endl;
    out << "total_neurons_fired: " << run_data.neurons_fired << std::endl;
}

std::ofstream sanafe::sim_trace_open_spike_trace(
        const std::filesystem::path &out_dir)
{
    const std::filesystem::path spike_path = out_dir / "spikes.csv";
    std::ofstream spike_file(spike_path);

    if (!spike_file.is_open())
    {
        throw std::runtime_error(
                "Error: Couldn't open trace file for writing.");
    }
    sim_trace_write_spike_header(spike_file);
    return spike_file;
}

std::ofstream sanafe::sim_trace_open_potential_trace(
        const std::filesystem::path &out_dir, const SpikingHardware &hw)
{
    const std::filesystem::path potential_path = out_dir / "potential.csv";
    std::ofstream potential_file(potential_path);

    if (!potential_file.is_open())
    {
        throw std::runtime_error(
                "Error: Couldn't open trace file for writing.");
    }
    sim_trace_write_potential_header(potential_file, hw);
    return potential_file;
}

std::ofstream sanafe::sim_trace_open_perf_trace(
        const std::filesystem::path &out_dir)
{
    const std::filesystem::path perf_path = out_dir / "perf.csv";
    std::ofstream perf_file(perf_path);
    if (!perf_file.is_open())
    {
        throw std::runtime_error(
                "Error: Couldn't open trace file for writing.");
    }
    sim_trace_write_perf_header(perf_file);

    return perf_file;
}

std::ofstream sanafe::sim_trace_open_message_trace(
        const std::filesystem::path &out_dir)
{
    const std::filesystem::path message_path = out_dir / "messages.csv";
    std::ofstream message_file(message_path);
    if (!message_file.is_open())
    {
        throw std::runtime_error(
                "Error: Couldn't open trace file for writing.");
    }
    sim_trace_write_message_header(message_file);
    return message_file;
}

void sanafe::sim_timestep(Timestep &ts, SpikingHardware &hw)
{
    Scheduler scheduler;

    // Start the next time-step, clear all buffers
    assert(hw.core_count > 0);
    ts = Timestep(ts.timestep, hw.core_count);
    sim_reset_measurements(hw);

    pipeline_process_neurons(ts, hw);
    pipeline_process_messages(ts, hw);

    scheduler.noc_width = hw.noc_width;
    scheduler.noc_height = hw.noc_height;
    scheduler.buffer_size = hw.noc_buffer_size;
    scheduler.core_count = hw.core_count;
    scheduler.max_cores_per_tile = hw.max_cores_per_tile;

    ts.sim_time = schedule_messages(ts.messages, scheduler);
    ts.energy = sim_calculate_energy(hw);

    for (auto &tile : hw.tiles)
    {
        ts.total_hops += tile.hops;
        for (auto &c : tile.cores)
        {
            for (const auto &syn : c.synapse)
            {
                ts.spike_count += syn->spikes_processed;
            }
            for (const auto &soma : c.soma)
            {
                ts.neurons_fired += soma->neurons_fired;
            }
            for (const auto &axon_out : c.axon_out_hw)
            {
                ts.packets_sent += axon_out.packets_out;
            }
        }
    }

    SIM_TRACE1("Spikes sent: %ld\n", ts.spike_count);
}

sanafe::Timestep::Timestep(const long int ts, const int core_count)
        : messages(std::vector<std::list<Message>>(core_count))
        , timestep(ts)
{
}

double sanafe::sim_estimate_network_costs(const Tile &src, Tile &dest)
{
    double network_latency;
    long int x_hops;
    long int y_hops;

    network_latency = 0.0;

    // Calculate the energy and time for sending spike packets
    x_hops = abs_diff(src.x, dest.x);
    y_hops = abs_diff(src.y, dest.y);
    // E-W hops

    if (src.x < dest.x)
    {
        dest.east_hops += x_hops;
        network_latency += (double) x_hops * src.latency_east_hop;
    }
    else
    {
        dest.west_hops += x_hops;
        network_latency += (double) x_hops * src.latency_west_hop;
    }

    // N-S hops
    if (src.y < dest.y)
    {
        dest.north_hops += y_hops;
        network_latency += (double) y_hops * src.latency_north_hop;
    }
    else
    {
        dest.south_hops += y_hops;
        network_latency += (double) y_hops * src.latency_south_hop;
    }

    dest.hops += (x_hops + y_hops);
    dest.messages_received++;
    SIM_TRACE1("xhops:%ld yhops%ld total hops:%ld latency:%e\n", x_hops, y_hops,
            x_hops + y_hops, network_latency);
    return network_latency;
}


// TODO: reimplement noise generation from file in C++
/*
double sanafe::sim_generate_noise(Neuron *n)
{
	assert(n != NULL);
	struct SomaUnit &soma_hw = *(n->soma_hw);
	int noise_val = 0;

	if (soma_hw.noise_type == NOISE_FILE_STREAM)
	{
		// With a noise stream, we have a file containing a series of
		//  random values. This is useful if we want to exactly
		//  replicate h/w without knowing how the stream is generated.
		//  We can record the random sequence and replicate it here
		std::string noise_str;
		char *str = &(noise_str[0]);
		// If we get to the end of the stream, by default reset it.
		//  However, it is unlikely the stream will be correct at this
		//  point
		if (feof(soma_hw.noise_stream))
		{
			INFO("Warning: At the end of the noise stream. "
			     "Random values are unlikely to be correct.\n");
			fseek(soma_hw.noise_stream, 0, SEEK_SET);
		}
		char *result = fgets(
			str, MAX_NOISE_FILE_ENTRY, soma_hw.noise_stream);
		if (result != NULL)
		{
			const int ret = sscanf(noise_str, "%d", &noise_val);
			TRACE2("noise val:%d\n", noise_val);

			if (ret < 1)
			{
				INFO("Error: invalid noise stream entry.\n");
			}
		}
	}

	// Get the number of noise bits required TODO: generalize
	int sign_bit = noise_val & 0x100;
	noise_val &= 0x7f; // TODO: hack, fixed for 8 bits
	if (sign_bit)
	{
		// Sign extend
		noise_val |= ~(0x7f);
	}

	return (double) noise_val;
}
*/

double sanafe::sim_calculate_energy(const SpikingHardware &hw)
{
    // Returns the total energy across the design, for this timestep
    double total_energy{0.0};
    double network_energy{0.0};
    double axon_in_energy{0.0};
    double synapse_energy{0.0};
    double soma_energy{0.0};
    double axon_out_energy{0.0};
    double model_simulated_energy{0.0};

    for (const auto &t : hw.tiles)
    {
        double total_hop_energy =
                (static_cast<double>(t.east_hops) * t.energy_east_hop);
        total_hop_energy +=
                (static_cast<double>(t.west_hops) * t.energy_west_hop);
        total_hop_energy +=
                (static_cast<double>(t.south_hops) * t.energy_south_hop);
        total_hop_energy +=
                (static_cast<double>(t.north_hops) * t.energy_north_hop);
        network_energy += total_hop_energy;
        TRACE1("east:%ld west:%ld north:%ld south:%ld\n", t.east_hops,
                t.west_hops, t.north_hops, t.south_hops);

        for (const auto &c : t.cores)
        {
            model_simulated_energy += c.energy;
            for (const auto &axon : c.axon_in_hw)
            {
                axon_in_energy += static_cast<double>(axon.spike_messages_in) *
                        axon.energy_spike_message;
                TRACE1("spikes in: %ld, energy:%e\n", axon.spike_messages_in,
                        axon.energy_spike_message);
            }
            // TODO: fix energy and latency for individual units
            /*
            for (const auto &syn : c.synapse)
            {
                synapse_energy += static_cast<double>(syn.spikes_processed) *
                        syn.energy_spike_op;
                TRACE1("synapse processed: %ld, energy:%e\n",
                        syn.spikes_processed, syn.energy_spike_op);
            }
            for (const auto &soma : c.soma)
            {
                soma_energy += static_cast<double>(soma.neuron_count) *
                        soma.energy_access_neuron;
                soma_energy += static_cast<double>(soma.neuron_updates) *
                        soma.energy_update_neuron;
                soma_energy += static_cast<double>(soma.neurons_fired) *
                        soma.energy_spiking;
                TRACE1("neurons:%ld updates:%ld, spiking:%ld\n",
                        soma.neuron_count, soma.neuron_updates,
                        soma.neurons_fired);
            }
            */
            for (const auto &axon : c.axon_out_hw)
            {
                axon_out_energy += static_cast<double>(axon.packets_out) *
                        axon.energy_access;
                TRACE1("packets: %ld, energy:%e\n", axon.packets_out,
                        axon.energy_access);
            }
        }
    }

    // TODO: clean up energy breakdown. Include model simulated energies
    //  in their respective units, and track energies of cores, tiles etc
    total_energy = model_simulated_energy + axon_in_energy + synapse_energy +
            soma_energy + axon_out_energy + network_energy;

    TRACE1("model_simulated_energy:%e\n", model_simulated_energy);
    TRACE1("axon_in_energy:%e\n", axon_in_energy);
    TRACE1("synapse_energy:%e\n", synapse_energy);
    TRACE1("soma_energy:%e\n", soma_energy);
    TRACE1("axon_out_energy:%e\n", axon_out_energy);
    TRACE1("network_energy:%e\n", network_energy);
    TRACE1("total:%e\n", total_energy);

    return total_energy;
}

void sanafe::sim_create_neuron_axons(MappedNeuron &pre_neuron)
{
    // Setup the connections between neurons and map them to hardware
    assert(pre_neuron.core != nullptr);

    // Figure out the unique set of cores that this neuron broadcasts to
    TRACE1("Counting connections for neuron nid:%s\n", pre_neuron.id.c_str());
    std::set<Core *> cores_out;
    for (MappedConnection &curr_connection : pre_neuron.connections_out)
    {
        TRACE1("Looking at connection id: %d\n", curr_connection.id);
        Core *dest_core = curr_connection.post_neuron->core;
        cores_out.insert(dest_core);
        TRACE1("Connected to dest core: %zu\n", dest_core->id);
    }

    TRACE1("Creating connections for neuron nid:%s to %zu core(s)\n",
            pre_neuron.id.c_str(), cores_out.size());
    for (Core *dest_core : cores_out)
    {
        // Create the axon, and add it to both the destination and
        //  source cores
        sim_allocate_axon(pre_neuron, *dest_core);
    }
    TRACE3("Counted all maps for nid:%s count: %d\n", pre_neuron.id.c_str());

    for (MappedConnection &curr_connection : pre_neuron.connections_out)
    {
        // Add every connection to the axon. Also link to the map in the
        //  post synaptic core / neuron
        Core &post_core = *(curr_connection.post_neuron->core);
        TRACE1("Adding connection:%d\n", curr_connection.id);
        sim_add_connection_to_axon(curr_connection, post_core);
    }
    TRACE1("Finished mapping connections to hardware for nid:%s.%zu.\n",
            pre_neuron.parent_group_id.c_str(), pre_neuron.id);
}

void sanafe::sim_add_connection_to_axon(MappedConnection &con, Core &post_core)
{
    // Add a given connection to the axon in the post-synaptic core
    TRACE3("Adding to connection to axon:%zu\n",
            post_core.axons_out.size() - 1);

    post_core.connections_in.push_back(&con);
    con.synapse_address = post_core.connections_in.size() - 1;

    // Access the most recently created axon in for the post-synaptic core
    AxonInModel &last_added_target_axon = post_core.axons_in.back();
    last_added_target_axon.synapse_addresses.push_back(con.synapse_address);
}

void sanafe::sim_print_axon_summary(SpikingHardware &hw)
{
    int in_count = 0;
    int out_count = 0;

    INFO("** Mapping summary **\n");
    for (Tile &tile : hw.tiles)
    {
        // For debug only, print the axon maps
        for (Core &core : tile.cores)
        {
            bool core_used = false;
            for (size_t k = 0; k < core.neurons.size(); k++)
            {
#ifdef DEBUG
                Neuron *n = c->neurons[k];
                TRACE2("\tnid:%s.%s ", n->group->id.c_str(), n->id.c_str());
                TRACE2("i:%d o:%d\n", n->maps_in_count, n->maps_out_count);
#endif
                core_used = true;
            }

            if (core_used)
            {
                // TODO: update these counts of axons
                //in_count += c.axon_in_hw.axons;
                //out_count += c.axon_out.map_count;
            }
        }
    }
    INFO("Total cores: %zu\n", hw.core_count);
    INFO("Average in map count: %lf\n",
            static_cast<double>(in_count) / hw.core_count);
    INFO("Average out map count: %lf\n",
            static_cast<double>(out_count) / hw.core_count);
}


void sanafe::sim_allocate_axon(MappedNeuron &pre_neuron, Core &post_core)
{
    // Create a new input axon at a receiving (destination) core
    //  Then create the output axon at the sending core. Finally
    //  update the presynaptic neuron and postsynaptic neuron

    Core &pre_core = *(pre_neuron.core);

    TRACE3("Adding connection to core.\n");
    // Allocate the axon and its connections at the post-synaptic core
    post_core.axons_in.emplace_back(AxonInModel());
    const size_t new_axon_in_address = post_core.axons_in.size() - 1;

    // Add the axon at the sending, pre-synaptic core
    TRACE1("Axon in address:%zu for core:%zu.%zu\n", new_axon_in_address,
            post_core.parent_tile_id, post_core.id);
    AxonOutModel out;
    out.dest_axon_id = new_axon_in_address;
    out.dest_core_offset = post_core.offset;
    out.dest_tile_id = post_core.parent_tile_id;
    out.src_neuron_id = pre_neuron.id;
    pre_core.axons_out.push_back(out);
    const size_t new_axon_out_address = pre_core.axons_out.size() - 1;

    // Then add the output axon to the sending pre-synaptic neuron
    pre_neuron.axon_out_addresses.push_back(new_axon_out_address);
    TRACE1("nid:%s.%zu cid:%zu.%zu added one output axon address %zu.\n",
            pre_neuron.parent_group_id.c_str(), pre_neuron.id.c_str(),
            pre_core.parent_tile_id, pre_core.offset, new_axon_out_address);
}

void sanafe::sim_reset_measurements(SpikingHardware &hw)
{
    // Reset any energy, time latency or other measurements of network
    //  hardware
    for (auto &t : hw.tiles)
    {
        // Reset tile
        t.energy = 0.0;

        t.hops = 0;
        t.east_hops = 0;
        t.west_hops = 0;
        t.south_hops = 0;
        t.north_hops = 0;
        t.messages_received = 0;
        for (auto &c : t.cores)
        {
            for (MappedNeuron &neuron : c.neurons)
            {
                neuron.status = sanafe::IDLE;
            }
            // Reset core
            c.energy = 0.0;
            c.next_message_generation_delay = 0.0;

            for (auto axon : c.axon_in_hw)
            {
                axon.spike_messages_in = 0L;
                axon.energy = 0.0;
                axon.time = 0;
            }

            for (auto &dendrite : c.dendrite)
            {
                dendrite->energy = 0.0;
                dendrite->time = 0.0;
            }

            for (auto &syn : c.synapse)
            {
                syn->energy = 0.0;
                syn->time = 0.0;
                syn->spikes_processed = 0;
            }

            for (auto &soma : c.soma)
            {
                soma->energy = 0.0;
                soma->time = 0.0;
                soma->neuron_updates = 0L;
                soma->neurons_fired = 0L;
            }

            for (auto &axon : c.axon_out_hw)
            {
                axon.energy = 0.0;
                axon.time = 0.0;
                axon.packets_out = 0;
            }

            // Reset the message buffer
            c.messages_in = std::vector<Message *>();
        }
    }
}

void sanafe::sim_trace_write_spike_header(std::ofstream &spike_trace_file)
{
    assert(spike_trace_file.is_open());
    spike_trace_file << "neuron,timestep" << std::endl;
}

void sanafe::sim_trace_write_potential_header(
        std::ofstream &potential_trace_file, const SpikingHardware &hw)
{
    // Write csv header for probe outputs - record which neurons have been
    //  probed
    assert(potential_trace_file.is_open());
    potential_trace_file << "timestep,";
    for (const auto &[group_name, group_neurons]: hw.mapped_neuron_groups)
    {
        for (const MappedNeuron *neuron : group_neurons)
        {
            if (neuron->log_potential)
            {
                potential_trace_file << "neuron " << group_name;
                potential_trace_file << "." << neuron->id << ",";
            }
        }
    }
    potential_trace_file << std::endl;
}

void sanafe::sim_trace_write_perf_header(std::ofstream &perf_trace_file)
{
    assert(perf_trace_file.is_open());
    perf_trace_file << "timestep,";
    perf_trace_file << "fired,";
    perf_trace_file << "packets,";
    perf_trace_file << "hops,";
    perf_trace_file << "sim_time,";
    perf_trace_file << "total_energy,";
    perf_trace_file << std::endl;
}

void sanafe::sim_trace_write_message_header(std::ofstream &message_trace_file)
{
    assert(message_trace_file.is_open());
    message_trace_file << "timestep,";
    message_trace_file << "src_neuron,";
    message_trace_file << "src_hw,";
    message_trace_file << "dest_hw,";
    message_trace_file << "hops,";
    message_trace_file << "spikes,";
    message_trace_file << "generation_latency,";
    message_trace_file << "network_latency,";
    message_trace_file << "processing_latency,";
    message_trace_file << "blocking_latency";
    message_trace_file << std::endl;
}

void sanafe::sim_trace_record_spikes(std::ofstream &spike_trace_file,
        const long int timestep, const SpikingHardware &hw)
{
    // A trace of all spikes that are generated
    assert(spike_trace_file.is_open());

    for (const auto &[group_name, group_neurons] : hw.mapped_neuron_groups)
    {
        for (const MappedNeuron *neuron : group_neurons)
        {
            if (neuron->log_spikes && (neuron->status == sanafe::FIRED))
            {
                spike_trace_file << neuron->parent_group_name << ".";
                spike_trace_file << neuron->id << ",";
                spike_trace_file << timestep;
                spike_trace_file << std::endl;
            }
        }
    }
}

void sanafe::sim_trace_record_potentials(std::ofstream &potential_trace_file,
        const long int timestep, const SpikingHardware &hw)
{
    // Each line of this csv file is the potential of all probed neurons for
    //  one time-step
    assert(potential_trace_file.is_open());
    SIM_TRACE1("Recording potential for timestep: %d\n", timestep);
    potential_trace_file << timestep << ",";

    long int potential_probe_count = 0;
    for (const auto &[group_name, group_neurons] : hw.mapped_neuron_groups)
    {
        for (const MappedNeuron *neuron : group_neurons)
        {
            if (neuron->log_potential)
            {
                potential_trace_file << neuron->soma_hw->get_potential(
                        neuron->mapped_address);
                potential_trace_file << ",";
                potential_probe_count++;
            }
        }
    }

    // Each timestep takes up a line in the respective csv file
    if (potential_probe_count > 0)
    {
        potential_trace_file << std::endl;
    }
}

void sanafe::sim_trace_perf_log_timestep(std::ofstream &out, const Timestep &ts)
{
    out << ts.timestep << ",";
    out << ts.neurons_fired << ",";
    out << ts.packets_sent << ",";
    out << ts.total_hops << ",";
    out << std::scientific << ts.sim_time << ",";
    out << std::scientific << ts.energy << ",";
    out << std::endl;
}

void sanafe::sim_trace_record_message(
        std::ofstream &message_trace_file, const Message &m)
{
    assert(message_trace_file.is_open());
    message_trace_file << m.timestep << ",";
    message_trace_file << m.src_neuron_group_id << ".";
    message_trace_file << m.src_neuron_id << ",";
    message_trace_file << m.src_tile_id << "." << m.src_core_offset << ",";

    if (m.placeholder)
    {
        message_trace_file << "x.x,";
    }
    else
    {
        message_trace_file << m.dest_tile_id << ".";
        message_trace_file << m.dest_core_offset << ",";
    }

    message_trace_file << m.hops << ",";
    message_trace_file << m.spikes << ",";
    message_trace_file << m.generation_delay << ",";
    message_trace_file << m.network_delay << ",";
    message_trace_file << m.receive_delay << ",";
    message_trace_file << m.blocked_delay << ",";
    message_trace_file << std::endl;
}

timespec sanafe::calculate_elapsed_time(
        const timespec &ts_start, const timespec &ts_end)
{
    // Calculate elapsed wall-clock time between ts_start and ts_end
    timespec ts_elapsed;

    ts_elapsed.tv_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
    ts_elapsed.tv_sec = ts_end.tv_sec - ts_start.tv_sec;
    if (ts_end.tv_nsec < ts_start.tv_nsec)
    {
        --ts_elapsed.tv_sec;
        constexpr double ns_in_second = 1000000000UL;
        ts_elapsed.tv_nsec += ns_in_second;
    }

    return ts_elapsed;
}

/*
double sim_calculate_time(const Architecture *const arch)
{
	Core *c = n->core;
	SIM_TRACE1("nid:%d sending spike(s).\n", n->id);
	int core_id = n->core->id;


	for (int k = 0; k < n->maps_out_count; k++)
	{
		Axon *dest_axon;
		Message *m;
		int message_index;

		dest_axon = n->maps_out[k];
		message_index = ts->message_queues[core_id].count;

		// Generate a spike message
		m = &(ts->messages[core_id][message_index]);
		arch_init_message(m);
		m->timestep = ts->timestep;
		m->src_neuron = n;
		m->spikes = dest_axon->connection_count;
		m->dest_neuron = dest_axon->connections[0]->post_neuron;
		// Add axon access cost to message latency and energy
		m->generation_latency =
			c->next_message.generation_delay +
			c->axon_out.latency_access;
		sim_message_fifo_push(&(ts->message_queues[core_id]), m);

		c->axon_out.packets_out++;

		// Record a spike message at all the connected cores (axons)
		dest_axon->spikes_received++;
		dest_axon->message = m;

		// Reset the next message in this core
		arch_init_message(&(c->next_message));
	}

	return 0.0;
}
*/

/*
int sim_poisson_input(const double firing_probability)
{
	// Simulate a single external input (as one neuron) for a timestep
	//  Return 1 if the input fires, 0 otherwise
	double rand_uniform;
	int input_fired;

	rand_uniform = (double) rand() / RAND_MAX;
	input_fired = rand_uniform < firing_probability;

	return input_fired;
}

int sim_rate_input(const double firing_rate, double *current)
{
	int input_fired;

	// Note: rate-based input (without randomization) is equivalent to a
	//  neuron with a fixed bias.
	TRACE2("rate input:%lf\n", firing_rate);
	*current += firing_rate;
	if (*current > 255.0)
	{
		*current = 0;
		input_fired = 1;
	}
	else
	{
		input_fired = 0;
	}

	TRACE2("input fired value:%d\n", input_fired);
	return input_fired;
}
*/
