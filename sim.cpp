// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  sim.cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <list>
#include <memory>
#include <optional>
#include <sstream>
#include <vector>

#include "arch.hpp"
#include "models.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"
#include "schedule.hpp"
#include "sim.hpp"

sanafe::Simulation::Simulation(Architecture &a, Network &n,
        const std::string &output_dir = ".", const bool record_spikes = false,
        const bool record_potentials = false, const bool record_perf = false,
        const bool record_messages = false)
        : arch(a)
        , net(n)
        , out_dir(output_dir)
{
    INFO("Initializing simulation.\n");
    total_energy = 0.0; // Joules
    total_sim_time = 0.0; // Seconds
    wall_time = 0.0; // Seconds
    total_timesteps = 0;
    total_spikes = 0;
    total_messages_sent = 0;
    total_neurons_fired = 0;

    // All logging disabled by default
    spike_trace_enabled = record_spikes;
    potential_trace_enabled = record_potentials;
    perf_trace_enabled = record_perf;
    message_trace_enabled = record_messages;

    // Do final mapping stages and sanity check
    arch_create_axons(arch);
    net.check_mapped();
}

sanafe::Simulation::~Simulation()
{
    // Close any open trace files
    spike_trace.close();
    potential_trace.close();
    perf_trace.close();
    message_trace.close();
}

sanafe::RunData::RunData(const long int start, const long int steps)
{
    timestep_start = start;
    timesteps_executed = steps;
    energy = 0.0;
    sim_time = 0.0;
    spikes = 0L;
    packets_sent = 0L;
    neurons_fired = 0L;
    return;
}

sanafe::RunData sanafe::Simulation::run(
        const long int timesteps, const long int heartbeat)
{
    RunData rd((total_timesteps + 1), timesteps);
    if (total_timesteps <= 0)
    {
        // If no timesteps have been simulated, open the trace files
        //  and simulate.
        // TODO: consider moving this back and initializing the
        //  simulation with a reference to an arch and network
        if (spike_trace_enabled)
        {
            spike_trace = sim_trace_open_spike_trace(out_dir);
        }
        if (potential_trace_enabled)
        {
            potential_trace = sim_trace_open_potential_trace(out_dir, net);
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
            // Print heart-beat every hundred timesteps
            INFO("*** Time-step %ld ***\n", timestep);
        }
        const Timestep ts = step();
        rd.energy += ts.energy;
        rd.sim_time += ts.sim_time;
        rd.spikes += ts.spike_count;
        rd.packets_sent += ts.packets_sent;
        rd.neurons_fired += ts.neurons_fired;
        rd.wall_time = wall_time;
    }

    return rd;
}

sanafe::Timestep sanafe::Simulation::step()
{
    // Run neuromorphic hardware simulation for one timestep
    //  Measure the CPU time it takes and accumulate the stats
    total_timesteps++;
    struct Timestep ts = Timestep(total_timesteps, arch.core_count);
    struct timespec ts_start, ts_end, ts_elapsed;

    // Run and measure the wall-clock time taken to run the simulation
    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    sim_timestep(ts, arch, net);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

    total_energy += ts.energy;
    total_sim_time += ts.sim_time;
    total_spikes += ts.spike_count;
    total_neurons_fired += ts.neurons_fired;
    total_messages_sent += ts.packets_sent;
    if (spike_trace_enabled)
    {
        sim_trace_record_spikes(spike_trace, total_timesteps, net);
    }
    if (potential_trace_enabled)
    {
        sim_trace_record_potentials(potential_trace, total_timesteps, net);
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
    wall_time += (double) ts_elapsed.tv_sec + (ts_elapsed.tv_nsec / 1.0e9);
    SIM_TRACE1("Time-step took: %fs.\n",
            static_cast<double>(
                    ts_elapsed.tv_sec + (ts_elapsed.tv_nsec / 1.0e9)));

    return ts;
}

double sanafe::Simulation::get_power()
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

sanafe::RunData sanafe::Simulation::get_run_summary()
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
        const std::filesystem::path &out_dir, const RunData &run_data)
{
    // Summarize and output the run data using a YAML format to the console
    sim_format_run_summary(std::cout, run_data);

    // Output the same YAML-formatted summary to the given output file
    const std::filesystem::path summary_filename("run_summary.yaml");
    const std::filesystem::path summary_path = out_dir / summary_filename;
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

    return;
}

std::ofstream sanafe::sim_trace_open_spike_trace(
        const std::filesystem::path &out_dir)
{
    // TODO: warning, this is specific to Linux
    // To be more portable, consider using the filesystem library in C++17
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
        const std::filesystem::path &out_dir, const Network &net)
{
    // TODO: warning, this is specific to Linux
    // To be more portable, consider using the filesystem library in C++17
    const std::filesystem::path potential_path = out_dir / "potential.csv";
    std::ofstream potential_file(potential_path);

    if (!potential_file.is_open())
    {
        throw std::runtime_error(
                "Error: Couldn't open trace file for writing.");
    }
    sim_trace_write_potential_header(potential_file, net);
    return potential_file;
}

std::ofstream sanafe::sim_trace_open_perf_trace(
        const std::filesystem::path &out_dir)
{
    // TODO: warning, this is specific to Linux
    // To be more portable, consider using the filesystem library in C++17
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

void sanafe::sim_timestep(Timestep &ts, Architecture &arch, Network &net)
{
    Scheduler scheduler;

    // Start the next time-step, clear all buffers
    ts = Timestep(ts.timestep, arch.core_count);
    sim_reset_measurements(net, arch);

    pipeline_process_neurons(ts, arch);
    pipeline_process_messages(ts, arch);

    scheduler.noc_width = arch.noc_width;
    scheduler.noc_height = arch.noc_height;
    scheduler.buffer_size = arch.noc_buffer_size;
    scheduler.core_count = arch.core_count;
    scheduler.max_cores_per_tile = arch.max_cores_per_tile;

    ts.sim_time = schedule_messages(ts.messages, scheduler);
    // Performance statistics for this time step
    ts.energy = sim_calculate_energy(arch);

    for (auto &tile : arch.tiles)
    {
        int tile_spike_count = 0;
        ts.total_hops += tile.hops;
        for (auto &c : tile.cores)
        {
            for (std::vector<SynapseUnit>::size_type k = 0;
                    k < c.synapse.size(); k++)
            {
                ts.spike_count += c.synapse[k].spikes_processed;
                tile_spike_count += c.synapse[k].spikes_processed;
            }
            for (std::vector<SomaUnit>::size_type k = 0; k < c.soma.size(); k++)
            {
                ts.neurons_fired += c.soma[k].neurons_fired;
            }
            for (std::vector<AxonOutUnit>::size_type k = 0;
                    k < c.axon_out_hw.size(); k++)
            {
                ts.packets_sent += c.axon_out_hw[k].packets_out;
            }
        }
    }

    SIM_TRACE1("Spikes sent: %ld\n", ts.spike_count);
    return;
}

sanafe::Timestep::Timestep(const long int ts, const int core_count)
{
    timestep = ts;
    spike_count = 0L;
    messages = std::vector<std::list<Message>>(core_count);
    neurons_fired = 0L;
    total_hops = 0L;
    energy = 0.0;
    sim_time = 0.0;
    packets_sent = 0L;
}

double sanafe::sim_estimate_network_costs(Tile &src, Tile &dest)
{
    double network_latency;
    long int x_hops, y_hops;

    network_latency = 0.0;

    // Calculate the energy and time for sending spike packets
    x_hops = abs(src.x - dest.x);
    y_hops = abs(src.y - dest.y);
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

double sanafe::sim_calculate_energy(const Architecture &arch)
{
    // Returns the total energy across the design, for this timestep
    double network_energy, synapse_energy, soma_energy, axon_out_energy;
    double axon_in_energy, total_energy;

    network_energy = 0.0;
    axon_in_energy = 0.0;
    synapse_energy = 0.0;
    soma_energy = 0.0;
    axon_out_energy = 0.0;

    for (auto &t : arch.tiles)
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

        for (auto &c : t.cores)
        {
            for (std::vector<AxonInUnit>::size_type k = 0;
                    k < c.axon_in_hw.size(); k++)
            {
                axon_in_energy += c.axon_in_hw[k].spike_messages_in *
                        c.axon_in_hw[k].energy_spike_message;
                TRACE1("spikes in: %ld, energy:%e\n",
                        c.axon_in_hw[k].spike_messages_in,
                        c.axon_in_hw[k].energy_spike_message);
            }
            for (std::vector<SynapseUnit>::size_type k = 0;
                    k < c.synapse.size(); k++)
            {
                synapse_energy += c.synapse[k].spikes_processed *
                        c.synapse[k].energy_spike_op;
                TRACE1("synapse processed: %ld, energy:%e\n",
                        c.synapse[k].spikes_processed,
                        c.synapse[k].energy_spike_op);
            }
            for (std::vector<SomaUnit>::size_type k = 0; k < c.soma.size(); k++)
            {
                soma_energy += static_cast<double>(c.soma[k].neuron_count) *
                        c.soma[k].energy_access_neuron;
                soma_energy += static_cast<double>(c.soma[k].neuron_updates) *
                        c.soma[k].energy_update_neuron;
                soma_energy += static_cast<double>(c.soma[k].neurons_fired) *
                        c.soma[k].energy_spiking;
                TRACE1("neurons:%ld updates:%ld, spiking:%ld\n",
                        c.soma[k].neuron_count, c.soma[k].neuron_updates,
                        c.soma[k].neurons_fired);
            }
            for (std::vector<AxonOutUnit>::size_type k = 0;
                    k < c.axon_out_hw.size(); k++)
            {
                axon_out_energy +=
                        static_cast<double>(c.axon_out_hw[k].packets_out) *
                        c.axon_out_hw[k].energy_access;
                TRACE1("packets: %ld, energy:%e\n",
                        c.axon_out_hw[k].packets_out,
                        c.axon_out_hw[k].energy_access);
            }
        }
    }

    total_energy = axon_in_energy + synapse_energy + soma_energy +
            axon_out_energy + network_energy;

    TRACE1("axon_in_energy:%e\n", axon_in_energy);
    TRACE1("synapse_energy:%e\n", synapse_energy);
    TRACE1("soma_energy:%e\n", soma_energy);
    TRACE1("axon_out_energy:%e\n", axon_out_energy);
    TRACE1("network_energy:%e\n", network_energy);
    TRACE1("total:%e\n", total_energy);

    return total_energy;
}

void sanafe::sim_reset_measurements(Network &net, Architecture &arch)
{
    for (auto &group : net.groups)
    {
        for (auto &n : group->neurons)
        {
            n.status = sanafe::IDLE;
        }
    }

    // Reset any energy, time latency or other measurements of network
    //  hardware
    for (auto &t : arch.tiles)
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
            // Reset core
            c.energy = 0.0;
            c.next_message_generation_delay = 0.0;

            for (std::vector<AxonInUnit>::size_type k = 0;
                    k < c.axon_in_hw.size(); k++)
            {
                c.axon_in_hw[k].spike_messages_in = 0L;
                c.axon_in_hw[k].energy = 0.0;
                c.axon_in_hw[k].time = 0;
            }

            for (std::vector<DendriteUnit>::size_type k = 0;
                    k < c.dendrite.size(); k++)
            {
                c.dendrite[k].energy = 0.0;
                c.dendrite[k].time = 0.0;
            }

            for (std::vector<SynapseUnit>::size_type k = 0;
                    k < c.synapse.size(); k++)
            {
                c.synapse[k].energy = 0.0;
                c.synapse[k].time = 0.0;
                c.synapse[k].spikes_processed = 0;
            }

            for (std::vector<SomaUnit>::size_type k = 0; k < c.soma.size(); k++)
            {
                c.soma[k].energy = 0.0;
                c.soma[k].time = 0.0;
                c.soma[k].neuron_updates = 0L;
                c.soma[k].neurons_fired = 0L;
            }

            for (std::vector<AxonOutUnit>::size_type k = 0;
                    k < c.axon_out_hw.size(); k++)
            {
                c.axon_out_hw[k].energy = 0.0;
                c.axon_out_hw[k].time = 0.0;
                c.axon_out_hw[k].packets_out = 0;
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

    return;
}

void sanafe::sim_trace_write_potential_header(
        std::ofstream &potential_trace_file, const Network &net)
{
    // Write csv header for probe outputs - record which neurons have been
    //  probed
    assert(potential_trace_file.is_open());
    potential_trace_file << "timestep,";
    for (auto &group : net.groups)
    {
        for (auto &n : group->neurons)
        {
            if (n.log_potential)
            {
                potential_trace_file << "neuron " << group->id;
                potential_trace_file << "." << n.id << ",";
            }
        }
    }
    potential_trace_file << std::endl;

    return;
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

    return;
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

    return;
}

void sanafe::sim_trace_record_spikes(
        std::ofstream &out, const long int timestep, const Network &net)
{
    // A trace of all spikes that are generated
    assert(out.is_open());

    for (auto &group : net.groups)
    {
        for (auto &n : group->neurons)
        {
            if (n.log_spikes && (n.status == sanafe::FIRED))
            {
                out << n.parent_group_id << "." << n.id << ",";
                out << timestep;
                out << std::endl;
            }
        }
    }

    return;
}

// TODO: should potential be a required value in the soma?
void sanafe::sim_trace_record_potentials(
        std::ofstream &out, const int timestep, const Network &net)
{
    // Each line of this csv file is the potential of all probed neurons for
    //  one time-step
    assert(out.is_open());
    SIM_TRACE1("Recording potential for timestep: %d\n", timestep);
    out << timestep << ",";

    long int potential_probe_count = 0;
    for (auto &group : net.groups)
    {
        for (auto &n : group->neurons)
        {
            if (n.log_potential)
            {
                out << n.soma_model->get_potential() << ",";
                potential_probe_count++;
            }
        }
    }

    // Each timestep takes up a line in the respective csv file
    if (potential_probe_count > 0)
    {
        out << std::endl;
    }

    return;
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

void sanafe::sim_trace_record_message(std::ofstream &out, const Message &m)
{
    assert(out.is_open());
    out << m.timestep << ",";
    out << m.src_neuron_group_id << "." << m.src_neuron_id << ",";
    out << m.src_tile_id << "." << m.src_core_offset << ",";

    if (m.placeholder)
    {
        out << "x.x,";
    }
    else
    {
        out << m.dest_tile_id << "." << m.dest_core_offset << ",";
    }

    out << m.hops << ",";
    out << m.spikes << ",";
    out << m.generation_delay << ",";
    out << m.network_delay << ",";
    out << m.receive_delay << ",";
    out << m.blocked_delay << ",";

    out << std::endl;

    return;
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
        ts_elapsed.tv_sec--;
        ts_elapsed.tv_nsec += 1000000000UL;
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