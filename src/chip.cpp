// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  chip.cpp
#include <algorithm>
#include <cassert>
#include <chrono>
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

#include <booksim_lib.hpp> // For cycle-accurate NoC simulation
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "arch.hpp"
#include "chip.hpp"
#include "core.hpp"
#include "message.hpp"
#include "models.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "plugins.hpp"
#include "print.hpp"
#include "utils.hpp" // For abs_diff
#include "schedule.hpp"
#include "tile.hpp"

// Track the number of SpikingChips current instantiated, only as a safety
//  mechanism for cycle-accurate runs (Booksim2): We throw an error
//  at launch for cycle-accurate runs when multiple SANA-FE simulations are
//  active, as Booksim2 stores a bunch of global state and weird things will
//  happen if we access the library across multiple threads or simulations.
//  If I get time to better encapsulate the Booksim2 library, this check
//  will no longer be needed...
std::atomic<int> sanafe::SpikingChip::chip_count = 0;

sanafe::SpikingChip::SpikingChip(const Architecture &arch,
        const std::filesystem::path &output_dir, const bool record_spikes,
        const bool record_potentials, const bool record_perf,
        const bool record_messages)
        : core_count(arch.core_count)
        , noc_width(arch.noc_width_in_tiles)
        , noc_height(arch.noc_height_in_tiles)
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
            for (const PipelineUnitConfiguration &unit_config :
                    core_config.pipeline_hw)
            {
                hardware_core.create_pipeline_unit(unit_config);
            }
            for (const AxonOutConfiguration &axon_config : core_config.axon_out)
            {
                hardware_core.create_axon_out(axon_config);
            }
        }
    }

    const std::vector<std::string> booksim_config_vec(
            std::begin(booksim_config_str), std::end(booksim_config_str));
    BookSimConfig new_config =
            booksim_load_config(std::move(booksim_config_vec));
    // Use a unique_ptr for the config so that we don't need to include Booksim
    //  library in the header (meaning that plugins using chip.hpp don't need to
    //  also include this library)
    booksim_config = std::make_unique<BookSimConfig>(new_config);
    chip_count++;
}

sanafe::SpikingChip::~SpikingChip()
{
    // Close any open trace files
    spike_trace.close();
    potential_trace.close();
    perf_trace.close();
    message_trace.close();

    chip_count--;
}

void sanafe::SpikingChip::load(const SpikingNetwork &net)
{
    map_neurons(net);
    map_connections(net);
}

void sanafe::SpikingChip::map_neurons(const SpikingNetwork &net)
{
    std::vector<const Neuron *> neurons_in_mapped_order;

    // Figure out which order to map the neurons in to hardware
    for (const auto &[name, group] : net.groups)
    {
        mapped_neuron_groups[name] =
                std::vector<MappedNeuron *>(group.neurons.size(), nullptr);
        for (const Neuron &neuron : group.neurons)
        {
            neurons_in_mapped_order.push_back(&neuron);
        }
    }
    INFO("Total neurons to map: %zu\n", neurons_in_mapped_order.size());

    std::sort(neurons_in_mapped_order.begin(), neurons_in_mapped_order.end(),
            [](const Neuron *const a, const Neuron *const b) {
                return a->mapping_order < b->mapping_order;
            });

    auto list_of_cores = cores();
    // Map all neurons in order
    for (const Neuron *neuron : neurons_in_mapped_order)
    {
        if (!neuron->core_address.has_value())
        {
            std::string error = "Neuron: " + neuron->parent_group_name + "." +
                    std::to_string(neuron->offset) + " not mapped.";
            INFO("%s", error.c_str());
            throw std::runtime_error(error);
        }
        TRACE1(CHIP, "Mapping neuron %s.%zu to core:%zu\n",
                neuron->parent_group_name.c_str(), neuron->offset,
                neuron->core_address.value().id);
        Core &mapped_core = list_of_cores[neuron->core_address.value().id];
        mapped_core.map_neuron(*neuron, total_neurons_mapped);
        ++total_neurons_mapped;
    }

    // Now that we mapped all neurons, index them using their group and neuron
    //  IDs so that we can easily connect neurons to each other later.
    //  Note we don't do this in the loop above because the vector will be
    //   changing sizes dynamically as we map.
    for (Core &core : list_of_cores)
    {
        for (MappedNeuron &mapped_neuron : core.neurons)
        {
            mapped_neuron_groups[mapped_neuron.parent_group_name]
                                [mapped_neuron.offset] = &mapped_neuron;
        }
    }
}

void sanafe::SpikingChip::map_connections(const SpikingNetwork &net)
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

    map_axons();

    // Set each connection attribute
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

                for (auto &[key, value] : con.synapse_attributes)
                {
                    if (value.forward_to_synapse)
                    {
                        mapped_con.synapse_hw->check_attribute(key);
                        mapped_con.synapse_hw->set_attribute_edge(
                                mapped_con.synapse_address, key, value);
                    }
                    if (value.forward_to_dendrite)
                    {
                        MappedNeuron &n = *(mapped_con.post_neuron);
                        n.dendrite_hw->check_attribute(key);
                        n.dendrite_hw->set_attribute_edge(
                                mapped_con.synapse_address, key, value);
                    }
                }
            }
        }
    }

    return;
}

sanafe::MappedConnection &sanafe::SpikingChip::map_connection(
        const Connection &con)
{
    auto list_of_cores = cores();

    auto &pre_group = mapped_neuron_groups.at(con.pre_neuron.group_name);
    MappedNeuron &pre_neuron =
            *(pre_group[con.pre_neuron.neuron_offset.value()]);

    auto &post_group = mapped_neuron_groups.at(con.post_neuron.group_name);
    MappedNeuron &post_neuron =
            *(post_group[con.post_neuron.neuron_offset.value()]);

    pre_neuron.connections_out.emplace_back();
    MappedConnection &mapped_con = pre_neuron.connections_out.back();
    mapped_con.pre_neuron = &pre_neuron;
    mapped_con.post_neuron = &post_neuron;

    // Map to synapse hardware unit
    Core &post_core = *(post_neuron.core);
    mapped_con.synapse_hw = post_core.pipeline_hw[0].get();

    bool choose_first_by_default = (con.synapse_hw_name.length() == 0);
    bool synapse_found = false;
    for (auto &hw : post_core.pipeline_hw)
    {
        if (hw->implements_synapse &&
                (choose_first_by_default || (con.synapse_hw_name == hw->name)))
        {
            mapped_con.synapse_hw = hw.get();
            synapse_found = true;
            break;
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

    mapped_con.build_message_processing_pipeline();
    return mapped_con;
}

void sanafe::SpikingChip::map_axons()
{
    TRACE1(CHIP, "Creating all connection maps.\n");
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

    TRACE1(CHIP, "Finished creating connection maps.\n");
    sim_print_axon_summary();
}

sanafe::RunData::RunData(const long int start, const long int steps)
        : timestep_start(start)
        , timesteps_executed(steps)
{
}

void sanafe::SpikingChip::retire_scheduled_messages(
        sanafe::RunData &rd, sanafe::Scheduler &scheduler)
{
    while (!scheduler.timesteps_to_write.empty())
    {
        Timestep ts;
        scheduler.timesteps_to_write.pop(ts);
        TRACE1(CHIP, "retiring ts:%ld\n", ts.timestep);
        if (perf_trace_enabled)
        {
            sim_trace_record_perf(perf_trace, ts);
        }
        if (message_trace_enabled)
        {
            for (auto &q : *(ts.messages))
            {
                for (Message &m : q)
                {
                    sim_trace_record_message(message_trace, m);
                }
            }
        }
        rd.total_energy += ts.total_energy;
        rd.synapse_energy += ts.synapse_energy;
        rd.dendrite_energy += ts.dendrite_energy;
        rd.soma_energy += ts.soma_energy;
        rd.network_energy += ts.network_energy;
        rd.sim_time += ts.sim_time;
        rd.spikes += ts.spike_count;
        rd.packets_sent += ts.packets_sent;
        rd.neurons_fired += ts.neurons_fired;
        rd.wall_time = wall_time;

        total_sim_time += ts.sim_time;
    }

    return;
}

sanafe::RunData sanafe::SpikingChip::sim(const long int timesteps,
        const long int heartbeat, const TimingModel timing_model,
        const int scheduler_threads)
{
    RunData rd((total_timesteps + 1), timesteps);
    Scheduler scheduler;

    scheduler.noc_width = noc_width;
    scheduler.noc_height = noc_height;
    scheduler.buffer_size = noc_buffer_size;
    scheduler.core_count = core_count;
    scheduler.max_cores_per_tile = max_cores_per_tile;
    scheduler.timing_model = timing_model;
    // TOOD: for now, hard code the number of threads. Then calculate the
    //  hard-coded thread count. The finally allow threads to dynamically change
    //  based on performance
    INFO("Creating %d scheduler threads\n", scheduler_threads);
    for (int tid = 0; tid < scheduler_threads; tid++)
    {
        TRACE1(CHIP, "Created scheduler thread:%d\n", tid);
        scheduler.scheduler_threads.emplace_back(
                &schedule_messages_thread, std::ref(scheduler), tid);
    }

    if (total_timesteps <= 0)
    {
        // If no timesteps have been simulated, open the trace files
        if (spike_trace_enabled)
        {
            spike_trace = sim_trace_open_spike_trace(out_dir);
        }
        if (potential_trace_enabled)
        {
            potential_trace = sim_trace_open_potential_trace(out_dir);
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
        Timestep ts = step(scheduler);
        // Retire and trace messages as we go along. Not strictly required, but
        //  avoids the message write queue growing too large
        retire_scheduled_messages(rd, scheduler);
    }

    schedule_stop_all_threads(scheduler, message_trace, rd);
    retire_scheduled_messages(rd, scheduler);

    return rd;
}

sanafe::Timestep sanafe::SpikingChip::step(Scheduler &scheduler)
{
    // Run neuromorphic hardware simulation for one timestep
    //  Measure the CPU time it takes and accumulate the stats
    std::chrono::high_resolution_clock timer;

    ++total_timesteps;
    Timestep ts = Timestep(total_timesteps, core_count);
    double ts_elapsed;

    // Run and measure the wall-clock time taken to run the simulation
    auto ts_start = timer.now();
    sim_hw_timestep(ts, scheduler);
    auto ts_end = timer.now();
    ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

    // Update global chip performance counters
    total_energy += ts.total_energy;
    synapse_energy += ts.synapse_energy;
    dendrite_energy += ts.dendrite_energy;
    soma_energy += ts.soma_energy;
    network_energy += ts.network_energy;

    total_spikes += ts.spike_count;
    total_neurons_fired += ts.neurons_fired;
    // The total_messages_sent is incremented during the simulation since it's
    //  used to calculate the mid, so nothing needs to be done here
    if (spike_trace_enabled)
    {
        sim_trace_record_spikes(spike_trace, total_timesteps);
    }
    if (potential_trace_enabled)
    {
        sim_trace_record_potentials(potential_trace, total_timesteps);
    }
    wall_time += ts_elapsed;

    return ts;
}

void sanafe::SpikingChip::reset()
{
    for (Tile tile : tiles)
    {
        for (Core core : tile.cores)
        {
            core.timestep_buffer.clear();
            for (auto hw : core.pipeline_hw)
            {
                hw->reset();
            }
        }
    }

    for (auto &[group_name, neurons] : mapped_neuron_groups)
    {
        for (MappedNeuron *neuron : neurons)
        {
            neuron->status = INVALID_NEURON_STATE;
        }
    }
}

double sanafe::SpikingChip::get_power() const
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

// Pipeline modeling
void sanafe::SpikingChip::process_neurons(Timestep &ts)
{
    auto core_list = cores();

    // Older versions of OpenMP don't support range-based for loops yet...
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    // codechecker_suppress [modernize-loop-convert]
    for (size_t idx = 0; idx < core_list.size(); idx++)
    {
        for (MappedNeuron &n : core_list[idx].get().neurons)
        {
            process_neuron(ts, n);
        }

        Core &core = core_list[idx];
        // Account for any remaining neuron processing
        bool placeholder_event = (core.next_message_generation_delay != 0.0);
        if (placeholder_event)
        {
            const MappedNeuron &last_neuron = core.neurons.back();
            Message placeholder =
                    Message(placeholder_mid, *this, last_neuron, ts.timestep);
            placeholder.generation_delay = core.next_message_generation_delay;
            // Create a dummy placeholder message
            auto &message_queue = *(ts.messages);
            message_queue[core.id].push_back(placeholder);
        }
    }
}

void sanafe::SpikingChip::process_messages(Timestep &ts)
{
    // Assign outgoing spike messages to their respective destination
    //  cores, and calculate network costs
    auto &message_queue = *(ts.messages);
    for (auto &q : message_queue)
    {
        for (auto &m : q)
        {
            if (!m.placeholder)
            {
                receive_message(m);
            }
        }
    }

    // Now process all messages at receiving cores
    auto core_list = cores();
    // Older versions of OpenMP don't support range-based for loops yet...
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    // codechecker_suppress [modernize-loop-convert]
    for (size_t idx = 0; idx < core_list.size(); idx++)
    {
        //INFO("omp thread:%d\n", omp_get_thread_num());
        Core &core = core_list[idx];
        TRACE1(CHIP, "Processing %zu message(s) for cid:%zu\n",
                core.messages_in.size(), core.id);
        for (auto &m : core.messages_in)
        {
            m.receive_delay += process_message(ts, core, m);
        }
    }
}

void sanafe::SpikingChip::receive_message(Message &m)
{
    assert(static_cast<size_t>(m.src_tile_id) < tiles.size());
    assert(static_cast<size_t>(m.dest_tile_id) < tiles.size());

    const Tile &src_tile = tiles[m.src_tile_id];
    Tile &dest_tile = tiles[m.dest_tile_id];

    m.network_delay = sim_estimate_network_costs(src_tile, dest_tile);
    m.hops = abs_diff(src_tile.x, dest_tile.x) +
            abs_diff(src_tile.y, dest_tile.y);

    Core &core = dest_tile.cores[m.dest_core_offset];
    core.messages_in.push_back(m);
}

void sanafe::SpikingChip::process_neuron(Timestep &ts, MappedNeuron &n)
{
    Core &c = *(n.core);
    bool simulate_buffer = (c.pipeline_config.buffer_position ==
                                   BUFFER_BEFORE_DENDRITE_UNIT) ||
            (c.pipeline_config.buffer_position == BUFFER_BEFORE_SOMA_UNIT);

    PipelineResult input{};
    if (simulate_buffer)
    {
        // Read then clear the pipeline buffer entry
        input = c.timestep_buffer[n.mapped_address];
        c.timestep_buffer[n.mapped_address] = PipelineResult{};
    }
    PipelineResult pipeline_output = execute_pipeline(
            n.neuron_processing_pipeline, ts, n, std::nullopt, input);
    n.core->next_message_generation_delay += pipeline_output.latency.value();
    pipeline_process_axon_out(ts, n);

    return;
}

double sanafe::SpikingChip::process_message(
        Timestep &ts, Core &core, Message &m)
{
    double message_processing_latency = pipeline_process_axon_in(core, m);

    assert(static_cast<size_t>(m.dest_axon_id) < core.axons_in.size());
    const AxonInModel &axon_in = core.axons_in[m.dest_axon_id];
    PipelineResult empty_input{}; // Empty/default struct used as dummy input

    for (const int synapse_address : axon_in.synapse_addresses)
    {
        MappedConnection &con = *(core.connections_in[synapse_address]);

        // In certain pipeline configurations, every synaptic lookup requires
        //  updates to the dendrite and/or soma units as well. Keep propagating
        //  outputs/inputs until we hit the time-step buffer, where outputs
        //  are stored as inputs ready for the next time-step
        MappedNeuron &n = *(con.post_neuron);
        PipelineResult pipeline_output = execute_pipeline(
                con.message_processing_pipeline, ts, n, &con, empty_input);
        core.timestep_buffer[n.mapped_address] = pipeline_output;
        m.receive_delay += pipeline_output.latency.value();
    }

    return message_processing_latency;
}

sanafe::PipelineResult sanafe::SpikingChip::execute_pipeline(
        const std::vector<PipelineUnit *> &pipeline, Timestep &ts,
        MappedNeuron &n, std::optional<MappedConnection *> con,
        const PipelineResult &input)
{
    double total_energy{0.0};
    double total_latency{0.0};

    PipelineResult output{input};
    for (auto unit : pipeline)
    {
        output = unit->process(ts, n, con, output);
        total_energy += output.energy.value();
        total_latency += output.latency.value();
        n.status = output.status;
    }

    output.energy = total_energy;
    output.latency = total_latency;
    return output;
}

sanafe::PipelineResult sanafe::PipelineUnit::process(Timestep &ts,
        MappedNeuron &n, std::optional<MappedConnection *> con,
        const PipelineResult &input)
{
    TRACE2(CHIP, "Updating nid:%zu (ts:%ld)\n", n.id, ts.timestep);
    set_time(ts.timestep);

    // Process inputs
    PipelineResult output = (this->*process_input_fn)(ts, n, con, input);

    // Post-processing on outputs
    (this->*process_output_fn)(n, con, output);

#ifndef NDEBUG
    check_outputs(n, output);
#endif
    energy += output.energy.value();

    return output;
}

double sanafe::SpikingChip::pipeline_process_axon_in(
        Core &core, const Message &m)
{
    assert(m.dest_axon_hw >= 0);
    assert(static_cast<size_t>(m.dest_axon_hw) < core.axon_in_hw.size());
    AxonInUnit &axon_unit = core.axon_in_hw[m.dest_axon_hw];
    axon_unit.spike_messages_in++;

    return axon_unit.latency_spike_message;
}

sanafe::PipelineResult sanafe::SpikingChip::pipeline_process_axon_out(
        Timestep &ts, MappedNeuron &n)
{
    PipelineResult axon_result{};
    axon_result.latency = 0.0;
    axon_result.energy = 0.0;

    if (!n.axon_out_input_spike)
    {
        return axon_result;
    }

    TRACE1(CHIP, "nid:%s.%zu sending spike message to %zu axons out\n",
            n.parent_group_name.c_str(), n.id, n.axon_out_addresses.size());
    for (const int axon_address : n.axon_out_addresses)
    {
        Message m(total_messages_sent, *this, n, ts.timestep, axon_address);
        // Add axon access cost to message latency and energy
        AxonOutUnit &axon_out_hw = *(n.axon_out_hw);
        axon_out_hw.energy += axon_out_hw.energy_access;

        m.generation_delay = n.core->next_message_generation_delay +
                axon_out_hw.latency_access;
        n.core->next_message_generation_delay = 0.0;

        auto &message_queue = *(ts.messages);
        message_queue[n.core->id].push_back(m);
        ++axon_out_hw.packets_out;
        ++total_messages_sent;

        axon_result.energy.value() += axon_out_hw.energy_access;
        axon_result.latency.value() += axon_out_hw.latency_access;
    }
    n.axon_out_input_spike = false;

    return axon_result;
}

std::vector<std::reference_wrapper<sanafe::Core>> sanafe::SpikingChip::cores()
{
    std::vector<std::reference_wrapper<Core>> all_cores_in_hw;

    for (Tile &tile : tiles)
    {
        std::copy(tile.cores.begin(), tile.cores.end(),
                std::back_inserter(all_cores_in_hw));
    }

    return all_cores_in_hw;
}

sanafe::AxonInUnit::AxonInUnit(const AxonInConfiguration &config)
        : name(config.name)
        , energy_spike_message(config.metrics.energy_message_in)
        , latency_spike_message(config.metrics.latency_message_in)
{
}

sanafe::AxonOutUnit::AxonOutUnit(const AxonOutConfiguration &config)
        : name(std::move(config.name))
        , energy_access(config.metrics.energy_message_out)
        , latency_access(config.metrics.latency_message_out)
{
}

sanafe::Core::Core(const CoreConfiguration &config)
        : pipeline_config(config.pipeline)
        , name(std::move(config.name))
        , id(config.address.id)
        , offset(config.address.offset_within_tile)
        , parent_tile_id(config.address.parent_tile_id)
        , log_energy(config.pipeline.log_energy)
        , log_latency(config.pipeline.log_latency)

{
    timestep_buffer.resize(pipeline_config.max_neurons_supported);
}

void sanafe::SpikingChip::sim_output_run_summary(
        const std::filesystem::path &output_dir, const RunData &run_data) const
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

#ifndef GIT_COMMIT
#define GIT_COMMIT "git-hash-unknown"
#endif

void sanafe::SpikingChip::sim_format_run_summary(
        std::ostream &out, const RunData &run_data) const
{
    out << "build_git_version: '" << GIT_COMMIT << "'\n";
    out << "timesteps_executed: " << run_data.timesteps_executed << "\n";
    out << "total_spikes: " << run_data.spikes << "\n";
    out << "total_messages_sent: " << run_data.packets_sent << "\n";
    out << "total_neurons_fired: " << run_data.neurons_fired << "\n";
    out << "sim_time: " << std::scientific << run_data.sim_time << "\n";
    // Give a more detailed energy breakdown, as this is often useful
    out << "energy:\n";
    out << "  synapse:" << std::scientific << run_data.synapse_energy << "\n";
    out << "  dendrite:" << std::scientific << run_data.dendrite_energy << "\n";
    out << "  soma:" << std::scientific << run_data.soma_energy << "\n";
    out << "  network: " << std::scientific << run_data.network_energy << "\n";
    out << "  total: " << std::scientific << run_data.total_energy << "\n";
    // Give a more detailed simulator walltime breakdown too
    out << "wall_time:\n";
    out << "  neuron_processing: " << std::fixed << neuron_processing_wall
        << "\n";
    out << "  message_processing: " << std::fixed << message_processing_wall
        << "\n";
    out << "  scheduler: " << std::fixed << scheduler_wall << "\n";
    out << "  energy: " << std::fixed << energy_stats_wall << "\n";
    double misc_wall = run_data.wall_time -
            (neuron_processing_wall + message_processing_wall + scheduler_wall +
                    energy_stats_wall);
    out << "  misc: " << std::fixed << misc_wall << "\n";
    out << "  total: " << std::fixed << run_data.wall_time << "\n";
}

std::ofstream sanafe::SpikingChip::sim_trace_open_spike_trace(
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

std::ofstream sanafe::SpikingChip::sim_trace_open_potential_trace(
        const std::filesystem::path &out_dir)
{
    const std::filesystem::path potential_path = out_dir / "potentials.csv";
    std::ofstream potential_file(potential_path);

    if (!potential_file.is_open())
    {
        throw std::runtime_error(
                "Error: Couldn't open trace file for writing.");
    }
    sim_trace_write_potential_header(potential_file);
    return potential_file;
}

std::ofstream sanafe::SpikingChip::sim_trace_open_perf_trace(
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

std::ofstream sanafe::SpikingChip::sim_trace_open_message_trace(
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

void sanafe::SpikingChip::forced_updates(const Timestep &ts)
{
    // You can optionally force a neuron to update its associated h/w
    //  every time-step, regardless of whether it received inputs or not.
    // Note that energy is accounted for, but latency is not considered here.
    auto core_list = cores();
#ifdef ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t idx = 0; idx < core_list.size(); idx++)
    {
        Core &core = core_list[idx];
        for (MappedNeuron &n : core.neurons)
        {
            if (n.force_synapse_update)
            {
                for (MappedConnection &con : n.connections_out)
                {
                    con.synapse_hw->set_time(ts.timestep);
                    PipelineResult result =
                            con.synapse_hw->update(con.synapse_address);
                    if (result.energy.has_value())
                    {
                        con.synapse_hw->energy += result.energy.value();
                    }
                    // Latency is not considered; as it isn't within
                    //  either neuron processing or message processing
                }
            }
            if (n.force_dendrite_update)
            {
                n.dendrite_hw->set_time(ts.timestep);
                sanafe::PipelineResult result = n.dendrite_hw->update(
                        n.mapped_address, std::nullopt, std::nullopt);
                if (result.energy.has_value())
                {
                    n.dendrite_hw->energy += result.energy.value();
                }
                // Latency is not considered; as it isn't within
                //  either neuron processing or message processing
            }
            // Note that soma updates will always be handed in the neuron
            //  processing loop and so don't need to be supported here
        }
    }

    return;
}

void sanafe::SpikingChip::sim_hw_timestep(Timestep &ts, Scheduler &scheduler)
{
    std::chrono::high_resolution_clock timer;

    // Start the next time-step, clear all buffers
    assert(core_count > 0);
    ts = Timestep(ts.timestep, core_count);
    sim_reset_measurements();

    auto neuron_processing_start_tm = timer.now();
    process_neurons(ts);
    auto neuron_processing_end_tm = timer.now();
    auto message_processing_start_tm = neuron_processing_end_tm;
    process_messages(ts);
    forced_updates(ts);
    auto message_processing_end_tm = timer.now();
    auto energy_calculation_start_tm = message_processing_end_tm;
    sim_calculate_energy(ts);
    for (auto &tile : tiles)
    {
        ts.total_hops += tile.hops;
        for (auto &c : tile.cores)
        {
            for (const auto &hw : c.pipeline_hw)
            {
                ts.spike_count += hw->spikes_processed;
                ts.neurons_fired += hw->neurons_fired;
            }
            for (const auto &axon_out : c.axon_out_hw)
            {
                ts.packets_sent += axon_out.packets_out;
            }
        }
    }
    TRACE1(CHIP, "Spikes sent: %ld\n", ts.spike_count);
    auto energy_calculation_end_tm = timer.now();
    auto scheduler_start_tm = energy_calculation_end_tm;

    if (scheduler.timing_model == TIMING_MODEL_SIMPLE)
    {
        TRACE1(CHIP, "Running simple timing model\n");
        schedule_messages_simple(ts, scheduler);
    }
    else if (scheduler.timing_model == TIMING_MODEL_DETAILED)
    {
        TRACE1(CHIP, "Running detailed timing model\n");
        schedule_messages_detailed(ts, scheduler);
    }
    else if (scheduler.timing_model == TIMING_MODEL_CYCLE_ACCURATE)
    {
        TRACE1(CHIP, "Running cycle-accurate timing model\n");
        if ((scheduler.scheduler_threads.size() > 1) || (chip_count > 1))
        {
            INFO("Error: Cannot run multiple simultaneous cycle-accurate "
                 "simulations. The Booksim2 library does not support "
                 "concurrent runs as it has a lot of global state. For now "
                 "it's simplest to just not allow for concurrent runs. If you "
                 "need to simulate in parallel, launch separate SANA-FE "
                 "processes.");
            throw std::runtime_error(
                    "Error: Cannot run multiple simultaneous cycle-accurate "
                    "simulations.");
        }
        schedule_messages_cycle_accurate(ts, *booksim_config, scheduler);
    }
    else
    {
        INFO("Error: Timing model:%d not recognized\n", scheduler.timing_model);
        throw std::invalid_argument("Timing model not recognized");
    }
    auto scheduler_end_tm = timer.now();

    // Calculate various simulator timings
    constexpr double ns_to_second = 1.0e-9;
    neuron_processing_wall +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                    neuron_processing_end_tm - neuron_processing_start_tm)
                    .count() *
            ns_to_second;
    message_processing_wall +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                    message_processing_end_tm - message_processing_start_tm)
                    .count() *
            ns_to_second;
    scheduler_wall += std::chrono::duration_cast<std::chrono::nanoseconds>(
                              scheduler_end_tm - scheduler_start_tm)
                              .count() *
            ns_to_second;
    energy_stats_wall +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                    energy_calculation_end_tm - energy_calculation_start_tm)
                    .count() *
            ns_to_second;
    TRACE1(CHIP, "neuron:%e message:%e scheduler:%e energy:%e\n",
            neuron_processing_wall, message_processing_wall, scheduler_wall,
            energy_stats_wall);
}

sanafe::Timestep::Timestep(const long int ts, const int core_count)
        : messages(
                  std::make_shared<std::vector<std::list<Message>>>(core_count))
        , timestep(ts)
{
}

double sanafe::SpikingChip::sim_estimate_network_costs(
        const Tile &src, Tile &dest)
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
    TRACE1(CHIP, "xhops:%ld yhops%ld total hops:%ld latency:%e\n", x_hops,
            y_hops, x_hops + y_hops, network_latency);
    return network_latency;
}

void sanafe::SpikingChip::sim_calculate_energy(Timestep &ts)
{
    // Returns the total energy across the design, for this timestep
    ts.network_energy = 0.0;
    ts.synapse_energy = 0.0;
    ts.dendrite_energy = 0.0;
    ts.soma_energy = 0.0;
    ts.total_energy = 0.0;

    double axon_in_energy{0.0};
    double axon_out_energy{0.0};
    double total_pipeline_energy{0.0};

    for (auto &t : tiles)
    {
        double total_hop_energy =
                (static_cast<double>(t.east_hops) * t.energy_east_hop);
        total_hop_energy +=
                (static_cast<double>(t.west_hops) * t.energy_west_hop);
        total_hop_energy +=
                (static_cast<double>(t.south_hops) * t.energy_south_hop);
        total_hop_energy +=
                (static_cast<double>(t.north_hops) * t.energy_north_hop);
        t.energy = total_hop_energy;
        ts.network_energy += total_hop_energy;
        TRACE1(CHIP, "east:%ld west:%ld north:%ld south:%ld\n", t.east_hops,
                t.west_hops, t.north_hops, t.south_hops);

        for (auto &c : t.cores)
        {
            for (const auto &axon : c.axon_in_hw)
            {
                axon_in_energy = static_cast<double>(axon.spike_messages_in) *
                        axon.energy_spike_message;
                TRACE1(CHIP, "spikes in: %ld, energy:%e\n",
                        axon.spike_messages_in, axon.energy_spike_message);
            }

            double pipeline_energy{0.0};
            for (auto &pipeline_unit : c.pipeline_hw)
            {
                // Separately track the total pipeline energy, as the same
                //  energy values may be added to multiple categories, i.e., if
                //  the pipeline h/w unit implements multiple functionality
                pipeline_energy += pipeline_unit->energy;

                if (pipeline_unit->implements_synapse)
                {
                    ts.synapse_energy += pipeline_unit->energy;
                }
                if (pipeline_unit->implements_dendrite)
                {
                    ts.dendrite_energy += pipeline_unit->energy;
                }
                if (pipeline_unit->implements_soma)
                {
                    ts.soma_energy += pipeline_unit->energy;
                }
            }

            for (const auto &axon : c.axon_out_hw)
            {
                axon_out_energy = axon.energy;
                TRACE1(CHIP, "packets: %ld, energy per packet:%e\n",
                        axon.packets_out, axon.energy_access);
            }

            c.energy = pipeline_energy + axon_in_energy + axon_out_energy;
            t.energy += c.energy;
            ts.network_energy += axon_in_energy;
            ts.network_energy += axon_out_energy;
            total_pipeline_energy += pipeline_energy;
        }
    }

    // TODO: should I keep these network units separate from the NoC costs?
    ts.total_energy = ts.network_energy + total_pipeline_energy;

    TRACE1(CHIP, "total_energy:%e\n", ts.total_energy);
    TRACE1(CHIP, "\tpipeline_energy:%e\n", total_pipeline_energy);
    TRACE1(CHIP, "\tnetwork_energy:%e\n", ts.network_energy);
    TRACE1(CHIP, "\taxon_in_energy:%e\n", axon_in_energy);
    TRACE1(CHIP, "\taxon_out_energy:%e\n", axon_out_energy);

    return;
}

void sanafe::SpikingChip::sim_create_neuron_axons(MappedNeuron &pre_neuron)
{
    // Setup the connections between neurons and map them to hardware
    assert(pre_neuron.core != nullptr);

    // Figure out the unique set of cores that this neuron broadcasts to
    TRACE1(CHIP, "Counting connections for neuron nid:%zu\n", pre_neuron.id);
    std::set<Core *> cores_out;
    for (MappedConnection &curr_connection : pre_neuron.connections_out)
    {
        Core *dest_core = curr_connection.post_neuron->core;
        cores_out.insert(dest_core);
        TRACE1(CHIP, "Connected to dest core: %zu\n", dest_core->id);
    }

    TRACE1(CHIP, "Creating connections for neuron nid:%zu to %zu core(s)\n",
            pre_neuron.id, cores_out.size());
    for (Core *dest_core : cores_out)
    {
        // Create the axon, and add it to both the destination and
        //  source cores
        sim_allocate_axon(pre_neuron, *dest_core);
    }
    TRACE3(CHIP, "Allocated all axons for nid:%zu count: %d\n", pre_neuron.id,
            pre_neuron.maps_out_count);

    for (MappedConnection &curr_connection : pre_neuron.connections_out)
    {
        // Add every connection to the axon. Also link to the map in the
        //  post synaptic core / neuron
        Core &post_core = *(curr_connection.post_neuron->core);
        //TRACE1(CHIP, "Adding connection:%d\n", curr_connection.id);
        sim_add_connection_to_axon(curr_connection, post_core);
    }
    TRACE1(CHIP, "Finished mapping connections to hardware for nid:%s.%zu.\n",
            pre_neuron.parent_group_name.c_str(), pre_neuron.id);
}

void sanafe::SpikingChip::sim_add_connection_to_axon(
        MappedConnection &con, Core &post_core)
{
    // Add a given connection to the axon in the post-synaptic core
    TRACE3(CHIP, "Adding to connection to axon:%zu\n",
            post_core.axons_out.size() - 1);

    // TODO: this is difficult because we don't want to duplicate all the stored
    //  information about synapses twice. However, the connection info may be
    //  useful to the synaptic model when creating the synapses. We can't really
    //  store pointers to the connections because this array will grow as we map
    //  connections and invalidate references.. we also have the challenge of
    //  referencing synapses with a single address value.. what does this mean
    //  with multiple synapse hw elements. unless maybe we leave the lines below
    //  here and go back to pointers again.
    //
    // TODO: these 3 lines should be moved to where we map the connection to h/w
    //  I think
    post_core.connections_in.push_back(&con);
    con.synapse_address = post_core.connections_in.size() - 1;
    con.synapse_hw->map_connection(con);

    // Access the most recently created axon in for the post-synaptic core
    AxonInModel &last_added_target_axon = post_core.axons_in.back();
    last_added_target_axon.synapse_addresses.push_back(con.synapse_address);
}

void sanafe::SpikingChip::sim_print_axon_summary()
{
    int in_count = 0;
    int out_count = 0;

    INFO("** Mapping summary **\n");
    for (Tile &tile : tiles)
    {
        // For debug only, print the axon maps
        for (Core &core : tile.cores)
        {
            bool core_used = false;
            for (size_t k = 0; k < core.neurons.size(); k++)
            {
                if (DEBUG_LEVEL_CHIP >= 2)
                {
                    MappedNeuron &n = core.neurons[k];
                    TRACE2(CHIP, "\tnid:%s.%zu ", n.parent_group_name.c_str(),
                            n.id);
                    TRACE2(CHIP, "i:%d o:%d\n", n.maps_in_count,
                            n.maps_out_count);
                }
                core_used = true;
            }

            if (core_used)
            {
                in_count += core.axons_in.size();
                out_count += core.axons_out.size();
            }
        }
    }
    INFO("Total cores: %zu\n", core_count);
    INFO("Average in map count: %lf\n",
            static_cast<double>(in_count) / core_count);
    INFO("Average out map count: %lf\n",
            static_cast<double>(out_count) / core_count);
}

void sanafe::SpikingChip::sim_allocate_axon(
        MappedNeuron &pre_neuron, Core &post_core)
{
    // Create a new input axon at a receiving (destination) core
    //  Then create the output axon at the sending core. Finally
    //  update the presynaptic neuron and postsynaptic neuron

    Core &pre_core = *(pre_neuron.core);

    TRACE3(CHIP, "Adding connection to core.\n");
    // Allocate the axon and its connections at the post-synaptic core
    post_core.axons_in.emplace_back(AxonInModel());
    const size_t new_axon_in_address = post_core.axons_in.size() - 1;

    // Add the axon at the sending, pre-synaptic core
    TRACE1(CHIP, "Axon in address:%zu for core:%zu.%zu\n", new_axon_in_address,
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
    TRACE1(CHIP, "nid:%s.%zu cid:%zu.%zu added one output axon address %zu.\n",
            pre_neuron.parent_group_name.c_str(), pre_neuron.id,
            pre_core.parent_tile_id, pre_core.offset, new_axon_out_address);
}

void sanafe::SpikingChip::sim_reset_measurements()
{
    // Reset any energy, time latency or other measurements of network
    //  hardware
    for (auto &t : tiles)
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

            for (auto axon : c.axon_in_hw)
            {
                axon.spike_messages_in = 0L;
                axon.energy = 0.0;
                axon.time = 0;
            }

            for (auto &hw : c.pipeline_hw)
            {
                hw->energy = 0.0;
                hw->time = 0.0;
                hw->spikes_processed = 0;
                hw->neuron_updates = 0L;
                hw->neurons_fired = 0L;
            }

            for (auto &axon : c.axon_out_hw)
            {
                axon.energy = 0.0;
                axon.time = 0.0;
                axon.packets_out = 0;
            }

            // Reset the message buffer
            c.messages_in = std::vector<Message>();
        }
    }

    return;
}

void sanafe::SpikingChip::sim_trace_write_spike_header(
        std::ofstream &spike_trace_file)
{
    assert(spike_trace_file.is_open());
    spike_trace_file << "neuron,timestep" << std::endl;
}

void sanafe::SpikingChip::sim_trace_write_potential_header(
        std::ofstream &potential_trace_file)
{
    // Write csv header for probe outputs - record which neurons have been
    //  probed
    assert(potential_trace_file.is_open());
    potential_trace_file << "timestep,";
    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
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
    potential_trace_file << "\n";
    potential_trace_file.flush();
}

std::map<std::string, double>
sanafe::SpikingChip::sim_trace_get_optional_traces()
{
    std::map<std::string, double> optional_perf_traces{};
    for (const Tile &t : tiles)
    {
        if (t.log_energy)
        {
            optional_perf_traces[t.name + ".energy"] = t.energy;
        }
        if (t.log_latency)
        {
            optional_perf_traces[t.name + ".latency"] = 0.0; // TODO
        }
        for (const Core &c : t.cores)
        {
            if (c.log_energy)
            {
                optional_perf_traces[t.name + "." + c.name + ".energy"] =
                        c.energy;
            }
            if (c.log_latency)
            {
                // TODO
                optional_perf_traces[t.name + "." + c.name + ".latency"] = 0.0;
            }

            for (const auto &hw : c.pipeline_hw)
            {
                if (hw->log_energy)
                {
                    optional_perf_traces[t.name + "." + c.name + "." +
                            hw->name + ".energy"] = hw->energy;
                }
                if (hw->log_latency)
                {
                    optional_perf_traces[t.name + "." + c.name + "." +
                            hw->name + ".latency"] = 0.0; // TODO
                }
            }
        }
    }

    return optional_perf_traces;
}

void sanafe::SpikingChip::sim_trace_write_perf_header(
        std::ofstream &perf_trace_file)
{
    assert(perf_trace_file.is_open());
    perf_trace_file << "timestep,";

    // Mandatory performance metrics
    perf_trace_file << "fired,";
    perf_trace_file << "packets,";
    perf_trace_file << "hops,";
    perf_trace_file << "spikes,";
    perf_trace_file << "sim_time,";
    perf_trace_file << "total_energy";

    // Optional performance metrics
    std::map<std::string, double> optional_perf_traces =
            sim_trace_get_optional_traces();
    for (auto &[name, trace] : optional_perf_traces)
    {
        perf_trace_file << "," << name;
    }
    perf_trace_file << "\n";
    perf_trace_file.flush();
}

void sanafe::SpikingChip::sim_trace_write_message_header(
        std::ofstream &message_trace_file)
{
    assert(message_trace_file.is_open());
    message_trace_file << "timestep,";
    message_trace_file << "mid,";
    message_trace_file << "src_neuron,";
    message_trace_file << "src_hw,";
    message_trace_file << "dest_hw,";
    message_trace_file << "hops,";
    message_trace_file << "spikes,";
    message_trace_file << "generation_latency,";
    message_trace_file << "network_latency,";
    message_trace_file << "processing_latency,";
    message_trace_file << "blocking_latency\n";
    message_trace_file.flush();
}

void sanafe::SpikingChip::sim_trace_record_spikes(
        std::ofstream &spike_trace_file, const long int timestep)
{
    // A trace of all spikes that are generated
    assert(spike_trace_file.is_open());

    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
    {
        for (const MappedNeuron *neuron : group_neurons)
        {
            if (neuron->log_spikes && (neuron->status == sanafe::FIRED))
            {
                spike_trace_file << neuron->parent_group_name << ".";
                spike_trace_file << neuron->id << ",";
                spike_trace_file << timestep;
                spike_trace_file << "\n";
                spike_trace_file.flush();
            }
        }
    }
}

void sanafe::SpikingChip::sim_trace_record_potentials(
        std::ofstream &potential_trace_file, const long int timestep)
{
    // Each line of this csv file is the potential of all probed neurons for
    //  one time-step
    assert(potential_trace_file.is_open());
    TRACE1(CHIP, "Recording potential for timestep: %ld\n", timestep);
    potential_trace_file << timestep << ",";

    long int potential_probe_count = 0;
    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
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
        potential_trace_file << "\n";
        potential_trace_file.flush();
    }
}

void sanafe::SpikingChip::sim_trace_record_perf(
        std::ofstream &perf_trace_file, const Timestep &ts)
{
    perf_trace_file << ts.timestep << ",";
    // Start with the mandatory counters and traces
    perf_trace_file << ts.neurons_fired << ",";
    perf_trace_file << ts.packets_sent << ",";
    perf_trace_file << ts.total_hops << ",";
    perf_trace_file << ts.spike_count << ",";
    perf_trace_file << std::scientific << ts.sim_time << ",";
    perf_trace_file << std::scientific << ts.total_energy;
    // Finish with the optional traces, which can be enabled or disabled per
    //  tile, core, or pipeline h/w unit
    const auto optional_perf_traces = sim_trace_get_optional_traces();
    for (auto &[name, trace] : optional_perf_traces)
    {
        perf_trace_file << "," << std::scientific << trace;
    }
    perf_trace_file << "\n";
    perf_trace_file.flush();
}

void sanafe::sim_trace_record_message(
        std::ofstream &message_trace_file, const Message &m)
{
    assert(message_trace_file.is_open());
    message_trace_file << m.timestep << ",";
    message_trace_file << m.mid << ",";
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
    message_trace_file << m.blocked_delay << "\n";
    message_trace_file.flush();
}


sanafe::BufferPosition sanafe::pipeline_parse_buffer_pos_str(
        const std::string &buffer_pos_str, const bool buffer_inside_unit)
{
    BufferPosition buffer_pos;
    // H/w either has SANA-FE insert and manage a time-step buffer on the h/w
    //  unit's inputs (before the unit), or can assume that the unit manages its
    //  own state that is internally buffered over consecutive time-steps e.g.,
    //  using a double buffer

    if (buffer_pos_str == "dendrite")
    {
        if (buffer_inside_unit)
        {
            buffer_pos = BUFFER_INSIDE_DENDRITE_UNIT;
        }
        else
        {
            buffer_pos = BUFFER_BEFORE_DENDRITE_UNIT;
        }
    }
    else if (buffer_pos_str == "soma")
    {
        if (buffer_inside_unit)
        {
            buffer_pos = BUFFER_INSIDE_SOMA_UNIT;
        }
        else
        {
            buffer_pos = BUFFER_BEFORE_SOMA_UNIT;
        }
    }
    else if (buffer_pos_str == "axon_out")
    {
        buffer_pos = BUFFER_BEFORE_AXON_OUT_UNIT;
    }
    else
    {
        INFO("Error: Buffer position %s not supported", buffer_pos_str.c_str());
        throw std::invalid_argument("Error: Buffer position not supported");
    }

    return buffer_pos;
}
