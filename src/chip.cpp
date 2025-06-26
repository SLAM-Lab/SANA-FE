// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  chip.cpp

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <booksim_lib.hpp> // For cycle-accurate NoC simulation
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "arch.hpp"
#include "chip.hpp"

#include "core.hpp"
#include "mapped.hpp"
#include "message.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"
#include "schedule.hpp"
#include "tile.hpp"
#include "timestep.hpp"
#include "utils.hpp" // For abs_diff

sanafe::RunData::RunData(const long int start)
        : timestep_start(start)
{
}

// Track the number of SpikingChips current instantiated, only as a safety
//  mechanism for cycle-accurate runs (Booksim2): We throw an error
//  at launch for cycle-accurate runs when multiple SANA-FE simulations are
//  active, as Booksim2 stores a bunch of global state and weird things will
//  happen if we access the library across multiple threads or simulations.
//  If I get time to better encapsulate the Booksim2 library, this check
//  will no longer be needed...
std::atomic<int> sanafe::SpikingChip::chip_count = 0;

sanafe::SpikingChip::SpikingChip(const Architecture &arch)
        : ts_sync_delay_table(arch.ts_sync_delay_table)
        , core_count(arch.core_count)
        , max_cores_per_tile(arch.max_cores_per_tile)
        , noc_width_in_tiles(arch.noc_width_in_tiles)
        , noc_height_in_tiles(arch.noc_height_in_tiles)
        , noc_buffer_size(arch.noc_buffer_size)
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

    std::vector<std::string> booksim_config_vec(
            std::begin(booksim_config_str), std::end(booksim_config_str));
    const BookSimConfig new_config =
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
    if (spike_trace.is_open())
    {
        spike_trace.close();
    }
    if (potential_trace.is_open())
    {
        potential_trace.close();
    }
    if (perf_trace.is_open())
    {
        perf_trace.close();
    }
    if (message_trace.is_open())
    {
        message_trace.close();
    }

    chip_count--;
}

void sanafe::SpikingChip::load(const SpikingNetwork &net)
{
    map_neurons(net);
    map_connections(net);
}

void sanafe::SpikingChip::map_neurons(const SpikingNetwork &net)
{
    std::vector<std::reference_wrapper<const Neuron>> neurons_in_mapped_order;

    // Figure out which order to map the neurons in to hardware
    for (const auto &[name, group] : net.groups)
    {
        mapped_neuron_groups[name].reserve(group.neurons.size());
        for (const Neuron &neuron : group.neurons)
        {
            neurons_in_mapped_order.emplace_back(neuron);
        }
    }
    INFO("Total neurons to map: %zu\n", neurons_in_mapped_order.size());
    std::sort(neurons_in_mapped_order.begin(), neurons_in_mapped_order.end(),
            [](const Neuron &a, const Neuron &b) {
                return a.mapping_order < b.mapping_order;
            });

    auto list_of_cores = cores();
    // Map all neurons in a given order, which may be non-obvious (e.g., if the
    //  user gives a specific nontrivial mapping ordering)
    for (const Neuron &neuron : neurons_in_mapped_order)
    {
        if (!neuron.core_address.has_value())
        {
            const std::string error = "Neuron: " + neuron.parent_group_name +
                    "." + std::to_string(neuron.offset) + " not mapped.";
            INFO("%s", error.c_str());
            throw std::runtime_error(error);
        }
        TRACE1(CHIP, "Mapping neuron %s.%zu to core:%zu\n",
                neuron.parent_group_name.c_str(), neuron.offset,
                neuron.core_address.value().id);
        Core &mapped_core = list_of_cores[neuron.core_address.value().id];
        mapped_core.map_neuron(neuron, total_neurons_mapped);
        ++total_neurons_mapped;
    }

    // Link mapped neurons to their neuron addresses (group.offset), making it
    //  possible to connect mapped neurons easily later
    track_mapped_neurons();
    // Record how many tiles and cores have been mapped in total
    track_mapped_tiles_and_cores();
}

void sanafe::SpikingChip::track_mapped_neurons()
{
    // Create a structure tracking each neuron address against the corresponding
    //  mapped neuron, object
    auto list_of_cores = cores();

    //  Push neurons first in their mapped order, causing this vector to
    //  initially be out of order i.e., with vector idx != correct offset.
    //  As a second pass we sort the vectors by offset order. Note we can't
    //  populate this vector earlier when we map neurons, because the core's
    //  neuron vector is dynamically populated. As the vector grows, previous
    //  references would be invalidated. Also, as we use reference_wrappers
    //  instead of pointers, we can't simply resize() the vector at the
    //  beginning and fill them as we go
    for (Core &core : list_of_cores)
    {
        for (MappedNeuron &mapped_neuron : core.neurons)
        {
            mapped_neuron_groups[mapped_neuron.parent_group_name].emplace_back(
                    std::ref(mapped_neuron));
        }
    }

    // Now for every neuron group, sort by offset so that we should get the
    //  vector in offset order
    for (auto &[group_name, neuron_refs] : mapped_neuron_groups)
    {
        std::sort(neuron_refs.begin(), neuron_refs.end(),
                [](const auto &a, const auto &b) {
                    return a.get().offset < b.get().offset;
                });

        // Validate that offset == index
        for (size_t i = 0; i < neuron_refs.size(); ++i)
        {
            if (neuron_refs[i].get().offset != i)
            {
                throw std::logic_error("Offset incorrect in group '" +
                        group_name + "' at index " + std::to_string(i) + "(" +
                        std::to_string(neuron_refs[i].get().offset) + ")");
            }
        }
    }
}

void sanafe::SpikingChip::track_mapped_tiles_and_cores() noexcept
{
    // Reset in case this gets called multiple times
    mapped_tiles = 0UL;
    mapped_cores = 0UL;

    for (const Tile &tile : tiles)
    {
        bool tile_used = false;
        for (const Core &core : tile.cores)
        {
            if (!core.neurons.empty())
            {
                tile_used = true;
                ++mapped_cores;
            }
        }

        if (tile_used)
        {
            ++mapped_tiles;
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
    forward_connection_attributes(net);
}

void sanafe::SpikingChip::forward_connection_attributes(
        const SpikingNetwork &net)
{
    // Set each connection attribute
    for (const auto &[name, group] : net.groups)
    {
        for (size_t nid = 0; nid < group.neurons.size(); ++nid)
        {
            const Neuron &pre_neuron = group.neurons[nid];
            MappedNeuron &mapped_neuron = mapped_neuron_groups[name][nid];
            for (size_t idx = 0; idx < pre_neuron.edges_out.size(); ++idx)
            {
                const Connection &con = pre_neuron.edges_out[idx];
                MappedConnection &mapped_con =
                        mapped_neuron.connections_out[idx];

                for (const auto &[key, value] : con.synapse_attributes)
                {
                    if (value.forward_to_synapse)
                    {
                        mapped_con.synapse_hw->check_attribute(key);
                        mapped_con.synapse_hw->set_attribute_edge(
                                mapped_con.synapse_address, key, value);
                    }
                    if (value.forward_to_dendrite)
                    {
                        MappedNeuron &n = mapped_con.post_neuron_ref;
                        n.dendrite_hw->check_attribute(key);
                        n.dendrite_hw->set_attribute_edge(
                                mapped_con.synapse_address, key, value);
                    }
                }
            }
        }
    }
}

sanafe::MappedConnection &sanafe::SpikingChip::map_connection(
        const Connection &con)
{
    if (!con.pre_neuron.neuron_offset.has_value())
    {
        throw std::invalid_argument("Pre neuron doesn't specify group offset");
    }
    if (!con.post_neuron.neuron_offset.has_value())
    {
        throw std::invalid_argument("Post neuron doesn't specify group offset");
    }

    auto &pre_group = mapped_neuron_groups.at(con.pre_neuron.group_name);
    MappedNeuron &pre_neuron = pre_group[con.pre_neuron.neuron_offset.value()];

    auto &post_group = mapped_neuron_groups.at(con.post_neuron.group_name);
    MappedNeuron &post_neuron =
            post_group[con.post_neuron.neuron_offset.value()];

    pre_neuron.connections_out.emplace_back(pre_neuron, post_neuron);
    MappedConnection &mapped_con = pre_neuron.connections_out.back();

    // Map to synapse hardware unit
    Core &post_core = *(post_neuron.core);
    mapped_con.synapse_hw = post_core.pipeline_hw[0].get();

    const bool choose_first_by_default = (con.synapse_hw_name.empty());
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

void sanafe::SpikingChip::flush_timestep_data(
        sanafe::RunData &rd, sanafe::Scheduler &scheduler)
{
    while (!scheduler.timesteps_to_write.empty())
    {
        TimestepHandle timestep_handle;
        const bool got_ts = scheduler.timesteps_to_write.pop(timestep_handle);
        Timestep &ts = timestep_handle.get();

        // Attempt to get the timestep. Note that if someone else is writing to
        //  the queue, we will have to wait until they are finished
        if (got_ts)
        {
            TRACE1(CHIP, "Retiring ts:%ld\n", ts.timestep);
            if (perf_trace.is_open())
            {
                sim_trace_record_perf(perf_trace, ts);
            }

            if (message_trace.is_open())
            {
                // Not crucial, but its nice to print messages in ID order.
                //  Copy all the messages into a single vector and sort
                std::vector<std::reference_wrapper<const Message>> all_messages;
                for (auto &q : ts.messages)
                {
                    for (const Message &m : q)
                    {
                        all_messages.emplace_back(m);
                    }
                }
                // Sort messages in message ID order (with placeholders last)
                std::sort(all_messages.begin(), all_messages.end(),
                        CompareMessagesByID{});
                // Save the messages in sorted order
                for (const Message &m : all_messages)
                {
                    sim_trace_record_message(message_trace, m);
                }
            }
            update_run_data(rd, ts);
            retire_timestep(ts);
        }
        // else someone else was writing, try again
    }
}

void sanafe::SpikingChip::update_run_data(
        sanafe::RunData &rd, const sanafe::Timestep &ts)
{
    rd.total_energy += ts.total_energy;
    rd.synapse_energy += ts.synapse_energy;
    rd.dendrite_energy += ts.dendrite_energy;
    rd.soma_energy += ts.soma_energy;
    rd.network_energy += ts.network_energy;
    rd.sim_time += ts.sim_time;
    rd.spikes += ts.spike_count;
    rd.packets_sent += ts.packets_sent;
    rd.neurons_updated += ts.neurons_updated;
    rd.neurons_fired += ts.neurons_fired;
}

sanafe::RunData sanafe::SpikingChip::sim(const long int timesteps,
        const TimingModel timing_model, const int scheduler_thread_count,
        const bool record_spikes,
        const bool record_potentials,
        const bool record_perf,
        const bool record_messages,
        std::string output_dir)
{
    RunData rd(total_timesteps + 1);
    rd.timesteps_executed += timesteps;

    Scheduler scheduler;
    scheduler.noc_width_in_tiles = noc_width_in_tiles;
    scheduler.noc_height_in_tiles = noc_height_in_tiles;
    scheduler.buffer_size = noc_buffer_size;
    scheduler.core_count = core_count;
    scheduler.max_cores_per_tile = max_cores_per_tile;
    scheduler.timing_model = timing_model;
    schedule_create_threads(scheduler, scheduler_thread_count);

    if (total_timesteps <= 0)
    {
        // If no timesteps have been simulated, open the trace files
        if (record_spikes)
        {
            spike_trace = sim_trace_open_spike_trace(output_dir);
        }
        if (record_potentials)
        {
            potential_trace = sim_trace_open_potential_trace(output_dir);
        }
        if (record_perf)
        {
            perf_trace = sim_trace_open_perf_trace(output_dir);
        }
        if (record_messages)
        {
            message_trace = sim_trace_open_message_trace(output_dir);
        }
    }

    for (long int timestep = 1; timestep <= timesteps; timestep++)
    {
        if ((timestep % heartbeat) == 0)
        {
            // Print a heart-beat message to show that the simulation is running
            INFO("*** Time-step %ld ***\n", timestep);
        }
        step(scheduler);
        flush_timestep_data(rd, scheduler);
    }

    schedule_stop_all_threads(scheduler);
    flush_timestep_data(rd, scheduler);

    return rd;
}

void sanafe::SpikingChip::sim_update_total_energy_and_counts(const Timestep &ts)
{
    // Update global chip performance counters
    total_energy += ts.total_energy;
    synapse_energy += ts.synapse_energy;
    dendrite_energy += ts.dendrite_energy;
    soma_energy += ts.soma_energy;
    network_energy += ts.network_energy;

    total_spikes += ts.spike_count;
    total_neurons_updated += ts.neurons_updated;
    total_neurons_fired += ts.neurons_fired;
}

void sanafe::SpikingChip::step(Scheduler &scheduler)
{
    // Run neuromorphic hardware simulation for one timestep
    //  Measure the CPU time it takes and accumulate the stats
    ++total_timesteps;
    //Timestep ts = Timestep(total_timesteps);
    //ts.set_cores(core_count);

    // Run and measure the wall-clock time taken to run the simulation
    const auto timestep_handle = sim_hw_timestep(total_timesteps, scheduler);
    sim_update_total_energy_and_counts(timestep_handle.get());
    // The total_messages_sent is incremented during the simulation since it's
    //  used to calculate the message id, so nothing needs to be done here
    if (spike_trace.is_open())
    {
        sim_trace_record_spikes(spike_trace, total_timesteps);
    }
    if (potential_trace.is_open())
    {
        sim_trace_record_potentials(potential_trace, total_timesteps);
    }

    return;
}

void sanafe::SpikingChip::sim_timestep_sync(Scheduler &scheduler) const
{
    // Simulate the end of timestep where all tiles synchronize. This
    //  could be using a fixed clock frequency or dynamically waiting until all
    //  cores finish processing. This can also account for other fixed costs
    //  on the chip e.g., initialization/monitoring by a housekeeping CPU
    //
    // TODO: in future we may want to have the option of simulating this
    //  dynamically e.g., by modeling synchronization messages in the NoC

    // Use a simple table based model for sync costs
    scheduler.timestep_sync_delay = ts_sync_delay_table.get(mapped_tiles);
}

void sanafe::SpikingChip::reset()
{
    for (Tile &tile : tiles)
    {
        for (Core &core : tile.cores)
        {
            core.timestep_buffer.clear();
            for (const auto &hw : core.pipeline_hw)
            {
                hw->reset();
            }
        }
    }

    for (auto &[group_name, neurons] : mapped_neuron_groups)
    {
        for (MappedNeuron &neuron : neurons)
        {
            neuron.status = invalid_neuron_state;
        }
    }
}

void sanafe::SpikingChip::retire_timestep(const Timestep &ts)
{
    total_sim_time += ts.sim_time;
}

double sanafe::SpikingChip::get_power() const noexcept
{
    double power = NAN; // Watts
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
    for (auto idx : core_list)
    {
        for (MappedNeuron &n : idx.get().neurons)
        {
            process_neuron(ts, n);
        }

        Core &core = idx;
        // Account for any remaining neuron processing
        const bool placeholder_event =
                (core.next_message_generation_delay != 0.0);
        if (placeholder_event)
        {
            const MappedNeuron &last_neuron = core.neurons.back();
            Message placeholder =
                    Message(placeholder_mid, *this, last_neuron, ts.timestep);
            placeholder.generation_delay = core.next_message_generation_delay;
            // Create a dummy placeholder message
            auto &message_queue = ts.messages;
            message_queue[core.id].push_back(std::move(placeholder));
        }
    }
}

void sanafe::SpikingChip::process_messages(Timestep &ts)
{
    // Assign outgoing spike messages to their respective destination
    //  cores, and calculate network costs
    auto &message_queue = ts.messages;
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
    // Older versions of OpenMP don't support range-based for loops yet
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (auto idx : core_list)
    {
#ifdef HAVE_OPENMP
        TRACE3(CHIP, "omp thread:%d\n", omp_get_thread_num());
#endif
        Core &core = idx;
        TRACE1(CHIP, "Processing %zu message(s) for cid:%zu\n",
                core.messages_in.size(), core.id);
        for (auto &m_ref : core.messages_in)
        {
            Message &m = m_ref;
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

    m.min_hop_delay = sim_estimate_network_costs(src_tile, dest_tile);
    m.hops = abs_diff(src_tile.x, dest_tile.x) +
            abs_diff(src_tile.y, dest_tile.y);

    Core &core = dest_tile.cores[m.dest_core_offset];
    core.messages_in.emplace_back(m);
}

void sanafe::SpikingChip::process_neuron(Timestep &ts, MappedNeuron &n)
{
    Core &c = *(n.core);
    const bool simulate_buffer = (c.pipeline_config.buffer_position ==
                                         buffer_before_dendrite_unit) ||
            (c.pipeline_config.buffer_position == buffer_before_soma_unit);

    PipelineResult input{};
    if (simulate_buffer)
    {
        // Read then clear the pipeline buffer entry
        input = c.timestep_buffer[n.mapped_address];
        c.timestep_buffer[n.mapped_address] = PipelineResult{};
    }
    const PipelineResult pipeline_output = execute_pipeline(
            n.neuron_processing_pipeline, ts, n, std::nullopt, input);
    n.core->next_message_generation_delay +=
            pipeline_output.latency.value_or(0.0);
    // INFO("nid:%s.%d +%e:%e\n", n.parent_group_name.c_str(), n.id,
    //         pipeline_output.latency.value_or(0.0),
    //         n.core->next_message_generation_delay);
    if (n.status == fired)
    {
        pipeline_process_axon_out(ts, n);
    }
}

double sanafe::SpikingChip::process_message(
        Timestep &ts, Core &core, Message &m)
{
    double message_processing_latency = pipeline_process_axon_in(core, m);

    assert(static_cast<size_t>(m.dest_axon_id) < core.axons_in.size());
    const AxonInModel &axon_in = core.axons_in[m.dest_axon_id];
    const PipelineResult
            empty_input{}; // Empty/default struct used as dummy input

    for (const size_t synapse_address : axon_in.synapse_addresses)
    {
        MappedConnection &con = *(core.connections_in[synapse_address]);

        // In certain pipeline configurations, every synaptic lookup requires
        //  updates to the dendrite and/or soma units as well. Keep propagating
        //  outputs/inputs until we hit the time-step buffer, where outputs
        //  are stored as inputs ready for the next time-step
        MappedNeuron &n = con.post_neuron_ref;
        const PipelineResult pipeline_output = execute_pipeline(
                con.message_processing_pipeline, ts, n, &con, empty_input);
        core.timestep_buffer[n.mapped_address] = pipeline_output;
        message_processing_latency += pipeline_output.latency.value_or(0.0);
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
    for (auto *unit : pipeline)
    {
        output = unit->process(ts, n, con, output);
        total_energy += output.energy.value_or(0.0);
        total_latency += output.latency.value_or(0.0);
        if (output.status != invalid_neuron_state)
        {
            n.status = output.status;
        }
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
    PipelineResult output = process_input_fn(ts, n, con, input);

    // Post-processing on outputs
    process_output_fn(n, con, output);

#ifndef NDEBUG
    check_outputs(n, output);
#endif
    energy += output.energy.value_or(0.0);

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

    TRACE1(CHIP, "nid:%s.%zu sending spike message to %zu axons out\n",
            n.parent_group_name.c_str(), n.id, n.axon_out_addresses.size());
    for (const size_t axon_address : n.axon_out_addresses)
    {
        // Save and increment atomically to ensure each message ID is unique,
        //  even when creating messages on multiple threads
        const long int id = total_messages_sent.fetch_add(1);
        Message m(id, *this, axon_address, n, ts.timestep);
        // Add axon access cost to message latency and energy
        AxonOutUnit &axon_out_hw = *(n.axon_out_hw);
        axon_out_hw.energy += axon_out_hw.energy_access;

        m.generation_delay = n.core->next_message_generation_delay +
                axon_out_hw.latency_access;
        n.core->next_message_generation_delay = 0.0;

        auto &message_queue = ts.messages;
        message_queue[n.core->id].push_back(std::move(m));
        ++axon_out_hw.packets_out;

        axon_result.energy.value() += axon_out_hw.energy_access;
        axon_result.latency.value() += axon_out_hw.latency_access;
    }

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
    out << "total_neurons_updated: " << run_data.neurons_updated << "\n";
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
    for (auto idx : core_list)
    {
        Core &core = idx;
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
}

void sanafe::SpikingChip::sim_update_ts_counters(Timestep &ts)
{
    for (auto &tile : tiles)
    {
        ts.total_hops += tile.hops;
        for (auto &c : tile.cores)
        {
            for (const auto &hw : c.pipeline_hw)
            {
                ts.spike_count += hw->spikes_processed;
                ts.neurons_updated += hw->neurons_updated;
                ts.neurons_fired += hw->neurons_fired;
            }
            for (const auto &axon_out : c.axon_out_hw)
            {
                ts.packets_sent += axon_out.packets_out;
            }
        }
    }
    TRACE1(CHIP, "Spikes sent: %ld\n", ts.spike_count);
}

sanafe::TimestepHandle sanafe::SpikingChip::sim_hw_timestep(
        const long int timestep, Scheduler &scheduler)
{
    // Start the next time-step, clear all buffers
    auto ts = TimestepHandle(timestep);
    Timestep &ts_data = ts.get(); // Get reference to timestep data
    ts_data.set_cores(core_count);
    sim_reset_measurements();

    auto neuron_processing_start_tm = std::chrono::high_resolution_clock::now();
    process_neurons(ts_data);
    auto neuron_processing_end_tm = std::chrono::high_resolution_clock::now();
    auto message_processing_start_tm = neuron_processing_end_tm;
    process_messages(ts_data);
    forced_updates(ts_data);
    // The timestep ends once cores are synchronized
    sim_timestep_sync(scheduler);
    auto message_processing_end_tm = std::chrono::high_resolution_clock::now();
    auto energy_calculation_start_tm = message_processing_end_tm;
    sim_calculate_ts_energy(ts_data);
    sim_update_ts_counters(ts_data);

    auto energy_calculation_end_tm = std::chrono::high_resolution_clock::now();
    auto scheduler_start_tm = energy_calculation_end_tm;
    if (scheduler.timing_model == timing_model_cycle_accurate)
    {
        check_booksim_compatibility(scheduler, chip_count);
    }
    schedule_messages(ts, scheduler, *booksim_config);
    auto scheduler_end_tm = std::chrono::high_resolution_clock::now();

    // Calculate various simulator timings
    constexpr double ns_to_second = 1.0e-9;
    const long int neuron_processing_cycles =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                    neuron_processing_end_tm - neuron_processing_start_tm)
                    .count();
    neuron_processing_wall +=
            static_cast<double>(neuron_processing_cycles) * ns_to_second;
    const long int message_processing_cycles =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                    message_processing_end_tm - message_processing_start_tm)
                    .count();
    message_processing_wall +=
            static_cast<double>(message_processing_cycles) * ns_to_second;
    const long int scheduler_cycles =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                    scheduler_end_tm - scheduler_start_tm)
                    .count();
    scheduler_wall += static_cast<double>(scheduler_cycles) * ns_to_second;
    const long int energy_stats_cycles =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                    energy_calculation_end_tm - energy_calculation_start_tm)
                    .count();
    energy_stats_wall +=
            static_cast<double>(energy_stats_cycles) * ns_to_second;
    TRACE1(CHIP, "neuron:%e message:%e scheduler:%e energy:%e\n",
            neuron_processing_wall, message_processing_wall, scheduler_wall,
            energy_stats_wall);

    return ts;
}

double sanafe::SpikingChip::sim_estimate_network_costs(
        const Tile &src, Tile &dest)
{
    double network_latency = NAN;
    size_t x_hops = 0;
    size_t y_hops = 0;

    network_latency = 0.0;

    // Calculate the energy and time for sending spike packets
    x_hops = abs_diff(src.x, dest.x);
    y_hops = abs_diff(src.y, dest.y);
    // E-W hops

    if (src.x < dest.x)
    {
        dest.east_hops += x_hops;
        network_latency += static_cast<double>(x_hops) * src.latency_east_hop;
    }
    else
    {
        dest.west_hops += x_hops;
        network_latency += static_cast<double>(x_hops) * src.latency_west_hop;
    }

    // N-S hops
    if (src.y < dest.y)
    {
        dest.north_hops += y_hops;
        network_latency += static_cast<double>(y_hops) * src.latency_north_hop;
    }
    else
    {
        dest.south_hops += y_hops;
        network_latency += static_cast<double>(y_hops) * src.latency_south_hop;
    }

    dest.hops += (x_hops + y_hops);
    dest.messages_received++;
    TRACE1(CHIP, "xhops:%ld yhops%ld total hops:%ld latency:%e\n", x_hops,
            y_hops, x_hops + y_hops, network_latency);
    return network_latency;
}

void sanafe::SpikingChip::sim_calculate_ts_energy(Timestep &ts)
{
    // Returns the total energy across the design, for this timestep
    ts.network_energy = 0.0;
    ts.synapse_energy = 0.0;
    ts.dendrite_energy = 0.0;
    ts.soma_energy = 0.0;
    ts.total_energy = 0.0;

    for (auto &tile : tiles)
    {
        ts.total_energy += sim_calculate_tile_energy(ts, tile);
    }

    TRACE1(CHIP, "total_energy:%e\n", ts.total_energy);
    TRACE1(CHIP, "\tnetwork_energy:%e\n", ts.network_energy);
}

double sanafe::SpikingChip::sim_calculate_tile_energy(Timestep &ts, Tile &tile)
{
    double total_hop_energy =
            (static_cast<double>(tile.east_hops) * tile.energy_east_hop);
    total_hop_energy +=
            (static_cast<double>(tile.west_hops) * tile.energy_west_hop);
    total_hop_energy +=
            (static_cast<double>(tile.south_hops) * tile.energy_south_hop);
    total_hop_energy +=
            (static_cast<double>(tile.north_hops) * tile.energy_north_hop);
    tile.energy = total_hop_energy;
    ts.network_energy += total_hop_energy;
    TRACE1(CHIP, "east:%ld west:%ld north:%ld south:%ld\n", tile.east_hops,
            tile.west_hops, tile.north_hops, tile.south_hops);

    for (auto &core : tile.cores)
    {
        tile.energy += sim_calculate_core_energy(ts, core);
    }

    return tile.energy;
}

double sanafe::SpikingChip::sim_calculate_core_energy(Timestep &ts, Core &core)
{
    double axon_in_energy{0.0};
    double axon_out_energy{0.0};
    for (const auto &axon : core.axon_in_hw)
    {
        axon_in_energy = static_cast<double>(axon.spike_messages_in) *
                axon.energy_spike_message;
        TRACE1(CHIP, "spikes in: %ld, energy:%e\n", axon.spike_messages_in,
                axon.energy_spike_message);
    }
    ts.network_energy += axon_in_energy;
    ts.network_energy += axon_out_energy;

    double pipeline_energy{0.0};
    for (auto &pipeline_unit : core.pipeline_hw)
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

    for (const auto &axon : core.axon_out_hw)
    {
        axon_out_energy = axon.energy;
        TRACE1(CHIP, "packets: %ld, energy per packet:%e\n", axon.packets_out,
                axon.energy_access);
    }

    core.energy = axon_in_energy;
    core.energy += pipeline_energy;
    core.energy += axon_out_energy;

    return core.energy;
}

void sanafe::SpikingChip::sim_create_neuron_axons(MappedNeuron &pre_neuron)
{
    // Setup the connections between neurons and map them to hardware
    assert(pre_neuron.core != nullptr);

    // Figure out the unique set of cores that this neuron broadcasts to
    TRACE1(CHIP, "Counting connections for neuron nid:%zu\n", pre_neuron.id);
    std::set<Core *> cores_out;
    for (const MappedConnection &curr_connection : pre_neuron.connections_out)
    {
        const MappedNeuron &post_neuron = curr_connection.post_neuron_ref;
        Core *dest_core = post_neuron.core;
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
        const MappedNeuron &post_neuron = curr_connection.post_neuron_ref;
        Core &post_core = *(post_neuron.core);
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
            post_core.axons_out.size() - 1UL);

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

void sanafe::SpikingChip::sim_print_axon_summary() const noexcept
{
    [[maybe_unused]] size_t in_count = 0UL;
    [[maybe_unused]] size_t out_count = 0UL;

    INFO("** Mapping summary **\n");
    for (const Tile &tile : tiles)
    {
        // For debug only, print the axon maps
        for (const Core &core : tile.cores)
        {
#if (DEBUG_LEVEL_CHIP == 2)
            for (const auto &n : core.neurons)
            {
                TRACE2(CHIP, "\tnid:%s.%zu ", n.parent_group_name.c_str(),
                        n.id);
                TRACE2(CHIP, "i:%d o:%d\n", n.maps_in_count, n.maps_out_count);
            }
#endif
            if (!core.neurons.empty())
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
    post_core.axons_in.emplace_back();
    const size_t new_axon_in_address = post_core.axons_in.size() - 1;

    // Add the axon at the sending, pre-synaptic core
    TRACE1(CHIP, "Axon in address:%zu for core:%zu.%zu\n", new_axon_in_address,
            post_core.parent_tile_id, post_core.id);
    AxonOutModel out;
    out.dest_axon_id = new_axon_in_address;
    out.dest_core_offset = post_core.offset;
    out.dest_tile_id = post_core.parent_tile_id;
    out.src_neuron_offset = pre_neuron.offset;
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
    // Reset any energy, time latency or other measurements of network hardware
    //  This is called every timestep
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

            for (auto &axon : c.axon_in_hw)
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
                hw->neurons_updated = 0L;
                hw->neurons_fired = 0L;
            }

            for (auto &axon : c.axon_out_hw)
            {
                axon.energy = 0.0;
                axon.time = 0.0;
                axon.packets_out = 0;
            }

            // Reset the message buffer
            c.messages_in = std::vector<std::reference_wrapper<Message>>();
        }
    }
}

void sanafe::SpikingChip::sim_trace_write_spike_header(
        std::ostream &spike_trace_file)
{
    assert(spike_trace_file.good());
    spike_trace_file << "neuron,timestep" << '\n';
}

void sanafe::SpikingChip::sim_trace_write_potential_header(
        std::ostream &potential_trace_file)
{
    // Write csv header for probe outputs - record which neurons have been
    //  probed
    assert(potential_trace_file.good());
    potential_trace_file << "timestep,";
    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
    {
        for (const MappedNeuron &neuron : group_neurons)
        {
            if (neuron.log_potential)
            {
                potential_trace_file << "neuron " << group_name;
                potential_trace_file << "." << neuron.offset << ",";
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
        std::ostream &perf_trace_file)
{
    assert(perf_trace_file.good());
    perf_trace_file << "timestep,";

    // Mandatory performance metrics
    perf_trace_file << "fired,";
    perf_trace_file << "updated,";
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
        std::ostream &message_trace_file)
{
    assert(message_trace_file.good());
    message_trace_file << "timestep,";
    message_trace_file << "mid,";
    message_trace_file << "src_neuron,";
    message_trace_file << "src_hw,";
    message_trace_file << "dest_hw,";
    message_trace_file << "hops,";
    message_trace_file << "spikes,";
    message_trace_file << "send_timestamp,";
    message_trace_file << "received_timestamp,";
    message_trace_file << "processed_timestamp,";
    message_trace_file << "generation_delay,";
    message_trace_file << "processing_delay,";
    message_trace_file << "network_delay,";
    message_trace_file << "blocking_delay,";
    message_trace_file << "min_hop_delay\n";
    message_trace_file.flush();
}

void sanafe::SpikingChip::sim_trace_record_spikes(
        std::ostream &spike_trace_file, const long int timestep)
{
    // A trace of all spikes that are generated
    assert(spike_trace_file.good());

    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
    {
        for (const MappedNeuron &neuron : group_neurons)
        {
            if (neuron.log_spikes && (neuron.status == sanafe::fired))
            {
                spike_trace_file << neuron.parent_group_name << ".";
                spike_trace_file << neuron.offset << ",";
                spike_trace_file << timestep;
                spike_trace_file << "\n";
                spike_trace_file.flush();
            }
        }
    }
}

void sanafe::SpikingChip::sim_trace_record_potentials(
        std::ostream &potential_trace_file, const long int timestep)
{
    // Each line of this csv file is the potential of all probed neurons for
    //  one time-step
    assert(potential_trace_file.good());
    TRACE1(CHIP, "Recording potential for timestep: %ld\n", timestep);
    potential_trace_file << timestep << ",";

    long int potential_probe_count = 0;
    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
    {
        for (const MappedNeuron &neuron : group_neurons)
        {
            if (neuron.log_potential)
            {
                potential_trace_file
                        << neuron.soma_hw->get_potential(neuron.mapped_address);
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
        std::ostream &perf_trace_file, const Timestep &ts)
{
    perf_trace_file << ts.timestep << ",";
    // Start with the mandatory counters and traces
    perf_trace_file << ts.neurons_fired << ",";
    perf_trace_file << ts.neurons_updated << ",";
    perf_trace_file << ts.packets_sent << ",";
    perf_trace_file << ts.total_hops << ",";
    perf_trace_file << ts.spike_count << ",";
    perf_trace_file << std::scientific << ts.sim_time << ",";
    perf_trace_file << std::scientific << ts.total_energy;
    // Finish with the optional traces, which can be enabled or disabled per
    //  tile, core, or pipeline h/w unit
    const auto optional_perf_traces = sim_trace_get_optional_traces();
    for (const auto &[name, trace] : optional_perf_traces)
    {
        perf_trace_file << "," << std::scientific << trace;
    }
    perf_trace_file << "\n";
    perf_trace_file.flush();
}

void sanafe::sim_trace_record_message(
        std::ostream &message_trace_file, const Message &m)
{
    assert(message_trace_file.good());
    message_trace_file << m.timestep << ",";
    message_trace_file << m.mid << ",";
    message_trace_file << m.src_neuron_group_id << ".";
    message_trace_file << m.src_neuron_offset << ",";
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
    message_trace_file << m.sent_timestamp << ",";
    message_trace_file << m.received_timestamp << ",";
    message_trace_file << m.processed_timestamp << ",";
    message_trace_file << m.generation_delay << ",";
    message_trace_file << m.receive_delay << ",";
    message_trace_file << m.network_delay << ",";
    message_trace_file << m.blocked_delay << ",";
    message_trace_file << m.min_hop_delay << "\n";
    message_trace_file.flush();
}

std::vector<sanafe::NeuronAddress> sanafe::SpikingChip::get_spikes() const
{
    std::vector<NeuronAddress> spikes;

    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
    {
        for (const MappedNeuron &neuron : group_neurons)
        {
            const NeuronAddress address{
                    neuron.parent_group_name, neuron.offset};
            if (neuron.log_spikes && (neuron.status == sanafe::fired))
            {
                spikes.emplace_back(address);
            }
        }
    }

    return spikes;
}

std::vector<double> sanafe::SpikingChip::get_potentials() const
{
    std::vector<double> potentials;

    for (const auto &[group_name, group_neurons] : mapped_neuron_groups)
    {
        for (const MappedNeuron &neuron : group_neurons)
        {
            const NeuronAddress address{
                    neuron.parent_group_name, neuron.offset};
            if (neuron.log_potential)
            {
                potentials.emplace_back(
                        neuron.soma_hw->get_potential(neuron.mapped_address));
            }
        }
    }

    return potentials;
}

void sanafe::SpikingChip::check_booksim_compatibility(
        const Scheduler &scheduler, const int /*sim_count*/)
{
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
}

sanafe::TimingModel sanafe::parse_timing_model(
        const std::string_view &timing_model_str)
{
    sanafe::TimingModel timing_model{
            sanafe::TimingModel::timing_model_detailed};

    if (timing_model_str == "simple")
    {
        timing_model = sanafe::timing_model_simple;
    }
    else if (timing_model_str == "detailed")
    {
        timing_model = sanafe::timing_model_detailed;
    }
    else if (timing_model_str == "cycle")
    {
        timing_model = sanafe::timing_model_cycle_accurate;
    }
    else
    {
        INFO("Error: Timing model %s not recognized, default is 'detailed'.\n",
                std::string(timing_model_str).c_str());
    }

    return timing_model;
}
