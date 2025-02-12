// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  chip.cpp
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
#include "chip.hpp"
#include "models.hpp"
#include "network.hpp"
#include "plugins.hpp"
#include "print.hpp"
#include "schedule.hpp"

sanafe::SpikingChip::SpikingChip(const Architecture &arch,
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

sanafe::SpikingChip::~SpikingChip()
{
    // Close any open trace files
    spike_trace.close();
    potential_trace.close();
    perf_trace.close();
    message_trace.close();
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
        if (!neuron->core_id.has_value())
        {
            std::string error = "Neuron: " + neuron->parent_group_id + "." +
                    std::to_string(neuron->id) + " not mapped.";
            INFO("%s", error.c_str());
            throw std::runtime_error(error);
        }
        TRACE1(CHIP, "Mapping neuron %s.%zu to core:%zu\n",
                neuron->parent_group_id.c_str(), neuron->id,
                neuron->core_id.value());
        Core &mapped_core = list_of_cores[neuron->core_id.value()];
        mapped_core.map_neuron(*neuron);
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
                                [mapped_neuron.id] = &mapped_neuron;
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

    // TODO: it makes sense to build the axons as we go, rather than
    //  running in another loop, refactor and restructure this
    // TODO: however, map_axons assumes that connections in are fixed and
    //  determined; if we resize these vectors, the references become
    //  invalidated
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

    // TODO: for debug, remove
    /*
    for (auto &t: tiles)
    {
        for (auto &c: t.cores)
        {
            INFO("**CORE:%zu**\n", c.id);
            for (size_t axon_address = 0; axon_address < c.axons_in.size(); ++axon_address)
            {
                auto &axon = c.axons_in[axon_address];
                INFO("AXON:%zu\n", axon_address);
                for (int address: axon.synapse_addresses)
                {
                    auto &syn = c.connections_in[address];
                    INFO("%d: %s.%zu->%s.%zu cx:%zu w:%d\n", address,
                            syn->pre_neuron->parent_group_name.c_str(),
                            syn->pre_neuron->id,
                            syn->post_neuron->parent_group_name.c_str(),
                            syn->post_neuron->id,
                            syn->post_neuron->mapped_address,
                            (int) syn->synapse_hw->weight(syn->synapse_address));
                }
            }
        }
    }
    */

    return;
}

sanafe::MappedConnection &sanafe::SpikingChip::map_connection(
        const Connection &con)
{
    auto list_of_cores = cores();

    auto &pre_group = mapped_neuron_groups.at(con.pre_neuron.group_name);
    MappedNeuron &pre_neuron = *(pre_group[con.pre_neuron.neuron_id.value()]);

    auto &post_group = mapped_neuron_groups.at(con.post_neuron.group_name);
    MappedNeuron &post_neuron =
            *(post_group[con.post_neuron.neuron_id.value()]);

    pre_neuron.connections_out.emplace_back(pre_neuron.connections_out.size());
    MappedConnection &mapped_con = pre_neuron.connections_out.back();
    mapped_con.pre_neuron = &pre_neuron;
    mapped_con.post_neuron = &post_neuron;
    mapped_con.dendrite_params = con.dendrite_params;

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
    sim_print_axon_summary(*this);
}

sanafe::RunData::RunData(const long int start, const long int steps)
        : timestep_start(start)
        , timesteps_executed(steps)
{
}

sanafe::RunData sanafe::SpikingChip::sim(
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
        TRACE1(CHIP, "Neurons fired in ts:%ld: %zu\n", timestep,
                ts.neurons_fired);
        rd.energy += ts.energy;
        rd.sim_time += ts.sim_time;
        rd.spikes += ts.spike_count;
        rd.packets_sent += ts.packets_sent;
        rd.neurons_fired += ts.neurons_fired;
        rd.wall_time = wall_time;
    }

    return rd;
}

sanafe::Timestep sanafe::SpikingChip::step()
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
    TRACE1(CHIP, "Time-step took: %fs.\n",
            static_cast<double>(
                    ts_elapsed.tv_sec + (ts_elapsed.tv_nsec / ns_in_second)));

    return ts;
}

void sanafe::SpikingChip::reset()
{
    for (auto &[group_name, neurons] : mapped_neuron_groups)
    {
        for (MappedNeuron *neuron : neurons)
        {
            // TODO: reset all h/w units, axon, synapse and dendrite too
            neuron->dendrite_hw->reset();
            neuron->soma_hw->reset();

            neuron->status = IDLE;
            neuron->dendrite_input_synapses.clear();
            neuron->soma_input_charge = 0.0;
            neuron->axon_out_input_spike = false;
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

// Pipeline models
void sanafe::pipeline_process_neurons(Timestep &ts, SpikingChip &hw)
{
    auto cores = hw.cores();

    // Older versions of OpenMP don't support range-based for loops yet...
#pragma omp parallel for schedule(dynamic)
    // codechecker_suppress [modernize-loop-convert]
    for (size_t idx = 0; idx < cores.size(); idx++)
    {
        Core &core = cores[idx];
        for (MappedNeuron &n : core.neurons)
        {
            pipeline_process_neuron(ts, hw, n);
        }

        if (core.next_message_generation_delay != 0.0)
        {
            // This message accounts for any remaining neuron processing
            const MappedNeuron &last_neuron = core.neurons.back();
            Message placeholder(hw, last_neuron, ts.timestep);
            placeholder.generation_delay = core.next_message_generation_delay;
            // Create a dummy placeholder message
            ts.messages[core.id].push_back(placeholder);
        }
    }
}

void sanafe::pipeline_process_messages(Timestep &ts, SpikingChip &hw)
{
    // Assign outgoing spike messages to their respective destination
    //  cores, and calculate network costs
    for (auto &q : ts.messages)
    {
        for (auto &m : q)
        {
            if (!m.placeholder)
            {
                pipeline_receive_message(hw, m);
            }
        }
    }

    // Now process all messages at receiving cores
    auto cores = hw.cores();
    // Older versions of OpenMP don't support range-based for loops yet...
#pragma omp parallel for schedule(dynamic)
    // codechecker_suppress [modernize-loop-convert]
    for (size_t idx = 0; idx < cores.size(); idx++)
    {
        Core &core = cores[idx];
        TRACE1(CHIP, "Processing %zu message(s) for cid:%zu\n",
                core.messages_in.size(), core.id);
        for (auto *m : core.messages_in)
        {
            m->receive_delay += pipeline_process_message(ts, core, *m);
        }
    }
}

void sanafe::pipeline_receive_message(SpikingChip &hw, Message &m)
{
    assert(static_cast<size_t>(m.src_tile_id) < arch.tiles.size());
    assert(static_cast<size_t>(m.dest_tile_id) < arch.tiles.size());
    const Tile &src_tile = hw.tiles[m.src_tile_id];
    Tile &dest_tile = hw.tiles[m.dest_tile_id];
    m.network_delay = sim_estimate_network_costs(src_tile, dest_tile);
    m.hops = abs_diff(src_tile.x, dest_tile.x) +
            abs_diff(src_tile.y, dest_tile.y);

    Core &core = dest_tile.cores[m.dest_core_offset];
    core.messages_in.push_back(&m);
}

void sanafe::pipeline_process_neuron(
        Timestep &ts, const SpikingChip &arch, MappedNeuron &n)
{
    TRACE1(CHIP, "Processing neuron: %s.%zu\n", n.parent_group_name.c_str(),
            n.id);
    double neuron_processing_latency = 0.0;

    // Update any H/W following the time-step buffer in pipeline order
    if (n.core->pipeline_config.buffer_position <=
                BUFFER_BEFORE_DENDRITE_UNIT)
    {
        neuron_processing_latency += pipeline_process_dendrite(ts, n);
    }
    if (n.core->pipeline_config.buffer_position <= BUFFER_BEFORE_SOMA_UNIT)
    {
        neuron_processing_latency += pipeline_process_soma(ts, n);
    }
    if (n.core->pipeline_config.buffer_position <= BUFFER_BEFORE_AXON_OUT_UNIT)
    {
        neuron_processing_latency += pipeline_process_axon_out(ts, arch, n);
    }

    // Hardware before the time-step buffer in the message processing pipeline
    //  may also need to be updated every time-step. For example, if a leak has
    //  to be applied every time-step. Do this after neuron processing has
    //  finished
    // TODO: move this into the message receiving loops somehow
    if ((n.core->pipeline_config.buffer_position >
                BUFFER_BEFORE_DENDRITE_UNIT) &&
            n.force_dendrite_update)
    {
        neuron_processing_latency += pipeline_process_dendrite(ts, n);
    }
    if (n.core->pipeline_config.buffer_position > BUFFER_BEFORE_SOMA_UNIT &&
            n.force_soma_update)
    {
        neuron_processing_latency += pipeline_process_soma(ts, n);
    }

    // Account for latencies and reset counters
    n.core->next_message_generation_delay += neuron_processing_latency;
    n.spike_count = 0;
}

double sanafe::pipeline_process_message(
        const Timestep &ts, Core &core, Message &m)
{
    // Simulate message m in the message processing pipeline. The message is
    //  sequentially handled by units up to the time-step buffer
    TRACE1(CHIP, "Receiving message for cid:%zu\n", core.id);
    double message_processing_latency = pipeline_process_axon_in(core, m);

    assert(static_cast<size_t>(m.dest_axon_id) < core.axons_in.size());
    const AxonInModel &axon_in = core.axons_in[m.dest_axon_id];
    if (core.id == 31)
    {
        //INFO("Accessing %zu synapses\n", axon_in.synapse_addresses.size());
    }

    for (const int synapse_address : axon_in.synapse_addresses)
    {
        MappedConnection &con = *(core.connections_in[synapse_address]);
        message_processing_latency += pipeline_process_synapse(ts, con);
        if (core.pipeline_config.buffer_position == BUFFER_BEFORE_DENDRITE_UNIT)
        {
            continue; // Process next synapse
        }
        // In certain pipeline configurations, every synaptic lookup requires
        //  updates to the dendrite and/or soma units as well
        MappedNeuron &n = *(con.post_neuron);
        message_processing_latency += pipeline_process_dendrite(ts, n);

        if (core.pipeline_config.buffer_position == BUFFER_BEFORE_SOMA_UNIT)
        {
            continue; // Process next synapse
        }
        message_processing_latency += pipeline_process_soma(ts, n);
        assert(core.pipeline_config.buffer_position ==
                BUFFER_BEFORE_AXON_OUT_UNIT);
    }

    return message_processing_latency;
}

double sanafe::pipeline_process_axon_in(Core &core, const Message &m)
{
    assert(m.dest_axon_hw >= 0);
    assert(static_cast<size_t>(m.dest_axon_hw) < core.axon_in_hw.size());
    AxonInUnit &axon_unit = core.axon_in_hw[m.dest_axon_hw];
    axon_unit.spike_messages_in++;

    return axon_unit.latency_spike_message;
}

double sanafe::pipeline_process_synapse(
        const Timestep &ts, MappedConnection &con)
{
    // Update all synapses to different neurons in one core. If a synaptic
    //  lookup, read and accumulate the synaptic weights. Otherwise, just
    //  update filtered current and any other connection properties
    TRACE1(CHIP, "Updating synapses for (cid:%zu)\n", con.pre_neuron->core->id);
    con.synapse_hw->set_time(ts.timestep);
    Core &dest_core = *(con.post_neuron->core);

    auto [synaptic_current, simulated_energy, simulated_latency] =
            con.synapse_hw->update(con.synapse_address, true);
    if (con.synapse_hw->default_energy_process_spike.has_value())
    {
        simulated_energy = con.synapse_hw->default_energy_process_spike.value();
    }
    if (con.synapse_hw->default_latency_process_spike.has_value())
    {
        simulated_latency =
                con.synapse_hw->default_latency_process_spike.value();
    }

    // Input and buffer synaptic info at the next hardware unit (dendrite unit)
    Synapse synapse_data = {synaptic_current, con};
    con.post_neuron->dendrite_input_synapses.push_back(synapse_data);
    con.post_neuron->spike_count++;
    assert(con.synapse_hw != nullptr);
    con.synapse_hw->spikes_processed++;
    TRACE1(CHIP, "(nid:%s.%zu->nid:%s.%zu) current:%lf\n",
            con.pre_neuron->parent_group_name.c_str(), con.pre_neuron->id,
            con.post_neuron->parent_group_name.c_str(), con.post_neuron->id,
            synaptic_current);

    // Check that both the energy and latency costs have been set, either within
    //  the synapse h/w model, or by default metrics
    if (!simulated_energy.has_value())
    {
        INFO("Error: No synapse energy model or metrics. Either return an "
             "energy value or set the attribute: 'energy_process_spike'\n");
        throw std::runtime_error("Error: No synapse energy model or metrics.");
    }
    if (!simulated_latency.has_value())
    {
        INFO("Error: No synapse latency model or metrics. Either return a "
             "latency value or set the attribute: 'latency_process_spike'\n");
        throw std::runtime_error("Error: No synapse latency model or metrics.");
    }

    dest_core.energy += simulated_energy.value();
    return simulated_latency.value();
}

std::pair<double, double> sanafe::pipeline_apply_default_dendrite_power_model(
        MappedNeuron &n, std::optional<double> energy,
        std::optional<double> latency)
{
    // Apply default energy and latency metrics if set
    if (n.dendrite_hw->default_energy_update.has_value())
    {
        energy = n.dendrite_hw->default_energy_update.value();
    }
    if (n.dendrite_hw->default_latency_update.has_value())
    {
        latency = n.dendrite_hw->default_latency_update.value();
    }

    // Check that both the energy and latency costs have been set, either within
    //  the dendrite h/w model, or by default metrics
    if (!energy.has_value())
    {
        INFO("Error: No dendrite energy model or metrics.\n");
        throw std::runtime_error(
                "Error: No dendrite energy model or metrics.\n");
    }
    if (!latency.has_value())
    {
        INFO("Error: No dendrite latency model or metrics.\n");
        throw std::runtime_error(
                "Error: No dendrite latency model or metrics.\n");
    }

    return {energy.value(), latency.value()};
}

double sanafe::pipeline_process_dendrite(const Timestep &ts, MappedNeuron &n)
{
    n.dendrite_hw->set_time(ts.timestep);
    TRACE2(CHIP, "Updating nid:%zu (ts:%ld)\n", n.id, ts.timestep);

    double total_latency{0.0};
    if (n.dendrite_input_synapses.empty())
    {
        // Update the dendrite h/w at least once, even if there are no inputs.
        //  This might be required e.g., if there is a leak to apply
        auto [current, model_energy, model_latency] =
                n.dendrite_hw->update(n.mapped_address, std::nullopt);
        n.soma_input_charge = current;
        auto [energy, latency] = pipeline_apply_default_dendrite_power_model(
                n, model_energy, model_latency);
        n.core->energy += energy;
        total_latency += latency;
    }
    else
    {
        // Update the dendrite h/w for all input synaptic currents
        for (const auto &synapse : n.dendrite_input_synapses)
        {
            auto [current, model_energy, model_latency] =
                    n.dendrite_hw->update(n.mapped_address, synapse);
            n.soma_input_charge = current;
            auto [energy, latency] =
                    pipeline_apply_default_dendrite_power_model(
                            n, model_energy, model_latency);
            n.core->energy += energy;
            total_latency += latency;
        }
        n.dendrite_input_synapses.clear();
    }

    // Finally, send dendritic current to the soma
    TRACE2(CHIP, "nid:%zu updating dendrite, soma_input_charge:%lf\n", n.id,
            n.soma_input_charge);
    return total_latency;
}

double sanafe::pipeline_process_soma(const Timestep &ts, MappedNeuron &n)
{
    TRACE1(CHIP, "nid:%s.%zu updating, current_in:%lf (ts:%lu)\n",
            n.parent_group_name.c_str(), n.id, n.soma_input_charge,
            ts.timestep);
    n.soma_hw->set_time(ts.timestep);

    std::optional<double> soma_current_in;
    if ((n.spike_count > 0) || (std::fabs(n.soma_input_charge) > 0.0))
    {
        soma_current_in = n.soma_input_charge;
        n.soma_input_charge = 0.0;
    }

    auto [neuron_status, simulated_energy, simulated_latency] =
            n.soma_hw->update(n.mapped_address, soma_current_in);

    if (n.soma_hw->default_energy_metrics.has_value())
    {
        simulated_energy =
                n.soma_hw->default_energy_metrics->energy_access_neuron;
    }
    if (n.soma_hw->default_latency_metrics.has_value())
    {
        simulated_latency =
                n.soma_hw->default_latency_metrics->latency_access_neuron;
    }

    if (neuron_status == INVALID_NEURON_STATE)
    {
        std::string error = "Soma model for nid: " + n.parent_group_name + "." +
                std::to_string(n.id) + " returned invalid state.\n";
        throw std::runtime_error(error);
    }
    else if ((neuron_status == sanafe::UPDATED) ||
            (neuron_status == sanafe::FIRED))
    {
        n.soma_hw->neuron_updates++;
        if (n.soma_hw->default_energy_metrics.has_value())
        {
            simulated_energy.value() +=
                    n.soma_hw->default_energy_metrics->energy_update_neuron;
        }
        if (n.soma_hw->default_latency_metrics.has_value())
        {
            simulated_latency.value() +=
                    n.soma_hw->default_latency_metrics->latency_update_neuron;
        }
    }

    if (neuron_status == sanafe::FIRED)
    {
        if (n.soma_hw->default_energy_metrics.has_value())
        {
            simulated_energy.value() +=
                    n.soma_hw->default_energy_metrics->energy_spike_out;
        }
        if (n.soma_hw->default_latency_metrics.has_value())
        {
            simulated_latency.value() +=
                    n.soma_hw->default_latency_metrics->latency_spike_out;
        }

        n.soma_hw->neurons_fired++;
        n.axon_out_input_spike = true;
        TRACE1(CHIP, "Neuron %s.%zu fired\n", n.parent_group_name.c_str(),
                n.id);
    }

    n.status = neuron_status;
    TRACE1(CHIP, "neuron status:%d\n", n.status);

    // Check that both the energy and latency costs have been set, either within
    //  the synapse h/w model, or by default metrics
    if (!simulated_energy.has_value())
    {
        INFO("Error: No soma energy model or metrics.\n");
        throw std::runtime_error("Error: No soma energy model or metrics.");
    }
    if (!simulated_latency.has_value())
    {
        INFO("Error: No soma latency model or metrics.\n");
        throw std::runtime_error("Error: No soma latency model or metrics.");
    }
    n.core->energy += simulated_energy.value();
    return simulated_latency.value();
}

double sanafe::pipeline_process_axon_out(
        Timestep &ts, const SpikingChip &hw, MappedNeuron &n)
{
    if (!n.axon_out_input_spike)
    {
        return 0.0;
    }

    TRACE1(CHIP, "nid:%s.%zu sending spike message to %zu axons out\n",
            n.parent_group_name.c_str(), n.id, n.axon_out_addresses.size());
    for (const int axon_address : n.axon_out_addresses)
    {
        Message m(hw, n, ts.timestep, axon_address);
        // Add axon access cost to message latency and energy
        AxonOutUnit &axon_out_hw = *(n.axon_out_hw);

        m.generation_delay = n.core->next_message_generation_delay;
        n.core->next_message_generation_delay = 0.0;
        m.generation_delay += axon_out_hw.latency_access;
        ts.messages[n.core->id].push_back(m);
        axon_out_hw.packets_out++;
    }
    n.axon_out_input_spike = false;

    return n.axon_out_hw->latency_access;
}

sanafe::BufferPosition sanafe::pipeline_parse_buffer_pos_str(
        const std::string &buffer_pos_str)
{
    BufferPosition buffer_pos;
    if (buffer_pos_str == "dendrite")
    {
        buffer_pos = BUFFER_BEFORE_DENDRITE_UNIT;
    }
    else if (buffer_pos_str == "soma")
    {
        buffer_pos = BUFFER_BEFORE_SOMA_UNIT;
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

size_t sanafe::abs_diff(const size_t a, const size_t b)
{
    // Returns the absolute difference between two unsigned (size_t) values
    return (a > b) ? (a - b) : (b - a);
}

sanafe::RunData sanafe::SpikingChip::get_run_summary() const
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

void sanafe::SynapseUnit::configure(
        std::string synapse_name, const ModelInfo &model, size_t core_id)
{
    model_parameters = model.model_parameters;
    plugin_lib = model.plugin_library_path;
    name = synapse_name;
    // TODO: figure a clean way of doing this
    host_core_id = core_id;

    if (model_parameters.find("energy_process_spike") != model_parameters.end())
    {
        default_energy_process_spike =
                static_cast<double>(model_parameters["energy_process_spike"]);
    }
    if (model_parameters.find("latency_process_spike") !=
            model_parameters.end())
    {
        default_latency_process_spike =
                static_cast<double>(model_parameters["latency_process_spike"]);
    }
}

void sanafe::SynapseUnit::add_connection(MappedConnection &con)
{
    // TODO: this is wasteful storing every synapse twice, once in the neuron
    //  and again in the synapse h/w unit
    mapped_connections_in.push_back(&con);
}

void sanafe::DendriteUnit::configure(
        std::string dendrite_name, const ModelInfo &model_details)
{
    model_parameters = model_details.model_parameters;
    plugin_lib = model_details.plugin_library_path;
    name = dendrite_name;
    model = model_details.name;
    if (model_parameters.find("energy_update") != model_parameters.end())
    {
        default_energy_update =
                static_cast<double>(model_parameters["energy_update"]);
    }
    if (model_parameters.find("latency_update") != model_parameters.end())
    {
        default_latency_update =
                static_cast<double>(model_parameters["latency_update"]);
    }
}

void sanafe::SomaUnit::configure(
        const std::string &soma_name, const ModelInfo &model_details)
{
    model_parameters = model_details.model_parameters;
    plugin_lib = model_details.plugin_library_path;
    name = soma_name;
    model = model_details.name;

    auto key_exists = [this](const std::string &key) {
        return model_parameters.find(key) != model_parameters.end();
    };

    const std::set<std::string> energy_metric_names{
            "energy_access_neuron", "energy_update_neuron", "energy_spike_out"};
    bool parse_energy_metrics = std::any_of(
            energy_metric_names.begin(), energy_metric_names.end(), key_exists);
    if (parse_energy_metrics)
    {
        for (const auto &metric : energy_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Error: Metric not defined: " + metric + "\n";
                INFO("%s", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaEnergyMetrics energy_metrics;
        energy_metrics.energy_access_neuron =
                static_cast<double>(model_parameters["energy_access_neuron"]);
        energy_metrics.energy_update_neuron =
                static_cast<double>(model_parameters["energy_update_neuron"]);
        energy_metrics.energy_spike_out =
                static_cast<double>(model_parameters["energy_spike_out"]);
        default_energy_metrics = energy_metrics;
    }

    const std::set<std::string> latency_metric_names{"latency_access_neuron",
            "latency_update_neuron", "latency_spike_out"};
    bool parse_latency_metrics = std::any_of(latency_metric_names.begin(),
            latency_metric_names.end(), key_exists);
    if (parse_latency_metrics)
    {
        for (const auto &metric : latency_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Error: Missing metric: " + metric + "\n";
                INFO("%s", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaLatencyMetrics latency_metrics;
        latency_metrics.latency_access_neuron =
                static_cast<double>(model_parameters["latency_access_neuron"]);
        latency_metrics.latency_update_neuron =
                static_cast<double>(model_parameters["latency_update_neuron"]);
        latency_metrics.latency_spike_out =
                static_cast<double>(model_parameters["latency_spike_out"]);
        default_latency_metrics = latency_metrics;
    }
}

sanafe::AxonOutUnit::AxonOutUnit(const AxonOutConfiguration &config)
        : name(std::move(config.name))
        , energy_access(config.metrics.energy_message_out)
        , latency_access(config.metrics.latency_message_out)
{
}

void sanafe::Core::map_neuron(const Neuron &neuron_to_map)
{
    TRACE1(CHIP, "Mapping nid:%s.%zu to core: %zu\n",
            neuron_to_map.parent_group_id.c_str(), neuron_to_map.id, id);

    if (neurons.size() >= pipeline_config.max_neurons_supported)
    {
        INFO("Error: Exceeded maximum neurons per core (%zu)",
                pipeline_config.max_neurons_supported);
        throw std::runtime_error("Error: Exceeded maximum neurons per core.");
    }

    // Map neuron model to dendrite and soma hardware units in this core.
    //  Search through all models implemented by this core and return the
    //  one that matches. If no dendrite / soma hardware is specified,
    //  default to the first one defined
    if (dendrite.empty())
    {
        INFO("Error: No dendrite units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No dendrite units defined");
    }
    DendriteUnit *mapped_dendrite = dendrite[0].get();
    if (neuron_to_map.dendrite_hw_name.length() > 0)
    {
        bool dendrite_found = false;
        for (auto &dendrite_hw : dendrite)
        {
            if (neuron_to_map.dendrite_hw_name == dendrite_hw->name)
            {
                mapped_dendrite = dendrite_hw.get();
                dendrite_found = true;
            }
        }
        if (!dendrite_found)
        {
            INFO("Error: Could not map neuron nid:%zu (hw:%s) "
                 "to any dendrite h/w.\n",
                    neuron_to_map.id, neuron_to_map.dendrite_hw_name.c_str());
            throw std::runtime_error(
                    "Error: Could not map neuron to dendrite h/w");
        }
    }

    if (soma.empty())
    {
        INFO("Error: No soma units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No soma units defined");
    }
    SomaUnit *mapped_soma = soma[0].get();
    if (neuron_to_map.soma_hw_name.length() > 0)
    {
        bool soma_found = false;
        for (auto &soma_hw : soma)
        {
            if (neuron_to_map.soma_hw_name == soma_hw->name)
            {
                mapped_soma = soma_hw.get();
                soma_found = true;
            }
        }
        if (!soma_found)
        {
            INFO("Error: Could not map neuron nid:%zu (hw:%s) "
                 "to any soma h/w.\n",
                    neuron_to_map.id, neuron_to_map.soma_hw_name.c_str());
            throw std::runtime_error("Error: Could not map neuron to soma h/w");
        }
    }
    mapped_soma->neuron_count++;

    if (axon_out_hw.empty())
    {
        INFO("Error: No axon out units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No axon out units defined");
    }
    AxonOutUnit *mapped_axon_out = &(axon_out_hw[0]);

    // Map the neuron to the core and its hardware units
    const size_t address = neurons.size();
    neurons.emplace_back(neuron_to_map, this, address, mapped_dendrite,
            mapped_soma, mapped_axon_out);

    return;
}

sanafe::AxonInUnit &sanafe::Core::create_axon_in(
        const AxonInConfiguration &config)
{
    axon_in_hw.emplace_back(config);
    TRACE1(CHIP, "New axon in h/w unit created (%zu.%zu)\n", parent_tile_id,
            id);

    return axon_in_hw.back();
}

sanafe::SynapseUnit &sanafe::Core::create_synapse(
        const SynapseConfiguration &config)
{
    // Create the synapse model
    if (config.model_info.plugin_library_path.has_value())
    {
        const std::filesystem::path plugin_lib_path =
                config.model_info.plugin_library_path.value();
        TRACE1(CHIP, "Creating synapse from plugin: %s.\n",
                plugin_lib_path.c_str());
        synapse.emplace_back(
                plugin_get_synapse(config.model_info.name, plugin_lib_path));
    }
    else
    {
        // Use built in models
        TRACE1(CHIP, "Creating synapse built-in model %s.\n",
                config.model_info.name.c_str());
        synapse.emplace_back(model_get_synapse(config.model_info.name));
    }

    auto &new_unit = synapse.back();
    new_unit->configure(config.name, config.model_info, id);
    TRACE1(CHIP, "New synapse h/w unit created\n");

    return *new_unit;
}

sanafe::DendriteUnit &sanafe::Core::create_dendrite(
        const DendriteConfiguration &config)
{
    TRACE1(CHIP, "New dendrite h/w unit created\n");

    if (config.model_info.plugin_library_path.has_value())
    {
        const std::filesystem::path &plugin_library_path =
                config.model_info.plugin_library_path.value();
        TRACE1(CHIP, "Creating dendrite from plugin %s.\n",
                plugin_library_path.c_str());
        dendrite.emplace_back(plugin_get_dendrite(
                config.model_info.name, plugin_library_path));
    }
    else
    {
        // Use built in models
        TRACE1(CHIP, "Creating dendrite built-in model %s.\n",
                config.model_info.name.c_str());
        dendrite.emplace_back(model_get_dendrite(config.model_info.name));
    }

    auto &unit = dendrite.back();
    unit->configure(config.name, config.model_info);
    return *unit;
}

sanafe::SomaUnit &sanafe::Core::create_soma(const SomaConfiguration &config)
{
    TRACE1(CHIP, "New soma h/w unit created (%s)\n", config.name.c_str());

    if (config.model_info.plugin_library_path.has_value())
    {
        // Use external plug-in
        auto &plugin_library_path =
                config.model_info.plugin_library_path.value();
        TRACE1(CHIP, "Creating soma from plugin %s.\n",
                plugin_library_path.c_str());
        soma.emplace_back(
                plugin_get_soma(config.model_info.name, plugin_library_path));
    }
    else
    {
        // Use built-in unit models
        TRACE1(CHIP, "Creating soma built-in model %s.\n",
                config.model_info.name.c_str());
        soma.emplace_back(model_get_soma(config.model_info.name));
    }
    auto &unit = soma.back();
    //INFO("Model created, now configuring\n");
    unit->configure(config.name, config.model_info);

    return *unit;
}

sanafe::AxonOutUnit &sanafe::Core::create_axon_out(
        const AxonOutConfiguration &config)
{
    axon_out_hw.emplace_back(config);
    TRACE1(CHIP, "New axon out h/w unit created: (%zu.%zu)\n", parent_tile_id,
            id);

    return axon_out_hw.back();
}

sanafe::Message::Message(
        const SpikingChip &hw, const MappedNeuron &n, const long int timestep)
        : timestep(timestep)
        , src_neuron_id(n.id)
        , src_neuron_group_id(n.parent_group_name)
{
    // If no axon was given create a message with no destination. By
    //  default, messages without destinations act as a placeholder for neuron
    //  processing
    const Core &src_core = *(n.core);
    const Tile &src_tile = hw.tiles[src_core.parent_tile_id];
    src_x = src_tile.x;
    src_y = src_tile.y;
    src_tile_id = src_tile.id;
    src_core_id = src_core.id;
    src_core_offset = src_core.offset;
}

sanafe::Message::Message(const SpikingChip &hw, const MappedNeuron &n,
        const long int timestep, const int axon_address)
        : Message(hw, n, timestep)
{
    const Core &src_core = *(n.core);
    const AxonOutModel &src_axon = src_core.axons_out[axon_address];
    const Tile &dest_tile = hw.tiles[src_axon.dest_tile_id];
    const Core &dest_core = dest_tile.cores[src_axon.dest_core_offset];
    const AxonInModel &dest_axon = dest_core.axons_in[src_axon.dest_axon_id];

    placeholder = false;
    spikes = dest_axon.synapse_addresses.size();
    dest_x = dest_tile.x;
    dest_y = dest_tile.y;
    dest_tile_id = dest_tile.id;
    dest_core_id = dest_core.id;
    dest_core_offset = dest_core.offset;
    dest_axon_id = src_axon.dest_axon_id;
    // TODO: support multiple axon output units, included in the synapse
    dest_axon_hw = 0;
}

sanafe::Tile::Tile(const TileConfiguration &config)
        : name(config.name)
        , energy_north_hop(config.power_metrics.energy_north_hop)
        , latency_north_hop(config.power_metrics.latency_north_hop)
        , energy_east_hop(config.power_metrics.energy_east_hop)
        , latency_east_hop(config.power_metrics.latency_east_hop)
        , energy_south_hop(config.power_metrics.energy_south_hop)
        , latency_south_hop(config.power_metrics.latency_south_hop)
        , energy_west_hop(config.power_metrics.energy_west_hop)
        , latency_west_hop(config.power_metrics.latency_west_hop)
        , id(config.id)
        , x(config.x)
        , y(config.y)
{
}

std::string sanafe::Tile::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Tile(tile=" << id << " cores=";
    ss << cores.size() << ")";

    return ss.str();
}

sanafe::Core::Core(const CoreConfiguration &config)
        : pipeline_config(config.pipeline)
        , name(std::move(config.name))
        , id(config.address.id)
        , offset(config.address.offset_within_tile)
        , parent_tile_id(config.address.parent_tile_id)
{
}

sanafe::MappedConnection::MappedConnection(const int connection_id)
        : post_neuron(nullptr)
        , pre_neuron(nullptr)
        , synapse_hw(nullptr)
        , id(connection_id)
{
}

sanafe::MappedNeuron::MappedNeuron(const Neuron &neuron_to_map,
        Core *mapped_core, const size_t address, DendriteUnit *mapped_dendrite,
        SomaUnit *mapped_soma, AxonOutUnit *mapped_axon_out)
        : parent_group_name{neuron_to_map.parent_group_id}
        , id(neuron_to_map.id)
        , core(mapped_core)
        , dendrite_hw(mapped_dendrite)
        , soma_hw(mapped_soma)
        , axon_out_hw(mapped_axon_out)
        , mapped_address(address)
        , mapping_order{neuron_to_map.mapping_order}
        , force_synapse_update{neuron_to_map.force_synapse_update}
        , force_dendrite_update{neuron_to_map.force_dendrite_update}
        , force_soma_update{neuron_to_map.force_soma_update}
        , log_spikes{neuron_to_map.log_spikes}
        , log_potential{neuron_to_map.log_potential}

{
    configure_models(neuron_to_map.model_parameters);
}

void sanafe::MappedNeuron::configure_models(
        const std::map<std::string, sanafe::ModelParam> &model_parameters)
{
    for (auto &[key, param] : model_parameters)
    {
        TRACE2(CHIP, "Forwarding param: %s (dendrite:%d soma:%d)\n",
                key.c_str(), param.forward_to_dendrite, param.forward_to_soma);
        if (param.forward_to_dendrite && (dendrite_hw != nullptr))
        {
            dendrite_hw->set_attribute(mapped_address, key, param);
        }
        if (param.forward_to_soma && (soma_hw != nullptr))
        {
            soma_hw->set_attribute(mapped_address, key, param);
        }
    }
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
        const std::filesystem::path &out_dir, const SpikingChip &hw)
{
    const std::filesystem::path potential_path = out_dir / "potentials.csv";
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

void sanafe::sim_timestep(Timestep &ts, SpikingChip &hw)
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
    //ts.sim_time = schedule_messages_simple(ts.messages, scheduler);
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

    TRACE1(CHIP, "Spikes sent: %ld\n", ts.spike_count);
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
    TRACE1(CHIP, "xhops:%ld yhops%ld total hops:%ld latency:%e\n", x_hops,
            y_hops, x_hops + y_hops, network_latency);
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
			TRACE2(CHIP, "noise val:%d\n", noise_val);

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

double sanafe::sim_calculate_energy(const SpikingChip &hw)
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
        TRACE1(CHIP, "east:%ld west:%ld north:%ld south:%ld\n", t.east_hops,
                t.west_hops, t.north_hops, t.south_hops);

        for (const auto &c : t.cores)
        {
            model_simulated_energy += c.energy;
            for (const auto &axon : c.axon_in_hw)
            {
                axon_in_energy += static_cast<double>(axon.spike_messages_in) *
                        axon.energy_spike_message;
                TRACE1(CHIP, "spikes in: %ld, energy:%e\n",
                        axon.spike_messages_in, axon.energy_spike_message);
            }
            // TODO: fix energy and latency for individual units
            /*
            for (const auto &syn : c.synapse)
            {
                synapse_energy += static_cast<double>(syn.spikes_processed) *
                        syn.energy_spike_op;
                TRACE1(CHIP, "synapse processed: %ld, energy:%e\n",
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
                TRACE1(CHIP, "neurons:%ld updates:%ld, spiking:%ld\n",
                        soma.neuron_count, soma.neuron_updates,
                        soma.neurons_fired);
            }
            */
            for (const auto &axon : c.axon_out_hw)
            {
                axon_out_energy += static_cast<double>(axon.packets_out) *
                        axon.energy_access;
                TRACE1(CHIP, "packets: %ld, energy:%e\n", axon.packets_out,
                        axon.energy_access);
            }
        }
    }

    // TODO: clean up energy breakdown. Include model simulated energies
    //  in their respective units, and track energies of cores, tiles etc
    total_energy = model_simulated_energy + axon_in_energy + synapse_energy +
            soma_energy + axon_out_energy + network_energy;

    TRACE1(CHIP, "model_simulated_energy:%e\n", model_simulated_energy);
    TRACE1(CHIP, "axon_in_energy:%e\n", axon_in_energy);
    TRACE1(CHIP, "synapse_energy:%e\n", synapse_energy);
    TRACE1(CHIP, "soma_energy:%e\n", soma_energy);
    TRACE1(CHIP, "axon_out_energy:%e\n", axon_out_energy);
    TRACE1(CHIP, "network_energy:%e\n", network_energy);
    TRACE1(CHIP, "total:%e\n", total_energy);

    return total_energy;
}

void sanafe::sim_create_neuron_axons(MappedNeuron &pre_neuron)
{
    // Setup the connections between neurons and map them to hardware
    assert(pre_neuron.core != nullptr);

    // Figure out the unique set of cores that this neuron broadcasts to
    TRACE1(CHIP, "Counting connections for neuron nid:%zu\n", pre_neuron.id);
    std::set<Core *> cores_out;
    for (MappedConnection &curr_connection : pre_neuron.connections_out)
    {
        TRACE1(CHIP, "Looking at connection id: %d\n", curr_connection.id);
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
        TRACE1(CHIP, "Adding connection:%d\n", curr_connection.id);
        sim_add_connection_to_axon(curr_connection, post_core);
    }
    TRACE1(CHIP, "Finished mapping connections to hardware for nid:%s.%zu.\n",
            pre_neuron.parent_group_name.c_str(), pre_neuron.id);
}

void sanafe::sim_add_connection_to_axon(MappedConnection &con, Core &post_core)
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
    con.synapse_hw->add_connection(con);

    // Access the most recently created axon in for the post-synaptic core
    AxonInModel &last_added_target_axon = post_core.axons_in.back();
    last_added_target_axon.synapse_addresses.push_back(con.synapse_address);
}

void sanafe::sim_print_axon_summary(SpikingChip &hw)
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

void sanafe::sim_reset_measurements(SpikingChip &hw)
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
        std::ofstream &potential_trace_file, const SpikingChip &hw)
{
    // Write csv header for probe outputs - record which neurons have been
    //  probed
    assert(potential_trace_file.is_open());
    potential_trace_file << "timestep,";
    for (const auto &[group_name, group_neurons] : hw.mapped_neuron_groups)
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
        const long int timestep, const SpikingChip &hw)
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
        const long int timestep, const SpikingChip &hw)
{
    // Each line of this csv file is the potential of all probed neurons for
    //  one time-step
    assert(potential_trace_file.is_open());
    TRACE1(CHIP, "Recording potential for timestep: %ld\n", timestep);
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
	TRACE1(CHIP, "nid:%d sending spike(s).\n", n->id);
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
