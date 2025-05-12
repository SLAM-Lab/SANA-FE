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

#include <booksim_lib.hpp>

#include "arch.hpp"
#include "chip.hpp"
#include "models.hpp"
#include "network.hpp"
#include "plugins.hpp"
#include "print.hpp"
#include "schedule.hpp"

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
    BookSimConfig new_config = booksim_load_config(booksim_config_vec);
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
        if (!neuron->core_id.has_value())
        {
            std::string error = "Neuron: " + neuron->parent_group_id + "." +
                    std::to_string(neuron->offset) + " not mapped.";
            INFO("%s", error.c_str());
            throw std::runtime_error(error);
        }
        TRACE1(CHIP, "Mapping neuron %s.%zu to core:%zu\n",
                neuron->parent_group_id.c_str(), neuron->offset,
                neuron->core_id.value());
        Core &mapped_core = list_of_cores[neuron->core_id.value()];
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

                for (auto &name_value_pair : con.synapse_params)
                {
                    if (name_value_pair.second.forward_to_synapse)
                    {
                        mapped_con.synapse_hw->set_attribute_edge(
                                mapped_con.synapse_address,
                                name_value_pair.first, name_value_pair.second);
                    }
                    if (name_value_pair.second.forward_to_dendrite)
                    {
                        MappedNeuron &n = *(mapped_con.post_neuron);
                        n.dendrite_hw->set_attribute_edge(
                                mapped_con.synapse_address,
                                name_value_pair.first, name_value_pair.second);
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

    pre_neuron.connections_out.emplace_back(pre_neuron.connections_out.size());
    MappedConnection &mapped_con = pre_neuron.connections_out.back();
    mapped_con.pre_neuron = &pre_neuron;
    mapped_con.post_neuron = &post_neuron;
    mapped_con.dendrite_params = con.dendrite_params;

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

sanafe::RunData sanafe::SpikingChip::sim(const long int timesteps,
        const long int heartbeat, const TimingModel timing_model)
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
        const Timestep ts = step(timing_model);
        TRACE1(CHIP, "Neurons fired in ts:%ld: %zu\n", timestep,
                ts.neurons_fired);
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
    }

    return rd;
}

sanafe::Timestep sanafe::SpikingChip::step(const TimingModel timing_model)
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
    sim_timestep(ts, timing_model);
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

    // Update global chip performance counters
    total_energy += ts.total_energy;
    synapse_energy += ts.synapse_energy;
    dendrite_energy += ts.dendrite_energy;
    soma_energy += ts.soma_energy;
    network_energy += ts.network_energy;

    total_sim_time += ts.sim_time;
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
    if (perf_trace_enabled)
    {
        sim_trace_record_perf(perf_trace, ts);
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
        Core &core = core_list[idx];
        if (core.pipeline_config.buffer_position < BUFFER_BEFORE_DENDRITE_UNIT)
        {
            throw std::logic_error("Error: Buffer must be after synaptic h/w");
        }
        for (MappedNeuron &n : core.neurons)
        {
            process_neuron(ts, n);
        }

        // Account for any remaining neuron processing
        bool placeholder_event = (core.next_message_generation_delay != 0.0);
        if (placeholder_event)
        {
            const MappedNeuron &last_neuron = core.neurons.back();
            Message placeholder(
                    placeholder_mid, *this, last_neuron, ts.timestep);
            placeholder.generation_delay = core.next_message_generation_delay;
            // Create a dummy placeholder message
            ts.messages[core.id].push_back(placeholder);
        }
    }
}

void sanafe::SpikingChip::process_messages(Timestep &ts)
{
    // Assign outgoing spike messages to their respective destination
    //  cores, and calculate network costs
    for (auto &q : ts.messages)
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
        Core &core = core_list[idx];
        TRACE1(CHIP, "Processing %zu message(s) for cid:%zu\n",
                core.messages_in.size(), core.id);
        for (auto *m : core.messages_in)
        {
            m->receive_delay += process_message(ts, core, *m);
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
    core.messages_in.push_back(&m);
}

void sanafe::MappedNeuron::build_neuron_processing_pipeline()
{
    if (core->pipeline_config.buffer_position <= BUFFER_INSIDE_DENDRITE_UNIT)
    {
        neuron_processing_pipeline.push_back(dendrite_hw);
    }
    if ((core->pipeline_config.buffer_position <= BUFFER_INSIDE_SOMA_UNIT) &&
            (soma_hw != dendrite_hw))
    {
        neuron_processing_pipeline.push_back(soma_hw);
    }
}

void sanafe::MappedConnection::build_message_processing_pipeline()
{
    MappedNeuron &n = *post_neuron;
    Core &mapped_core = *(n.core);

    // We don't support putting the buffer inside or before the synapse unit, so
    //  unconditionally push the synapse h/w. This is because putting the buffer
    //  here could cause a spike sent that shouldn't be
    message_processing_pipeline.push_back(synapse_hw);
    if ((mapped_core.pipeline_config.buffer_position >
                BUFFER_BEFORE_DENDRITE_UNIT) &&
            (n.dendrite_hw != synapse_hw))
    {
        message_processing_pipeline.push_back(n.dendrite_hw);
    }
    if ((mapped_core.pipeline_config.buffer_position >
                BUFFER_BEFORE_SOMA_UNIT) &&
            (n.soma_hw != n.dendrite_hw))
    {
        message_processing_pipeline.push_back(n.soma_hw);
    }

    return;
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

void sanafe::PipelineUnit::check_outputs(
        const MappedNeuron &n, const PipelineResult &result)
{
    // Check the hw returns a valid value for the next unit or network to
    //  process
    if (implements_soma && result.status == INVALID_NEURON_STATE)
    {
        throw std::runtime_error("Soma output; should return valid "
                                 "neuron state.\n");
    }
    else if (!implements_soma && (implements_synapse || implements_dendrite) &&
            !result.current.has_value())
    {
        throw std::runtime_error("Synaptic or dendritic output; should return "
                                 "synaptic/dendritic current\n");
    }

    return;
}

double sanafe::SpikingChip::process_message(
        Timestep &ts, Core &core, Message &m)
{
    double message_processing_latency = pipeline_process_axon_in(core, m);

    PipelineResult pipeline_output;

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
        output = unit->process_input(ts, n, con, output);
        total_energy += output.energy.value();
        total_latency += output.latency.value();
        if (output.status != INVALID_NEURON_STATE)
        {
            n.status = output.status;
        }
    }

    output.energy = total_energy;
    output.latency = total_latency;
    return output;
}

sanafe::PipelineResult sanafe::PipelineUnit::process_input(Timestep &ts,
        MappedNeuron &n, std::optional<MappedConnection *> con,
        const PipelineResult &input)
{
    TRACE2(CHIP, "Updating nid:%zu (ts:%ld)\n", n.id, ts.timestep);
    set_time(ts.timestep);

    PipelineResult output{};
    if (implements_synapse) // Synapse is input interface
    {
        if (!con.has_value())
        {
            throw std::logic_error(
                    "Error: Pipeline error, didn't receive "
                    "synaptic connection info. Check that no h/w unit is being "
                    "invoked before this one in the pipeline.");
        }
        output = update(con.value()->synapse_address, true);
        ++spikes_processed;
    }
    else if (implements_dendrite) // Dendrite is input interface
    {
        std::optional<size_t> synapse_address{std::nullopt};
        if (con.has_value() && (con.value() != nullptr))
        {
            synapse_address = con.value()->synapse_address;
        }
        output = update(n.mapped_address, input.current, synapse_address);
    }
    else if (implements_soma) // Soma is input interface
    {
        output = update(n.mapped_address, input.current);
    }

    if (implements_soma) // Soma is output interface
    {
        output = calculate_soma_default_energy_latency(n, output);
        update_soma_activity(n, output);
    }
    else if (implements_dendrite) // Dendrite is output interface
    {
        output = calculate_dendrite_default_energy_latency(n, output);
    }
    else if (implements_synapse) // Synapse is output interface
    {
        output = calculate_synapse_default_energy_latency(
                *(con.value()), output);
    }
    check_outputs(n, output);

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

sanafe::PipelineResult
sanafe::PipelineUnit::calculate_synapse_default_energy_latency(
        MappedConnection &con, const PipelineResult &simulation_result)
{
    PipelineResult updated_result{simulation_result};

    bool energy_simulated = simulation_result.energy.has_value();
    bool latency_simulated = simulation_result.latency.has_value();

    bool default_synapse_energy_metrics_set =
            con.synapse_hw->default_energy_process_spike.has_value();
    if (energy_simulated && default_synapse_energy_metrics_set)
    {
        std::string error("Error: Synapse unit simulates energy and also has "
                          "default energy metrics set.");
        throw std::logic_error(error);
    }
    if (default_synapse_energy_metrics_set)
    {
        updated_result.energy =
                con.synapse_hw->default_energy_process_spike.value();
    }

    bool default_synapse_latency_metrics_set =
            con.synapse_hw->default_latency_process_spike.has_value();
    if (latency_simulated && default_synapse_latency_metrics_set)
    {
        std::string error("Error: Synapse unit simulates latency and also has "
                          "default latency metrics set. Remove the default "
                          "metric from the architecture description.");
        throw std::logic_error(error);
    }

    if (default_synapse_latency_metrics_set)
    {
        if (simulation_result.latency.has_value())
        {
            std::string error(
                    "Error: Synapse unit simulates latency and also has "
                    "default latency metrics set. Remove the default "
                    "metric from the architecture description.");
            throw std::logic_error(error);
        }
        updated_result.latency =
                con.synapse_hw->default_latency_process_spike.value();
    }

    if (!updated_result.energy.has_value())
    {
        std::string error(
                "Error: Synapse unit does not simulate energy or provide "
                "a default energy cost in the architecture "
                "description.");
        throw std::logic_error(error);
    }
    if (!updated_result.latency.has_value())
    {
        std::string error("Error: Synapse unit does not simulate latency or "
                          "provide a default latency cost in the architecture "
                          "description.");
        throw std::logic_error(error);
    }

    return updated_result;
}

sanafe::PipelineResult
sanafe::PipelineUnit::calculate_dendrite_default_energy_latency(
        MappedNeuron &n, const PipelineResult &simulation_result)
{
    PipelineResult updated_result{simulation_result};

    bool energy_simulated = simulation_result.energy.has_value();
    bool latency_simulated = simulation_result.latency.has_value();

    bool default_dendrite_energy_metrics_set =
            n.dendrite_hw->default_energy_update.has_value();
    if (energy_simulated && default_dendrite_energy_metrics_set)
    {
        std::string error("Error: Dendrite unit simulates energy and also has "
                          "default energy metrics set.");
        throw std::logic_error(error);
    }
    if (default_dendrite_energy_metrics_set)
    {
        updated_result.energy = n.dendrite_hw->default_energy_update.value();
    }

    bool default_dendrite_latency_metrics_set =
            n.dendrite_hw->default_latency_update.has_value();
    if (latency_simulated && default_dendrite_latency_metrics_set)
    {
        std::string error("Error: Dendrite unit simulates latency and also has "
                          "default energy metrics set.");
        throw std::logic_error(error);
    }

    if (default_dendrite_latency_metrics_set)
    {
        if (simulation_result.latency.has_value())
        {
            std::string error(
                    "Error: Dendrite unit simulates energy and also has default energy metrics set.");
            throw std::logic_error(error);
        }
        updated_result.latency = n.dendrite_hw->default_latency_update.value();
    }

    if (!updated_result.energy.has_value())
    {
        std::string error(
                "Error: Dendrite unit does not simulate energy or provide "
                "a default energy cost in the architecture "
                "description.");
        throw std::logic_error(error);
    }
    if (!updated_result.latency.has_value())
    {
        std::string error("Error: Dendrite unit does not simulate latency or "
                          "provide a default latency cost in the architecture "
                          "description.");
        throw std::logic_error(error);
    }

    return updated_result;
}

sanafe::PipelineResult
sanafe::PipelineUnit::calculate_soma_default_energy_latency(
        MappedNeuron &n, const PipelineResult &simulation_result)
{
    PipelineResult updated_result{simulation_result};

    bool energy_simulated = simulation_result.energy.has_value();
    bool latency_simulated = simulation_result.latency.has_value();

    bool soma_energy_metrics_set =
            n.soma_hw->default_soma_energy_metrics.has_value();
    if (energy_simulated && soma_energy_metrics_set)
    {
        std::string error(
                "Error: Soma unit simulates energy and also has "
                "default energy metrics set. Remove the default energy metrics "
                "from the architecture description.");
        throw std::logic_error(error);
    }
    if (soma_energy_metrics_set)
    {
        updated_result.energy =
                n.soma_hw->default_soma_energy_metrics->energy_access_neuron;
    }

    bool soma_latency_metrics_set =
            n.soma_hw->default_soma_energy_metrics.has_value();
    if (latency_simulated && soma_latency_metrics_set)
    {
        std::string error(
                "Error: Soma unit simulates latency and also has "
                "default energy costs set. Remove the default latency metrics "
                "from the architecture description");
        throw std::logic_error(error);
    }
    if (soma_latency_metrics_set)
    {
        if (simulation_result.latency.has_value())
        {
            std::string error(
                    "Error: Soma unit simulates energy and also has "
                    "default energy costs set. Remove default energy costs from the "
                    "architecture description.");
            throw std::logic_error(error);
        }
        updated_result.latency =
                n.soma_hw->default_soma_latency_metrics->latency_access_neuron;
    }

    if ((simulation_result.status == sanafe::UPDATED) ||
            (simulation_result.status == sanafe::FIRED))
    {
        if (n.soma_hw->default_soma_energy_metrics.has_value())
        {
            updated_result.energy.value() +=
                    n.soma_hw->default_soma_energy_metrics->energy_update_neuron;
        }
        if (n.soma_hw->default_soma_latency_metrics.has_value())
        {
            updated_result.latency.value() +=
                    n.soma_hw->default_soma_latency_metrics
                            ->latency_update_neuron;
        }
    }
    if (simulation_result.status == sanafe::FIRED)
    {
        if (n.soma_hw->default_soma_energy_metrics.has_value())
        {
            updated_result.energy.value() +=
                    n.soma_hw->default_soma_energy_metrics->energy_spike_out;
        }
        if (n.soma_hw->default_soma_latency_metrics.has_value())
        {
            updated_result.latency.value() +=
                    n.soma_hw->default_soma_latency_metrics->latency_spike_out;
        }
    }

    if (!updated_result.energy.has_value())
    {
        std::string error("Error: Soma unit does not simulate energy or "
                          "provide default energy costs in the architecture "
                          "description.");
        throw std::logic_error(error);
    }
    if (!updated_result.latency.has_value())
    {
        std::string error("Error: Soma unit does not simulate latency or "
                          "provide default latency costs in the architecture "
                          "description.");
        throw std::logic_error(error);
    }

    return updated_result;
}

void sanafe::PipelineUnit::update_soma_activity(
        MappedNeuron &n, const PipelineResult &simulation_result)
{
    if ((simulation_result.status == sanafe::UPDATED) ||
            (simulation_result.status == sanafe::FIRED))
    {
        n.soma_hw->neuron_updates++;

        if (simulation_result.status == sanafe::FIRED)
        {
            n.soma_hw->neurons_fired++;
            n.axon_out_input_spike = true;
            TRACE1(CHIP, "Neuron %s.%zu fired\n", n.parent_group_name.c_str(),
                    n.id);
        }
    }
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

        ts.messages[n.core->id].push_back(m);
        ++axon_out_hw.packets_out;
        ++total_messages_sent;

        axon_result.energy.value() += axon_out_hw.energy_access;
        axon_result.latency.value() += axon_out_hw.latency_access;
    }
    n.axon_out_input_spike = false;

    return axon_result;
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

size_t sanafe::abs_diff(const size_t a, const size_t b)
{
    // Returns the absolute difference between two unsigned (size_t) values
    return (a > b) ? (a - b) : (b - a);
}

sanafe::RunData sanafe::SpikingChip::get_run_summary() const
{
    // Store the summary data in a string to string mapping
    RunData run_data(0, total_timesteps);

    run_data.spikes = total_spikes;
    run_data.packets_sent = total_messages_sent;
    run_data.wall_time = wall_time;
    run_data.neurons_fired = total_neurons_fired;

    run_data.total_energy = total_energy;
    run_data.synapse_energy = synapse_energy;
    run_data.dendrite_energy = dendrite_energy;
    run_data.soma_energy = soma_energy;
    run_data.network_energy = network_energy;

    run_data.sim_time = total_sim_time;

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

void sanafe::PipelineUnit::configure(
        std::string unit_name, const ModelInfo &model)
{
    model_parameters = model.model_parameters;
    plugin_lib = model.plugin_library_path;
    name = unit_name;
    log_energy = model.log_energy;
    log_latency = model.log_latency;

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
        default_soma_energy_metrics = energy_metrics;
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
        default_soma_latency_metrics = latency_metrics;
    }

    // Finally, forward all parameters from the architecture description to the
    //  model. This might be useful if you want to define any additional
    //  model-specific attributes here, e.g., fault-rate or maximum memory size.
    for (auto &[key, param] : model_parameters)
    {
        set_attribute_hw(key, param);
    }
}

sanafe::AxonOutUnit::AxonOutUnit(const AxonOutConfiguration &config)
        : name(std::move(config.name))
        , energy_access(config.metrics.energy_message_out)
        , latency_access(config.metrics.latency_message_out)
{
}

void sanafe::Core::map_neuron(
        const Neuron &neuron_to_map, const size_t neuron_id)
{
    TRACE1(CHIP, "Mapping nid:%s.%zu to core: %zu\n",
            neuron_to_map.parent_group_id.c_str(), neuron_to_map.offset, id);

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
    if (pipeline_hw.empty())
    {
        INFO("Error: No pipeline units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No units defined");
    }
    PipelineUnit *mapped_dendrite;

    bool choose_first_dendrite_by_default =
            (neuron_to_map.dendrite_hw_name.length() == 0);
    bool dendrite_found = false;
    for (auto &hw : pipeline_hw)
    {
        if (hw->implements_dendrite &&
                (choose_first_dendrite_by_default ||
                        neuron_to_map.dendrite_hw_name == hw->name))
        {
            mapped_dendrite = hw.get();
            dendrite_found = true;
            break;
        }
    }
    if (!dendrite_found)
    {
        INFO("Error: Could not map neuron nid:%zu (hw:%s) "
             "to any dendrite h/w.\n",
                neuron_to_map.offset, neuron_to_map.dendrite_hw_name.c_str());
        throw std::runtime_error("Error: Could not map neuron to dendrite h/w");
    }

    PipelineUnit *mapped_soma;
    bool choose_first_soma_by_default =
            (neuron_to_map.soma_hw_name.length() == 0);
    bool soma_found = false;
    for (auto &hw : pipeline_hw)
    {
        if (hw->implements_soma &&
                (choose_first_soma_by_default ||
                        neuron_to_map.soma_hw_name == hw->name))
        {
            mapped_soma = hw.get();
            soma_found = true;
            break;
        }
    }
    if (!soma_found)
    {
        INFO("Error: Could not map neuron nid:%zu (hw:%s) "
             "to any soma h/w.\n",
                neuron_to_map.offset, neuron_to_map.soma_hw_name.c_str());
        throw std::runtime_error("Error: Could not map neuron to soma h/w");
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
    neurons.emplace_back(neuron_to_map, this, neuron_id, address,
            mapped_dendrite, mapped_soma, mapped_axon_out);

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

sanafe::PipelineUnit &sanafe::Core::create_pipeline_unit(
        const PipelineUnitConfiguration &config)
{
    // Create the synapse model
    if (config.model_info.plugin_library_path.has_value())
    {
        const std::filesystem::path plugin_lib_path =
                config.model_info.plugin_library_path.value();
        TRACE1(CHIP, "Creating unit from plugin: %s.\n",
                plugin_lib_path.c_str());
        pipeline_hw.emplace_back(
                plugin_get_hw(config.model_info.name, plugin_lib_path));
    }
    else
    {
        // Use built in models
        TRACE1(CHIP, "Creating built-in model %s.\n",
                config.model_info.name.c_str());
        pipeline_hw.emplace_back(
                model_get_pipeline_unit(config.model_info.name));
    }

    auto &new_unit = pipeline_hw.back();
    new_unit->implements_synapse = config.implements_synapse;
    new_unit->implements_dendrite = config.implements_dendrite;
    new_unit->implements_soma = config.implements_soma;
    TRACE1(CHIP, "implements synapse:%d dendrite:%d soma:%d\n",
            new_unit->implements_synapse, new_unit->implements_dendrite,
            new_unit->implements_soma);
    new_unit->configure(config.name, config.model_info);
    TRACE1(CHIP, "New h/w unit created (%s) in core:%zu\n", config.name.c_str(),
            id);

    return *new_unit;
}

sanafe::AxonOutUnit &sanafe::Core::create_axon_out(
        const AxonOutConfiguration &config)
{
    axon_out_hw.emplace_back(config);
    TRACE1(CHIP, "New axon out h/w unit created: (%zu.%zu)\n", parent_tile_id,
            id);

    return axon_out_hw.back();
}

sanafe::Message::Message(const long int id, const SpikingChip &hw,
        const MappedNeuron &n, const long int timestep)
        : timestep(timestep)
        , mid(id)
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

sanafe::Message::Message(const long int id, const SpikingChip &hw,
        const MappedNeuron &n, const long int timestep, const int axon_address)
        : Message(id, hw, n, timestep)
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
        , log_energy(config.power_metrics.log_energy)
        , log_latency(config.power_metrics.log_latency)
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
        , log_energy(config.pipeline.log_energy)
        , log_latency(config.pipeline.log_latency)

{
    timestep_buffer.resize(pipeline_config.max_neurons_supported);
}

sanafe::MappedConnection::MappedConnection(const int connection_id)
        : post_neuron(nullptr)
        , pre_neuron(nullptr)
        , synapse_hw(nullptr)
        , id(connection_id)
{
}

sanafe::MappedNeuron::MappedNeuron(const Neuron &neuron_to_map,
        Core *mapped_core, const size_t nid, const size_t mapped_address,
        PipelineUnit *mapped_dendrite, PipelineUnit *mapped_soma,
        AxonOutUnit *mapped_axon_out)
        : parent_group_name(neuron_to_map.parent_group_id)
        , offset(neuron_to_map.offset)
        , id(nid)
        , core(mapped_core)
        , dendrite_hw(mapped_dendrite)
        , soma_hw(mapped_soma)
        , axon_out_hw(mapped_axon_out)
        , mapped_address(mapped_address)
        , mapping_order(neuron_to_map.mapping_order)
        , force_synapse_update(neuron_to_map.force_synapse_update)
        , force_dendrite_update(neuron_to_map.force_dendrite_update)
        , force_soma_update(neuron_to_map.force_soma_update)
        , log_spikes(neuron_to_map.log_spikes)
        , log_potential(neuron_to_map.log_potential)

{
    set_model_attributes(neuron_to_map.model_parameters);
    build_neuron_processing_pipeline();
}

void sanafe::MappedNeuron::set_model_attributes(
        const std::map<std::string, sanafe::ModelParam> &model_parameters)
{
    for (auto &[key, param] : model_parameters)
    {
        TRACE2(CHIP, "Forwarding param: %s (dendrite:%d soma:%d)\n",
                key.c_str(), param.forward_to_dendrite, param.forward_to_soma);
        if (param.forward_to_dendrite && (dendrite_hw != nullptr))
        {
            dendrite_hw->set_attribute_neuron(mapped_address, key, param);
        }
        if (param.forward_to_soma && (soma_hw != nullptr))
        {
            soma_hw->set_attribute_neuron(mapped_address, key, param);
        }
    }
}

void sanafe::SpikingChip::sim_output_run_summary(
        const std::filesystem::path &output_dir) const
{
    const RunData run_data = get_run_summary();
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
    out << "wall_time: " << std::fixed << run_data.wall_time << "\n";
    out << "total_neurons_fired: " << run_data.neurons_fired << "\n";
    out << "sim_time: " << std::scientific << run_data.sim_time << "\n";
    // Give a more detailed energy breakdown, as this is often useful
    out << "energy:\n";
    out << "  synapse:" << std::scientific << run_data.synapse_energy << "\n";
    out << "  dendrite:" << std::scientific << run_data.dendrite_energy << "\n";
    out << "  soma:" << std::scientific << run_data.soma_energy << "\n";
    out << "  network: " << std::scientific << run_data.network_energy << "\n";
    out << "  total: " << std::scientific << run_data.total_energy << "\n";
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

void sanafe::SpikingChip::sim_timestep(
        Timestep &ts, const TimingModel timing_model)
{
    Scheduler scheduler;

    // Start the next time-step, clear all buffers
    assert(core_count > 0);
    ts = Timestep(ts.timestep, core_count);
    sim_reset_measurements();

    process_neurons(ts);
    process_messages(ts);
    forced_updates(ts);

    scheduler.noc_width = noc_width;
    scheduler.noc_height = noc_height;
    scheduler.buffer_size = noc_buffer_size;
    scheduler.core_count = core_count;
    scheduler.max_cores_per_tile = max_cores_per_tile;

    if (timing_model == TIMING_MODEL_SIMPLE)
    {
        TRACE1(CHIP, "Running simple timing model\n");
        ts.sim_time = schedule_messages_simple(ts.messages, scheduler);
    }
    else if (timing_model == TIMING_MODEL_DETAILED)
    {
        TRACE1(CHIP, "Running detailed timing model\n");
        ts.sim_time = schedule_messages(ts.messages, scheduler);
    }
    else if (timing_model == TIMING_MODEL_CYCLE_ACCURATE)
    {
        TRACE1(CHIP, "Running cycle-accurate timing model\n");
        if (chip_count > 1)
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
        ts.sim_time =
                schedule_messages_cycle_accurate(ts.messages, *booksim_config);
    }
    else
    {
        INFO("Error: Timing model:%d not recognized\n", timing_model);
        throw std::invalid_argument("Timing model not recognized");
    }

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
}

sanafe::Timestep::Timestep(const long int ts, const int core_count)
        : messages(std::vector<std::list<Message>>(core_count))
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
            c.messages_in = std::vector<Message *>();
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
    potential_trace_file << std::endl;
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

            for (const auto hw : c.pipeline_hw)
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
        perf_trace_file <<  "," << name;
    }
    perf_trace_file << "\n";
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
    message_trace_file << "blocking_latency";
    message_trace_file << std::endl;
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
                spike_trace_file << std::endl;
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
        potential_trace_file << std::endl;
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
}

void sanafe::SpikingChip::sim_trace_record_message(
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
    message_trace_file << m.blocked_delay;
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
