// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  pipeline.cpp
#include <cassert>
#include <cmath>

#include <omp.h>

#include "arch.hpp"
#include "models.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"
#include "sim.hpp"

void sanafe::pipeline_process_neurons(Timestep &ts, Architecture &arch)
{
    auto cores = arch.cores();

    // Older versions of OpenMP don't support range-based for loops yet...
    // TODO: Figure a way so that the timestep and Architecture structs aren't
    //  shared between threads (dangerous)
#pragma omp parallel for schedule(dynamic) default(none) shared(cores, ts, arch)
    // codechecker_suppress [modernize-loop-convert]
    for (size_t idx = 0; idx < cores.size(); idx++)
    {
        Core &core = cores[idx];
        for (Neuron *n : core.neurons)
        {
            pipeline_process_neuron(ts, arch, *n);
        }

        if (core.next_message_generation_delay != 0.0)
        {
            // This message accounts for any remaining neuron processing
            const Neuron &last_neuron = *(core.neurons.back());
            Message placeholder(arch, last_neuron, ts.timestep);
            placeholder.generation_delay = core.next_message_generation_delay;
            ts.messages[core.id].push_back(placeholder);
        }
    }
}

void sanafe::pipeline_process_messages(Timestep &ts, Architecture &arch)
{
    // Assign outgoing spike messages to their respective destination
    //  cores, and calculate network costs
    for (auto &q : ts.messages)
    {
        for (auto &m : q)
        {
            if (!m.placeholder)
            {
                pipeline_receive_message(arch, m);
            }
        }
    }

    // Now process all messages at receiving cores
    auto cores = arch.cores();
    // Older versions of OpenMP don't support range-based for loops yet...
#pragma omp parallel for schedule(dynamic) default(none) shared(cores, ts, arch)
    // codechecker_suppress [modernize-loop-convert]
    for (size_t idx = 0; idx < cores.size(); idx++)
    {
        Core &core = cores[idx];
        TRACE1("Processing %zu message(s) for cid:%zu\n",
                core.messages_in.size(), core.id);
        for (auto *m : core.messages_in)
        {
            m->receive_delay += pipeline_process_message(ts, core, *m);
        }
    }
}

void sanafe::pipeline_receive_message(Architecture &arch, Message &m)
{
    assert(static_cast<size_t>(m.src_tile_id) < arch.tiles.size());
    assert(static_cast<size_t>(m.dest_tile_id) < arch.tiles.size());
    const Tile &src_tile = arch.tiles[m.src_tile_id];
    Tile &dest_tile = arch.tiles[m.dest_tile_id];
    m.network_delay = sim_estimate_network_costs(src_tile, dest_tile);
    m.hops = abs_diff(src_tile.x, dest_tile.x) +
            abs_diff(src_tile.y, dest_tile.y);

    Core &core = dest_tile.cores[m.dest_core_offset];
    core.messages_in.push_back(&m);
}

void sanafe::pipeline_process_neuron(
        Timestep &ts, const Architecture &arch, Neuron &n)
{
    SIM_TRACE1("Processing neuron: %d.%d\n", n.id, n.parent_group_id);
    double neuron_processing_latency = 0.0;

    // Update any H/W following the time-step buffer in pipeline order
    if ((n.core->pipeline_config.buffer_position <=
                BUFFER_BEFORE_DENDRITE_UNIT))
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

    // Finally see if any H/W units have been forced to update every time-step.
    //  This may be needed if the H/W has some state e.g., a leak, that has
    //  to be constantly updated every time-step for enabled neurons
    if (n.dendrite_hw->force_update)
    {
        neuron_processing_latency += pipeline_process_dendrite(ts, n);
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
    SIM_TRACE1("Receiving message for cid:%d\n", c.id);
    double message_processing_latency = pipeline_process_axon_in(core, m);

    assert(static_cast<size_t>(m.dest_axon_id) < core.axons_in.size());
    const AxonInModel &axon_in = core.axons_in[m.dest_axon_id];
    for (const int synapse_address : axon_in.synapse_addresses)
    {
        Connection &con = *(core.connections_in[synapse_address]);
        message_processing_latency +=
                pipeline_process_synapse(ts, con, synapse_address);
        if (core.pipeline_config.buffer_position ==
                BUFFER_BEFORE_DENDRITE_UNIT)
        {
            continue; // Process next synapse
        }
        // In certain pipeline configurations, every synaptic lookup requires
        //  updates to the dendrite and/or soma units as well
        Neuron &n = *(con.post_neuron);
        message_processing_latency += pipeline_process_dendrite(ts, n);

        if (core.pipeline_config.buffer_position ==
                BUFFER_BEFORE_SOMA_UNIT)
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
        const Timestep &ts, Connection &con, const int synapse_address)
{
    // Update all synapses to different neurons in one core. If a synaptic
    //  lookup, read and accumulate the synaptic weights. Otherwise, just
    //  update filtered current and any other connection properties
    SIM_TRACE1("Updating synapses for (cid:%d)\n", c.id);
    while (con.last_updated < ts.timestep)
    {
        SIM_TRACE1("Updating synaptic current (last updated:%d, ts:%ld)\n",
                con.last_updated, ts.timestep);
        con.synapse_model->update();
        con.last_updated++;
    }
    Synapse synapse_data = {0.0, con};
    synapse_data.current = con.synapse_model->update(synapse_address, false);

    // Input and buffer synaptic info at the next hardware unit (dendrite unit)
    con.post_neuron->dendrite_input_synapses.push_back(synapse_data);
    con.post_neuron->spike_count++;
    assert(con.synapse_hw != nullptr);
    con.synapse_hw->spikes_processed++;
    SIM_TRACE1("(nid:%s.%s->nid:%s.%s) con->current:%lf\n",
            con.pre_neuron->parent_group_id.c_str(), con.pre_neuron->id.c_str(),
            con.post_neuron->parent_group_id.c_str(),
            con.post_neuron->id.c_str(), con.current);

    return con.synapse_hw->latency_spike_op;
}

double sanafe::pipeline_process_dendrite(const Timestep &ts, Neuron &n)
{
    while (n.dendrite_last_updated < ts.timestep)
    {
        TRACE2("Updating nid:%s dendritic current "
               "(last_updated:%d, ts:%ld)\n",
                n.id.c_str(), n.dendrite_last_updated, ts.timestep);
        n.soma_input_charge = n.dendrite_model->update();
        n.dendrite_last_updated++;
    }
    for (const auto &synapse : n.dendrite_input_synapses)
    {
        n.soma_input_charge = n.dendrite_model->update(synapse, false);
    }
    n.dendrite_input_synapses.clear();

    // Finally, send dendritic current to the soma
    TRACE2("nid:%s updating dendrite, soma_input_charge:%lf\n", n.id.c_str(),
            n.soma_input_charge);

    return 0.0;
}

double sanafe::pipeline_process_soma(const Timestep &ts, Neuron &n)
{
    TRACE1("nid:%s updating, current_in:%lf\n", n.id.c_str(),
            n.soma_input_charge);
    double soma_processing_latency = 0.0;
    while (n.soma_last_updated < ts.timestep)
    {
        std::optional<double> soma_current_in;
        if ((n.spike_count > 0) || (std::fabs(n.soma_input_charge) > 0.0))
        {
            soma_current_in = n.soma_input_charge;
            n.soma_input_charge = 0.0;
        }
        n.status = n.soma_model->update(soma_current_in);
        if (n.forced_spikes > 0)
        {
            n.status = sanafe::FIRED;
            n.forced_spikes--;
        }

        soma_processing_latency += n.soma_hw->latency_access_neuron;
        if ((n.status == sanafe::UPDATED) || (n.status == sanafe::FIRED))
        {
            soma_processing_latency += n.soma_hw->latency_update_neuron;
            n.soma_hw->neuron_updates++;
        }
        if (n.status == sanafe::FIRED)
        {
            soma_processing_latency += n.soma_hw->latency_spiking;
            n.soma_hw->neurons_fired++;
            n.axon_out_input_spike = true;
            SIM_TRACE1("Neuron %d.%d fired\n", n.parent_group_id, n.id);
        }

        n.soma_last_updated++;
    }

    SIM_TRACE1("neuron status:%d\n", n.status);

    return soma_processing_latency;
}

double sanafe::pipeline_process_axon_out(
        Timestep &ts, const Architecture &arch, Neuron &n)
{
    if (!n.axon_out_input_spike)
    {
        return 0.0;
    }

    TRACE1("nid:%s.%s sending spike message to %zu axons out\n",
            n.parent_group_id.c_str(), n.id.c_str(),
            n.axon_out_addresses.size());
    for (const int axon_address : n.axon_out_addresses)
    {
        Message m(arch, n, ts.timestep, axon_address);
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
