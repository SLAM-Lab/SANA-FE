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
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < arch.cores_vec.size(); i++)
    {
        Core &c = arch.cores_vec[i];
        for (Neuron *n : c.neurons)
        {
            pipeline_process_neuron(ts, arch, *n);
        }

        if (c.neurons.size() > 0)
        {
            // Add a dummy message to account for neuron
            //  processing that does not result in any sent
            //  messages. To do this, set the dest neuron
            //  set as invalid with a 0 receiving latency)
            Message dummy_message = c.next_message;
            dummy_message.dummy_message = true;
            dummy_message.src_neuron = c.neurons.back();
            dummy_message.receive_delay = 0.0;
            dummy_message.network_delay = 0.0;
            dummy_message.hops = 0;
            dummy_message.timestep = ts.timestep;
            assert(c.id >= 0);
            ts.messages[c.id].push_back(dummy_message);
        }
    }
}

void sanafe::pipeline_receive_messages(Timestep &ts, Architecture &arch)
{
    // Assign outgoing spike messages to their respective destination
    //  cores, and calculate network costs
    for (auto &q : ts.messages)
    {
        for (auto &m : q)
        {
            if (!m.dummy_message)
            {
                const size_t src_tile_id = m.src_neuron->core->parent_tile_id;
                assert(src_tile_id < arch.tiles_vec.size());
                assert(m.dest_tile_id >= 0);
                assert(static_cast<size_t>(m.dest_tile_id) <
                        arch.tiles_vec.size());
                Tile &src_tile = arch.tiles_vec[src_tile_id];
                Tile &dest_tile = arch.tiles_vec[m.dest_tile_id];
                m.network_delay =
                        sim_estimate_network_costs(src_tile, dest_tile);
                m.hops = abs(src_tile.x - dest_tile.x) +
                        abs(src_tile.y - dest_tile.y);

                Core &core = dest_tile.cores_vec[m.dest_core_offset];
                core.messages_in.push_back(&m);
            }
        }
    }

    // Now process all messages at receiving cores
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < arch.cores_vec.size(); i++)
    {
        Core &core = arch.cores_vec[i];
        TRACE1("Processing %lu message(s) for cid:%d\n",
                core.messages_in.size(), core.id);
        for (auto m : core.messages_in)
        {
            m->receive_delay += pipeline_process_message(ts, arch, core, *m);
        }
    }
}

void sanafe::pipeline_process_neuron(
        Timestep &ts, Architecture &arch, Neuron &n)
{
    SIM_TRACE1("Processing neuron: %d.%d\n", n.id, n.parent_group_id);
    double neuron_processing_latency = 0.0;
    if (n.core->timestep_buffer_position == BUFFER_BEFORE_DENDRITE_UNIT)
    {
        neuron_processing_latency += pipeline_process_dendrite(ts, arch, n);
    }
    if ((n.core->timestep_buffer_position == BUFFER_BEFORE_DENDRITE_UNIT) ||
            (n.core->timestep_buffer_position == BUFFER_BEFORE_SOMA_UNIT))
    {
        neuron_processing_latency += pipeline_process_soma(ts, arch, n);
    }
    if ((n.core->timestep_buffer_position == BUFFER_BEFORE_DENDRITE_UNIT) ||
            (n.core->timestep_buffer_position == BUFFER_BEFORE_SOMA_UNIT) ||
            (n.core->timestep_buffer_position == BUFFER_BEFORE_AXON_OUT_UNIT))
    {
        neuron_processing_latency += pipeline_process_axon_out(ts, arch, n);
    }

    n.core->next_message.generation_delay += neuron_processing_latency;
    n.spike_count = 0;

    return;
}

double sanafe::pipeline_process_message(
        Timestep &ts, Architecture &arch, Core &core, Message &m)
{
    // Simulate message m in the message processing pipeline. The message is
    //  sequentially handled by units up to the time-step buffer
    SIM_TRACE1("Receiving message for cid:%d\n", c.id);
    double message_processing_latency = pipeline_process_axon_in(core, m);

    assert(static_cast<size_t>(m.dest_axon_id) < core.axons_in.size());
    AxonInModel &axon_in = core.axons_in[m.dest_axon_id];
    for (const int synapse_address : axon_in.synapse_addresses)
    {
        Connection &con = *(core.connections_in[synapse_address]);
        message_processing_latency +=
                pipeline_process_synapse(ts, arch, con, synapse_address);
        if (core.timestep_buffer_position == BUFFER_BEFORE_DENDRITE_UNIT)
        {
            continue; // Process next synapse
        }
        // In certain pipeline configurations, every synaptic lookup requires
        //  updates to the dendrite and/or soma units as well
        Neuron &n = *(con.post_neuron);
        message_processing_latency += pipeline_process_dendrite(ts, arch, n);

        if (core.timestep_buffer_position == BUFFER_BEFORE_SOMA_UNIT)
        {
            continue; // Process next synapse
        }
        message_processing_latency += pipeline_process_soma(ts, arch, n);
        assert(core.timestep_buffer_position == BUFFER_BEFORE_AXON_OUT_UNIT);
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

double sanafe::pipeline_process_synapse(Timestep &ts, Architecture &arch,
        Connection &con, const int synapse_address)
{
    // Update all synapses to different neurons in one core. If a synaptic
    //  lookup, read and accumulate the synaptic weights. Otherwise, just
    //  update filtered current and any other connection properties
    SIM_TRACE1("Updating synapses for (cid:%d)\n", c.id);
    while (con.last_updated < ts.timestep)
    {
        SIM_TRACE1("Updating synaptic current (last updated:%ld, ts:%ld)\n",
                con.last_updated, ts.timestep);
        con.synapse_model->update();
        con.last_updated++;
    }
    Synapse synapse_data = {0.0, con.dest_compartment};
    synapse_data.current = con.synapse_model->update(synapse_address, false);

    // Input and buffer synaptic info at the next hardware unit (dendrite unit)
    con.post_neuron->dendrite_input_synapses.push_back(synapse_data);
    con.post_neuron->spike_count++;
    assert(con.synapse_hw != nullptr);
    con.synapse_hw->spikes_processed++;
    SIM_TRACE1("(nid:%d.%d->nid:%d.%d) con->current:%lf\n",
            con.pre_neuron->parent_group_id, con.pre_neuron->id,
            con.post_neuron->parent_group_id, con.post_neuron->id, con.current);

    return con.synapse_hw->latency_spike_op;
}

double sanafe::pipeline_process_dendrite(
        Timestep &ts, Architecture &arch, Neuron &n)
{
    double latency;
    latency = 0.0;

    while (n.dendrite_last_updated < ts.timestep)
    {
        TRACE2("Updating nid:%d dendritic current "
               "(last_updated:%d, ts:%ld)\n",
                n.id, n.dendrite_last_updated, ts.timestep);
        n.soma_input_charge = n.dendrite_model->update();
        n.dendrite_last_updated++;
    }
    for (const auto &synapse : n.dendrite_input_synapses)
    {
        n.soma_input_charge = n.dendrite_model->update(
                synapse.current, synapse.dest_compartment, false);
    }
    n.dendrite_input_synapses.clear();

    // Finally, send dendritic current to the soma
    TRACE2("nid:%d updating dendrite, soma_input_charge:%lf\n", n.id,
            n.soma_input_charge);

    return latency;
}

double sanafe::pipeline_process_soma(
        Timestep &ts, Architecture &arch, Neuron &n)
{
    SomaUnit *const soma = n.soma_hw;

    SIM_TRACE1("nid:%d updating, current_in:%lf\n", n.id, n.soma_input_charge);
    double soma_processing_latency = 0.0;
    while (n.soma_last_updated < ts.timestep)
    {
        std::optional<double> soma_current_in = std::nullopt;
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
            soma->neuron_updates++;

            if (n.status == sanafe::FIRED)
            {
                soma_processing_latency += n.soma_hw->latency_spiking;
                soma->neurons_fired++;
                n.axon_out_input_fired = true;
                SIM_TRACE1("Neuron %d.%d fired\n", n.parent_group_id, n.id);
            }
        }

        n.soma_last_updated++;
    }

    SIM_TRACE1("neuron status:%d\n", n.status);

    return soma_processing_latency;
}

double sanafe::pipeline_process_axon_out(
        Timestep &ts, Architecture &arch, Neuron &n)
{
    if (!n.axon_out_input_fired)
    {
        return 0.0;
    }

    TRACE1("nid:%d.%d sending spike message to %lu axons out\n",
            n.parent_group_id, n.id, n.axon_out_addresses.size());
    for (int address : n.axon_out_addresses)
    {
        Core &src_core = *(n.core);
        const Tile &src_tile = arch.tiles_vec[src_core.parent_tile_id];
        AxonOutModel &src_axon = src_core.axons_out[address];

        // Generate a spike message
        Tile &dest_tile = arch.tiles_vec[src_axon.dest_tile_id];
        Core &dest_core = dest_tile.cores_vec[src_axon.dest_core_offset];
        AxonInModel &dest_axon = dest_core.axons_in[src_axon.dest_axon_id];

        Message m;
        m.timestep = ts.timestep;
        m.src_neuron = &n;
        m.dest_axon_id = src_axon.dest_axon_id;
        m.spikes = dest_axon.synapse_addresses.size();
        m.dummy_message = false;
        m.dest_tile_id = dest_tile.id;
        m.dest_core_id = dest_core.id;
        m.dest_core_offset = dest_core.offset;
        m.src_x = src_tile.x;
        m.dest_x = dest_tile.x;
        m.src_y = src_tile.y;
        m.dest_y = dest_tile.y;
        // TODO: support multiple axon output units, included in the synapse
        m.dest_axon_hw = 0;

        // Add axon access cost to message latency and energy
        AxonOutUnit &axon_out_hw = *(n.axon_out_hw);
        m.generation_delay += src_core.next_message.generation_delay +
                axon_out_hw.latency_access;
        m.network_delay = 0.0;
        m.receive_delay = 0.0;
        axon_out_hw.packets_out++;

        ts.messages[src_core.id].push_back(m);
        // Reset the next message in this core
        src_core.next_message = Message();
    }
    n.axon_out_input_fired = false;

    return n.axon_out_hw->latency_access;
}
