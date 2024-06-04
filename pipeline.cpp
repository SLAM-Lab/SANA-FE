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
        for (std::vector<Neuron>::size_type k = 0; k < c.neurons.size(); k++)
        {
            Neuron &n = *(c.neurons[k]);
            pipeline_process_neuron(ts, arch, n);
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

// TODO: replace ts with messages
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
        Core &c = arch.cores_vec[i];
        TRACE1("Processing %lu message(s) for cid:%d\n", c.messages_in.size(),
                c.id);
        for (auto m : c.messages_in)
        {
            m->receive_delay += pipeline_process_message(ts, arch, c, *m);
        }
    }
}

void sanafe::pipeline_process_neuron(Timestep &ts, Architecture &arch,
        Neuron &n)
{
    Core &c = *(n.core);
    n.processing_latency = 0.0;

    SIM_TRACE1("Processing neuron: %d.%d\n", n.id, n.parent_group_id);

    if (c.buffer_pos == BUFFER_BEFORE_SYNAPSE)
    {
        INFO("Error: Not implemented\n");
        throw std::invalid_argument("Error: Buffer pos not implemented\n");
    }
    else if (c.buffer_pos == BUFFER_BEFORE_DENDRITE)
    {
        n.processing_latency = pipeline_update_dendrite(ts, arch, n);
    }
    else if (c.buffer_pos == BUFFER_BEFORE_SOMA)
    {
        n.processing_latency = pipeline_update_soma(ts, arch, n);
    }
    else if (c.buffer_pos == BUFFER_BEFORE_AXON_OUT)
    {
        if (n.neuron_status == sanafe::FIRED)
        {
            n.processing_latency = n.soma_hw->latency_spiking;
            pipeline_neuron_send_spike_messages(ts, arch, n);
        }
    }
    SIM_TRACE1("Updating neuron %d.%d.\n", n.parent_group_id, n.id);

    c.next_message.generation_delay += n.processing_latency;
    n.spike_count = 0;
}

double sanafe::pipeline_process_message(
        Timestep &ts, Architecture &arch, Core &c, Message &m)
{
    // We receive a spike and process up to the time-step buffer
    double message_processing_latency = 0.0;

    SIM_TRACE1("Receiving messages for cid:%d\n", c.id);
    if (c.buffer_pos >= BUFFER_BEFORE_SYNAPSE)
    {
        assert(m.dest_axon_id >= 0);
        assert(static_cast<size_t>(m.dest_axon_id) < c.axons_in.size());
        AxonInModel &a = c.axons_in[m.dest_axon_id];

        assert(m.dest_axon_hw >= 0);
        assert(static_cast<size_t>(m.dest_axon_hw) < c.axon_in_hw.size());

        AxonInUnit &hw = c.axon_in_hw[m.dest_axon_hw];
        message_processing_latency += hw.latency_spike_message;
        hw.spike_messages_in++;
        for (int s : a.synapse_addresses)
        {
            message_processing_latency += pipeline_update_synapse(ts, arch, c, s);
        }
    }

    return message_processing_latency;
}

void sanafe::pipeline_neuron_send_spike_messages(
        Timestep &ts, Architecture &arch, Neuron &n)
{
    SIM_TRACE1("nid:%d.%d sending spike message to %lu axons out\n",
            n.parent_group_id, n.id, n.axon_out_addresses.size());
    for (int address : n.axon_out_addresses)
    {
        Core &src_core = *(n.core);
        const Tile &src_tile = arch.tiles_vec[src_core.parent_tile_id];
        AxonOutModel &src_axon = src_core.axons_out[address];
        const int dest_address = src_axon.dest_axon_id;

        // Generate a spike message
        Tile &dest_tile = arch.tiles_vec[src_axon.dest_tile_id];
        Core &dest_core = dest_tile.cores_vec[src_axon.dest_core_offset];
        AxonInModel &dest_axon = dest_core.axons_in[dest_address];

        // TODO: figure some constructor for the message i.e., required
        //  fields?
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

        // TODO: support multiple axon output units
        //  I would need to figure how we map a synapse to a specific
        //  axon out
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

    return;
}

double sanafe::pipeline_update_synapse(
        Timestep &ts, Architecture &arch, Core &c, const int synapse_address)
{
    // Update all synapses to different neurons in one core. If a synaptic
    //  lookup, read and accumulate the synaptic weights. Otherwise, just
    //  update filtered current and any other connection properties
    double latency;
    latency = 0.0;
    Connection &con = *(c.connections_in[synapse_address]);
    Neuron &post_neuron = *(con.post_neuron);

    assert(con.synapse_hw != NULL);
    SIM_TRACE1("Updating synapses for (cid:%d)\n", c.id);
    Synapse synapse_data;
    synapse_data.dest_compartment = con.dest_compartment;
    while (con.last_updated < ts.timestep)
    {
        SIM_TRACE1("Updating synaptic current (last updated:%ld, ts:%ld)\n",
                con.last_updated, ts.timestep);
        con.last_updated++;
        synapse_data.current = con.synapse_model->update();
    }
    // Store the synapse in a buffer belonging to the post-synaptic neuron
    synapse_data.current = con.synapse_model->update(synapse_address, false);

    post_neuron.dendrite_input_synapses.push_back(synapse_data);
    post_neuron.spike_count++;
    con.synapse_hw->spikes_processed++;
    latency += con.synapse_hw->latency_spike_op;
    TRACE2("Sending spike to nid:%d, current:%lf\n", post_neuron->id,
            syn.current);

    SIM_TRACE1("(nid:%d.%d->nid:%d.%d) con->current:%lf\n",
            con.pre_neuron->parent_group_id, con.pre_neuron->id,
            con.post_neuron->parent_group_id, con.post_neuron->id, con.current);

    if (c.buffer_pos != BUFFER_BEFORE_DENDRITE)
    {
        latency += pipeline_update_dendrite(ts, arch, post_neuron);
    }

    return latency;
}

double sanafe::pipeline_update_dendrite(Timestep &ts, Architecture &arch, Neuron &n)
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
    if (n.core->buffer_pos != BUFFER_BEFORE_SOMA)
    {
        latency += pipeline_update_soma(ts, arch, n);
    }

    return latency;
}

double sanafe::pipeline_update_soma(Timestep &ts, Architecture &arch, Neuron &n)
{
    SomaUnit *const soma = n.soma_hw;

    SIM_TRACE1("nid:%d updating, current_in:%lf\n", n.id, n.soma_input_charge);
    double latency = 0.0;
    while (n.soma_last_updated < ts.timestep)
    {
        std::optional<double> soma_current_in = std::nullopt;
        if ((n.spike_count > 0) || (std::fabs(n.soma_input_charge) > 0.0))
        {
            soma_current_in = n.soma_input_charge;
            n.soma_input_charge = 0.0;
        }
        n.neuron_status = n.soma_model->update(soma_current_in);
        if (n.forced_spikes > 0)
        {
            n.neuron_status = sanafe::FIRED;
            n.forced_spikes--;
        }

        latency += n.soma_hw->latency_access_neuron;
        if ((n.neuron_status == sanafe::UPDATED) ||
                (n.neuron_status == sanafe::FIRED))
        {
            latency += n.soma_hw->latency_update_neuron;
            soma->neuron_updates++;
            if (n.neuron_status == sanafe::FIRED)
            {
                SIM_TRACE1("Neuron %d.%d fired\n", n.parent_group_id, n.id);
                pipeline_neuron_send_spike_messages(ts, arch, n);
                latency += n.soma_hw->latency_spiking;
                soma->neurons_fired++;
            }
        }

        n.soma_last_updated++;
    }

    SIM_TRACE1("neuron status:%d\n", n.neuron_status);

    return latency;
}
