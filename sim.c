// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  sim.c

// TODO:
// j* ust have a single set of message lists for each core and write to them
//  in the sim routine. Instead of storing a message in each axon and copying
//  the data
// * figure out about this timestep struct. maybe the timestep struct goes instead
//  of the whole simulation struct? The simulation struct is for the higher
//  level routines? When we simulate a timestep, we should only return details
//  about the timestep. The bigger simulator loop should accumulate results.

// * account for final dummy message which adds all remaining neuron processing
//    time. Need to support 1 extra message per core in theory
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "print.h"
#include "sim.h"
#include "network.h"
#include "arch.h"

void sim_init_fifo(struct message_fifo *f)
{
	f->count = 0;
	f->head = NULL;
	f->tail = NULL;
	f->next = NULL;
}

void arch_init_message_scheduler(struct message_scheduler *const scheduler,
	struct architecture *const arch)
{
	sim_init_fifo(&(scheduler->network_fifo));

	scheduler->tile_count = arch->tile_count;
	scheduler->core_count = 0;
	// Count the total number of cores in the design
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		scheduler->core_count += t->core_count;
		scheduler->spikes_buffered[i] = 0;
	}
	for (int i = 0; i < scheduler->core_count; i++)
	{
		sim_init_fifo(&(scheduler->send_fifo[i]));
		sim_init_fifo(&(scheduler->receive_fifo[i]));
	}

	return;
}

struct timestep sim_timestep(struct simulation *const sim,
	struct network *const net, struct architecture *const arch,
	struct message_scheduler *const scheduler)
{
	int core_count;

	// Start the next time-step
	sim_init_timestep(&(sim->ts));
	sim_reset_measurements(net, arch);

	sim_process_neurons(sim, net, arch);
	sim_input_spikes(net);
	sim_receive_messages(sim, arch);

	core_count = 0;

	// TODO: it might make sense to just write to this buffer in the main
	//  loop? What sharing will actually happen? Will there be any
	//  serious false sharing? there shouldn't be
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			for (struct message *m =
				sim_message_fifo_pop(&(c->messages_sent));
				m != NULL;
				m = sim_message_fifo_pop(&(c->messages_sent)))
			{
				// Copy all message into global buffers from the
				//  core structs
				int message_count =
					sim->ts.messages_sent_count[core_count];
				sim->ts.messages[core_count][message_count] =*m;
				m = &(sim->ts.messages[core_count][message_count]);
				(sim->ts.messages_sent_count[core_count])++;
				sim_message_fifo_push(
					&(scheduler->send_fifo[core_count]), m);
			}
			core_count++;
		}
	}

	sim->ts.sim_time = sim_schedule_messages(sim, scheduler);
	// Performance statistics for this time step
	sim->ts.energy = sim_calculate_energy(arch);

	core_count = 0;
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);

			for (int k = 0; k < ARCH_MAX_UNITS; k++)
			{
				sim->ts.spike_count +=
					c->synapse[k].spikes_processed;
				sim->ts.total_neurons_fired +=
					c->soma[k].neurons_fired;
			}
			sim->total_messages_sent +=
				sim->ts.messages_sent_count[core_count];
			core_count++;
		}
	}
	if (sim->perf_fp)
	{
		sim_perf_log_timestep(&(sim->ts), sim->perf_fp);
	}

	sim->timesteps++;
	sim->total_energy += sim->ts.energy;
	sim->total_sim_time += sim->ts.sim_time;
	sim->total_spikes += sim->ts.spike_count;
	sim->total_neurons_fired += sim->ts.total_neurons_fired;

	TRACE1("Spikes sent: %ld\n", sim->total_spikes);
	return sim->ts;
}

void sim_init_sim(struct simulation *sim)
{
	sim->total_energy = 0.0;   // Joules
	sim->total_sim_time = 0.0; // Seconds
	sim->wall_time = 0.0;	   // Seconds
	sim->timesteps = 0;
	sim->total_spikes = 0;
	sim->total_messages_sent = 0;
	sim->total_neurons_fired = 0;

	// All logging disabled by default
	sim->log_perf = 0;
	sim->log_potential = 0;
	sim->log_spikes = 0;
	sim->log_messages = 0;

	sim->potential_trace_fp = NULL;
	sim->spike_trace_fp = NULL;
	sim->perf_fp = NULL;
	sim->message_trace_fp = NULL;
	sim->stats_fp = NULL;

	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		sim->ts.messages_sent_count[i] = 0;
	}
}

void sim_init_timestep(struct timestep *const ts)
{
	ts->timestep = -1;
	ts->spike_count = 0;
	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		ts->messages_sent_count[i] = 0;
	}
	ts->total_neurons_fired = 0;
	ts->spikes = 0;
	ts->total_hops = 0;
	ts->energy = 0.0;
	ts->sim_time = 0.0;
}

void sim_process_neurons(struct simulation *const sim, struct network *net,
	struct architecture *arch)
{
#pragma omp parallel for
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			for (int k = 0; k < c->neuron_count; k++)
			{
				struct neuron *n = c->neurons[k];
				assert(n != NULL);
				sim_process_neuron(sim, n);
			}

			// Add a dummy message to account for neuron processing
			// TODO: make this cleaner, when we just write all
			//  messages into the main buffer...
			struct message *m = &(c->neuron_processing_latency);
			sim_message_fifo_push(&(c->messages_sent), m);
		}
	}
	if (sim->log_spikes)
	{
		sim_trace_record_spikes(sim, net);
	}
	if (sim->log_potential)
	{
		sim_trace_record_potentials(sim, net);
	}
}

void sim_receive_messages(struct simulation *const sim,
	struct architecture *arch)
{
#pragma omp parallel for
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			for (int k = 0; k < c->axon_in.map_count; k++)
			{
				struct connection_map *axon =
					&(c->axon_in.map[k]);
				if (axon->spikes_received > 0)
				{
					struct neuron *pre_neuron =
						axon->pre_neuron;
					assert(pre_neuron != NULL);
					struct core *pre_core =
						pre_neuron->core;
					assert(pre_core != NULL);
					struct tile *pre_tile = pre_core->t;
					assert(pre_tile != NULL);
					axon->message.network_latency =
						sim_estimate_network_costs(
							pre_tile, t);
					axon->message.receive_latency =
						sim_pipeline_receive(
							sim, c, axon);
				}
			}
		}
	}
}

double sim_estimate_network_costs(struct tile *const src,
	struct tile *const dest)
{
	double network_latency;
	long int x_hops, y_hops;

	assert(src != NULL);
	assert(dest != NULL);

	network_latency = src->time_spike_within_tile;

	// Calculate the energy and time for sending spike packets
	x_hops = abs(src->x - dest->x);
	y_hops = abs(src->y - dest->y);
	// E-W hops

	if (src->x < dest->x)
	{
		dest->east_hops += x_hops;
		network_latency += (double) x_hops * src->latency_east_hop;
	}
	else
	{
		dest->west_hops += x_hops;
		network_latency += (double) x_hops * src->latency_west_hop;
	}

	// N-S hops
	if (src->y < dest->y)
	{
		dest->north_hops += y_hops;
		network_latency += (double) y_hops * src->latency_north_hop;
	}
	else
	{
		dest->south_hops += y_hops;
		network_latency += (double) y_hops * src->latency_south_hop;
	}

	dest->hops += (x_hops + y_hops);
	dest->messages_received++;
	TRACE1("xhops:%ld yhops%ld total hops:%ld latency:%e\n", x_hops, y_hops,
		t->hops, network_latency);
	return network_latency;
}

struct message *sim_get_next_message(struct core *c, const int timestep)
{
	/*
	struct message *m;
	double neuron_processing_latency;

	// Get the next message
	if (c->messages_left == 0)
	{
		struct neuron *src_neuron;

		neuron_processing_latency = 0.0;
		while (c->neurons_left > 0)
		{
			assert(c->curr_neuron < c->neuron_count);
			src_neuron = c->neurons[c->curr_neuron];
			assert(src_neuron != NULL);

			neuron_processing_latency +=
				src_neuron->processing_latency;
			if (!(src_neuron->is_init) || !(src_neuron->fired) ||
				(src_neuron->maps_out_count == 0))
			{
				// Neuron isn't initialized, didn't fire or
				//  has no outgoing connections, so skip it.
				//  This neuron can't generate any messages.
				c->curr_neuron++;
				//INFO("\t curr_neuron:%d neuron_processing_latency:%e\n", c->curr_neuron, neuron_processing_latency);
				c->neurons_left--;
				c->curr_axon = 0;
			}
			else if (src_neuron->fired)
			{
				TRACE2("(cid:%d.%d) nid:%d.%d fired, "
					"time:%e\n",
					c->t->id, c->id,
					src_neuron->group->id,
					src_neuron->id, c->time);
				c->messages_left = src_neuron->maps_out_count;
				c->curr_axon = 0;
				break;
			}
		}
	}
	else
	{
		neuron_processing_latency = 0.0;
	}

	if (c->messages_left > 0)
	{
		m = &(c->curr_message);
		arch_init_message(m);

		m->timestep = timestep;
		m->src_neuron = c->neurons[c->curr_neuron];
		m->dest_axon = m->src_neuron->maps_out[c->curr_axon];
		m->spikes = m->dest_axon->connection_count;
		m->dest_neuron = m->dest_axon->connections[0]->post_neuron;

		// Add axon access cost to message latency and energy
		m->generation_latency = neuron_processing_latency +
			c->axon_out.time_access;
		m->network_latency = m->dest_axon->network_latency;
		m->receive_latency = m->dest_axon->receive_latency;
		c->axon_out.packets_out++;

		assert(m->src_neuron->maps_out != NULL);
		assert(m->src_neuron->maps_out[0] != NULL);
		assert(m->dest_axon != NULL);
		assert(m->dest_axon->connections != NULL);
		assert(m->dest_axon->spikes_received > 0);

		c->curr_axon++;
		c->messages_left--;
		if (c->messages_left == 0)
		{
			// If we processed the last message for this neuron,
			//  go to the next neuron
			c->curr_neuron++;
			c->neurons_left--;
			c->curr_axon = 0;
		}
	}
	else
	{
		assert(c->neurons_left == 0);
		// No messages left to process, but track any additional core
		//  processing time after the last message. Assume we've
		//  accounted for all mapped neurons
		m = NULL;
		c->latency_after_last_message = neuron_processing_latency;
	}
	*/
	return NULL;
}

struct message *sim_message_fifo_pop(struct message_fifo *queue)
{
	struct message *m;

	assert(queue->count >= 0);
	if (queue->count == 0)
	{
		m = NULL;
	}
	else
	{
		assert(queue->tail != NULL);
		assert(queue->head != NULL);
		queue->count--;
		m = queue->tail;
		queue->tail = queue->tail->next;
		if (queue->count <= 0)
		{
			queue->head = NULL;
		}
	}

	return m;
}

void sim_message_fifo_push(struct message_fifo *queue, struct message *m)
{
	assert(queue->count >= 0);
	assert(m != NULL);

	m->next = NULL;
	if (queue->count == 0)
	{
		queue->tail = m;
	}
	else
	{
		queue->head->next = m;
	}
	queue->head = m;
	queue->count++;
}

double sim_schedule_messages(const struct simulation *const sim,
	struct message_scheduler *scheduler)
{
	struct message_fifo *priority_queue;
	priority_queue = sim_init_timing_priority(scheduler);
	double last_timestamp, t;

	last_timestamp = 0.0;
	// Setup timing counters
	TRACE1("Scheduling global order of messages.\n");

	// While queue isn't empty
	while (priority_queue != NULL)
	{
		// Get the core with the earliest simulation time
		struct message_fifo *q =
			sim_pop_priority_queue(&priority_queue);
		struct message *m = sim_message_fifo_pop(q);
		last_timestamp = fmax(last_timestamp, m->generation_latency);

		if (m->dest_neuron != NULL)
		{
			if (m->dest_neuron->core->t->is_blocking)
			{
				m->blocked_latency = fmax(m->blocked_latency,
					m->dest_neuron->core->t->blocked_until -
					m->sent_timestamp);
				// Update the core global time, blocking until
				//  the receiving tile is free
				m->sent_timestamp = fmax(m->sent_timestamp,
					m->dest_neuron->core->t->blocked_until);
			}
			if (m->dest_neuron->core->is_blocking)
			{
				// Track how long the message is blocked for
				m->blocked_latency = fmax(m->blocked_latency,
					m->dest_neuron->core->blocked_until -
					m->sent_timestamp);
				// Update the core global time, blocking until
				//  the receiving core is free
				m->sent_timestamp = fmax(m->sent_timestamp,
					m->dest_neuron->core->blocked_until);

				if (m->sent_timestamp <
					m->dest_neuron->core->blocked_until)
				{
					// If we were trying to send a spike to
					//  a blocked core, also block the tile
					//  for this duration as well
					m->dest_neuron->core->t->blocked_until =
					m->dest_neuron->core->blocked_until;
					// Update the core global time, blocking
					//  until the receiving core is free
					m->sent_timestamp =
					m->dest_neuron->core->blocked_until;
				}
			}

			// Set time-stamps, calculating when the receiving H/W will be
			//  busy until
			m->sent_timestamp += m->network_latency;
			last_timestamp =
				fmax(last_timestamp, m->sent_timestamp);
			// TODO: for some reason this seems quite important for DVS
			//  gesture accuracy. The core is busy until the message is
			//  delivered by the network
			m->dest_neuron->core->blocked_until = fmax(
				(m->dest_neuron->core->blocked_until +
				m->network_latency + m->receive_latency),
				(m->sent_timestamp + m->receive_latency));
			m->processed_timestamp =
				m->dest_neuron->core->blocked_until;
			last_timestamp =
				fmax(last_timestamp, m->processed_timestamp);

			TRACE2("\t(cid:%d.%d) synapse at %d.%d busy until %e\n",
				c->t->id, c->id, m->dest_core->t->id,
				m->dest_core->id, m->dest_core->blocked_until);
		}

		if (sim->log_messages)
		{
			sim_trace_record_message(sim, m);
		}

		// The time that the last message sent is the time that we
		//  start the next message's processing
		t = m->sent_timestamp;
		// Get the next message, neuron or core
		m = q->tail;
		// Regardless of whether we are sending a message, add the
		//  processing time
		if (m != NULL)
		{
			m->sent_timestamp = t + m->generation_latency;
			last_timestamp =
				fmax(last_timestamp, m->sent_timestamp);
			sim_insert_priority_queue(&priority_queue, q);
		}
		else
		{
			TRACE2("\t(cid:%d.%d) finished simulating\n", c->t->id,
				c->id);
		}

		if (priority_queue != NULL)
		{
			TRACE2("\t(cid:%d.%d) time:%e\n",
				(*priority_queue)->t->id,
				(*priority_queue)->id, (*priority_queue)->time);
		}
	}
	TRACE1("Neurons fired: %ld\n", ts->total_neurons_fired);

	return last_timestamp;
}

/*
int sim_schedule_messages_old(const struct simulation *const sim,
	struct timestep *const ts, struct architecture *const arch)
{
	struct core *priority_queue;
	priority_queue = sim_init_timing_priority(arch, ts);

	// Setup timing counters
	TRACE1("Scheduling global order of messages.\n");

	// While queue isn't empty
	while (priority_queue != NULL)
	{
		// Get the core with the earliest simulation time
		struct core *c = sim_pop_priority_queue(&priority_queue);
		struct message *m = &(c->curr_message);

		if (m->dest_neuron->core->t->is_blocking)
		{
			m->blocked_latency = fmax(m->blocked_latency,
				m->dest_neuron->core->t->blocked_until -
				c->time);
			// Update the core global time, blocking until
			//  the receiving tile is free
			c->time = fmax(c->time,
				m->dest_neuron->core->t->blocked_until);
		}
		if (m->dest_neuron->core->is_blocking)
		{
			// Track how long the message is blocked for
			m->blocked_latency = fmax(m->blocked_latency,
				m->dest_neuron->core->blocked_until - c->time);
			// Update the core global time, blocking until
			//  the receiving core is free
			c->time = fmax(c->time,
				m->dest_neuron->core->blocked_until);

			// Calculate when the receiving h/w will be busy
			//  until
			if (c->time < m->dest_neuron->core->blocked_until)
			{
				// If we were trying to send a spike to
				//  a blocked core, also block the tile
				//  for this duration as well
				m->dest_neuron->core->t->blocked_until =
					m->dest_neuron->core->blocked_until;
				// Update the core global time, blocking
				//  until the receiving core is free
				c->time = m->dest_neuron->core->blocked_until;
			}
		}

		// Set time-stamps, calculating when the receiving H/W will be
		//  busy until
		c->time += m->network_latency;
		// TODO: for some reason this seems quite important for DVS
		//  gesture accuracy. The core is busy until the message is
		//  delivered by the network

		m->dest_neuron->core->blocked_until = fmax(
			(m->dest_neuron->core->blocked_until +
			m->network_latency + m->receive_latency),
			(c->time + m->receive_latency));

		TRACE2("\t(cid:%d.%d) synapse at %d.%d busy until %e\n",
			c->t->id, c->id, m->dest_core->t->id,
			m->dest_core->id, m->dest_core->blocked_until);

		TRACE2("cid:%d axon:%d (%d messages left)\n", c->id,
			c->curr_axon, c->messages_left);
		if (sim->log_messages)
		{
			sim_trace_record_message(sim, m);
		}

		// Get the next message, neuron or core
		m = sim_get_next_message(c, ts->timestep);
		// Regardless of whether we are sending a message, add the
		//  processing time
		if (m != NULL)
		{
			c->time += m->generation_latency;
			sim_insert_priority_queue(&priority_queue, c);
		}
		else
		{
			TRACE2("\t(cid:%d.%d) finished simulating\n", c->t->id,
				c->id);
		}

		if (priority_queue != NULL)
		{
			TRACE2("\t(cid:%d.%d) time:%e\n",
				(*priority_queue)->t->id,
				(*priority_queue)->id, (*priority_queue)->time);
		}
	}
	TRACE1("Neurons fired: %ld\n", ts->total_neurons_fired);

	// Add the processing time of neurons after sending the last message
	//INFO("** Final times **\n");
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			//INFO("c->time:%e, c->latency_after_last_message:%e\n", c->time, c->latency_after_last_message);
			c->time += c->latency_after_last_message;
		}
	}

	return 0;
}
*/

void sim_process_neuron(struct simulation *const sim, struct neuron *n)
{
	if (!n->is_init)
	{
		return;
	}

	struct core *c = n->core;
	n->processing_latency = 0.0;

	// TODO: I think there's a flaw here. We update synapses for the core,
	//  not for the neuron. Then each neuron does have a single dendritic
	//  tree. So it's just this, maybe we make it "process neurons"
	if (c->buffer_pos == BUFFER_SYNAPSE)
	{
		for (int i = 0; i < n->maps_in_count; i++)
		{
			struct connection_map *axon = &(n->maps_in[i]);
			n->processing_latency +=
				sim_update_synapse(sim, axon, 1);
		}
	}
	else if (c->buffer_pos == BUFFER_DENDRITE)
	{
		// Go through all synapses connected to this neuron and update
		//  all synaptic currents into the dendrite
		for (int i = 0; i < n->maps_in_count; i++)
		{
			struct connection_map *a = &(n->maps_in[i]);
			for (int j = 0; j < a->connection_count; j++)
			{
				struct connection *con = a->connections[j];
				n->processing_latency += sim_update_dendrite(
					sim, n, con->current);
			}
		}
	}
	else if (c->buffer_pos == BUFFER_SOMA)
	{
		n->processing_latency = sim_update_soma(sim, n, n->charge);
	}
	else if (c->buffer_pos == BUFFER_AXON_OUT)
	{
		if (n->fired)
		{
			struct soma_processor *soma = n->soma_hw;
			n->processing_latency = soma->latency_spiking;
			sim_neuron_send_spike_message(sim, n);
		}
	}
	TRACE1("Updating neuron %d.%d.\n", n->group->id, n->id);

	c->neuron_processing_latency.generation_latency +=
		n->processing_latency;
	n->update_needed = 0;
	n->spike_count = 0;
}

double sim_pipeline_receive(
	struct simulation *const sim, struct core *c,
	struct connection_map *axon)
{
	// We receive a spike and process up to the time-step buffer
	double message_processing_latency = 0.0;

	TRACE1("Receiving messages for cid:%d\n", c->id);
	if (c->buffer_pos >= BUFFER_SYNAPSE)
	{
		int synaptic_lookup = 1;
		message_processing_latency =
			sim_update_synapse(sim, axon, synaptic_lookup);
	}

	return message_processing_latency;
}

struct message_fifo *sim_init_timing_priority(
	struct message_scheduler *scheduler)
{
	struct message_fifo *priority_queue;
	priority_queue = NULL;

	TRACE1("Initializing priority queue.\n");
	for (int i = 0; i < scheduler->core_count; i++)
	{
		struct message *m = scheduler->send_fifo[i].tail;
		if (m != NULL) // messages
		{
			m->sent_timestamp = m->generation_latency;
			sim_insert_priority_queue(&priority_queue,
				&(scheduler->send_fifo[i]));
		}
		else
		{
			TRACE1("No messages for core %d\n", i);
		}
	}

#ifdef DEBUG2
	int i = 0;
	for (struct core *curr = priority_queue; curr != NULL; curr = curr->next)
	{
		// TODO
	}
#endif

	return priority_queue;
}

struct message_fifo *sim_pop_priority_queue(
	struct message_fifo **priority_queue)
{
	struct message_fifo *curr;

	// Pop the first element from the priority queue
	curr = *priority_queue;
	*priority_queue = (*priority_queue)->next;

	// For safety, remove current element from queue and unlink
	curr->next = NULL;
	return curr;
}

void sim_insert_priority_queue(struct message_fifo **priority_queue,
	struct message_fifo *core_message_fifo)
{
	struct message_fifo *next;
	int list_depth = 0;

	// TODO: implement heap-based priority queue rather than list-based.
	//  Will achieve O(lg N) insertion time rather than O(N)

	assert(priority_queue != NULL);
	assert(core_message_fifo != NULL);

	//INFO("Inserting into priority queue.\n");
	if (core_message_fifo->tail == NULL)
	{
		INFO("error?\n");
	}
	if (((*priority_queue) == NULL) ||
		(core_message_fifo->tail->sent_timestamp <=
		(*priority_queue)->tail->sent_timestamp))
	{
		// Queue is empty or this is the earliest time (highest
		//  priority), make this core the head of the queue
		core_message_fifo->next = (*priority_queue);
		*priority_queue = core_message_fifo;
	}
	else
	{
		struct message_fifo *curr = *priority_queue;
		next = curr->next;

		// Reinsert core into the correct place in the priority list
		while (next != NULL)
		{
			if (core_message_fifo->tail->sent_timestamp <
				next->tail->sent_timestamp)
			{
				break;
			}
			curr = next;
			next = curr->next;
			list_depth++;
		}
		curr->next = core_message_fifo;
		core_message_fifo->next = next;
	}

#ifdef DEBUG
	TRACE3("*** Priority queue ***\n");
	for (struct core *tmp = *priority_queue; tmp != NULL;
		tmp = tmp->next_timing)
	{
		// TRACE3
		TRACE3("tmp->time:%e (id:%d)\n", tmp->time, tmp->id);
		assert((tmp->next_timing == NULL) ||
			tmp->time <= tmp->next_timing->time);
	}
	// TRACE2
	TRACE3("List depth = %d\n", list_depth+1);
#endif

	return;
}

int sim_input_spikes(struct network *net)
{
	int input_spike_count = 0;

	// Seed all externally input spikes in the network for this timestep
	for (int i = 0; i < net->external_input_count; i++)
	{
		struct input *in = &(net->external_inputs[i]);

		if (in == NULL)
		{
			break;
		}

		if (in->type == INPUT_EVENT)
		{
			in->send_spike = (in->spike_val > 0.0);
		}
		else if (in->type == INPUT_POISSON)
		{
			in->send_spike = sim_poisson_input(in->rate);
		}
		else // INPUT_RATE)
		{
			in->send_spike =
				sim_rate_input(in->rate, &(in->spike_val));
		}

		if (!in->send_spike)
		{
			TRACE3("Not sending spike\n");
			continue;
		}
		for (int j = 0; j < in->post_connection_count; j++)
		{
			// Send a spike to all neurons connected to this input
			//  Normally, we would have a number of input dimensions
			//  for a given network
			struct neuron *post_neuron;
			struct connection *con;

			con = &(in->connections[j]);
			assert(con);

			post_neuron = con->post_neuron;
			TRACE3("nid:%d Energy before: %lf\n", post_neuron->id,
				post_neuron->current);
			if (post_neuron->core->buffer_pos == BUFFER_SOMA)
			{
				post_neuron->charge += con->weight;
			}
			else
			{
				post_neuron->current += con->weight;
			}
			TRACE3("nid:%d Energy after: %lf\n", post_neuron->id,
				post_neuron->current);

			con->synapse_hw->time += con->synapse_hw->latency_spike_op;

			post_neuron->update_needed = 1;
			post_neuron->spike_count++;
			input_spike_count++;
		}
		TRACE1("Sent spikes to %d connections\n",
			in->post_connection_count);

		// If inputting sets of events, then reset the spike after
		//  it's processed i.e. inputs are one-shot. For poisson and
		//  rate inputs, their values stay unchanged until the user
		//  sets them again (e.g. until the next input is presented)
		if (in->type == INPUT_EVENT)
		{
			in->spike_val = 0.0;
		}
	}

	return input_spike_count;
}

double sim_update_synapse(struct simulation *const sim,
	struct connection_map *axon, const int synaptic_lookup)
{
	// Update all synapses to different neurons in one core. If a synaptic
	//  lookup, read and accumulate the synaptic weights. Otherwise, just
	//  update filtered current and any other connection properties
	struct core *post_core;
	double latency, min_synaptic_resolution;
	//int spike_ops;

	latency = 0.0;
	post_core = axon->connections[0]->post_neuron->core;
	//spike_ops = axon->connection_count;

	TRACE1("Updating synapses for (cid:%d)\n", axon->pre_neuron->id);
	while (axon->last_updated <= sim->timesteps)
	{
		TRACE1("Updating synaptic current (last updated:%ld, ts:%ld)\n",
			axon->last_updated, sim->timesteps);
		if (axon->active_synapses > 0)
		{
			for (int i = 0; i < axon->connection_count; i++)
			{
				struct neuron *post_neuron;
				struct connection *con = axon->connections[i];

				post_neuron = con->post_neuron;
				con->current *= con->synaptic_current_decay;

				// "Turn off" synapses that have basically no
				//  synaptic current left to decay (based on
				//  the weight resolution)
				min_synaptic_resolution = (1.0 /
					con->synapse_hw->weight_bits);
				if (fabs(con->current) <
					min_synaptic_resolution)
				{
					con->current = 0.0;
					axon->active_synapses--;
				}

				TRACE2("(nid:%d->nid:%d) con->current:%lf\n",
					con->pre_neuron->id,
					con->post_neuron->id, con->current);
				if (post_core->buffer_pos != BUFFER_DENDRITE)
				{
					latency += sim_update_dendrite(
						sim, post_neuron, con->current);
				}
			}
		}
		axon->last_updated++;
	}

	if (synaptic_lookup)
	{
		if (axon->connection_count > 0)
		{
			latency += post_core->axon_in.latency_spike_message;
			post_core->axon_in.spike_messages_in++;
		}
		axon->active_synapses = axon->connection_count;

		for (int i = 0; i < axon->connection_count; i++)
		{
			struct neuron *post_neuron;
			struct connection *con = axon->connections[i];

			con->current += con->weight;
			post_neuron = con->post_neuron;
			post_neuron->update_needed = 1;
			post_neuron->spike_count++;

			assert(con->synapse_hw != NULL);
			con->synapse_hw->spikes_processed++;
			TRACE2("Sending spike to nid:%d, current:%lf\n",
				post_neuron->id, con->current);
			latency += con->synapse_hw->latency_spike_op;

			if (post_core->buffer_pos != BUFFER_DENDRITE)
			{
				latency += sim_update_dendrite(
					sim, post_neuron, con->current);
			}
		}
	}

	return latency;
}

double sim_update_dendrite(
	struct simulation *const sim, struct neuron *n, const double charge)
{
	// TODO: Support dendritic operations, combining the current in
	//  different neurons in some way, and writing the result to an output
	double dendritic_current, latency;
	latency = 0.0;

	dendritic_current = 0.0;
	while (n->dendrite_last_updated <= sim->timesteps)
	{
		TRACE3("Updating dendritic current (last_updated:%d, ts:%ld)\n",
			n->dendrite_last_updated, sim->timesteps);
		n->charge *= n->dendritic_current_decay;
		n->dendrite_last_updated++;
		dendritic_current = n->charge;
		TRACE2("nid:%d charge:%lf\n", n->id, n->charge);
	}

	// Update dendritic tap currents
	// TODO: implement multi-tap models
	TRACE2("Charge:%lf\n", charge);
	dendritic_current += charge;
	n->charge += charge;

	// Finally, send dendritic current to the soma
	TRACE2("nid:%d updating dendrite, charge:%lf\n", n->id, n->charge);
	if (n->core->buffer_pos != BUFFER_SOMA)
	{
		latency += sim_update_soma(sim, n, dendritic_current);
	}

	return latency;
}

double sim_update_soma(
	struct simulation *const sim, struct neuron *n, const double current_in)
{
	double latency = 0.0;
	struct soma_processor *soma = n->soma_hw;

	TRACE1("nid:%d updating, current_in:%lf\n", n->id, current_in);
	if ((soma->model == NEURON_LIF) ||
		(soma->model == NEURON_STOCHASTIC_LIF))
	{
		latency += sim_update_soma_lif(sim, n, current_in);
	}
	else if (soma->model == NEURON_TRUENORTH)
	{
		latency += sim_update_soma_truenorth(sim, n, current_in);
	}
	else
	{
		INFO("Neuron model not recognised: %d", soma->model);
		assert(0);
	}

	return latency;
}

double sim_generate_noise(struct neuron *n)
{
	assert(n != NULL);
	struct soma_processor *soma_hw = n->soma_hw;
	int noise_val = 0;
	int ret;

	if (soma_hw->noise_type == NOISE_FILE_STREAM)
	{
		// With a noise stream, we have a file containing a series of
		//  random values. This is useful if we want to exactly
		//  replicate h/w without knowing how the stream is generated.
		//  We can record the random sequence and replicate it here
		char noise_str[MAX_NOISE_FILE_ENTRY];
		// If we get to the end of the stream, by default reset it.
		//  However, it is unlikely the stream will be correct at this
		//  point
		if (feof(soma_hw->noise_stream))
		{
			INFO("Warning: At the end of the noise stream. "
			     "Random values are unlikely to be correct.\n");
			fseek(soma_hw->noise_stream, 0, SEEK_SET);
		}
		fgets(noise_str, MAX_NOISE_FILE_ENTRY, soma_hw->noise_stream);
		ret = sscanf(noise_str, "%d", &noise_val);
		TRACE2("noise val:%d\n", noise_val);
		if (ret < 1)
		{
			INFO("Error: invalid noise stream entry.\n");
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

double sim_update_soma_lif(
	struct simulation *const sim, struct neuron *n, const double current_in)

{
	struct soma_processor *soma = n->soma_hw;
	double random_potential, latency = 0.0;

	// Calculate the change in potential since the last update e.g.
	//  integate inputs and apply any potential leak
	TRACE1("Updating potential, before:%f\n", n->potential);

	if ((soma->model == NEURON_LIF) ||
		(soma->model == NEURON_STOCHASTIC_LIF))
	{
		while (n->soma_last_updated <= sim->timesteps)
		{
			n->potential *= n->leak_decay;
			n->soma_last_updated++;
		}
	}

	// Add randomized noise to potential if enabled
	if ((soma->model == NEURON_STOCHASTIC_LIF) &&
		(soma->noise_type == NOISE_FILE_STREAM))
	{
		random_potential = sim_generate_noise(n);
		n->potential += random_potential;
	}

	// Add the synaptic / dendrite current to the potential
	//printf("n->bias:%lf n->potential before:%lf current_in:%lf\n", n->bias, n->potential, current_in);
	n->potential += current_in + n->bias;
	n->charge = 0.0;
	//printf("n->bias:%lf n->potential after:%lf\n", n->bias, n->potential);

	TRACE1("Updating potential, after:%f\n", n->potential);

	// Check against threshold potential (for spiking)
	if (((n->bias != 0.0) && (n->potential > n->threshold)) ||
		((n->bias == 0.0) && (n->potential >= n->threshold)))
	{
		struct soma_processor *soma = n->soma_hw;

		if (n->group->reset_mode == NEURON_RESET_HARD)
		{
			n->potential = n->reset;
		}
		else if (n->group->reset_mode == NEURON_RESET_SOFT)
		{
			n->potential -= n->threshold;
		}
		n->fired = 1;
		soma->neurons_fired++;
		latency += soma->latency_spiking;
		if (n->core->buffer_pos != BUFFER_AXON_OUT)
		{
			sim_neuron_send_spike_message(sim, n);
		}
	}

	// Check against reverse threshold
	if (n->potential < n->reverse_threshold)
	{
		if (n->group->reverse_reset_mode == NEURON_RESET_SOFT)
		{
			n->potential -= n->reverse_threshold;
		}
		else if (n->group->reverse_reset_mode == NEURON_RESET_HARD)
		{
			n->potential = n->reverse_reset;
		}
		else if (n->group->reverse_reset_mode == NEURON_RESET_SATURATE)
		{
			n->potential = n->reverse_threshold;
		}
	}

	// Update soma, if there are any received spikes, there is a non-zero
	//  bias or we force the neuron to update every time-step
	if ((fabs(n->potential) > 0.0) || n->spike_count ||
		(fabs(n->bias) > 0.0) || n->force_update)
	{
		latency += n->soma_hw->latency_update_neuron;
		soma->neuron_updates++;
	}

	latency += n->soma_hw->latency_access_neuron;

	return latency;
}

double sim_update_soma_truenorth(
	struct simulation *const sim, struct neuron *n, const double current_in)
{
	struct soma_processor *soma = n->soma_hw;
	double v, latency = 0.0;
	int randomize_threshold;

	// Apply leak
	while (n->soma_last_updated <= sim->timesteps)
	{
		// Linear leak
		if (soma->leak_towards_zero)
		{
			// TODO: what happens if we're above zero but by less
			//  than the leak amount (for convergent), will we
			//  oscillate between the two? Does it matter
			if (n->potential > 0.0)
			{
				n->potential -= n->leak_bias;
			}
			else if (n->potential < 0.0)
			{
				n->potential += n->leak_bias;
			}
			// else equals zero, so no leak is applied
		}
		else
		{
			n->potential += n->leak_decay;
		}
		n->soma_last_updated++;
	}

	// Add the synaptic currents, processed by the dendrite
	n->potential += current_in + n->bias;
	n->current = 0.0;
	n->charge = 0.0;

	// Apply thresholding and reset
	v = n->potential;
	randomize_threshold = (n->random_range_mask != 0);
	if (randomize_threshold)
	{
		unsigned int r = rand() & n->random_range_mask;
		v += (double) r;
	}

	TRACE2("v:%lf +vth:%lf mode:%d -vth:%lf mode:%d\n", v, n->threshold,
		n->group->reset_mode, n->reverse_threshold,
		n->group->reverse_reset_mode);
	if (v >= n->threshold)
	{
		struct soma_processor *soma = n->soma_hw;
		int reset_mode = n->group->reset_mode;

		if (reset_mode == NEURON_RESET_HARD)
		{
			n->potential = n->reset;
		}
		else if (reset_mode == NEURON_RESET_SOFT)
		{
			n->potential -= n->threshold;
		}
		else if (reset_mode == NEURON_RESET_SATURATE)
		{
			n->potential = n->threshold;
		}
		n->fired = 1;
		soma->neurons_fired++;
		latency += soma->latency_spiking;
		if (n->core->buffer_pos != BUFFER_AXON_OUT)
		{
			sim_neuron_send_spike_message(sim, n);
		}
	}
	else if (v <= n->reverse_threshold)
	{
		int reset_mode = n->group->reverse_reset_mode;
		if (reset_mode == NEURON_RESET_HARD)
		{
			n->potential = n->reverse_reset;
		}
		else if (reset_mode == NEURON_RESET_SOFT)
		{
			n->potential += n->reverse_threshold;
		}
		else if (reset_mode == NEURON_RESET_SATURATE)
		{
			n->potential = n->reverse_threshold;
		}
		// No spike is generated
	}
	TRACE2("potential:%lf threshold %lf\n", n->potential, n->threshold);

	return latency;
}

void sim_neuron_send_spike_message(const struct simulation *const sim,
	struct neuron *n)
{
	struct core *c = n->core;
	TRACE1("nid:%d sending spike(s).\n", n->id);

	for (int k = 0; k < n->maps_out_count; k++)
	{
		struct connection_map *dest_axon = n->maps_out[k];
		// TODO: allocate the message from some external message buffer
		struct message *m = &(dest_axon->message);

		// Generate a spike message
		arch_init_message(m);
		// TODO: just copy the message over
		m->timestep = sim->timesteps;
		m->src_neuron = n;
		m->spikes = dest_axon->connection_count;
		m->dest_neuron = dest_axon->connections[0]->post_neuron;
		// Add axon access cost to message latency and energy
		m->generation_latency =
			c->neuron_processing_latency.generation_latency +
			c->axon_out.latency_access;
		arch_init_message(&(c->neuron_processing_latency));
		c->neuron_processing_latency.generation_latency = 0.0;
		// Mark these values as invalid, they need to be set
		//  later
		m->network_latency = NAN;
		m->receive_latency = NAN;
		assert(m->src_neuron->maps_out != NULL);
		assert(m->src_neuron->maps_out[0] != NULL);

		c->axon_out.packets_out++;
		sim_message_fifo_push(&(c->messages_sent), m);

		// Record a spike message at all the connected cores (axons)
		dest_axon->spikes_received++;
	}

	return;
}

double sim_calculate_energy(const struct architecture *const arch)
{
	// Returns the total energy across the design, for this timestep
	double network_energy, synapse_energy, soma_energy, axon_out_energy;
	double axon_in_energy, total_energy;

	network_energy = 0.0;
	axon_in_energy = 0.0;
	synapse_energy = 0.0;
	soma_energy = 0.0;
	axon_out_energy = 0.0;

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		double total_hop_energy = t->east_hops * t->energy_east_hop;

		total_hop_energy += t->west_hops * t->energy_west_hop;
		total_hop_energy += t->south_hops * t->energy_south_hop;
		total_hop_energy +=  t->north_hops * t->energy_north_hop;
		network_energy += total_hop_energy;

		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);

			axon_in_energy += c->axon_in.spike_messages_in *
				c->axon_in.energy_spike_message;
			for (int k = 0; k < ARCH_MAX_UNITS; k++)
			{
				synapse_energy +=
					c->synapse[k].spikes_processed *
					c->synapse[k].energy_spike_op;
				soma_energy += c->soma[k].neuron_count *
					c->soma[k].energy_access_neuron;
				soma_energy += c->soma[k].neuron_updates *
					c->soma[k].energy_update_neuron;
				soma_energy += c->soma[k].neurons_fired *
					c->soma[k].energy_spiking;
			}

			axon_out_energy += c->axon_out.packets_out *
				c->axon_out.energy_access;
		}
	}

	total_energy = axon_in_energy + synapse_energy + soma_energy +
		axon_out_energy + network_energy;

	return total_energy;
}

void sim_reset_measurements(struct network *net, struct architecture *arch)
{
	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);
			// Neurons can be manually forced to update, for example
			//  if they have a constant input bias
			n->update_needed |=
				(n->force_update || (n->bias != 0.0));
			n->processing_latency = 0.0;
			n->fired = 0;

			for (int k = 0; k < n->maps_out_count; k++)
			{
				struct connection_map *a = n->maps_out[k];
				a->spikes_received = 0;
				// Just to be safe, zero out and clear the
				//  message struct
				arch_init_message(&(a->message));
			}
		}
	}

	// Reset any energy, time latency or other measurements of network
	//  hardware
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		// Reset tile
		t->energy = 0.0;
		t->blocked_until = 0.0;

		t->hops = 0;
		t->east_hops = 0;
		t->west_hops = 0;
		t->south_hops = 0;
		t->north_hops = 0;
		t->messages_received = 0;
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			// Reset core
			c->energy = 0.0;
			c->blocked_until = 0.0;
			arch_init_message(&(c->neuron_processing_latency));
			c->neuron_processing_latency.generation_latency = 0.0;

			c->axon_in.energy = 0.0;
			c->axon_in.time = 0;
			c->dendrite.energy = 0.0;
			c->dendrite.time = 0.0;

			for (int k = 0; k < c->synapse_count; k++)
			{
				c->synapse[k].energy = 0.0;
				c->synapse[k].time = 0.0;
				c->synapse[k].spikes_processed = 0;
			}

			for (int k = 0; k < c->soma_count; k++)
			{
				c->soma[k].energy = 0.0;
				c->soma[k].time = 0.0;
				c->soma[k].neuron_updates = 0;
				c->soma[k].neurons_fired = 0;
			}

			c->axon_out.energy = 0.0;
			c->axon_out.time = 0.0;
			c->axon_out.packets_out = 0;
		}
	}
}

void sim_perf_write_header(FILE *fp)
{
	fprintf(fp, "time,");
	fprintf(fp, "fired,");
	fprintf(fp, "packets,");
	fprintf(fp, "hops,");
	fprintf(fp, "total_energy,");
	fprintf(fp, "\n");
}

void sim_perf_log_timestep(const struct timestep *const ts, FILE *fp)
{
	fprintf(fp, "%le,", ts->sim_time);
	fprintf(fp, "%ld,", ts->total_neurons_fired);
	// TODO: have total sum again
	//fprintf(fp, "%ld,", ts->messages_sent_count);
	fprintf(fp, "%ld,", ts->total_hops);
	fprintf(fp, "%le,", ts->energy);
	fprintf(fp, "\n");
}

void sim_write_summary(FILE *fp, const struct simulation *sim)
{
	// Write the simulation summary to file
	fprintf(fp, "git_version: %s\n", GIT_COMMIT);
	fprintf(fp, "energy: %e\n", sim->total_energy);
	fprintf(fp, "time: %e\n", sim->total_sim_time);
	fprintf(fp, "total_spikes: %ld\n", sim->total_spikes);
	fprintf(fp, "total_packets: %ld\n", sim->total_messages_sent);
	fprintf(fp, "total_neurons_fired: %ld\n", sim->total_neurons_fired);
	fprintf(fp, "wall_time: %lf\n", sim->wall_time);
	fprintf(fp, "timesteps: %ld\n", sim->timesteps);
}

void sim_spike_trace_write_header(const struct simulation *const sim)
{
	assert(sim->spike_trace_fp != NULL);
	fprintf(sim->spike_trace_fp, "gid.nid,timestep\n");

	return;
}

void sim_potential_trace_write_header(
	const struct simulation *const sim, const struct network *net)
{
	// Write csv header for probe outputs - record which neurons have been
	//  probed
	for (int i = 0; i < net->external_input_count; i++)
	{
		const struct input *in = &(net->external_inputs[i]);

		if (sim->potential_trace_fp && sim->log_potential)
		{
			fprintf(sim->potential_trace_fp, "i.%d,", in->id);
		}
	}

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		const struct neuron_group *group = &(net->groups[i]);
		for (int j = 0; j < group->neuron_count; j++)
		{
			const struct neuron *n = &(group->neurons[j]);
			if (sim->potential_trace_fp && sim->log_potential &&
				n->log_potential)
			{
				fprintf(sim->potential_trace_fp, "%d.%d,",
					group->id, n->id);
			}
		}
	}

	if (sim->potential_trace_fp)
	{
		fputc('\n', sim->potential_trace_fp);
	}

	return;
}

void sim_message_trace_write_header(const struct simulation *const sim)
{
	assert(sim->message_trace_fp);
	fprintf(sim->message_trace_fp, "timestep,src_neuron,");
	fprintf(sim->message_trace_fp, "src_hw,dest_hw,hops,spikes,");
	fprintf(sim->message_trace_fp, "generation_latency,network_latency,");
	fprintf(sim->message_trace_fp, "processing_latency,blocking_latency\n");
}

void sim_trace_record_spikes(
	const struct simulation *const sim, const struct network *net)
{
	// A trace of all spikes that are generated
	int spike_probe_count = 0;
	assert(sim->spike_trace_fp != NULL);

	for (int i = 0; i < net->external_input_count; i++)
	{
		const struct input *in = &(net->external_inputs[i]);
		if (in->send_spike)
		{
			fprintf(sim->spike_trace_fp, "i.%d,%d,", in->id,
				in->send_spike);
		}
	}

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		const struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			const struct neuron *n = &(group->neurons[j]);
			if (n->log_spikes && n->fired)
			{
				fprintf(sim->spike_trace_fp, "%d.%d,%ld\n",
					n->group->id, n->id, sim->timesteps);
				spike_probe_count++;
			}
		}
	}

	return;
}

void sim_trace_record_potentials(
	const struct simulation *const sim, const struct network *net)
{
	// Each line of this csv file is the potential of all probed neurons for
	//  one time-step
	int potential_probe_count = 0;

	for (int i = 0; i < net->external_input_count; i++)
	{
		const struct input *in = &(net->external_inputs[i]);
		if (sim->potential_trace_fp)
		{
			fprintf(sim->potential_trace_fp, "%lf,", in->spike_val);
		}
	}

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		const struct neuron_group *group = &(net->groups[i]);
		for (int j = 0; j < group->neuron_count; j++)
		{
			const struct neuron *n = &(group->neurons[j]);
			if (sim->potential_trace_fp && n->log_potential)
			{
				fprintf(sim->potential_trace_fp, "%lf,",
					n->potential);
				potential_probe_count++;
			}
		}
	}

	// Each timestep takes up a line in the respective csv file
	if (sim->potential_trace_fp && (potential_probe_count > 0))
	{
		fputc('\n', sim->potential_trace_fp);
	}

	return;
}

void sim_trace_record_message(
	const struct simulation *const sim, const struct message *const m)
{
	fprintf(sim->message_trace_fp, "%ld,", m->timestep);
	assert(m->src_neuron != NULL);
	fprintf(sim->message_trace_fp, "%d.%d,", m->src_neuron->group->id,
		m->src_neuron->id);
	assert(m->src_neuron->core != NULL);
	assert(m->src_neuron->core->t != NULL);
	fprintf(sim->message_trace_fp, "%d.%d,", m->src_neuron->core->t->id,
		m->src_neuron->core->id);

	assert(m->dest_neuron->core != NULL);
	assert(m->dest_neuron->core->t != NULL);
	fprintf(sim->message_trace_fp, "%d.%d,", m->dest_neuron->core->t->id,
		m->dest_neuron->core->id);
	fprintf(sim->message_trace_fp, "%d,", m->hops);
	fprintf(sim->message_trace_fp, "%d,", m->spikes);
	fprintf(sim->message_trace_fp, "%le,", m->generation_latency);
	fprintf(sim->message_trace_fp, "%le,", m->network_latency);
	fprintf(sim->message_trace_fp, "%le,", m->receive_latency);
	fprintf(sim->message_trace_fp, "%le\n", m->blocked_latency);
	fprintf(sim->message_trace_fp, "%le\n", m->sent_timestamp);
	fprintf(sim->message_trace_fp, "%le\n", m->processed_timestamp);

	return;
}

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
