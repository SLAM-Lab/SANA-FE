// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  sim.c

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "print.h"
#include "sim.h"
#include "network.h"
#include "arch.h"

#define MAX_MESSAGES_PER_HOP 12

void sim_init_fifo(struct message_fifo *f)
{
	f->count = 0;
	f->head = NULL;
	f->tail = NULL;
	f->next = NULL;
}

void sim_timestep(struct timestep *const ts,
	struct network *const net, struct architecture *const arch)
{
	// Start the next time-step
	sim_init_timestep(ts);
	sim_reset_measurements(net, arch);

	// TODO: reimplement user input spikes again
	//sim_input_spikes(net);
	sim_process_neurons(ts, net, arch);
	sim_receive_messages(ts, arch);

	ts->sim_time = sim_schedule_messages(ts->message_queues,
		arch->noc_width, arch->noc_height);
	//ts->sim_time = sim_schedule_messages_simplest(ts->message_queues);
	//ts->sim_time = sim_schedule_messages_blocking(ts->message_queues);
	// Performance statistics for this time step
	ts->energy = sim_calculate_energy(arch);

	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);

			for (int k = 0; k < c->synapse_count; k++)
			{
				ts->spike_count +=
					c->synapse[k].spikes_processed;
			}
			for (int k = 0; k < c->soma_count; k++)
			{
				ts->total_neurons_fired +=
					c->soma[k].neurons_fired;
			}
			ts->packets_sent += c->axon_out.packets_out;
		}
		ts->total_hops += t->hops;

	}

	TRACE1("Spikes sent: %ld\n", sim->total_spikes);
	return;
}

struct simulation *sim_init_sim(void)
{
	struct simulation *sim =
		(struct simulation *) malloc(sizeof(struct simulation));
	if (sim == NULL)
	{
		INFO("Error: Couldn't create sim struct.\n");
		exit(1);
	}

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
		sim_init_fifo(&(sim->ts.message_queues[i]));
	}

	return sim;
}

void sim_init_timestep(struct timestep *const ts)
{
	ts->spike_count = 0L;
	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		sim_init_fifo(&(ts->message_queues[i]));
	}
	ts->total_neurons_fired = 0L;
	ts->spikes = 0L;
	ts->total_hops = 0L;
	ts->energy = 0.0;
	ts->sim_time = 0.0;
	ts->packets_sent = 0L;
}

void sim_process_neurons(struct timestep *const ts, struct network *net,
	struct architecture *arch)
{
#pragma omp parallel for
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			struct message *dummy_message;
			int message_count;

			for (int k = 0; k < c->neuron_count; k++)
			{
				struct neuron *n = c->neurons[k];
				assert(n != NULL);
				sim_process_neuron(ts, n);
			}

			if (c->neuron_count > 0)
			{
				message_count = ts->message_queues[c->id].count;
				// Add a dummy message to account for neuron processing
				//  that does not result in any sent messages. To do
				//  this, set the dest neuron set as invalid with a 0
				//  receiving latency)
				dummy_message = &(ts->messages[c->id][message_count]);
				*dummy_message = c->next_message;
				dummy_message->src_neuron =
						c->neurons[c->neuron_count-1];

				dummy_message->dest_neuron = NULL;
				dummy_message->receive_delay = 0.0;
				dummy_message->network_delay = 0.0;
				dummy_message->hops = 0;
				dummy_message->timestep = ts->timestep;

				sim_message_fifo_push(&(ts->message_queues[c->id]),
					dummy_message);
			}
		}
	}
}

void sim_receive_messages(struct timestep *const ts,
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
					axon->message->network_delay =
						sim_estimate_network_costs(
							pre_tile, t);
					axon->message->hops =
						abs(pre_tile->x - t->x) +
						abs(pre_tile->y - t->y);
					axon->message->receive_delay =
						sim_pipeline_receive(
							ts, c, axon);
				}
			}
		}
	}
}

double sim_estimate_network_costs(struct tile *const src,
	struct tile *const dest)
{
	double network_delay;
	long int x_hops, y_hops;

	assert(src != NULL);
	assert(dest != NULL);

	network_delay = 0.0;

	// Calculate the energy and time for sending spike packets
	x_hops = abs(src->x - dest->x);
	y_hops = abs(src->y - dest->y);
	// E-W hops

	if (src->x < dest->x)
	{
		dest->east_hops += x_hops;
		network_delay += (double) x_hops * src->latency_east_hop;
	}
	else
	{
		dest->west_hops += x_hops;
		network_delay += (double) x_hops * src->latency_west_hop;
	}

	// N-S hops
	if (src->y < dest->y)
	{
		dest->north_hops += y_hops;
		network_delay += (double) y_hops * src->latency_north_hop;
	}
	else
	{
		dest->south_hops += y_hops;
		network_delay += (double) y_hops * src->latency_south_hop;
	}

	dest->hops += (x_hops + y_hops);
	dest->messages_received++;
	TRACE1("xhops:%ld yhops%ld total hops:%ld latency:%e\n", x_hops, y_hops,
		t->hops, network_delay);
	return network_delay;
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
	assert(queue != NULL);
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

void sim_update_noc_message_counts(
	const struct message *m, struct noc_timings *noc, const int message_in)
{
	int src_x, src_y, dest_x, dest_y;

	assert(m != NULL);
	assert(m->src_neuron != NULL);
	assert(m->dest_neuron != NULL);
	src_x = m->src_neuron->core->t->x;
	src_y = m->src_neuron->core->t->y;
	dest_x = m->dest_neuron->core->t->x;
	dest_y = m->dest_neuron->core->t->y;

	// Go along x path, then y path (dimension order routing), and increment
	//  or decrement counter depending on the operation
	int x_increment, y_increment;
	//double adjust_original = (1.0 / (1.0+m->hops));
	double adjust = (1.0 / (2.0+m->hops));
	//double adjust = (1.0 / (1.0+m->hops));

	if (!message_in)
	{
		adjust *= -1;
		//adjust_original *= -1;
	}

	if (src_x < dest_x)
	{
		x_increment = 1;
	}
	else
	{
		x_increment = -1;
	}
	if (src_y < dest_y)
	{
		y_increment = 1;
	}
	else
	{
		y_increment = -1;
	}
	int prev_direction = 4 + (m->src_neuron->core->id % ARCH_MAX_CORES_PER_TILE);
	for (int x = src_x; x != dest_x; x += x_increment)
	{
		//noc->messages_in_flight[x][src_y] += adjust_original;

		int direction;
		if (x_increment > 0)
		{
			direction = EAST;
		}
		else
		{
			direction = WEST;
		}
		if (x == src_x)
		{
			noc->noc_messages_in[x][src_y][4 + (m->src_neuron->core->id % ARCH_MAX_CORES_PER_TILE)] += adjust;
			//noc->noc_messages_in[x][src_y][m->src_neuron->core->id % 4] += adjust / 2;
		}
		else
		{
			noc->noc_messages_in[x][src_y][direction] += adjust;
		}
		//noc->noc_messages_out[x][src_y][direction] += adjust;
		prev_direction = direction;
	}
	for (int y = src_y; y != dest_y; y += y_increment)
	{
		int direction;
		if (src_x == dest_x && y == src_y)
		{
			noc->noc_messages_in[dest_x][y][4 + (m->src_neuron->core->id % ARCH_MAX_CORES_PER_TILE)] += adjust;
			//noc->noc_messages_in[dest_x][y][m->src_neuron->core->id % 4] += adjust / 2;
		}
		else
		{
			noc->noc_messages_in[dest_x][y][prev_direction] += adjust;
		}

		if (y_increment > 0)
		{
			direction = NORTH;
		}
		else
		{
			direction = SOUTH;
		}
		noc->noc_messages_out[dest_x][y][direction] += adjust;
		//noc->messages_in_flight[dest_x][y] += adjust_original;
		prev_direction = direction;
	}

	if ((src_x == dest_x) && (src_y == dest_y))
	{
		//noc->noc_messages_in[dest_x][dest_y][m->src_neuron->core->id % 4] += adjust / 2;
		noc->noc_messages_in[dest_x][dest_y][4 + (m->src_neuron->core->id % ARCH_MAX_CORES_PER_TILE)] += adjust;
	}
	else
	{
		noc->noc_messages_in[dest_x][dest_y][prev_direction] += adjust;
	}
	noc->noc_messages_out[dest_x][dest_y][4 + (m->dest_neuron->core->id % ARCH_MAX_CORES_PER_TILE)] += adjust;
	//noc->noc_messages_out[dest_x][dest_y][m->dest_neuron->core->id % 4] += adjust;
	//noc->messages_in_flight[dest_x][dest_y] += adjust_original;

	// Update rolling averages and message counters
	if (message_in)  // Message entering NoC
	{
		noc->mean_in_flight_receive_delay +=
			(m->receive_delay - noc->mean_in_flight_receive_delay) /
				(noc->messages_in_noc + 1);
		noc->mean_hops += (m->hops - noc->mean_hops) / (noc->messages_in_noc + 1);
		noc->messages_in_noc++;
	}
	else // Message leaving the NoC
	{
		//printf("Message leaving NoC\n");
		if ((noc->messages_in_noc) > 1)
		{
			noc->mean_in_flight_receive_delay +=
			(noc->mean_in_flight_receive_delay - m->receive_delay) /
				(noc->messages_in_noc - 1);
			noc->mean_hops += (noc->mean_hops - m->hops) /
				(noc->messages_in_noc - 1);
		}
		else
		{
			noc->mean_in_flight_receive_delay = 0.0;
		}

		noc->messages_in_noc--;
	}
	/*
	long int noc_capacity = 0;
	long int total_links_in, links_in, total_links_out, links_out;
	total_links_out = 0;
	total_links_in = 0;

	for (int x = 0; x < noc->noc_width; x++)
	{
		for (int y = 0; y < noc->noc_height; y++)
		{
			links_out = 0;
			links_in = 0;
			double total_in, total_out;
			total_in = 0.0;
			total_out = 0.0;
			for (int link = 0; link < 4+ARCH_MAX_CORES_PER_TILE; link++)
			{
				if (noc->noc_messages_in[x][y][link] > 0.01)
				{
					noc_capacity +=
						MAX_MESSAGES_PER_HOP;
					links_in++;
					total_links_in++;
				}

				total_in += noc->noc_messages_in[x][y][link];
				total_out += noc->noc_messages_out[x][y][link];
				if (noc->noc_messages_out[x][y][link] > 0.01)
				{
					links_out++;
					total_links_out++;
				}
			}
			// Links out
			//
			for (int link = 4; link < 4+ARCH_MAX_CORES_PER_TILE; link++)
			{
				if (noc->noc_messages_out[x][y][link] > 0.01)
				{
					noc_capacity +=
						MAX_MESSAGES_PER_HOP;
				}
			}
			//INFO("x:%d y:%d links in:%ld links out:%ld total_in:%lf total out:%lf\n", x, y, links_in, links_out, total_in, total_out);
			//
		}
	}
	//INFO("total links in:%ld total links out:%ld\n", total_links_in, total_links_out);
	*/

	/*
	noc->noc_capacity = MAX(noc->noc_capacity, noc_capacity);
	noc->noc_capacity = 6144;
	if (noc->noc_capacity <= 0.0)
	{
		noc->utilization = 1.0;
	}
	else
	{
		noc->utilization =
			(double) noc->messages_in_noc / noc->noc_capacity;
		// Clip between [0,1]
		noc->utilization = fmin(1.0, noc->utilization);
		noc->utilization = fmax(0.0, noc->utilization);
	}

	TRACE1("total messages:%ld total capacity:%ld ratio:%lf\n",
		noc->messages_in_noc, noc->noc_capacity, noc->utilization);
	if (m->dest_neuron != NULL)
	{
		TRACE1("m:%d->%d in:%d\n", m->src_neuron->id, m->dest_neuron->id, message_in);
	}
	*/

	return;
}

double sim_calculate_messages_along_route(struct message *m,
	struct noc_timings *noc)
{
	int src_x, dest_x, src_y, dest_y;
	int x_increment, y_increment;
	double messages_along_route, messages_along_route_original;

	assert(m != NULL);
	src_x = m->src_neuron->core->t->x;
	dest_x = m->dest_neuron->core->t->x;
	src_y = m->src_neuron->core->t->y;
	dest_y = m->dest_neuron->core->t->y;
	messages_along_route = 0;
	messages_along_route_original = 0;
	int direction;

	if (src_x < dest_x)
	{
		x_increment = 1;
	}
	else
	{
		x_increment = -1;
	}
	if (src_y < dest_y)
	{
		y_increment = 1;
	}
	else
	{
		y_increment = -1;
	}
	//messages_along_route_original +=
	//		noc->messages_in_flight[src_x][src_y];
	int prev_direction;
	for (int x = src_x; x != dest_x; x += x_increment)
	{
		int direction = 0;
		if (x_increment > 0)
		{
			direction = EAST;
		}
		else
		{
			direction = WEST;
		}
		messages_along_route_original +=
			noc->messages_in_flight[x][src_y];
		if (x == src_x)
		{
			messages_along_route += noc->noc_messages_in[x][src_y][
				4+m->src_neuron->core->id % ARCH_MAX_CORES_PER_TILE];
		}
		else
		{
			messages_along_route += noc->noc_messages_in[x][src_y][direction];
		}
		//messages_along_route +=
		//	noc->noc_messages_out[x][src_y][direction];
		prev_direction = direction;
	}

	for (int y = src_y; y != dest_y; y += y_increment)
	{
		if (src_x == dest_x && y == src_y)
		{
			messages_along_route += noc->noc_messages_in[dest_x][y][
				4+m->src_neuron->core->id % ARCH_MAX_CORES_PER_TILE];
		}
		else
		{
			messages_along_route +=
				noc->noc_messages_in[dest_x][y][prev_direction];
		}

		if (y_increment > 0)
		{
			direction = NORTH;
		}
		else
		{
			direction = SOUTH;
		}
		messages_along_route_original +=
			noc->messages_in_flight[dest_x][y];
		//messages_along_route += noc->noc_messages_out[dest_x][y][direction];
		prev_direction = direction;
	}
	messages_along_route_original +=
		noc->messages_in_flight[dest_x][dest_y];
	// Handle the last tile
	if (src_x == dest_x && src_y == dest_y)
	{
		messages_along_route +=
			noc->noc_messages_in[dest_x][dest_y][4+m->src_neuron->core->id % ARCH_MAX_CORES_PER_TILE];
	}
	else
	{
		messages_along_route +=
			noc->noc_messages_in[dest_x][dest_y][prev_direction];
	}
	//messages_along_route +=
	//	noc->noc_messages_out[dest_x][dest_y][m->dest_neuron->core->id % 4];

	//printf("original:%f new:%f\n", messages_along_route_original, messages_along_route);
	assert(messages_along_route >= -0.1);
	//assert(fabs(messages_along_route - messages_along_route_original) < 0.01);
	return (int) messages_along_route;
}

void sim_update_noc(const double t, struct noc_timings *noc)
{
	// TODO: pass noc dimensions in
	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		struct message_fifo *q;
		struct message *m;

		q = &(noc->messages_received[i]);
		m = q->tail;
		while (m != NULL)
		{
			if (m->in_noc && (t >= m->received_timestamp))
			{
				// Mark the message as not in the NoC, moving it
				//  from the network to the receiving core
				m->in_noc = 0;
				sim_message_fifo_pop(q);
				// Go along the message path and decrement tile
				//  message counters
				sim_update_noc_message_counts(m, noc, 0);
			}
			m = m->next;
		}
	}

	return;
}

double sim_schedule_messages(struct message_fifo *const messages_sent,
				const int noc_width, const int noc_height)
{
	struct message_fifo *priority_queue;
	struct noc_timings noc;
	double last_timestamp;

	// TODO: we need more information about the NoC dimensions here,
	//  known from the architecture description
	noc.messages_in_noc = 0;
	noc.mean_in_flight_receive_delay = 0.0;
	noc.mean_hops = 0.0;
	noc.noc_capacity = 0;
	noc.noc_width = noc_width;
	noc.noc_height = noc_height;
	for (int x = 0; x < noc.noc_width; x++)
	{
		for (int y = 0; y < noc.noc_height; y++)
		{
			noc.messages_in_flight[x][y] = 0.0;
			for (int k = 0; k < 4+ARCH_MAX_CORES_PER_TILE; k++)
			{
				noc.noc_messages_in[x][y][k] = 0.0;
				noc.noc_messages_out[x][y][k] = 0.0;
			}
		}
	}
	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		sim_init_fifo(&(noc.messages_received[i]));
	}

	priority_queue = sim_init_timing_priority(messages_sent);
	last_timestamp = 0.0;
	// Setup timing counters
	TRACE1("Scheduling global order of messages.\n");

	// Each core has a queue of received messages.
	// A structure tracks how many in-flight
	//  messages are in the NoC and occupy each tile. We need to
	//  track the number of messages passing through each tile at
	//  the point of sending, and the average processing delay of
	//  all of those messages. When a message is added to the NoC
	//  we update the average, when a message is removed from the
	//  NoC we update the average. If a route becomes congested to
	//  the point it exceeds the capacity (based on the spikes
	//  buffered per hop) then we delay based on the average time
	//  to process network messages
	while (priority_queue != NULL)  // Queue isn't empty
	{
		struct message_fifo *q;
		struct message *next_message, *m;

		// Get the core's queue with the earliest simulation time
		q = sim_pop_priority_queue(&priority_queue);
		m = sim_message_fifo_pop(q);
		last_timestamp = fmax(last_timestamp, m->sent_timestamp);

		// Update the Network-on-Chip state
		sim_update_noc(m->sent_timestamp, &noc);

		// Messages without a destination (neuron) are dummy messages,
		//  that account for processing time that does not result in any
		//  spike messages. Normal messages are sent from a src neuron to
		//  a dest neuron
		if (m->dest_neuron != NULL)
		{
			const int dest_core = m->dest_neuron->core->id;
			// First, figure whether we are able to send a message
			//  into the network i.e., is the route to the dest core
			//  saturated and likely to block? Sum along the route
			//  and see if we have saturated route.
			double earliest_received_time;
			int messages_along_route;
			int path_capacity;

			messages_along_route = (int)
				sim_calculate_messages_along_route(m, &noc);

			// TODO: parameterize average path capacity
			//path_capacity = (m->hops+2) * MAX_MESSAGES_PER_HOP;
			path_capacity = (m->hops+1) * MAX_MESSAGES_PER_HOP;
			if (messages_along_route > path_capacity)
			{
				m->sent_timestamp +=
					(messages_along_route - path_capacity) * noc.mean_in_flight_receive_delay;
				//m->sent_timestamp += noc.mean_in_flight_receive_delay;
			}

			// Now, push the message into the right receiving queue
			//  Calculate the network delay and when the message
			//  is received
			m->in_noc = 1;
			sim_message_fifo_push(
				&(noc.messages_received[dest_core]), m);
			// Update the rolling average for message
			//  receiving times in-flight in the network
			sim_update_noc_message_counts(m, &noc, 1);

			//double network_delay = fmin(path_capacity, messages_along_route) * noc.mean_in_flight_receive_delay / (m->hops+2);
			//double network_delay = messages_along_route * noc.mean_in_flight_receive_delay / (m->hops+2);
			double network_delay = messages_along_route * noc.mean_in_flight_receive_delay / (m->hops+1);
			//double network_delay = messages_along_route * noc.mean_in_flight_receive_delay / (noc.mean_hops+1);
			//double network_delay = fmin(messages_along_route, path_capacity) * noc.mean_in_flight_receive_delay; // * noc.utilization;
			TRACE1("path capacity:%d messages:%d\n", path_capacity, messages_along_route);
			//double network_delay = messages_along_route * noc.mean_in_flight_receive_delay / 5.0;
			//double network_delay = m->network_delay;

			//printf("messages along:%d\n", messages_along_route);
			earliest_received_time =
				m->sent_timestamp + fmax(m->network_delay, network_delay);
			m->received_timestamp =
				fmax(m->dest_neuron->core->blocked_until,
					earliest_received_time);
			m->dest_neuron->core->blocked_until = fmax(
				(m->dest_neuron->core->blocked_until +
				m->receive_delay),
				(earliest_received_time + m->receive_delay));
			m->processed_timestamp =
				m->dest_neuron->core->blocked_until;
			last_timestamp =
				fmax(last_timestamp, m->processed_timestamp);
		}

		// Get the next message for this core
		next_message = q->tail;
		if (next_message != NULL)
		{
			// If applicable, schedule this next message immediately
			//  after the current message finishes sending
			next_message->sent_timestamp = m->sent_timestamp +
				next_message->generation_delay;
			last_timestamp = fmax(last_timestamp,
				next_message->sent_timestamp);
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
double sim_schedule_messages_blocking(struct message_fifo *const messages_sent)
{
	// TODO: this code doesn't really seem to work at all... it overly
	//  aggressively blocks and so doesn't really add any value
	struct message_fifo *priority_queue;
	double last_timestamp, t;

	priority_queue = sim_init_timing_priority(messages_sent);
	last_timestamp = 0.0;
	// Setup timing counters
	TRACE1("Scheduling global order of messages.\n");

	double average_receiving_delay[128][128];
	double average_send_time[128][128];
	int path_counts[128][128];
	int path_hops[128][128];
	int messages_until_path_full[128][128];

	for (int i = 0; i < 128; i++)
	{
		if ((messages_sent[i].count > 0))
		{
			struct message *m = (messages_sent[i]).tail;

			for (int j = 0; j < 128; j++)
			{
				path_counts[i][j] = 0;
				average_receiving_delay[i][j] = 0.0;
				average_send_time[i][j] = 0.0;
			}

			while ((m != NULL) && (m->dest_neuron != NULL))
			{
				assert(i == m->src_neuron->core->id);
				assert(m->dest_neuron->core->id >= 0 &&
					m->dest_neuron->core->id < ARCH_MAX_CORES);
				average_send_time[i][m->dest_neuron->core->id] +=
					m->generation_delay;
				average_receiving_delay[i][m->dest_neuron->core->id] +=
					m->receive_delay;
				path_counts[i][m->dest_neuron->core->id]++;
				path_hops[i][m->dest_neuron->core->id] = m->hops;
				m = m->next;
			}
		}
	}
	for (int i = 0; i < 128; i++)
	{
		for (int j = 0; j < 128; j++)
		{
			if (path_counts[i][j] > 0)
			{
				average_receiving_delay[i][j] /= path_counts[i][j];
				average_send_time[i][j] /= path_counts[i][j];

				double send_rate = 1.0 / average_send_time[i][j];
				double receive_rate = 1.0 / average_receiving_delay[i][j];
				double fill_rate = send_rate / (send_rate - receive_rate);
				int path_capacity = 32 + path_hops[i][j] * 16;
				messages_until_path_full[i][j] = (int) path_capacity * fill_rate;
			}
			else
			{
				average_receiving_delay[i][j] = 0.0;
			}
		}
	}

	// While queue isn't empty
	while (priority_queue != NULL)
	{
		// Get the core with the earliest simulation time
		struct message_fifo *q =
			sim_pop_priority_queue(&priority_queue);
		struct message *m = sim_message_fifo_pop(q);
		last_timestamp = fmax(last_timestamp, m->generation_delay);

		if (m->dest_neuron != NULL)
		{
			// Set time-stamps, calculating when the receiving H/W
			//  will be busy until
			last_timestamp =
				fmax(last_timestamp, m->sent_timestamp);

			if (messages_until_path_full[m->src_neuron->core->id][m->dest_neuron->core->id] <= 0)
			{
				printf("adding:%e\n", average_receiving_delay[m->src_neuron->core->id][m->dest_neuron->core->id]);
				m->sent_timestamp +=
					average_receiving_delay[m->src_neuron->core->id][m->dest_neuron->core->id];
			}
			else
			{
				messages_until_path_full[m->src_neuron->core->id][m->dest_neuron->core->id]--;
			}

			m->dest_neuron->core->blocked_until = fmax(
				(m->dest_neuron->core->blocked_until + m->receive_delay),
				(m->sent_timestamp + m->network_delay + m->receive_delay));
			m->processed_timestamp =
				m->dest_neuron->core->blocked_until;

			last_timestamp =
				fmax(last_timestamp, m->processed_timestamp);

			INFO("\t(cid:%d.%d) synapse at %d.%d busy until %e\n",
				m->src_neuron->core->t->id, m->src_neuron->core->id, m->dest_neuron->core->t->id,
				m->dest_neuron->core->id, m->dest_neuron->core->blocked_until);
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
			m->sent_timestamp = t + m->generation_delay;
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
*/

double sim_schedule_messages_simplest(struct message_fifo *const messages_sent)
{
	struct message_fifo *priority_queue;
	double last_timestamp, t;

	priority_queue = sim_init_timing_priority(messages_sent);
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
		last_timestamp = fmax(last_timestamp, m->generation_delay);

		if (m->dest_neuron != NULL)
		{
			// Set time-stamps, calculating when the receiving H/W
			//  will be busy until
			last_timestamp =
				fmax(last_timestamp, m->sent_timestamp);
			m->dest_neuron->core->blocked_until = fmax(
				(m->dest_neuron->core->blocked_until +
				m->network_delay + m->receive_delay),
				(m->sent_timestamp + m->receive_delay));
			m->processed_timestamp =
				m->dest_neuron->core->blocked_until;
			last_timestamp =
				fmax(last_timestamp, m->processed_timestamp);

			TRACE2("\t(cid:%d.%d) synapse at %d.%d busy until %e\n",
				c->t->id, c->id, m->dest_core->t->id,
				m->dest_core->id, m->dest_core->blocked_until);
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
			m->sent_timestamp = t + m->generation_delay;
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

void sim_process_neuron(struct timestep *const ts, struct neuron *n)
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
				sim_update_synapse(ts, axon, 1);
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
					ts, n, con->current);
			}
		}
	}
	else if (c->buffer_pos == BUFFER_SOMA)
	{
		n->processing_latency = sim_update_soma(ts, n, n->charge);
	}
	else if (c->buffer_pos == BUFFER_AXON_OUT)
	{
		if (n->fired)
		{
			struct soma_processor *soma = n->soma_hw;
			n->processing_latency = soma->latency_spiking;
			sim_neuron_send_spike_message(ts, n);
		}
	}
	TRACE1("Updating neuron %d.%d.\n", n->group->id, n->id);

	c->next_message.generation_delay += n->processing_latency;
	n->update_needed = 0;
	n->spike_count = 0;
}

double sim_pipeline_receive(
	struct timestep *const ts, struct core *c,
	struct connection_map *axon)
{
	// We receive a spike and process up to the time-step buffer
	double message_processing_latency = 0.0;

	TRACE1("Receiving messages for cid:%d\n", c->id);
	if (c->buffer_pos >= BUFFER_SYNAPSE)
	{
		int synaptic_lookup = 1;
		message_processing_latency =
			sim_update_synapse(ts, axon, synaptic_lookup);
	}

	return message_processing_latency;
}

struct message_fifo *sim_init_timing_priority(
	struct message_fifo *const message_queues)
{
	struct message_fifo *priority_queue;
	priority_queue = NULL;

	TRACE1("Initializing priority queue.\n");
	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		if ((message_queues[i]).count > 0) // messages
		{
			struct message *m = (message_queues[i]).tail;
			assert(m != NULL);
			m->sent_timestamp = m->generation_delay;
			sim_insert_priority_queue(&priority_queue,
				&(message_queues[i]));
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
	assert(core_message_fifo->tail != NULL);

	//INFO("Inserting into priority queue.\n");
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

double sim_update_synapse(struct timestep *const ts,
	struct connection_map *axon, const int synaptic_lookup)
{
	// Update all synapses to different neurons in one core. If a synaptic
	//  lookup, read and accumulate the synaptic weights. Otherwise, just
	//  update filtered current and any other connection properties
	struct core *post_core;
	double latency, min_synaptic_resolution;

	latency = 0.0;
	post_core = axon->connections[0]->post_neuron->core;

	TRACE1("Updating synapses for (cid:%d)\n", axon->pre_neuron->id);
	while (axon->last_updated <= ts->timestep)
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
						ts, post_neuron, con->current);
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
					ts, post_neuron, con->current);
			}
		}
	}

	return latency;
}

double sim_update_dendrite(
	struct timestep *const ts, struct neuron *n, const double charge)
{
	// TODO: Support dendritic operations, combining the current in
	//  different neurons in some way, and writing the result to an output
	double dendritic_current, latency;
	latency = 0.0;

	dendritic_current = 0.0;
	while (n->dendrite_last_updated <= ts->timestep)
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
		latency += sim_update_soma(ts, n, dendritic_current);
	}

	return latency;
}

double sim_update_soma(
	struct timestep *const ts, struct neuron *n, const double current_in)
{
	double latency = 0.0;
	struct soma_processor *soma = n->soma_hw;

	TRACE1("nid:%d updating, current_in:%lf\n", n->id, current_in);
	if ((soma->model == NEURON_LIF) ||
		(soma->model == NEURON_STOCHASTIC_LIF))
	{
		latency += sim_update_soma_lif(ts, n, current_in);
	}
	else if (soma->model == NEURON_TRUENORTH)
	{
		latency += sim_update_soma_truenorth(ts, n, current_in);
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
	struct timestep *const ts, struct neuron *n, const double current_in)

{
	struct soma_processor *soma = n->soma_hw;
	double random_potential, latency = 0.0;

	// Calculate the change in potential since the last update e.g.
	//  integate inputs and apply any potential leak
	TRACE1("Updating potential, before:%f\n", n->potential);

	while (n->soma_last_updated <= ts->timestep)
	{
		n->potential *= n->leak_decay;
		n->soma_last_updated++;
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
			sim_neuron_send_spike_message(ts, n);
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
		// Neuron is turned on and potential write
		latency += n->soma_hw->latency_update_neuron;
		soma->neuron_updates++;
	}

	latency += n->soma_hw->latency_access_neuron;

	return latency;
}

double sim_update_soma_truenorth(
	struct timestep *const ts, struct neuron *n, const double current_in)
{
	struct soma_processor *soma = n->soma_hw;
	double v, latency = 0.0;
	int randomize_threshold;

	// Apply leak
	while (n->soma_last_updated <= ts->timestep)
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
			sim_neuron_send_spike_message(ts, n);
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

void sim_neuron_send_spike_message(struct timestep *const ts,
	struct neuron *n)
{
	struct core *c = n->core;
	TRACE1("nid:%d sending spike(s).\n", n->id);
	int core_id = n->core->id;


	for (int k = 0; k < n->maps_out_count; k++)
	{
		struct connection_map *dest_axon;
		struct message *m;
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
		m->generation_delay =
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
			for (int k = 0; k < c->synapse_count; k++)
			{
				synapse_energy +=
					c->synapse[k].spikes_processed *
					c->synapse[k].energy_spike_op;
			}
			for (int k = 0; k < c->soma_count; k++)
			{

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
			arch_init_message(&(c->next_message));

			c->axon_in.spike_messages_in = 0L;
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
				c->soma[k].neuron_updates = 0L;
				c->soma[k].neurons_fired = 0L;
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
	fprintf(fp, "%ld,", ts->packets_sent);
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
	fprintf(sim->message_trace_fp, "generation_delay,network_delay,");
	fprintf(sim->message_trace_fp, "processing_latency,blocking_latency,");
	fprintf(sim->message_trace_fp, "sent_timestamp,processed_timestamp\n");
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
	if (m->dest_neuron != NULL)
	{
		assert(m->dest_neuron->core != NULL);
		assert(m->dest_neuron->core->t != NULL);
		fprintf(sim->message_trace_fp, "%d.%d,", m->dest_neuron->core->t->id,
			m->dest_neuron->core->id);
	}
	else
	{
		// Dummy message
		fprintf(sim->message_trace_fp, "x.x,");
	}
	fprintf(sim->message_trace_fp, "%d,", m->hops);
	fprintf(sim->message_trace_fp, "%d,", m->spikes);
	fprintf(sim->message_trace_fp, "%le,", m->generation_delay);
	fprintf(sim->message_trace_fp, "%le,", m->network_delay);
	fprintf(sim->message_trace_fp, "%le,", m->receive_delay);
	fprintf(sim->message_trace_fp, "%le,", m->blocked_latency);
	fprintf(sim->message_trace_fp, "%le,", m->sent_timestamp);
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
