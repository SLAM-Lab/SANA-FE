// sim.c
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "print.h"
#include "sim.h"
#include "network.h"
#include "arch.h"

void sim_copy_buffers(struct neuron *n)
{
	TRACE("Current before: %lf.\n", n->current);
	// TODO: remove hack
	if (n->core->synapse->buffer_in)
	{
		// Copy over all the synapses addresses to access
		n->charge = n->charge_buffer;
	}
	if (n->core->dendrite->buffer_in)
	{
		// Copy all synapse addresses
		n->d_currents[0] = n->d_currents_buffer[0];
	}
	if (n->core->soma->buffer_in)
	{
		n->current = n->current_buffer;
	}
	if (n->core->axon_out->buffer_in)
	{
		// Copy packets
		n->fired = n->fired_buffer;
	}
	//n->current *= n->current_decay;
	n->d_currents_buffer[0] = 0.0;
	n->current_buffer = 0.0;
	n->charge_buffer = 0.0;
	n->fired_buffer = 0;
	TRACE("Current after: %lf.\n", n->current);
}

struct sim_stats sim_timestep(struct network *const net,
				struct architecture *const arch,
				FILE *probe_spike_fp,
				FILE *probe_potential_fp,
				FILE *perf_fp)
{
	struct sim_stats stats;
	long int spikes_sent;

	sim_reset_measurements(net, arch);
#if 0
	// TODO: remove this hack, just to simulate inputs
	for (int i = 0; i < net->external_input_count; i++)
	{
		// About 1% of neurons spiking
		net->external_inputs[i].val = 0.005;
	}
	INFO("Seeded %d inputs\n", net->external_input_count);
#endif
	sim_update(net);
	sim_probe_log_timestep(probe_spike_fp, probe_potential_fp, net);

	spikes_sent = sim_input_spikes(net);
	//INFO("Input spikes sent: %ld\n", spikes_sent);
	spikes_sent += sim_route_spikes(arch, net);

	// Performance statistics for this time step
	stats.time_steps = 1;
	stats.total_sim_time = sim_calculate_time(arch, &(stats.network_time));
	stats.total_sim_time = sim_calculate_time_old(arch);
	stats.total_energy = sim_calculate_energy(arch, stats.total_sim_time);
	stats.total_packets_sent = sim_calculate_packets(arch);
	stats.total_spikes = spikes_sent;

	//INFO("Spikes sent: %ld\n", spikes_sent);
	if (perf_fp)
	{
		sim_perf_log_timestep(perf_fp, arch, net, &stats);
	}
	return stats;
}

void sim_update(struct network *net)
{
	//#pragma omp parallel for
	TRACE("Updating %d group(s).\n", net->neuron_group_count);
	net->total_neurons_fired = 0;
	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);
			sim_copy_buffers(n);
		}
	}

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		TRACE("Updating %d neuron(s) in gid:%d.\n",
						group->neuron_count, group->id);
		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);
			if (!n->is_init)
			{
				continue;
			}
			sim_pipeline_send(net, n);
			n->update_needed = 0;
			n->spike_count = 0;
		}
	}
}

struct core *sim_init_timing_priority(struct architecture *arch)
{
	struct core *next = NULL;

	// Initialize in reverse order, but assuming all cores start time
	//  synchronized (time ==), this is arbitrary
	for (int i = (arch->tile_count-1); i >= 0; i--)
	//for (int i = 11; i >= 0; i--)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = (t->core_count-1); j >= 0; j--)
		{
			struct core *c = &(t->cores[j]);
			assert(c->time == 0.0);
			c->next_timing = next;
			next = c;
		}
	}

	for (struct core *curr = next; curr != NULL; curr = curr->next_timing)
	{
		TRACE("curr:%d.%d\n", curr->t->id, curr->id);\
	}

	return next;
}

struct core *sim_update_timing_queue(struct architecture *arch,
						struct core *top_priority)
{
	struct core *next, *tmp;

	if ((top_priority == NULL) || (top_priority->next_timing == NULL) ||
			(top_priority->time <= top_priority->next_timing->time))
	{
		// Do nothing
		return top_priority;
	}
	else
	{
		struct core *curr = top_priority;
		next = top_priority->next_timing;

		tmp = top_priority;
		top_priority = next;
		// reinsert core into the correct place in the priority list
		while (next != NULL)
		{
			if (tmp->time < next->time)
			{
				break;
			}
			curr = next;
			next = curr->next_timing;
		}
		curr->next_timing = tmp;
		tmp->next_timing = next;
	}

	for (struct core *tmp = top_priority; tmp != NULL;
							tmp = tmp->next_timing)
	{
		//INFO("tmp->time:%e (id:%d)\n", tmp->time, tmp->id);
		assert((tmp->next_timing == NULL) ||
					tmp->time <= tmp->next_timing->time);
	}

	return top_priority;
}

int sim_route_spikes_new(struct network *net, struct neuron *n)
{
	int total_spike_count;

	total_spike_count = 0;
	// First functionally simulate the effect of spikes on the neuron
	//  Update the charge at the receiving neurons, based on the synaptic
	//  weight
	if (!n->is_init || !n->fired)
	{
		return 0;
	}
	TRACE("Routing for nid:%d post_connection_count:%d\n",
			n->id, n->post_connection_count);
	// Send all spikes to receiving neurons, updating
	//  their synaptic current
	for (int k = 0; k < n->post_connection_count; k++)
	{
		struct neuron *post_neuron;
		struct connection *connection_ptr;

		connection_ptr = &(n->connections[k]);

		post_neuron = connection_ptr->post_neuron;
		post_neuron->update_needed = 1;
		post_neuron->spike_count++;
		total_spike_count++;
		assert(connection_ptr->pre_neuron->id == n->id);
		sim_pipeline_receive(post_neuron, connection_ptr);
	}

	//TRACE("Neurons fired: %ld\n", net->total_neurons_fired);

	return total_spike_count;
}

double sim_pipeline_receive(struct neuron *n, struct connection *connection_ptr)
{
	// We receive a spike and process up to the first buffer encountered

	if (n->core->axon_in->buffer_in)
	{
		// TODO
		return n->core->axon_in->time;
	}
	else
	{
		// TODO
	}

	if (n->core->dendrite->buffer_in)
	{
		n->charge_buffer += connection_ptr->weight;
		return n->core->dendrite->time;
	}
	else
	{
		// Perform synaptic read
		n->charge += connection_ptr->weight;
		sim_update_synapse(n);
	}

	// Copy the charge to the soma current
	if (n->soma_hw->buffer_in)
	{
		n->current_buffer = n->charge;
		return n->soma_hw->time;
	}
	else
	{
		assert(0);
		n->current = n->charge;
		sim_update_soma(n);
	}

	return n->axon_out->time;
}

int sim_route_spikes(struct architecture *arch, struct network *net)
{
	struct core *top_priority;
	const int core_count = 128; // TODO calculate
	int total_spike_count, cores_left;

	total_spike_count = 0;
	//net->total_neurons_fired = 0;

	top_priority = sim_init_timing_priority(arch);

	// now we go through all neurons to see which ones spiked
	//  for all neurons that spike, we send a spike to the corresponding axon output
	//  do we care about the axon output, surely only care about sending packets to other axon_inputs
	//   but we will store the packet tracking struct inside the axon output struct
	//  that axon_output sends to all relevant axon_inputs (we look up the dest neuron which must store a link to that)
	//  the axon_output and input must store the relevant router if it applies i.e. the next link in the chain

	// Functional vs hardware modelling. Should we separate the two. I
	//  thought it would be nice to be able to simulate the network
	//  without hardware. That means
	//  for all neurons:
	//      if neuron has fired:
	//          send spike to their connections (update charge for all receiving neurons)
	//          mark neuron as not fired? (although we could just reset in reset measurements)
	// The timing model instead iterates through hardware

	// First functionally simulate the effect of spikes on the neuron
	//  Update the charge at the receiving neurons, based on the synaptic
	//  weight
	/*
	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);
		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);

			if (!n->is_init || !n->fired)
			{
				continue;
			}
			net->total_neurons_fired++;
			TRACE("Routing for nid:%d post_connection_count:%d\n",
					n->id, n->post_connection_count);
			// Send all spikes to receiving neurons, updating
			//  their synaptic current
			for (int k = 0; k < n->post_connection_count; k++)
			{
				struct neuron *post_neuron;
				struct connection *connection_ptr;

				connection_ptr = &(n->connections[k]);
				assert(connection_ptr->pre_neuron->id == n->id);

				post_neuron = connection_ptr->post_neuron;
				post_neuron->update_needed = 1;
				if (post_neuron->core->dendrite->buffer_in)
				{
					post_neuron->charge_buffer +=
						connection_ptr->weight;
				}
				else
				{
					post_neuron->charge +=
							connection_ptr->weight;
				}

				// Copy the charge to the soma current
				if (post_neuron->soma_hw->buffer_in)
				{
					post_neuron->current_buffer =
							post_neuron->charge;
				}
				else
				{
					post_neuron->current =
							post_neuron->charge;
				}

				post_neuron->spike_count++;
				total_spike_count++;
			}
		}
	}
	*/
	// Setup timing counters
	cores_left = core_count;
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			c->neurons_left = c->neuron_count;
			c->curr_neuron = 0;
		}
	}

	while(cores_left && (top_priority != NULL))
	{
		// Get the core with the earliest simulation time
		struct core *c = top_priority;
		// Process the next neuron for this core
		int nid = c->curr_neuron;
		struct neuron *n = c->neurons[nid];
		// Get the current update status
		TRACE("curr_neuron:%d neuron_count:%d\n", c->curr_neuron, c->neuron_count);
		fflush(stdout);
		if ((c->curr_neuron >= c->neuron_count) || (!n->is_init))
		{
			TRACE("\t(cid:%d.%d): Neuron %d not used, get next\n",
					c->t->id, c->id, c->curr_neuron);
			c->status = NEURON_FINISHED;
		}
		// Hack to ignore time measurements for the input layer
		/*
		else if (c->neurons[c->curr_neuron]->group->id == 0) // input layer
		{
			c->status = NEURON_FINISHED;
		}
		*/

		if (c->status == UPDATE_NEURON)
		{
			// Model effect of stage prior to sending message
			//  This should just be the last stage latency
			c->time += n->latency;
			TRACE("(cid:%d.%d) nid:%d.%d update, time:%e\n",
				c->t->id, c->id, n->group->id, n->id, c->time);

			if (n->fired)
			{
				for (int k = 0; k < core_count; k++)
				{
					c->spikes_sent_per_core[k] = 0;
					c->axon_map[k] = NULL;
				}
				for (int k = 0; k < n->post_connection_count; k++)
				{
					struct connection *conn =
							&(n->connections[k]);
					struct neuron *post_neuron =
								conn->post_neuron;
					struct core *post_core = post_neuron->core;
					// TODO: need a function to figure this out
					//  or we store an absolute id in the core
					//  OR we use a 2D array
					int axon_id = (post_core->t->id * 4) + post_core->id;
					c->spikes_sent_per_core[axon_id]++;
					c->axon_map[axon_id] = post_core;
					// TODO: actually simulate the axon
					//  at the input. We send a message and
					//  not a single connection. The input
					//  axon does the reverse lookup, checking
					//  all the synapses from this core
					sim_pipeline_receive(post_neuron, conn);
				}
				TRACE("\t(cid:%d.%d) [", c->t->id, c->id);
				for (int k = 0; k < core_count; k++)
				{
					//printf("%d,",
					//	c->spikes_sent_per_core[k]);
				}
				TRACE("]\n");

				c->status = SEND_SPIKES;
				c->curr_axon = 0;
			}
			else
			{
				TRACE("\t(cid:%d.%d) Neuron didn't fire\n",
							c->t->id, c->id);
				c->status = NEURON_FINISHED;
			}
		}
		else if (c->status == SEND_SPIKES)
		{
			TRACE("(cid:%d.%d) send spikes\n", c->t->id, c->id);
			// Figure the cost of sending each packet, assuming that
			//  each packet can be sent in parallel update the
			//  timing of the receiving synapse.
			// Get next packet to send
			while ((c->curr_axon < 128) &&
				c->spikes_sent_per_core[c->curr_axon] == 0)
			{
				c->curr_axon++;
			}
			if (c->curr_axon >= 128)
			{
				// No cores left to send to
				c->status = NEURON_FINISHED;
			}
			else if (c->spikes_sent_per_core[c->curr_axon])
			{
				struct core *post_core = c->axon_map[c->curr_axon];
				double spike_processing_time = 0.0;
				int memory_accesses, spike_ops;

				TRACE("\t(cid:%d.%d) Sending %d spikes to %d.%d\n",
					c->t->id, c->id,
					c->spikes_sent_per_core[post_core->id],
					post_core->t->id, post_core->id);

				if (post_core->synapse[0].time > c->time)
				{
					TRACE("\t(cid:%d.%d) blocked until %e\n",
						c->t->id, c->id,
					post_core->synapse[0].time);
				}

				// Modeling tile as occupied while the network
				//  is used
				// Cost of accessing this axon entry
				c->time += 7.6e-9;
				c->axon_out[0].time += 7.6e-9;

				// Block the neuron pipeline until the
				//  receiving synapse is ready
				struct tile *tile_pre = c->t;
				struct tile *tile_post = post_core->t;
				//c->time = fmax(c->time, tile_post->busy_until);
				c->time = fmax(c->time,
					post_core->synapse[0].time);

				assert(tile_pre != NULL);
				assert(tile_post != NULL);
				// Calculate the energy and time for
				//  sending spike packets
				int x_hops, y_hops;

				// Axon access cost (currently assuming
				//  all axon access are in parallel,
				//  this might not be valid but
				//  could be good enough)
				c->axon_out[0].total_packets_sent++;
				c->axon_out[0].energy += 40.8e-12;
				c->time += c->axon_out[0].time_spike_within_tile;
				c->axon_out[0].time +=
					c->axon_out[0].time_spike_within_tile;
				tile_post->busy_until =
					fmax(tile_post->busy_until, c->time);
				//spike_processing_time +=
				//	c->axon_out[0].time_spike_within_tile;
				x_hops = abs(tile_pre->x - tile_post->x);
				y_hops = abs(tile_pre->y - tile_post->y);
				// E-W hops
				// TODO: record the energy of each hop at the
				//  corresponding router, means we need to
				//  iterate through all the x and y routers
				// Maybe it makes sense to have this as two
				//  options?
				tile_pre->energy += x_hops *
					tile_pre->energy_east_west_hop;
				//spike_processing_time += x_hops *
				//	tile_pre->time_east_west_hop;
				c->time += x_hops *
					tile_pre->time_east_west_hop;
				c->axon_out[0].time += x_hops *
					tile_pre->time_east_west_hop;
				//printf("xhops:%d 1 hop:%e total:%e\n", x_hops, tile_pre->time_east_west_hop, x_hops * tile_pre->time_east_west_hop);

				// N-S hops
				tile_pre->energy += y_hops *
					tile_pre->energy_north_south_hop;
				//spike_processing_time += y_hops *
				//	tile_pre->time_north_south_hop;
				c->time += y_hops *
					tile_pre->time_north_south_hop;
				c->axon_out[0].time += y_hops *
					tile_pre->time_north_south_hop;

				arch->total_hops += x_hops;
				arch->total_hops += y_hops;

				// Read word from memory, this is a very
				//  simplified model
				memory_accesses =
				(c->spikes_sent_per_core[c->curr_axon] + 7) / 8;
				post_core->synapse[0].energy += memory_accesses * 55.1e-12;
				//spike_processing_time += memory_accesses * 28.0e-9;
				spike_processing_time += memory_accesses * 23.0e-9;
				post_core->synapse[0].energy += c->spikes_sent_per_core[c->curr_axon] * post_core->synapse[0].energy_spike_op;
				//spike_ops = spikes_sent_per_core[post_core->id];
				spike_ops =
					(c->spikes_sent_per_core[c->curr_axon]
								+ 3) / 4;
				post_core->synapse[0].time += spike_ops * 4.0e-9;

				post_core->synapse[0].time =
					c->time + spike_processing_time;
					TRACE("\t(cid:%d.%d) synapse at %d.%d busy until %e\n",
					c->t->id, c->id,
					post_core->t->id, post_core->id,
					post_core->synapse[0].time);
				c->curr_axon++;
			}
		}
		else if (c->status == NEURON_FINISHED)
		{
			TRACE("(cid:%d.%d) neuron finished\n",
				c->t->id, c->id);
			if (n != NULL)
			{
				// Reset the neuron for the next time step
				n->fired = 0;
			}
			// Get the next neuron
			c->neurons_left--;
			c->curr_neuron++;
			if (c->neurons_left <= 0)
			{
				TRACE("\t(cid:%d.%d) finished simulating\n",
							c->t->id, c->id);
				// Remove core from queue, core has finished
				cores_left--;
				assert(cores_left >= 0);
				top_priority = top_priority->next_timing;
			}
			c->status = UPDATE_NEURON;
		}

		top_priority = sim_update_timing_queue(arch, top_priority);
		if (top_priority != NULL)
		{
			TRACE("\t(cid:%d.%d) time:%e\n",
				top_priority->t->id, top_priority->id,
				top_priority->time);
		}
	}
	TRACE("Neurons fired: %ld\n", net->total_neurons_fired);

	return total_spike_count;
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
			TRACE("not sending spike\n");
			continue;
		}
		for (int j = 0; j < in->post_connection_count; j++)
		{
			// Send a spike to all neurons connected to this input
			//  Normally, we would have a number of input dimensions
			//  for a given network
			struct neuron *post_neuron;
			struct connection *connection_ptr;

			connection_ptr = &(in->connections[j]);
			assert(connection_ptr);

			post_neuron = connection_ptr->post_neuron;
			TRACE("nid:%d Energy before: %lf\n",
					post_neuron->id, post_neuron->current);
			if (post_neuron->core->soma->buffer_in)
			{
				post_neuron->current_buffer += connection_ptr->weight;
			}
			else
			{
				post_neuron->current += connection_ptr->weight;
			}
			TRACE("nid:%d Energy after: %lf\n",
					post_neuron->id, post_neuron->current);

			post_neuron->synapse_hw->energy +=
					post_neuron->synapse_hw->energy_spike_op;
			post_neuron->synapse_hw->time +=
					post_neuron->synapse_hw->time_spike_op;

			post_neuron->update_needed = 1;
			post_neuron->spike_count++;
			input_spike_count++;
		}
		TRACE("Sent spikes to %d connections\n",
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
	//TRACE("Processed %d inputs.\n", arch->max_external_inputs);

	return input_spike_count;
}

void sim_pipeline_send(struct network *net, struct neuron *n)
{
	// The neuron (state update) contains four main components
	// 1) synapse updates
	// 2) dendrite updates
	// 3) potential (soma) updates
	// 4) Axon updates

	// Depending on the first buffered stage, start processing from the
	//  first buffered value. Then compute for the rest of the pipeline
	//  regardless of the rest.
	int first_buffered = -1;

	if (n->core->synapse->buffer_in)
	{
		first_buffered = 0;
		sim_update_synapse(n);
	}
	else if (first_buffered >= 0 || n->core->dendrite->buffer_in)
	{
		first_buffered = 1;
		sim_update_dendrite(n);
	}
	else if (first_buffered >= 0 || n->core->soma->buffer_in)
	{
		first_buffered = 2;
		sim_update_soma(n);
		// route spikes so go to axon_in, synapse
	}
	else if (first_buffered >= 0 || n->core->axon_out->buffer_in)
	{
		first_buffered = 3;
		sim_update_axon(n);
		// route spikes, axon_in, synapse
	}

	TRACE("Updating neuron %d.%d.\n", n->group->id, n->id);
	if (n->spike_count || n->force_update)
	{
		n->latency += n->soma_hw->time_active_neuron_update;
		n->soma_hw->energy +=
			n->soma_hw->energy_active_neuron_update;
	}
	else
	{
		// inactive neuron update cost
		n->soma_hw->energy += 47.5e-12;
		n->soma_hw->time += 6.1e-9;
		n->latency += 6.1e-9;
	}

	if (n->fired)
	{
		// In update axon?
		net->total_neurons_fired++;
		sim_route_spikes_new(net, n);
	}
}

void sim_update_synapse(struct neuron *n)
{
	// Current based (CUBA) LIF neuron model as implemented by Loihi
	//  Rather than iterate over all synapses we can simplify and just
	//  track total current. This is what nengo-loihi did, I'm not sure if
	//  this would ha ve to be changed if we had a more complicated synapse
	//  model


	return;
}

void sim_update_dendrite(struct neuron *n)
{
	// TODO: Support dendritic operations, combining the current in
	//  different neurons in some way, and writing the result to
	//  an output neuron
	return;
}

void sim_update_soma(struct neuron *n)
{
	// Calculate the change in potential since the last update e.g.
	//  integate inputs and apply any potential leak
	TRACE("Updating potential, before:%f\n", n->potential);
	//n->potential *= n->potential_decay;

	// Add the spike potential
	n->potential += n->current + n->bias;
	n->current = 0.0;
	n->charge = 0.0;

	// Clamp min potential
	//n->potential = (n->potential < n->reset) ?
	//				n->reset : n->potential;
	TRACE("Updating potential, after:%f\n", n->potential);

	if ((n->force_update && n->potential > n->threshold) ||
		(!n->force_update && n->potential >= n->threshold))
	{
		struct axon_output *out = n->axon_out;
		if (out->buffer_in)
		{
			n->fired_buffer = 1;
		}
		else
		{
			n->fired = 1;
		}

		n->potential = n->reset;

		// The spiking update time is much larger for the neuron than
		//  the normal update, in addition to the within-tile latency
		// TODO: add to the soma processor not the axon
		out->energy += 60.0e-12;
		n->latency += 30.0e-9;
		n->soma_hw->time += 30.0e-9;
		// TODO: separate h/w perf stuff from the general network
		//  calculations
		assert(out != NULL);
		out->energy += out->energy_spike_within_tile;
		TRACE("out->time=%e\n", out->time);
		TRACE("nid %d fired.\n", n->id);
	}
	else
	{
		struct axon_output *out = n->axon_out;
		if (out->buffer_in)
		{
			n->fired_buffer = 0;
		}
		else
		{
			n->fired = 0;
		}
	}

	return;
}

void sim_update_axon(struct neuron *n)
{
	// TODO
	return;
}

double sim_calculate_time(const struct architecture *const arch,
							double *network_time)
{
	// Returns the simulation time of the current timestep.
	//  This is calculated by finding the simulation time of each core,
	//  and simply taking the maximum of this.
	double max_time = 0.0;

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		double max_core_time = 0.0;
		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);
			double this_core_time = 0.0;
			this_core_time =
					fmax(c->time, c->synapse[0].time);
			max_core_time = fmax(max_core_time, this_core_time);
		}

		max_time = fmax(max_time, max_core_time);
	}

	// Add the mesh-wide barrier sync time (assuming worst case of 32 tiles)
	max_time += arch->time_barrier;
	//INFO("Simulated time for step is:%es\n", max_time);

	return max_time;
}

double sim_calculate_time_old(const struct architecture *const arch)
{
	// Returns the simulation time of the current timestep.
	//  This is calculated by finding the simulation time of each core,
	//  and simply taking the maximum of this.
	double tile_time, max_time = 0.0;

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		double max_core_time = 0.0;
		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);
			double unit_time, this_core_time = 0.0;

			unit_time = 0.0;
			for (int k = 0; k < c->axon_in_count; k++)
			{
				unit_time = fmax(unit_time, c->axon_in[k].time);
			}
			this_core_time += unit_time;

			unit_time = 0.0;
			for (int k = 0; k < c->synapse_count; k++)
			{
				unit_time = fmax(unit_time, c->synapse[k].time);
			}
			this_core_time += unit_time;

			unit_time = 0.0;
			for (int k = 0; k < c->dendrite_count; k++)
			{
				unit_time =
					fmax(unit_time, c->dendrite[k].time);
			}
			this_core_time += unit_time;

			unit_time = 0.0;
			for (int k = 0; k < c->soma_count; k++)
			{
				unit_time = fmax(unit_time, c->soma[k].time);
			}
			this_core_time += unit_time;

			unit_time = 0.0;
			for (int k = 0; k < c->axon_out_count; k++)
			{
				unit_time =
					fmax(unit_time, c->axon_out[k].time);
			}
			this_core_time += unit_time;

			max_core_time = fmax(max_core_time, this_core_time);
		}
		tile_time = max_core_time + t->time;
		max_time = fmax(max_time, tile_time);
	}

	// Add the mesh-wide barrier sync time (assuming worst case of 32 tiles)
	max_time += arch->time_barrier;
	INFO("Simulated time for step is:%es\n", max_time);

	return max_time;
}

double sim_calculate_energy(const struct architecture *const arch, const double time)
{
	// Returns the total energy across the design, for this timestep
	double total_energy, dendrite_energy, soma_energy, synapse_energy;
	double axon_energy, network_energy;

	synapse_energy = 0.0;
	dendrite_energy = 0.0;
	soma_energy = 0.0;
	axon_energy = 0.0;
	network_energy = 0.0;

	// TODO: fix this part, just get energy from soma and synapse units
	/*
	for (int i = 0; i < arch->max_neurons; i++)
	{
		const struct neuron *n = &(arch->neurons[i]);
		if (n->neuron_used)
		{
			neuron_energy += n->energy;
			for (int j = 0; j < n->post_connection_count; j++)
			{
				const struct synapse *s = &(n->synapses[j]);
				synapse_energy += s->energy;
			}
		}
	}
	*/

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		network_energy += t->energy;

		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);

			for (int k = 0; k < c->synapse_count; k++)
			{
				const struct synapse_processor *s =
							&(c->synapse[k]);
				synapse_energy += s->energy;
			}

			for (int k = 0; k < c->soma_count; k++)
			{
				const struct soma_processor *s =
							&(c->soma[k]);
				soma_energy += s->energy;
			}

			for (int k = 0; k < c->axon_out_count; k++)
			{
				const struct axon_output *out =
							&(c->axon_out[k]);
				axon_energy += out->energy;
			}
		}
	}

	total_energy = synapse_energy + dendrite_energy + soma_energy +
						axon_energy + network_energy;
	// TODO: formalize leakage amount
	total_energy += time * 0.05;
	total_energy += 1000.0;

	return total_energy;
}

long int sim_calculate_packets(const struct architecture *arch)
{
	// Calculate the total number of packets sent this timestep
	long int total_packets = 0;

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);
			for (int k = 0; k < ARCH_MAX_PROCESSORS; k++)
			{
				const struct axon_output *out =
							&(c->axon_out[k]);
				total_packets += out->total_packets_sent;
			}
		}
	}

	return total_packets;
}

void sim_reset_measurements(struct network *net, struct architecture *arch)
{
	// TODO: refactor in network resets vs architecture design resets?
	// Reset any neuron related measurements or settings for the beginning
	//  of the timestep
	arch->total_hops = 0;
	net->total_neurons_fired = 0;

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);
			n->update_needed |= n->force_update;
			n->latency = 0.0;
		}
	}

	// Reset any energy, time latency or other measurements of network
	//  hardware
	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);
		// Reset tile
		t->energy = 0.0;
		t->time = 0.0;
		t->busy_until = 0.0;
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			// Reset core
			c->time = 0.0;
			c->energy = 0.0;
			c->status = UPDATE_NEURON;

			for (int k = 0; k < ARCH_MAX_PROCESSORS; k++)
			{
				struct axon_input *in =
						&(c->axon_in[k]);
				in->energy = 0.0;
				in->time = 0;
			}

			for (int k = 0; k < ARCH_MAX_PROCESSORS; k++)
			{
				struct synapse_processor *s =
						&(c->synapse[k]);
				s->energy = 0.0;
				s->time = 0.0;
				s->busy_until = 0.0;
			}

			for (int k = 0; k < ARCH_MAX_PROCESSORS; k++)
			{
				struct dendrite_processor *d =
						&(c->dendrite[k]);
				d->energy = 0.0;
				d->time = 0.0;
			}

			for (int k = 0; k < ARCH_MAX_PROCESSORS; k++)
			{
				struct soma_processor *s =
						&(c->soma[k]);
				s->energy = 0.0;
				s->time = 0.0;
			}

			for (int k = 0; k < ARCH_MAX_PROCESSORS; k++)
			{
				struct axon_output *out =
						&(c->axon_out[k]);
				out->energy = 0.0;
				out->time = 0.0;
				out->total_packets_sent = 0;
			}
		}
	}
}

void sim_perf_write_header(FILE *fp, const struct architecture *arch)
{
	fprintf(fp, "time,");
	fprintf(fp, "fired,");
	fprintf(fp, "packets,");
	fprintf(fp, "hops,");
	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);

			for (int k = 0; k < c->axon_out_count; k++)
			{
				const struct axon_output *out =
							&(c->axon_out[k]);
				fprintf(fp, "o[%d.%d.%d].energy,",
							t->id, c->id, out->id);
				fprintf(fp, "o[%d.%d.%d].time,",
							t->id, c->id, out->id);
			}

			for (int k = 0; k < c->synapse_count; k++)
			{
				const struct synapse_processor *s =
							&(c->synapse[k]);

				fprintf(fp, "s[%d.%d.%d].energy,",
							t->id, c->id, s->id);
				fprintf(fp, "s[%d.%d.%d].time,",
							t->id, c->id, s->id);
			}

			for (int k = 0; k < c->soma_count; k++)
			{
				const struct soma_processor *s = &(c->soma[k]);

				fprintf(fp, "+[%d.%d.%d].energy,",
							t->id, c->id, s->id);
				fprintf(fp, "+[%d.%d.%d].time,",
							t->id, c->id, s->id);
			}
		}
	}

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		fprintf(fp, "t[%d].energy,", t->id);
		fprintf(fp, "t[%d].time,", t->id);
	}

	fprintf(fp, "\n");
}

void sim_perf_log_timestep(FILE *fp, const struct architecture *arch,
				const struct network *net, const struct sim_stats *stats)
{
	// Log the energy and time simulated at every neuron, dump it to
	//  a big csv file. Then we can post process it to pull out the parallel
	//  time. Time doesn't make sense per neuron, only per parallel
	//  block. Pull out energy for synapses, routers, axons and neurons
	fprintf(fp, "%e,", stats->total_sim_time);
	fprintf(fp, "%ld,", net->total_neurons_fired);
	fprintf(fp, "%ld,", stats->total_packets_sent);
	fprintf(fp, "%ld,", arch->total_hops);

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);

			for (int k = 0; k < c->axon_out_count; k++)
			{
				const struct axon_output *out =
							&(c->axon_out[k]);
				fprintf(fp, "%e,", out->energy);
				fprintf(fp, "%e,", out->time);
			}

			for (int k = 0; k < c->synapse_count; k++)
			{
				const struct synapse_processor *s =
							&(c->synapse[k]);

				fprintf(fp, "%e,", s->energy);
				fprintf(fp, "%e,", s->time);
			}

			for (int k = 0; k < c->soma_count; k++)
			{
				const struct soma_processor *s = &(c->soma[k]);

				fprintf(fp, "%e,", s->energy);
				fprintf(fp, "%e,", s->time);
			}
		}
	}

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		fprintf(fp, "%e,", t->energy);
		fprintf(fp, "%e,", t->time);
	}

	fprintf(fp, "\n");

	/*
	for (int i = 0; i < arch->max_neurons; i++)
	{
		double synapse_energy;
		struct neuron *n = &(arch->neurons[i]);

		// Calculate the total energy at all synapses
		synapse_energy = 0.0;
		//INFO("cid:%d energy:%e", n->id, n->energy);
		for (int j = 0; j < n->post_connection_count; j++)
		{
			struct synapse *s = &(n->synapses[j]);
			synapse_energy += s->energy;
		}

		if (n->neuron_used)
		{
			fprintf(fp, "%e,%e,", n->energy, synapse_energy);
		}
	}

	for (int i = 0; i < arch->max_axon_outputs; i++)
	{
		struct axon_output *out = &(arch->axon_outputs[i]);
		fprintf(fp, "%e,", out->energy);
	}

	for (int i = 0; i < arch->max_routers; i++)
	{
		struct router *r = &(arch->routers[i]);
		fprintf(fp, "%e,", r->energy);
	}

	fprintf(fp, "\n");
	*/
}

void sim_write_summary(FILE *fp, const struct sim_stats *stats)
{
	// Write the simulation result to file
	fprintf(fp, "git_version: %s\n", GIT_COMMIT);
	fprintf(fp, "energy: %e\n", stats->total_energy);
	fprintf(fp, "time: %e\n", stats->total_sim_time);
	fprintf(fp, "total_spikes: %ld\n", stats->total_spikes);
	fprintf(fp, "total_packets: %ld\n", stats->total_packets_sent);
	//fprintf(fp, "network_time: %e\n", stats->network_time);
}

void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp,
						const struct network *net)
{
	// Write csv header for probe outputs - record which neurons have been
	//  probed
	for (int i = 0; i < net->external_input_count; i++)
	{
		const struct input *in = &(net->external_inputs[i]);

		if (spike_fp)
		{
			fprintf(spike_fp, "i.%d,", in->id);
		}
		if (potential_fp)
		{
			fprintf(potential_fp, "i.%d,", in->id);
		}
	}

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		const struct neuron_group *group = &(net->groups[i]);
		for (int j = 0; j < group->neuron_count; j++)
		{
			const struct neuron *n = &(group->neurons[j]);

			if (spike_fp && n->log_spikes)
			{
				fprintf(spike_fp, "%d.%d,", group->id, n->id);
			}
			if (potential_fp && n->log_voltage)
			{
				fprintf(potential_fp, "%d.%d,",
							group->id, n->id);
			}
		}
	}

	if (spike_fp)
	{
		fputc('\n', spike_fp);
	}
	if (potential_fp)
	{
		fputc('\n', potential_fp);
	}

	return;
}

void sim_probe_log_timestep(FILE *spike_fp, FILE *potential_fp,
				const struct network *net)
{
	// Each line of this csv file is the output of each probed neuron
	//  Currently we only probe the voltage / potential and the spikes
	//  (e.g. for a spike raster plot).
	//
	// In the future, we might want to generalise this to probe any of the
	//  variable parameters
	int spike_probe_count, potential_probe_count;

	spike_probe_count = 0;
	potential_probe_count = 0;

	for (int i = 0; i < net->external_input_count; i++)
	{
		const struct input *in = &(net->external_inputs[i]);

		if (spike_fp)
		{
			fprintf(spike_fp, "%d,", in->send_spike);
		}

		if (potential_fp)
		{
			fprintf(potential_fp, "%lf,", in->spike_val);
		}
	}

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		const struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			const struct neuron *n = &(group->neurons[j]);

			if (spike_fp && n->log_spikes)
			{
				fprintf(spike_fp, "%d,", n->fired);
				spike_probe_count++;
			}
			if (potential_fp && n->log_voltage)
			{
				fprintf(potential_fp, "%lf,", n->potential);
				potential_probe_count++;
			}
		}
	}

	// Each timestep takes up a line in the respective csv file
	if (spike_fp && (spike_probe_count > 0))
	{
		fputc('\n', spike_fp);
	}
	if (potential_fp && (potential_probe_count > 0))
	{
		fputc('\n', potential_fp);
	}

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
	TRACE("rate input:%lf\n", firing_rate);
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

	TRACE("input fired value:%d\n", input_fired);
	return input_fired;
}
