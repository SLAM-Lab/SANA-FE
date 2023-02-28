// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  sim.c
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "print.h"
#include "sim.h"
#include "network.h"
#include "arch.h"

struct sim_stats sim_timestep(struct network *const net,
				struct architecture *const arch,
				const int timestep,
				FILE *probe_spike_fp,
				FILE *probe_potential_fp,
				FILE *perf_fp)
{
	struct sim_stats stats;
	long int spikes_sent;

	// Start the next time-step
	sim_reset_measurements(net, arch);

	TRACE("Updating %d group(s).\n", net->neuron_group_count);
	net->total_neurons_fired = 0;

	for (int i = 0; i < arch->tile_count; i++)
	{
		struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			for (int k = 0; k < c->neuron_count; k++)
			{
				if (c->neurons[k] != NULL)
				{
					sim_start_next_timestep(c->neurons[k]);
					sim_process_neuron(net, c->neurons[k],
								timestep);
				}
			}
		}
	}
	sim_probe_log_timestep(probe_spike_fp, probe_potential_fp, net);

	// Send all spike messages
	spikes_sent = sim_input_spikes(net);
	spikes_sent += sim_send_messages(arch, net, timestep);

	// Performance statistics for this time step
	stats.time_steps = 1;
	stats.total_sim_time = sim_calculate_time(arch, &(stats.network_time));
	stats.total_energy = sim_calculate_energy(arch, stats.total_sim_time);
	stats.total_packets_sent = sim_calculate_packets(arch);
	stats.total_spikes = spikes_sent;

	INFO("Spikes sent: %ld\n", spikes_sent);
	if (perf_fp)
	{
		sim_perf_log_timestep(perf_fp, arch, net, &stats);
	}
	return stats;
}

void sim_process_neuron(struct network *net, struct neuron *n,
							const int timestep)
{
	// The neuron (state update) contains four main components
	// 1) synapse updates
	// 2) dendrite updates
	// 3) potential (soma) updates
	// 4) Axon updates

	// Depending on the first buffered stage, start processing from the
	//  first buffered value. Then compute for the rest of the pipeline
	//  regardless of the rest.
	if (!n->is_init)
	{
		return;
	}

	struct core *c = n->core;
	if (c->buffer_pos == BUFFER_SYNAPSE)
	{
		sim_update_synapse(n);
	}
	else if (c->buffer_pos == BUFFER_DENDRITE)
	{
		sim_update_dendrite(n);
	}
	else if (c->buffer_pos == BUFFER_SOMA)
	{
		double current_in = n->current_buffer;
		n->current_buffer = 0.0;
		n->processing_time = sim_update_soma(n, current_in, timestep);
	}
	TRACE("Updating neuron %d.%d.\n", n->group->id, n->id);

	n->update_needed = 0;
	n->spike_count = 0;
}

void sim_start_next_timestep(struct neuron *n)
{
	if (n == NULL)
	{
		return;
	}
	struct core *c = n->core;
	if (c->buffer_pos == BUFFER_SYNAPSE)
	{
		// Copy over all the synapses addresses to access
		n->charge = n->charge_buffer;
		n->charge_buffer = 0.0;
	}
	else if (c->buffer_pos == BUFFER_DENDRITE)
	{
		// Copy all synapse addresses
		n->d_currents[0] = n->d_currents_buffer[0];
		n->d_currents_buffer[0] = 0.0;
	}
	else if (c->buffer_pos == BUFFER_AXON_OUT)
	{
		// Copy packets
		n->fired = n->fired_buffer;
		n->fired_buffer = 0;
	}

	//n->current *= n->current_decay;
}

double sim_pipeline_receive(struct neuron *n, struct connection *connection_ptr,
							const int timestep)
{
	// We receive a spike and process up to the time-step buffer

	struct core *c = n->core;
	if (c->buffer_pos == BUFFER_AXON_IN)
	{
		// TODO
		return n->core->axon_in.time;
	}
	else
	{
		// TODO
	}

	if (c->buffer_pos == BUFFER_DENDRITE)
	{
		n->charge_buffer += connection_ptr->weight;
		return n->core->dendrite.time;
	}
	else
	{
		// Perform synaptic read
		n->charge += connection_ptr->weight;
		sim_update_synapse(n);
	}

	// Copy the charge to the soma current
	if (c->buffer_pos == BUFFER_SOMA)
	{
		n->current_buffer = n->charge;
		return n->soma_hw->time;
	}
	else
	{
		assert(0);
		return sim_update_soma(n, n->charge, timestep);
	}
}

int sim_send_messages(struct architecture *arch, struct network *net,
							const int timestep)
{
	struct core *top_priority;
	const int core_count = ARCH_MAX_CORES * ARCH_MAX_TILES; // TODO calculate
	int total_spike_count, cores_left;

	total_spike_count = 0;
	top_priority = sim_init_timing_priority(arch);

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
		int spike_ops;
		int nid = c->curr_neuron;
		struct neuron *n = c->neurons[nid];
		// Get the current update status
		if ((c->curr_neuron >= c->neuron_count) ||
						(n == NULL) || (!n->is_init))
		{
			TRACE("\t(cid:%d.%d): Neuron %d not used, get next\n",
					c->t->id, c->id, c->curr_neuron);
			c->status = NEURON_FINISHED;
		}

		if (c->status == UPDATE_NEURON)
		{
			// Model effect of stage prior to sending message
			//  This should just be the last stage latency
			c->time += n->processing_time;
			TRACE("(cid:%d.%d) nid:%d.%d update, time:%e\n",
				c->t->id, c->id, n->group->id, n->id, c->time);

			if (n->fired)
			{
				net->total_neurons_fired++;
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
					assert(post_neuron != NULL);
					assert(post_core != NULL);

					// TODO: need a function to figure this out
					//  or we store an absolute id in the core
					//  OR we use a 2D array
					int axon_id = (post_core->t->id * 4) + post_core->id;
					c->spikes_sent_per_core[axon_id]++;
					total_spike_count++;
					c->axon_map[axon_id] = post_core;
					sim_pipeline_receive(post_neuron, conn,
								timestep);
				}
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
			while ((c->curr_axon < core_count) &&
				c->spikes_sent_per_core[c->curr_axon] == 0)
			{
				c->curr_axon++;
			}
			if (c->curr_axon >= core_count)
			{
				// No cores left to send to
				c->status = NEURON_FINISHED;
			}
			else if (c->spikes_sent_per_core[c->curr_axon])
			{
				struct core *post_core = c->axon_map[c->curr_axon];
				double network_latency = 0.0;
				double spike_processing_time = 0.0;
				int memory_accesses;

				TRACE("\t(cid:%d.%d) Sending %d spikes to %d.%d\n",
					c->t->id, c->id,
					c->spikes_sent_per_core[post_core->id],
					post_core->t->id, post_core->id);

				struct tile *tile_pre = c->t;
				struct tile *tile_post = post_core->t;

				assert(tile_pre != NULL);
				assert(tile_post != NULL);

				// Calculate the energy and time for
				//  sending spike packets
				int x_hops, y_hops;

				// Axon access cost (currently assuming
				//  all axon access are in parallel,
				//  this might not be valid but
				//  could be good enough)
				c->axon_out.packets_out++;
				c->time += c->axon_out.time_access;

				c->axon_out.energy +=
					c->axon_out.energy_access;
				c->axon_out.energy +=
					c->t->energy_spike_within_tile;
				c->time += c->t->time_spike_within_tile;

				x_hops = abs(tile_pre->x - tile_post->x);
				y_hops = abs(tile_pre->y - tile_post->y);
				// E-W hops
				// Maybe it makes sense to have this as two
				//  options?
				tile_pre->energy += x_hops *
					tile_pre->energy_east_west_hop;
				network_latency += x_hops *
					tile_pre->time_east_west_hop;

				// N-S hops
				tile_pre->energy += y_hops *
					tile_pre->energy_north_south_hop;
				network_latency += y_hops *
					tile_pre->time_north_south_hop;

				if (post_core->is_blocking)
				{
					c->time = fmax(c->time,
						post_core->blocked_until);
				}
				if (tile_post->is_blocking)
				{
					c->time = fmax(c->time,
						tile_post->blocked_until);
				}
				arch->total_hops += x_hops;
				arch->total_hops += y_hops;
				tile_post->blocked_until = c->time + network_latency;

				c->time += network_latency;
				// Read word from memory, this is a very
				//  simplified model
				spike_ops = c->spikes_sent_per_core[c->curr_axon];
				memory_accesses =
					(spike_ops + (c->synapse.weights_per_word - 1)) /
					(c->synapse.weights_per_word);
				post_core->synapse.energy += memory_accesses * c->synapse.energy_memory_access;
				spike_processing_time += memory_accesses * c->synapse.time_memory_access;

				post_core->synapse.energy += spike_ops *
					post_core->synapse.energy_spike_op;
				spike_processing_time += spike_ops *
					post_core->synapse.time_spike_op;

				if (post_core->buffer_pos == BUFFER_SYNAPSE)
				{

				}
				else
				{
					post_core->synapse.total_spikes +=
						spike_ops;
					post_core->synapse.memory_reads +=
						memory_accesses;
				}

				if (post_core->is_blocking)
				{
					post_core->blocked_until =
						c->time + spike_processing_time;
				}
				else
				{
					post_core->blocked_until +=
						spike_processing_time;
				}

				TRACE("\t(cid:%d.%d) synapse at %d.%d busy until %e\n",
					c->t->id, c->id,
					post_core->t->id, post_core->id,
					post_core->synapse.busy_until);
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

struct core *sim_init_timing_priority(struct architecture *arch)
{
	struct core *next = NULL;

	// Initialize in reverse order, but assuming all cores start time
	//  synchronized (time ==), this is arbitrary
	for (int i = (arch->tile_count-1); i >= 0; i--)
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
			if (post_neuron->core->buffer_pos == BUFFER_SOMA)
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

double sim_update_soma_lif(struct neuron *n, const double current_in, const int timestep);
double sim_update_soma_truenorth(struct neuron *n, const double current_in, const int timestep);


double sim_update_soma(struct neuron *n, const double current_in,
							const int timestep)
{
	struct soma_processor *soma = n->soma_hw;

	if ((soma->model == NEURON_LIF) || (soma->model == NEURON_IF))
	{
		return sim_update_soma_lif(n, current_in, timestep);
	}
	else if (soma->model == NEURON_TRUENORTH)
	{
		return sim_update_soma_truenorth(n, current_in, timestep);
	}
	else
	{
		INFO("Neuron model not recognised: %d", soma->model);
		assert(0);
		return 0.0;
	}
}

double sim_update_soma_lif(struct neuron *n, const double current_in,
							const int timestep)

{
	struct soma_processor *soma = n->soma_hw;
	double latency = 0.0;

	// Calculate the change in potential since the last update e.g.
	//  integate inputs and apply any potential leak
	TRACE("Updating potential, before:%f\n", n->potential);

	if (soma->model == NEURON_LIF)
	{
		while (n->soma_last_updated <= timestep)
		{
			n->potential *= n->potential_decay;
			n->soma_last_updated++;
		}
	}

	// Add the spike potential
	n->potential += current_in + n->bias;
	n->current = 0.0;
	n->charge = 0.0;

#if 0
	// Clamp min potential
	n->potential = (n->potential < n->reset) ?
					n->reset : n->potential;
#endif
	TRACE("Updating potential, after:%f\n", n->potential);

	if ((n->force_update && (n->potential > n->threshold)) ||
		(!n->force_update && n->potential >= n->threshold))
	{
		n->potential = n->reset;
		sim_neuron_send_spike(n);
		TRACE("nid %d fired.\n", n->id);

		latency += soma->time_spiking;
		soma->energy += soma->energy_spiking;
		soma->spikes_sent++;
	}

	if (n->spike_count || n->force_update)
	{
		latency += n->soma_hw->time_active_neuron_update;
		soma->updates++;
		n->soma_hw->energy +=
			n->soma_hw->energy_active_neuron_update;
	}
	else
	{
		// inactive neuron update cost
		n->soma_hw->energy += n->soma_hw->energy_inactive_neuron_update;
		latency += n->soma_hw->time_inactive_neuron_update;
	}

	return latency;
}

double sim_update_soma_truenorth(struct neuron *n, const double current_in,
							const int timestep)
{
	struct soma_processor *soma = n->soma_hw;
	double v, latency = 0.0;
	int randomize_threshold;

	// Apply leak
	while (n->soma_last_updated <= timestep)
	{
		// Linear leak
		if (soma->leak_towards_zero)
		{
			// TODO: what happens if we're above zero but by less
			//  than the leak amount (for convergent), will we
			//  oscillate between the two? Does it matter
			if (n->potential > 0.0)
			{
				n->potential -= n->potential_decay;
			}
			else if (n->potential < 0.0)
			{
				n->potential += n->potential_decay;
			}
			// else equals zero, so no leak is applied
		}
		else
		{
			n->potential += n->potential_decay;
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

	INFO("v:%lf +vth:%lf mode:%d -vth:%lf mode:%d\n",
		v, n->threshold, n->group->reset_mode, n->reverse_threshold,
		n->group->reverse_reset_mode);
	if (v >= n->threshold)
	{
		int reset_mode = n->group->reset_mode;
		INFO("pos reset:%d\n", reset_mode);
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
		latency += sim_neuron_send_spike(n);
	}
	else if (v <= n->reverse_threshold)
	{
		int reset_mode = n->group->reverse_reset_mode;
		INFO("neg reset:%d\n", reset_mode);
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
	INFO("potential:%lf threshold %lf\n", n->potential, n->threshold);

	return latency;
}

double sim_neuron_send_spike(struct neuron *n)
{
	struct soma_processor *soma = n->soma_hw;
	double latency = soma->time_spiking;

	soma->energy += soma->energy_spiking;
	soma->spikes_sent++;

	if (n->core->buffer_pos == BUFFER_AXON_OUT)
	{
		n->fired_buffer = 1;
	}
	else
	{
		n->fired = 1;
	}

	return latency;
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
					fmax(c->time, c->blocked_until);
			max_core_time = fmax(max_core_time, this_core_time);
		}

		max_time = fmax(max_time, max_core_time);
	}

	// Add the mesh-wide barrier sync time (assuming worst case of 32 tiles)
	max_time += arch->time_barrier;
	//INFO("Simulated time for step is:%es\n", max_time);

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

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		network_energy += t->energy;

		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);
			synapse_energy += c->synapse.energy;
			soma_energy += c->soma.energy;
			axon_energy += c->axon_out.energy;
		}
	}

	total_energy = synapse_energy + dendrite_energy + soma_energy +
						axon_energy + network_energy;
	// Apply leakage factor
	//total_energy += time * 0.65;

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
			total_packets += c->axon_out.packets_out;
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
			n->processing_time = 0.0;
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
		t->blocked_until = 0.0;
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			// Reset core
			c->time = 0.0;
			c->energy = 0.0;
			c->status = UPDATE_NEURON;
			c->blocked_until = 0.0;

			c->axon_in.energy = 0.0;
			c->axon_in.time = 0;

			c->synapse.energy = 0.0;
			c->synapse.time = 0.0;

			c->dendrite.energy = 0.0;
			c->dendrite.time = 0.0;

			c->soma.energy = 0.0;
			c->soma.time = 0.0;

			c->axon_out.energy = 0.0;
			c->axon_out.time = 0.0;
			c->axon_out.packets_out = 0;
		}
	}
}

void sim_perf_write_header(FILE *fp, const struct architecture *arch)
{
	fprintf(fp, "time,");
	fprintf(fp, "fired,");
	fprintf(fp, "packets,");
	fprintf(fp, "hops,");
	fprintf(fp, "total_energy,");
	fprintf(fp, "wall_time,");
	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);
			fprintf(fp, "o[%d.%d].energy,",
						t->id, c->id);
			fprintf(fp, "o[%d.%d].time,",
						t->id, c->id);


			fprintf(fp, "s[%d.%d].energy,",
						t->id, c->id);
			fprintf(fp, "s[%d.%d].time,",
						t->id, c->id);

			fprintf(fp, "+[%d.%d].energy,",
						t->id, c->id);
			fprintf(fp, "+[%d.%d].time,",
						t->id, c->id);
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
	fprintf(fp, "%lf,", stats->total_energy);
	fprintf(fp, "%lf,", stats->wall_time);

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);

		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);

			fprintf(fp, "%e,", c->axon_out.energy);
			fprintf(fp, "%e,", c->axon_out.time);

			fprintf(fp, "%e,", c->synapse.energy);
			fprintf(fp, "%e,", c->synapse.time);

			fprintf(fp, "%e,", c->soma.energy);
			fprintf(fp, "%e,", c->soma.time);
		}
	}

	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		fprintf(fp, "%e,", t->energy);
		fprintf(fp, "%e,", t->time);
	}

	fprintf(fp, "\n");
}

void sim_write_summary(FILE *fp, const struct sim_stats *stats)
{
	// Write the simulation result to file
	fprintf(fp, "git_version: %s\n", GIT_COMMIT);
	fprintf(fp, "energy: %e\n", stats->total_energy);
	fprintf(fp, "time: %e\n", stats->total_sim_time);
	fprintf(fp, "total_spikes: %ld\n", stats->total_spikes);
	fprintf(fp, "total_packets: %ld\n", stats->total_packets_sent);
	fprintf(fp, "wall_time: %lf\n", stats->wall_time);
	fprintf(fp, "time_steps: %d\n", stats->time_steps);
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
