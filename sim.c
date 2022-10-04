// sim.c
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
	spikes_sent = sim_input_spikes(net);
	INFO("Input spikes sent: %ld\n", spikes_sent);
	spikes_sent += sim_route_spikes(net);
	sim_update(net);

	// Performance statistics for this time step
	stats.time_steps = 1;
	stats.total_sim_time = sim_calculate_time(arch, &(stats.network_time));
	stats.total_energy = sim_calculate_energy(arch, stats.total_sim_time);
	stats.total_packets_sent = sim_calculate_packets(arch);
	stats.total_spikes = spikes_sent;

	INFO("Spikes sent: %ld\n", spikes_sent);
	sim_probe_log_timestep(probe_spike_fp, probe_potential_fp, net);
	if (perf_fp)
	{
		sim_perf_log_timestep(perf_fp, arch);
	}
	return stats;
}

void sim_update(struct network *net)
{
	//#pragma omp parallel for
	TRACE("Updating %d group(s).\n", net->neuron_group_count);
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

			if (n->update_needed)
			{
				sim_update_neuron(n);
			}
			else
			{
				TRACE("Skipping neuron %d.%d\n",
							n->group->id, n->id);
			}
		}
	}
}

int sim_route_spikes(struct network *net)
{
	int total_spike_count, total_neurons_fired;

	total_spike_count = 0;
	total_neurons_fired = 0;

	// now we go through all neurons to see which ones spiked
	//  for all neurons that spike, we send a spike to the corresponding axon output
	//  do we care about the axon output, surely only care about sending packets to other axon_inputs
	//   but we will store the packet tracking struct inside the axon output struct
	//  that axon_output sends to all relevant axon_inputs (we look up the dest neuron which must store a link to that)
	//  the axon_output and input must store the relevant router if it applies i.e. the next link in the chain

	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);
		for (int j = 0; j < group->neuron_count; j++)
		{
			struct axon_output *pre_axon;
			struct neuron *n = &(group->neurons[j]);
			int prev_core_id;

			if (!n->is_init || !n->fired)
			{
				continue;
			}
			total_neurons_fired++;

			pre_axon = n->axon_out;
			assert(pre_axon != NULL); // TODO for now require h/w mapping

			// Generate all the spikes for a spiking neuron
			//  Only generate one spike packet per core, that the
			//  neuron is broadcasting to
			assert(n->post_connection_count >= 0);

			// First reset tracking for packets that this neuron
			//  should sent in the timestep
			TRACE("Routing for nid:%d post_connection_count:%d\n",
					n->id, n->post_connection_count);
			for (int k = 0; k < n->post_connection_count; k++)
			{
				struct neuron *post_neuron;
				struct connection *connection_ptr;
				struct axon_input *axon_in;

				connection_ptr = &(n->connections[k]);
				assert(connection_ptr != NULL);
				post_neuron = connection_ptr->post_neuron;
				assert(post_neuron != NULL);
				axon_in = post_neuron->axon_in;

				axon_in->packet_size = 0;
			}

			// Next send all spikes to receiving neurons, updating
			//  their synaptic current and accounting for
			//  performance costs of those updates
			prev_core_id = -1;
			for (int k = 0; k < n->post_connection_count; k++)
			{
				struct neuron *post_neuron;
				struct connection *connection_ptr;
				struct axon_input *post_axon;

				connection_ptr = &(n->connections[k]);
				if(!connection_ptr->pre_neuron)
				{
					printf("group: %d\n", i);
					printf("neuron: %d\n", j);
					printf("id: %d\n", n->id);
					printf("connection: %d\n", k);
					printf("count: %d\n",
						n->post_connection_count);
					assert(0);
				}
				assert(connection_ptr->pre_neuron->id == n->id);

				post_neuron = connection_ptr->post_neuron;

				post_neuron->update_needed = 1;
				post_neuron->charge += connection_ptr->weight;
				post_neuron->spike_count++;
				total_spike_count++;

				// Hardware specific updates
				post_axon = post_neuron->axon_in;

				// TODO: refactor to function call
				if ((post_neuron->core->id != prev_core_id)  ||
								(k % 8) == 0)
				{
					post_neuron->synapse_hw->energy +=
								55.1e-12;
					post_neuron->synapse_hw->time +=
								15.0e-9;
				}
				prev_core_id = post_neuron->core->id;
				/*
				post_neuron->synapse_hw->energy +=
					post_neuron->synapse_hw->energy_spike_op;
				post_neuron->synapse_hw->time +=
					post_neuron->synapse_hw->time_spike_op;
				*/

				post_neuron->synapse_hw->energy +=
					post_neuron->synapse_hw->energy_spike_op;
				if ((k % 4) == 0)
				{
					// Four way parallelism
					post_neuron->synapse_hw->time += 4.0e-9;
				}

				// TODO: support AER and other representations
				//  Loihi sends a destination axon index
				//  So we need to store this with the neuron
				//  mapped, i.e. all the indexes
				post_axon->packet_size = 1;
			}

			// Finally model the network behavior and the
			//  performance cost of sending packets through the NoC
			for (int k = 0; k < n->post_connection_count; k++)
			{
				struct neuron *post_neuron;
				struct axon_input *post_axon;

				post_neuron = n->connections[k].post_neuron;
				post_axon = post_neuron->axon_in;

				// TODO: I think the code below should be a
				//  function and sim_send_spike can be inlined
				//  again?
				if ((pre_axon != NULL) && (post_axon != NULL) &&
						(post_axon->packet_size > 0))
				{
					struct tile *tile_pre = pre_axon->t;
					struct tile *tile_post = post_axon->t;

					assert(tile_pre != NULL);
					assert(tile_post != NULL);
					// Calculate the energy and time for
					//  sending spike packets
					int x_hops, y_hops;

					pre_axon->energy += 40.8e-12;
					pre_axon->time += 7.6e-9;

					x_hops =
						abs(tile_pre->x - tile_post->x);
					y_hops =
						abs(tile_pre->y - tile_post->y);
					// E-W hops
					// TODO: record the energy of each hop at the
					//  corresponding router, means we need to
					//  iterate through all the x and y routers
					// Maybe it makes sense to have this as two
					//  options?
					tile_pre->energy += x_hops *
						tile_pre->energy_east_west_hop;
					tile_pre->time += x_hops *
						tile_pre->time_east_west_hop;
					pre_axon->time += x_hops * tile_pre->time_east_west_hop;
					//printf("xhops:%d 1 hop:%e total:%e\n", x_hops, tile_pre->time_east_west_hop, x_hops * tile_pre->time_east_west_hop);

					// N-S hops
					tile_pre->energy += y_hops *
						tile_pre->energy_north_south_hop;
					tile_pre->time += y_hops *
						tile_pre->time_north_south_hop;
					pre_axon->time += y_hops * tile_pre->time_north_south_hop;
				}
				pre_axon->total_packets_sent++;
				post_axon->packet_size = 0;
			}
			n->fired = 0; // Reset the neuron for the next time step
		}

	}
	INFO("Neurons fired: %d\n", total_neurons_fired);

	return total_spike_count;
}

int sim_input_spikes(struct network *net)
{
	int input_spike_count = 0;

	// Seed all externally input spikes in the network for this timestep
	for (int i = 0; i < net->external_input_count; i++)
	{
		int send_spike;
		struct input *in = &(net->external_inputs[i]);

		if (in == NULL)
		{
			break;
		}

		if (in->type == INPUT_EVENT)
		{
			send_spike = (in->spike_val > 0.0);
		}
		else if (in->type == INPUT_POISSON)
		{
			send_spike = sim_poisson_input(in->rate);
		}
		else // INPUT_RATE)
		{
			send_spike = sim_rate_input(in->rate, &(in->spike_val));
		}

		if (!send_spike)
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
			post_neuron->charge += connection_ptr->weight;
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

void sim_update_neuron(struct neuron *n)
{
	// The neuron (state update) contains four main components
	// 1) synapse updates
	// 2) dendrite updates
	// 3) potential (soma) updates
	// 4) Axon updates
	sim_update_synapse_cuba(n);
	sim_update_dendrite(n);
	sim_update_potential(n);
	sim_update_axon(n);

	TRACE("Updating neuron %d.%d.\n", n->group->id, n->id);
	if (n->spike_count || n->force_update)
	{
		n->soma_hw->energy +=
			n->soma_hw->energy_active_neuron_update;
		n->soma_hw->time +=
			n->soma_hw->time_active_neuron_update;
	}
	else
	{
		n->soma_hw->energy += 47.5e-12;
		n->soma_hw->time += 6.1-9;
	}
	/*
	else if (n->force_update)
	{
		n->soma_hw->energy +=
			n->soma_hw->energy_inactive_neuron_update;
		n->soma_hw->time +=
			n->soma_hw->time_inactive_neuron_update;
	}
	*/
}

void sim_update_synapse_cuba(struct neuron *n)
{
	// Current based (CUBA) LIF neuron model as implemented by Loihi
	//  Rather than iterate over all synapses we can simplify and just
	//  track total current. This is what nengo-loihi did, I'm not sure if
	//  this would ha ve to be changed if we had a more complicated synapse
	//  model
	TRACE("Current before: %lf.\n", n->current);
	// TODO: remove hack
	n->current = n->charge;
	//n->current *= n->current_decay;
	n->charge = 0.0;
	TRACE("Current after: %lf.\n", n->current);

	return;
}

void sim_update_dendrite(struct neuron *n)
{
	// TODO: Support dendritic operations, combining the current in
	//  different neurons in some way, and writing the result to
	//  an output neuron
	return;
}

void sim_update_potential(struct neuron *n)
{
	// Calculate the change in potential since the last update e.g.
	//  integate inputs and apply any potential leak
	TRACE("Updating potential, before:%f\n", n->potential);
	//n->potential *= n->potential_decay;

	// Add the spike potential
	n->potential += n->current + n->bias;

	// Clamp min potential
	//n->potential = (n->potential < n->reset) ?
	//				n->reset : n->potential;

	TRACE("Updating potential, after:%f\n", n->potential);

	if (n->potential > n->threshold)
	{
		struct axon_output *out = n->axon_out;
		n->fired = 1;
		n->potential = n->reset;

		// The spiking update time is much larger for the neuron than
		//  the normal update, in addition to the within-tile latency
		// TODO: add to the soma processor not the axon
		out->energy += 60.0e-12;
		out->time += 30.0e-9;

		// Add the "within-tile" spike energy, this is the minimum cost
		//  of sending a spike
		// TODO: separate h/w perf stuff from the general network
		//  calculations
		assert(out != NULL);
		out->energy += out->energy_spike_within_tile;
		out->time += out->time_spike_within_tile;
		TRACE("out->time=%e\n", out->time);
		TRACE("nid %d fired.\n", n->id);
	}

	return;
}

void sim_update_axon(struct neuron *n)
{
	// TODO
	return;
}

double sim_calculate_time(const struct architecture *const arch, double *network_time)
{
	// Returns the simulation time of the current timestep.
	//  This is calculated by finding the simulation time of each core,
	//  and simply taking the maximum of this.
	double max_time = 0.0;

	double max_time_minus_network = 0.0;
	for (int i = 0; i < arch->tile_count; i++)
	{
		const struct tile *t = &(arch->tiles[i]);
		double max_core_time = 0.0;
		double max_core_time_minus_network = 0.0;
		for (int j = 0; j < t->core_count; j++)
		{
			const struct core *c = &(t->cores[j]);
			double unit_time, this_core_time = 0.0;
			double spike_time, neuron_time;
			double core_time_minus_network;


			unit_time = 0.0;
			spike_time = 0.0;
			neuron_time = 0.0;
			for (int k = 0; k < c->axon_in_count; k++)
			{
				unit_time = fmax(unit_time, c->axon_in[k].time);
			}
			spike_time += unit_time;

			unit_time = 0.0;
			for (int k = 0; k < c->synapse_count; k++)
			{
				unit_time = fmax(unit_time, c->synapse[k].time);
			}
			spike_time += unit_time;

			unit_time = 0.0;
			for (int k = 0; k < c->dendrite_count; k++)
			{
				unit_time =
					fmax(unit_time, c->dendrite[k].time);
			}
			neuron_time += unit_time;

			unit_time = 0.0;
			for (int k = 0; k < c->soma_count; k++)
			{
				unit_time = fmax(unit_time, c->soma[k].time);
			}
			neuron_time += unit_time;

			core_time_minus_network = fmax(spike_time, neuron_time);
			max_core_time_minus_network =
				fmax(max_core_time_minus_network, core_time_minus_network);

			unit_time = 0.0;
			for (int k = 0; k < c->axon_out_count; k++)
			{
				unit_time =
					fmax(unit_time, c->axon_out[k].time);
			}
			neuron_time += unit_time;

			this_core_time = fmax(spike_time, neuron_time);
			max_core_time = fmax(max_core_time, this_core_time);
		}

		max_time = fmax(max_time, max_core_time);
		max_time_minus_network = fmax(max_time_minus_network, max_core_time_minus_network);
	}

	//max_time += max_network_time;
	// Add the mesh-wide barrier sync time (assuming worst case of 32 tiles)
	*network_time = max_time - max_time_minus_network;
	max_time += arch->time_barrier;
	INFO("Simulated time for step is:%es\n", max_time);
	INFO("Simulated network time for step is:%es\n", *network_time);
	INFO("Network time is %lf%% of time\n", (*network_time / max_time) * 100.0);

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
	//total_energy += time * 0.05;

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
	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);
			n->update_needed = n->force_update;
			n->spike_count = 0;
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
		for (int j = 0; j < t->core_count; j++)
		{
			struct core *c = &(t->cores[j]);
			// Reset core
			c->time = 0.0;
			c->energy = 0.0;

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

void sim_perf_log_timestep(FILE *fp, const struct architecture *arch)
{
	// Log the energy and time simulated at every neuron, dump it to
	//  a big csv file. Then we can post process it to pull out the parallel
	//  time. Time doesn't make sense per neuron, only per parallel
	//  block. Pull out energy for synapses, routers, axons and neurons
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
	fprintf(fp, "network_time: %e\n", stats->network_time);
}

void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp,
						const struct network *net)
{
	// Write csv header for probe outputs - record which neurons have been
	//  probed
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
	input_fired = (rand_uniform < firing_probability);

	return input_fired;
}

int sim_rate_input(const double firing_rate, double *current)
{
	int input_fired;

	// Note: rate-based input (without randomization) is equivalent to a
	//  neuron with a fixed bias.
	TRACE("rate input:%lf\n", firing_rate);

	*current += firing_rate;
	if (*current >= 1.0)
	{
		*current = 0.0;
		input_fired = 1;
	}
	else
	{
		input_fired = 0;
	}

	TRACE("input fired value:%d\n", input_fired);
	return input_fired;
}
