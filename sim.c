// sim.c
// Neuromorphic simulator kernel
//
// Time-step based simulation, based on loop:
// 1) seed any input spikes
// 2) route spikes
// 3) update neurons and check firing

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "sim.h"
#include "tech.h"
#include "network.h"
#include "arch.h"

struct sim_results sim_timestep(const struct technology *tech,
				struct architecture *arch,
				FILE *probe_spike_fp,
				FILE *probe_potential_fp)
{
	struct sim_results results;
	long int spikes_sent;

	spikes_sent = sim_input_spikes(tech, arch);
	INFO("Input spikes sent: %ld\n", spikes_sent);
	sim_reset_measurements(arch);
	spikes_sent += sim_route_spikes(tech, arch);
	sim_update(tech, arch);

	// Return results for this timestep
	results.total_energy = sim_calculate_energy(arch);
	results.total_sim_time = sim_calculate_time(tech, arch);
	results.total_spikes = spikes_sent;
	INFO("Spikes sent: %ld\n", spikes_sent);

	sim_probe_log_timestep(probe_spike_fp, probe_potential_fp, arch);
	return results;
}

void sim_update(const struct technology *tech, struct architecture *arch)
{
	//#pragma omp parallel for
	for (int i = 0; i < arch->max_compartments; i++)
	{
		struct compartment *c = &(arch->compartments[i]);
		if (!c->compartment_used)
		{
			continue;
		}

		if (c->update_needed)
		{
			sim_update_compartment(tech, c);
		}
		else
		{
			double *timer = c->time;
			c->energy += tech->energy_inactive_neuron_update;
			(*timer) += tech->time_inactive_neuron_update;
		}
	}
}

int sim_route_spikes(const struct technology *tech, struct architecture *arch)
{
	//const long int max_spikes = tech->max_compartments * tech->fan_out;
	int total_spike_count, total_neurons_fired;

	total_spike_count = 0;
	total_neurons_fired = 0;

	//#pragma omp parallel for

	// now we go through all neurons to see which ones spiked
	//  for all neurons that spike, we send a spike to the corresponding axon output
	//  do we care about the axon output, surely only care about sending packets to other axon_inputs
	//   but we will store the packet tracking struct inside the axon output struct
	//  that axon_output sends to all relevant axon_inputs (we look up the dest neuron which must store a link to that)
	//  the axon_output and input must store the relevant router if it applies i.e. the next link in the chain

	for (int i = 0; i < arch->max_compartments; i++)
	{
		struct compartment *c = &(arch->compartments[i]);
		struct axon_output *axon_out = c->axon_out;

		if (!c->compartment_used || !c->fired)
		{
			continue;
		}

		assert(axon_out != NULL);
		for (int j = 0; j < arch->max_axon_inputs; j++)
		{
			axon_out->packets_sent[j] = 0;
		}
		total_neurons_fired++;

		// Generate all the spikes for a spiking neuron
		//  Only generate one spike packet per core, that the
		//  neuron is broadcasting to
		assert(c->post_connection_count >= 0);
		//TRACE("n:%d post_connection_count:%d\n", n->id,
		//				n->post_connection_count);

		for (int j = 0; j < c->post_connection_count; j++)
		{
			struct compartment *post_neuron;
			struct synapse *synapse_ptr;
			struct axon_input *axon_in;
			struct axon_output *axon_out;

			synapse_ptr = &(c->synapses[j]);
			axon_out = c->axon_out;
			assert(c);
			assert(synapse_ptr);
			if(!synapse_ptr->pre_neuron)
			{
				printf("%d\n", i);
				printf("%d\n", j);
				printf("%d\n", c->id);
				printf("%d\n", j);
				printf("%d\n", c->post_connection_count);
				assert(0);
			}
			assert(synapse_ptr->pre_neuron->id == c->id);

			post_neuron = synapse_ptr->post_neuron;
			axon_in = post_neuron->axon_in;

			post_neuron->update_needed = 1;
			post_neuron->current += synapse_ptr->weight;
			post_neuron->energy += tech->energy_spike_op;
			*(post_neuron->time) += tech->time_spike_op;

			post_neuron->spike_count++;
			total_spike_count++;

			// Mark a packet as sent to the
			//  core. We will only send one packet per core
			axon_out->packets_sent[axon_in->id]++;
		}

		// Estimate the energy needed to send spike packets
		for (int j = 0; j < arch->max_axon_inputs; j++)
		{
			// TODO: I think the code below should be a
			//  function and sim_send_spike can be inlined
			//  again?
			if (axon_out->packets_sent[j])
			{
				struct axon_input *axon_in =
						&(arch->axon_inputs[j]);
				struct router *router_pre = axon_out->r;
				struct router *router_post = axon_in->r;

				// Calculate the energy and time for
				//  sending spike packets
				int x_hops, y_hops;
				double *time = c->time;

				x_hops = abs(router_pre->x - router_post->x);
				y_hops = abs(router_pre->y - router_post->y);
				// E-W hops
				c->energy += x_hops *
					tech->energy_east_west_hop;
				*time += x_hops *
					tech->time_east_west_hop;
				// N-S hops
				c->energy += y_hops *
					tech->energy_north_south_hop;
				*time += y_hops *
					tech->time_north_south_hop;
			}
		}
		c->fired = 0; // Reset the neuron for the next time step
	}
	INFO("Neurons fired: %d\n", total_neurons_fired);

	return total_spike_count;
}

int sim_input_spikes(const struct technology *tech, struct architecture *arch)
{
	int input_spike_count = 0;

	// Seed all externally input spikes in the network for this timestep
	for (int i = 0; i < arch->max_external_inputs; i++)
	{
		struct input *in = &(arch->external_inputs[i]);

		if ((in == NULL) || (!in->send_spike))
		{
			continue;
		}

		for (int j = 0; j < in->post_connection_count; j++)
		{
			// Send a spike to all neurons connected to this input
			//  Normally, we would have a number of input dimensions
			//  for a given network
			struct compartment *post_neuron;
			struct synapse *synapse_ptr;

			synapse_ptr = &(in->synapses[j]);
			assert(synapse_ptr);

			post_neuron = synapse_ptr->post_neuron;
			post_neuron->current += synapse_ptr->weight;
			post_neuron->energy += tech->energy_spike_op;
			*(post_neuron->time) += tech->time_spike_op;
			post_neuron->update_needed = 1;

			post_neuron->spike_count++;
			input_spike_count++;
		}

		// Reset the input ready for the next timestep
		in->send_spike = 0;
	}
	TRACE("Processed %d inputs.\n", arch->max_external_inputs);

	return input_spike_count;
}

void sim_update_compartment(const struct technology *tech,
							struct compartment *c)
{
	// The neuron (state update) contains four main components
	// 1) synapse updates
	// 2) dendrite updates
	// 3) potential (soma) updates
	// 4) Axon updates
	sim_update_synapse_cuba(tech, c);
	sim_update_dendrite(tech, c);
	sim_update_potential(tech, c);
	sim_update_axon(tech, c);

	if (c->spike_count)
	{
		if (c->bias)
		{
			c->energy += tech->energy_active_neuron_update;
			*(c->time) += tech->time_active_neuron_update;
		}
		else
		{
			c->energy += tech->energy_inactive_neuron_update;
			*(c->time) += tech->time_inactive_neuron_update;
		}
	}
}

void sim_update_synapse_cuba(const struct technology *tech,
							struct compartment *c)
{
	// Current based (CUBA) LIF neuron model as implemented by Loihi
	//  Rather than iterate over all synapses we can simplify and just
	//  track total current. This is what nengo-loihi did, I'm not sure if
	//  this would ha ve to be changed if we had a more complicated synapse
	//  model
	TRACE("Current before: %lf.\n", c->current);
	c->current *= c->current_decay;
	TRACE("Current after: %lf.\n", c->current);

	return;
}

void sim_update_dendrite(const struct technology *tech, struct compartment *n)
{
	// TODO: Support dendritic operations, combining the current in
	//  different compartments in some way, and writing the result to
	//  an output compartment
	return;
}

void sim_update_potential(const struct technology *tech, struct compartment *c)
{
	// Calculate the change in potential since the last update e.g.
	//  integate inputs and apply any potential leak
	TRACE("Updating potential, before:%f\n", c->potential);
	c->potential *= c->potential_decay;
	// Add the spike potential
	c->potential += c->current + c->bias;
	// Clamp min potential
	c->potential = (c->potential < c->reset) ?
					c->reset : c->potential;

	TRACE("Updating potential, after:%f\n", c->potential);

	if (c->potential > c->threshold)
	{
		c->fired = 1;
		c->potential = c->reset;
		// Add the "within-tile" spike energy, this is the minimum cost
		//  of sending a spike
		c->energy += tech->energy_spike_within_tile;
		*(c->time) += tech->time_spike_within_tile;
		TRACE("nid %d fired.\n", c->id);
	}

	return;
}

void sim_update_axon(const struct technology *tech, struct compartment *c)
{
	// TODO
	return;
}

double sim_calculate_time(const struct technology *tech,
						struct architecture *arch)
{
	// Returns the simulation time of the current timestep.
	//  This is calculated by finding the simulation time of each core,
	//  and simply taking the maximum of this.
	double max_time;

	max_time = 0.0; // s
	for (int i = 0; i < arch->max_timers; i++)
	{
		//printf("%lf,", arch->timers[i]);
		max_time = fmax(max_time, arch->timers[i]);
	}

	// Add the mesh-wide barrier sync time (assuming worst case of 32 tiles)
	max_time += tech->time_mesh_barrier;
	INFO("Simulated time for step is:%es\n", max_time);

	return max_time;
}

double sim_calculate_energy(struct architecture *arch)
{
	double total_energy = 0.0;

	for (int i = 0; i < arch->max_compartments; i++)
	{
		struct compartment *n = &(arch->compartments[i]);
		if (n->compartment_used)
		{
			total_energy += n->energy;
		}
	}

	return total_energy;
}

void sim_reset_measurements(struct architecture *arch)
{
	// Reset energy and spike count measurements to 0
	for (int i = 0; i < arch->max_compartments; i++)
	{
		struct compartment *c = &(arch->compartments[i]);

		c->energy = 0;
		c->update_needed = c->force_update;
		c->spike_count = 0;
	}

	// Reset all timers
	for (int i = 0; i < arch->max_timers; i++)
	{
		arch->timers[i] = 0.0;
	}
}

void sim_perf_write_header(FILE *fp)
{
	fprintf(fp, "git_version: %s\n", GIT_COMMIT);
}

void sim_log_perf(FILE *fp, const struct architecture *arch)
{
	// Log the energy and time simulated at every compartment, dump it to
	//  a big csv file. Then we can post process it to pull out the parallel
	//  time. Time doesn't make sense per compartment, only per parallel
	//  block. Pull out energy for synapses, routers, axons and compartments
}

void sim_write_summary(FILE *fp, const struct sim_results *results)
{
	// Write the simulation result to file
	fprintf(fp, "energy: %e\n", results->total_energy);
	fprintf(fp, "time: %e\n", results->total_sim_time);
	fprintf(fp, "total_spikes: %ld\n", results->total_spikes);
}

void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp,
						const struct architecture *arch)
{
	// Write csv header for probe outputs - record which neurons have been
	//  probed
	for (int i = 0; i < arch->max_compartments; i++)
	{
		const struct compartment *c = &(arch->compartments[i]);

		if (spike_fp && c->log_spikes)
		{
			fprintf(spike_fp, "%d,", c->id);
		}
		if (potential_fp && c->log_voltage)
		{
			fprintf(potential_fp, "%d,", c->id);
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
				const struct architecture *arch)
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

	for (int i = 0; i < arch->max_compartments; i++)
	{
		const struct compartment *c = &(arch->compartments[i]);

		if (spike_fp && c->log_spikes)
		{
			fprintf(spike_fp, "%d,", c->fired);
			spike_probe_count++;
		}
		if (potential_fp && c->log_voltage)
		{
			fprintf(potential_fp, "%lf,", c->potential);
			potential_probe_count++;
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

int sim_input(const double firing_probability)
{
	// Simulate a single external input (as one neuron) for a timestep
	//  Return 1 if the input fires, 0 otherwise
	double rand_uniform;
	int input_fired;

	rand_uniform = (double) rand() / RAND_MAX;
	input_fired = (rand_uniform < firing_probability);

	return input_fired;
}
