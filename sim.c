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
	sim_reset_measurements(arch);
	spikes_sent += sim_route_spikes(tech, arch);
	sim_update_state(tech, arch);

	// Return results for this timestep
	results.total_energy = sim_calculate_energy(arch);
	results.total_sim_time = sim_calculate_time(tech, arch);
	results.total_spikes = spikes_sent;
	INFO("Spikes sent: %ld\n", spikes_sent);

	sim_probe_log_timestep(probe_spike_fp, probe_potential_fp, arch);
	return results;
}

void sim_update_state(const struct technology *tech,
						struct architecture *arch)
{
	//#pragma omp parallel for
	for (int i = 0; i < arch->max_neurons; i++)
	{
		//TRACE("Processing %d spikes.\n", c->spike_count);
		struct neuron *n = &(arch->neurons[i]);

		if (n->active)
		{
			sim_update_active(tech, n);
		}
		else
		{
#if 0
			sim_update_inactive(tech, n);
#endif
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

	for (int i = 0; i < arch->max_neurons; i++)
	{
		struct neuron *n = &(arch->neurons[i]);
		struct axon_output *axon_out = n->axon_out;

		if (!n->active || !n->fired)
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
		assert(n->post_connection_count >= 0);
		//TRACE("n:%d post_connection_count:%d\n", n->id,
		//				n->post_connection_count);

		for (int j = 0; j < n->post_connection_count; j++)
		{
			struct neuron *post_neuron;
			struct synapse *synapse_ptr;
			struct axon_input *axon_in;
			struct axon_output *axon_out;

			synapse_ptr = &(n->synapses[j]);
			axon_out = n->axon_out;
			assert(n);
			assert(synapse_ptr);
			if(!synapse_ptr->pre_neuron)
			{
				printf("%d\n", i);
				printf("%d\n", j);
				printf("%d\n", n->id);
				printf("%d\n", j);
				printf("%d\n", n->post_connection_count);
				assert(0);
			}
			assert(synapse_ptr->pre_neuron->id == n->id);

			post_neuron = synapse_ptr->post_neuron;
			axon_in = post_neuron->axon_in;

			post_neuron->current += synapse_ptr->weight;
			post_neuron->energy += tech->energy_spike_op;
			post_neuron->time += tech->time_spike_op;

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

				x_hops = abs(router_pre->x - router_post->x);
				y_hops = abs(router_pre->y - router_post->y);
				// E-W hops
				n->energy += x_hops *
					tech->energy_east_west_hop;
				n->time += x_hops *
					tech->time_east_west_hop;
				// N-S hops
				n->energy += y_hops *
					tech->energy_north_south_hop;
				n->time += y_hops *
					tech->time_north_south_hop;
			}
		}
		n->fired = 0; // Reset the neuron for the next time step
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
			struct neuron *post_neuron;
			struct synapse *synapse_ptr;

			synapse_ptr = &(in->synapses[j]);
			assert(synapse_ptr);

			post_neuron = synapse_ptr->post_neuron;
			post_neuron->current += synapse_ptr->weight;
			post_neuron->energy += tech->energy_spike_op;
			post_neuron->time += tech->time_spike_op;

			//post_core->spike_count++;
			input_spike_count++;
		}

		// Reset the input ready for the next timestep
		in->send_spike = 0;
	}
	TRACE("Processed %d inputs.\n", arch->max_external_inputs);

	return input_spike_count;
}

void sim_update_active(const struct technology *tech, struct neuron *n)
{
	// The neuron (state update) contains four main components
	// 1) synapse updates
	// 2) dendrite updates
	// 3) LIF (soma) updates
	// 4) Axon updates
	sim_update_synapse_cuba(tech, n);
	sim_update_dendrite(tech, n);
	sim_update_lif(tech, n);
	sim_update_axon(tech, n);

	n->energy += tech->energy_active_neuron_update;
	n->time += tech->time_active_neuron_update;
}

void sim_update_inactive(const struct technology *tech,
						struct neuron *n)
{
	// TODO: figure what inactive neuron update even means...
	n->energy += tech->energy_inactive_neuron_update;
	n->time += tech->time_inactive_neuron_update;
}

void sim_update_synapse_cuba(const struct technology *tech, struct neuron *n)
{
	// Current based (CUBA) LIF neuron model as implemented by Loihi
	//  Rather than iterate over all synapses we can simplify and just
	//  track total current. This is what nengo-loihi did, I'm not sure if
	//  this would ha ve to be changed if we had a more complicated synapse
	//  model
	TRACE("Current before: %lf.\n", n->current);
	n->current *= n->current_decay;
	TRACE("Current after: %lf.\n", n->current);

	return;
}

void sim_update_dendrite(const struct technology *tech, struct neuron *n)
{
	// TODO
	return;
}

void sim_update_lif(const struct technology *tech, struct neuron *n)
{
	// Calculate the decay in potential since the last update i.e. the
	//  leak
	TRACE("Updating potential, before:%f\n", n->potential);
	n->potential *= n->potential_decay;
	// Add the spike potential
	n->potential += n->current + n->bias;
	// Clamp min potential
	n->potential = (n->potential < n->reset) ?
					n->reset : n->potential;

	TRACE("Updating potential, after:%f\n", n->potential);

	if (n->potential > n->threshold)
	{
		n->fired = 1;
		n->potential = n->reset;
		// Add the "within-tile" spike energy, this is
		//  the minimum cost of sending a spike
		n->energy += tech->energy_spike_within_tile;
		n->time += tech->time_spike_within_tile;
		TRACE("nid %d fired.\n", n->id);
	}

	return;
}

void sim_update_axon(const struct technology *tech, struct neuron *n)
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
	// TODO: adapt to general arch
	//  Use the timing tree to calculate, based on parallel, time multiplexed
	//  and globally synchronized blocks
	double max_time, core_time;

	max_time = 0.0; // s
	core_time = 0.0;

	// TODO: remove this hack, it's just here to make Loihi sims work until
	//  I have a more general timing structure
	for (int i = 0; i < arch->max_neurons; i++)
	{
		struct neuron *n = &(arch->neurons[i]);
		if ((i % 1024) == 0)
		{
			max_time = fmax(max_time, core_time);
			core_time = 0.0;
		}
		core_time += n->time;
	}
	max_time = fmax(max_time, core_time);

	// Add the mesh-wide barrier sync time (assuming worst case of 32 tiles)
	max_time += tech->time_mesh_barrier;
	INFO("Simulated time for step is:%es\n", max_time);

	return max_time;
}

double sim_calculate_energy(struct architecture *arch)
{
	double total_energy = 0.0;

	for (int i = 0; i < arch->max_neurons; i++)
	{
		struct neuron *n = &(arch->neurons[i]);
		total_energy += n->energy;
	}

	return total_energy;
}

void sim_reset_measurements(struct architecture *arch)
{
	// Reset time, energy and spike count measurements to 0
	for (int i = 0; i < arch->max_neurons; i++)
	{
		struct neuron *n = &(arch->neurons[i]);

		//n->spike_count = 0;
		n->energy = 0;
		n->time = 0;
	}
}

void sim_write_results(FILE *fp, const struct sim_results *results)
{
	// Write the simulation result to file
	fprintf(fp, "energy: %e\n", results->total_energy);
	fprintf(fp, "time: %e\n", results->total_sim_time);
	fprintf(fp, "total_spikes: %ld\n", results->total_spikes);
	fprintf(fp, "git_version: %s\n", GIT_COMMIT);
}

void sim_probe_write_header(FILE *spike_fp, FILE *potential_fp,
						const struct architecture *arch)
{
	// Write csv header for probe outputs - record which neurons have been
	//  probed
	for (int i = 0; i < arch->max_neurons; i++)
	{
		const struct neuron *n = &(arch->neurons[i]);

		if (spike_fp && n->log_spikes)
		{
			fprintf(spike_fp, "%d,", n->id);
		}
		if (potential_fp && n->log_voltage)
		{
			fprintf(potential_fp, "%d,", n->id);
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

	for (int i = 0; i < arch->max_neurons; i++)
	{
		const struct neuron *n = &(arch->neurons[i]);

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
