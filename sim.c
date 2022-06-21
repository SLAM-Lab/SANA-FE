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
#include <time.h>
#include <math.h>
#include <omp.h>

#include "sim.h"
#include "tech.h"
#include "network.h"


// TODO: move this out to main
//  Inside this routine is probably where we either generate our spike train
//   or read in a spike train from a file and pass it into the model (assuming
//   we almost always just do this a step at a time). The only thing worth
//   having a wrapper for is timing the step. I guess some cases may not need
//   any inputs at all (just empty inputs vector).
void sim_run(const int timesteps, struct technology *tech, struct core *cores,
			struct sim_results *results, struct input *inputs)
{
	// Run neuromorphic hardware simulation
	for (int i = 0; i < timesteps; i++)
	{
		struct timespec ts_start, ts_end, ts_elapsed;
		struct sim_results timestep_results;

		INFO("*** Time-step %d ***\n", i+1);
		// Measure the wall-clock time taken to run the simulation
		//  on the host machine
		clock_gettime(CLOCK_MONOTONIC, &ts_start);

		timestep_results = sim_timestep(tech, cores, inputs);
		// Accumulate totals for the entire simulation
		// TODO: make a function
		results->total_energy += timestep_results.total_energy;
		results->total_sim_time += timestep_results.total_sim_time;
		results->total_spikes += timestep_results.total_spikes;

		// Calculate elapsed time
		clock_gettime(CLOCK_MONOTONIC, &ts_end);
		ts_elapsed = sim_calculate_elapsed_time(ts_start, ts_end);
		results->wall_time +=
			(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9);
		INFO("Time-step took: %fs.\n",
			(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9));
	}
}


// TODO: make inputs not const and then we can reset them after sending the spike?
struct sim_results sim_timestep(const struct technology *tech,
				struct core *cores, const struct input *inputs)
{
	struct sim_results results;
	long int spikes_sent;

	spikes_sent = sim_input_spikes(tech, cores, inputs);
	sim_reset_measurements(cores, tech->max_cores);
	spikes_sent += sim_route_spikes(tech, cores);

	sim_update_neurons(tech, cores);

	results.total_energy = sim_calculate_energy(cores, tech->max_cores);
	results.total_sim_time = sim_calculate_time(tech, cores);
	INFO("Spikes sent: %ld\n", spikes_sent);
	results.total_spikes = spikes_sent;

	return results;
}

void sim_update_neurons(const struct technology *tech, struct core *cores)
{
	#pragma omp parallel for
	for (int i = 0; i < tech->max_cores; i++)
	{
		struct core *c = &(cores[i]);
		TRACE("Processing %d spikes.\n", c->spike_count);

		// Update all active neurons first
		for (int j = 0; j < c->compartments; j++)
		{
			struct neuron *n = &(c->neurons[j]);
			// Update simulation state, updating synaps current,
			//  neuron potential and TODO dendrite / axon
			//  currents
			sim_update_potential(tech, n, c);
		}

		// TODO: figure what inactive neuron update even means...
		// Uncomment this in some form
		// Update the remainder of the neurons in the core (inactive)
		//for (int j = c->compartments; j < MAX_COMPARTMENTS; j++)
		//{
		//	struct neuron *n = &(c->neurons[j]);
		//	n->energy += tech->energy_inactive_neuron_update;
		//	n->time += tech->time_inactive_neuron_update;
		//}
	}
}

int sim_route_spikes(const struct technology *tech, struct core *cores)
{
	const long int max_spikes = tech->max_compartments * tech->fan_out;
	int total_spike_count, total_neurons_fired;

	total_spike_count = 0;
	total_neurons_fired = 0;

	#pragma omp parallel for
	for (int i = 0; i < tech->max_cores; i++)
	{
		struct core *c = &(cores[i]);

		for (int j = 0; j < c->compartments; j++)
		{
			struct neuron *n;

			for (int k = 0; k < tech->max_cores; k++)
			{
				c->packets_sent[k] = 0;
			}

			n = &(c->neurons[j]);
			if (n->fired == 0)
			{
				continue;
			}

			#pragma omp atomic
			total_neurons_fired++;

			// Generate all the spikes for a spiking neuron
			//  Only generate one spike packet per core, that the
			//  neuron is broadcasting to
			assert(n->post_connection_count >= 0);
			assert(n->post_connection_count <= tech->fan_out);

			for (int k = 0; k < n->post_connection_count; k++)
			{
				struct core *post_core;
				struct neuron *post_neuron;
				struct synapse *synapse_ptr;

				synapse_ptr = &(c->synapses[j][k]);
				assert(n);
				assert(synapse_ptr);
				if(!synapse_ptr->pre_neuron)
				{
					printf("%d\n", i);
					printf("%d\n", j);
					printf("%d\n", n->id);
					printf("%d\n", k);
					printf("%d\n", n->post_connection_count);
					assert(0);
				}
				assert(synapse_ptr->pre_neuron->id == n->id);

				post_neuron = synapse_ptr->post_neuron;
				post_core = &(cores[post_neuron->core_id]);

				#pragma omp atomic
				post_neuron->current += synapse_ptr->weight;
				#pragma omp atomic
				post_core->energy += tech->energy_spike_op;
				#pragma omp atomic
				post_core->time += tech->time_spike_op;

				#pragma omp atomic
				post_core->spike_count++;
				assert(post_core->spike_count <= max_spikes);
				#pragma omp atomic
				total_spike_count++;

				// Mark a packet as sent to the
				//  core. We will only send one packet per core
				#pragma omp atomic
				c->packets_sent[post_core->id]++;
			}

			// Estimate the energy needed to send spike packets
			for (int k = 0; k < tech->max_cores; k++)
			{
				// TODO: I think the code below should be a
				//  function and sim_send_spike can be inlined
				//  again?
				if (c->packets_sent[k])
				{
					// Calculate the energy and time for
					//  sending spike packets
					struct core *post_core = &(cores[k]);
					int x_hops, y_hops;

					x_hops = abs(post_core->x - c->x);
					y_hops = abs(post_core->y - c->y);
					// E-W hops
					#pragma omp atomic
					post_core->energy += x_hops *
						tech->energy_east_west_hop;
					#pragma omp atomic
					post_core->time += x_hops *
						tech->time_east_west_hop;
					// N-S hops
					#pragma omp atomic
					post_core->energy += y_hops *
						tech->energy_north_south_hop;
					#pragma omp atomic
					post_core->time += y_hops *
						tech->time_north_south_hop;
				}
			}

			n->fired = 0; // Reset the neuron for the next time step
		}
	}
	INFO("Neurons fired: %d\n", total_neurons_fired);

	return total_spike_count;
}

//void sim_send_spike(struct synapse *s)
//{
//
//}

int sim_input_spikes(const struct technology *tech, struct core *cores,
			const struct input *inputs)
{
	int input_spike_count = 0;

	// Seed all externally input spikes in the network for this timestep
	for (int i = 0; i < tech->max_inputs; i++)
	{
		const struct input *in = &(inputs[i]);

		if ((in == NULL) || !in->send_spike)
		{
			continue;
		}

		for (int j = 0; j < in->post_connection_count; j++)
		{
			// Send a spike to all neurons connected to this input
			//  Normally, we would have a number of input dimensions
			//  for a given network
			struct core *post_core;
			struct neuron *post_neuron;
			struct synapse *synapse_ptr;

			synapse_ptr = &(in->synapses[j]);
			assert(synapse_ptr);

			post_neuron = synapse_ptr->post_neuron;
			post_core = &(cores[post_neuron->core_id]);

			post_neuron->current += synapse_ptr->weight;
			post_core->energy += tech->energy_spike_op;
			post_core->time += tech->time_spike_op;

			post_core->spike_count++;
			input_spike_count++;
		}
	}
	TRACE("Processed %d inputs.\n", total_inputs);

	return input_spike_count;
}

void sim_update_potential(const struct technology *tech, struct neuron *n,
                          					struct core *c)
{
	// The neuron (state update) contains four main components
	// 1) synapse updates
	// 2) dendrite updates
	// 3) LIF (soma) updates
	// 4) Axon updates
	sim_update_synapse_cuba(tech, c, n);
	sim_update_dendrite(tech, c, n);
	sim_update_lif(tech, c, n);
	sim_update_axon(tech, c, n);

	c->energy += tech->energy_active_neuron_update;
	c->time += tech->time_active_neuron_update;
}

void sim_update_synapse_cuba(const struct technology *tech, struct core *c,
							struct neuron *n)
{
	// Current based (CUBA) LIF neuron model as implemented by Loihi
	//  Rather than iterate over all synapses we can simplify and just
	//  track total current. This is what nengo-loihi did, I'm not sure if
	//  this would have to be changed if we had a more complicated synapse
	//  model
	n->current *= n->current_decay;

	return;
}

void sim_update_dendrite(const struct technology *tech, struct core *c,
							struct neuron *n)
{
	// TODO
	return;
}

void sim_update_lif(const struct technology *tech, struct core *c,
							struct neuron *n)
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

	// TODO: do fire here as well
	if (n->potential > n->threshold)
	{
		n->fired = 1;
		n->potential = n->reset;
		// Add the "within-tile" spike energy, this is
		//  the minimum cost of sending a spike
		c->energy += tech->energy_spike_within_tile;
		c->time += tech->time_spike_within_tile;
		TRACE("nid %d fired.\n", n->id);
	}

	return;
}

void sim_update_axon(const struct technology *tech, struct core *c,
							struct neuron *n)
{
	// TODO
	return;
}

double sim_calculate_time(const struct technology *tech, struct core *cores)
{
	// Returns the simulation time of the current timestep.
	//  This is calculated by finding the simulation time of each core,
	//  and simply taking the maximum of this.  We then add the sync
	//  time, which I can't quite figure out from the Loihi paper but
	//  I'll just add an upper bound..
	double max_time = 0; // s

	for (int i = 0; i < tech->max_cores; i++)
	{
		struct core *c = &(cores[i]);
		max_time = fmax(max_time, c->time);
	}

	// Add the mesh-wide barrier sync time (assuming worst case of 32 tiles)
	max_time += tech->time_mesh_barrier;
	INFO("Simulated time for step is:%es\n", max_time);

	return max_time;
}

double sim_calculate_energy(struct core *cores, const int max_cores)
{
	double total_energy = 0.0;

	for (int i = 0; i < max_cores; i++)
	{
		struct core *c = &(cores[i]);
		total_energy += c->energy;
	}

	return total_energy;
}

void sim_reset_measurements(struct core *cores, const int max_cores)
{
	// Reset time, energy and spike count measurements to 0
	for (int i = 0; i < max_cores; i++)
	{
		struct core *c = &(cores[i]);
		c->spike_count = 0;
		c->energy = 0;
		c->time = 0;
	}
}

struct timespec sim_calculate_elapsed_time(struct timespec ts_start,
							struct timespec ts_end)
{
	// Calculate elapsed wall-clock time between ts_start and ts_end
	struct timespec ts_elapsed;

	ts_elapsed.tv_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	ts_elapsed.tv_sec = ts_end.tv_sec - ts_start.tv_sec;
	if (ts_end.tv_nsec < ts_start.tv_nsec)
	{
		ts_elapsed.tv_sec--;
		ts_elapsed.tv_nsec += 1000000000UL;
	}

	return ts_elapsed;
}

void sim_write_results(FILE *fp, struct sim_results *results)
{
	// Write the simulation result to file
	fprintf(fp, "energy: %e\n", results->total_energy);
	fprintf(fp, "time: %e\n", results->total_sim_time);
	fprintf(fp, "total_spikes: %ld\n", results->total_spikes);
	fprintf(fp, "git_version: %s\n", GIT_COMMIT);
}
