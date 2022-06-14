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

void sim_run(const int timesteps, struct technology *tech, struct core *cores,
                				struct sim_results *results)
{
	// Run neuromorphic hardware simulation
	for (int i = 0; i < timesteps; i++)
	{
		struct timespec ts_start, ts_end, ts_elapsed;

		INFO("*** Time-step %d ***\n", i+1);
		// Measure the wall-clock time taken to run the simulation
		//  on the host machine
		clock_gettime(CLOCK_MONOTONIC, &ts_start);

		sim_timestep(tech, results, cores);

		// Calculate elapsed time
		clock_gettime(CLOCK_MONOTONIC, &ts_end);
		ts_elapsed = sim_calculate_elapsed_time(ts_start, ts_end);
		results->wall_time +=
			(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9);
		INFO("Time-step took: %fs.\n",
			(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9));
	}
}

void sim_timestep(const struct technology *tech, struct sim_results *results,
							struct core *cores)
{
	long int spikes_sent;

	sim_reset_measurements(cores, tech->max_cores);
	sim_seed_input_spikes(cores, tech->max_cores);

	spikes_sent = sim_route_spikes(tech, cores);
	results->total_spikes += spikes_sent;
	INFO("Spikes sent: %ld\n", spikes_sent);
	sim_update_neurons(tech, cores);

	results->total_energy += sim_calculate_energy(cores, tech->max_cores);
	results->total_sim_time += sim_calculate_time(tech, cores);
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
			// Add synaptic current and apply leak to all neurons
			sim_update_potential(tech, n, c);
		}

		// Update the remainder of the neurons in the core (inactive)
		//for (int j = c->compartments; j < MAX_COMPARTMENTS; j++)
		//{
		//	struct neuron *n = &(c->neurons[j]);
		//	n->energy += tech->energy_inactive_neuron_update;
		//	n->time += tech->time_inactive_neuron_update;
		//}

		// Check to see which neurons have fired
		for (int j = 0; j < c->compartments; j++)
		{
			struct neuron *n = &(c->neurons[j]);
			if (n->potential > n->threshold)
			{
				n->fired = 1;
				n->potential = n->reset;
				// Add the "within-tile" spike energy, this is
				//  the minimum cost of sending a spike
				c->energy += tech->energy_spike_within_tile;
				c->time += tech->time_spike_within_tile;
			}
		}
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
		struct neuron *n, *post_neuron;
		struct core *post_core, *c;
		struct synapse *synapse_ptr;

		// TODO: this is kind of horrible, maybe we should enforce a
		//  MAX_CORES define, and the technology file must be <= this
		//  limit. To avoid quite as much dynamic allocation
		unsigned int *spike_packets = (unsigned int *)
				malloc(tech->max_cores * sizeof(unsigned int));
		if (spike_packets == NULL)
		{
			INFO("Error: Failed to allocate memory.\n");
			exit(1);
		}

		c = &(cores[i]);
		for (int j = 0; j < c->compartments; j++)
		{
			for (int k = 0; k < tech->max_cores; k++)
			{
				spike_packets[k] = 0;
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
				spike_packets[post_core->id]++;
			}

			// Estimate the energy needed to send spike packets
			for (int k = 0; k < tech->max_cores; k++)
			{
				int x_hops, y_hops;

				// TODO: I think the code below should be a
				//  function and sim_send_spike can be inlined
				//  again?
				if (spike_packets[k])
				{
					// Calculate the energy and time for
					//  sending spike packets
					post_core = &(cores[k]);
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
		free(spike_packets);
	}
	INFO("Total neurons firing: %d\n", total_neurons_fired);

	return total_spike_count;
}

void sim_send_spike(struct synapse *s)
{

}

void sim_seed_input_spikes(struct core *cores, const int max_cores)
{
	// Seed all externally input spikes in the network for this timestep
	#pragma omp parallel for
	for (int i = 0; i < max_cores; i++)
	{
		struct core *c = &(cores[i]);

		for (int j = 0; j < c->compartments; j++)
		{
			struct neuron *n = &(c->neurons[j]);
			int is_input = (n->input_rate > 0);

			if (is_input)
			{
				n->fired |= sim_input(n->input_rate);
			}
		}
	}
}

void sim_update_potential(const struct technology *tech, struct neuron *n,
                          					struct core *c)
{
	// Current based (CUBA) LIF neuron model as implemented by Loihi
	n->current *= n->current_decay;

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

	c->energy += tech->energy_active_neuron_update;
	c->time += tech->time_active_neuron_update;
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
