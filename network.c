// network.c
// Utility functions for creating user defined spiking networks
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "network.h"
#include "sim.h"

// Each hardware timestep corresponds to a simulation of the spiking network for
//  dt seconds. This relates to the LIF time constant.
const double dt = 1.0e-3; // Seconds

void network_read_csv(FILE *fp, struct neuron **neuron_ptrs, struct core *cores,
			const int max_cores, const struct technology *tech)
{
	// Build arbitrary spiking network from a csv file
	//
	// SNN is defined in a single csv file, with one row per neuron.
	//  To build the network, make two passes over the file.
	//  1) Create all the neurons in the network, map these to hardware
	//  2) Add the synapses and connect neurons to form a network
	//  We can't connect the neurons until we have made an initial pass
	//  to define them all.  This is why a two pass approach was needed.
	//
	// The intention is probably to have the csv machine generated.
	//  There are a number of neuron fields, followed by up to FAN_OUT
	//  synapses each with a smaller number of fields.
	//
	// See network.h to see all the fields and what they mean
	struct core *c;
	struct neuron *n, *src, *dest;
	struct synapse *s;
	char *token, *line;
	float weight;
        int neuron_count, core_id, field_count, ret, neuron_id, dest_id;

	const int max_fields = NEURON_FIELDS + (tech->fan_out * SYNAPSE_FIELDS);
	const int max_neurons = tech->max_compartments * tech->max_cores;

	// Allocate memory to hold the string fields for each neuron description
	//  In the file, one line describes one neuron
	char **neuron_fields = (char **) malloc(max_fields * sizeof(char *));
	if (neuron_fields == NULL)
	{
		INFO("Error: Failed to allocate memory for network inputs.\n");
		exit(1);
	}
	for (int i = 0; i < max_fields; i++)
	{
		neuron_fields[i] = (char *) malloc(MAX_FIELD_LEN * sizeof(char));
		if (neuron_fields[i] == NULL)
		{
			INFO("Error: Failed to allocate memory of network inputs.\n");
			exit(1);
		}
	}

	network_init(cores, max_cores, tech);
	line = (char *) malloc(sizeof(char) * MAX_CSV_LINE);
	if (line == NULL)
	{
		INFO("Error: Couldn't allocate memory for network inputs.\n");
		exit(1);
	}

	neuron_count = 0;

	while (fgets(line, MAX_CSV_LINE, fp))
	{
		if (neuron_count >= max_neurons)
		{
			INFO("Error: inputting too many neurons, max is %d",
								max_neurons);
			exit(1);
		}

		// Zero initialize all fields, each entry is variable length
		//  The neuron fields are fixed, but this is follow by a
		//  variable number of synapses (0-MAX_FANOUT).  The synaptic
		//  fields are null delimited
		for (int i = 0; i < max_fields; i++)
		{
			neuron_fields[i][0] = '\0';
		}

		field_count = 0;
		token = strtok(line, ",");
		while ((token != NULL) && (field_count < NEURON_FIELDS))
		{
			// Only parse the neuron fields, the rest of the fields
			//  contain the synaptic data which we'll read next pass
			strncpy(neuron_fields[field_count], token,
								MAX_FIELD_LEN);
			token = strtok(NULL, ",");
			field_count++;
		}
		if (field_count < NEURON_FIELDS)
		{
			TRACE("nid:%d Number of fields read < %d, "
				"ignoring line.\n", neuron_count,
								NEURON_FIELDS);
			continue;
		}

		for (int i = 0; i < field_count; i++)
		{
			TRACE("nid:%d Parsed field: %s\n", neuron_count,
							neuron_fields[i]);
		}

		ret = sscanf(neuron_fields[NEURON_ID], "%d", &neuron_id);
		// The first field in the CSV is mostly just a label for
		//  readability.  It won't have any effect on the input
		if (ret <= 0)
		{
			// Couldn't parse the neuron id field
			if (neuron_count == 0)
			{
				// Some csv files might have a header on the
				//  first line(s), skip these
				TRACE("Header detected, skipping.\n");
				continue;
			}
			else
			{
				INFO("Error: couldn't parse %s.\n",
						neuron_fields[NEURON_ID]);
				exit(1);
			}
			continue;
		}
		else if (neuron_id != neuron_count)
		{
			INFO("Error: #line (%d) != #neurons (%d).\n",
						neuron_id, neuron_count);
			exit(1);
		}

		sscanf(neuron_fields[CORE_ID], "%d", &core_id);
		c = &(cores[core_id]);
		assert(c->id == core_id);

		// Add neuron to specified neuromorphic core
		n = &(c->neurons[c->compartments]);

		c->compartments++;
		if (c->compartments > tech->max_compartments)
		{
			INFO("Error: For core %d, #compartments (%d) > %d.\n",
				c->id, c->compartments, tech->max_compartments);
		}

		// Parse related neuron parameters
		sscanf(neuron_fields[THRESHOLD_VOLTAGE], "%lf",
							&(n->threshold));
		sscanf(neuron_fields[RESET_VOLTAGE], "%lf", &(n->reset));
		sscanf(neuron_fields[INPUT_RATE], "%lf", &(n->input_rate));
		sscanf(neuron_fields[RECORD_SPIKES], "%d", &(n->log_spikes));
		sscanf(neuron_fields[RECORD_VOLTAGE], "%d", &(n->log_voltage));

		// Keep track of which CSV line corresponds to which physical
		//  neuron in the hardware.  This info will be important for
		//  linking the synapse data to neuron compartments
		neuron_ptrs[neuron_count] = n;

		TRACE("Added nid:%d cid:%d vt:%lf r%lf in:%d log_s:%d "
                        "log_v:%d\n", n->id, n->core_id, n->threshold, n->reset,
			n->is_input, n->log_spikes, n->log_voltage);
		neuron_count++;
	}
	INFO("Created %d neurons.\n", neuron_count);

	// Next parse the whole file again, but this time read the synapse data
	fseek(fp, 0, SEEK_SET);
	while (fgets(line, MAX_CSV_LINE, fp))
	{
		for (int i = 0; i < max_fields; i++)
		{
			neuron_fields[i][0] = '\0';
		}

		// Read all csv fields into a buffer
		field_count = 0;
		token = strtok(line, ",");
		while (token != NULL)
		{
			// This time read all the fields in the line, we're
			//  interested in the synapse data
			strncpy(neuron_fields[field_count], token,
								MAX_FIELD_LEN);
			token = strtok(NULL, ",");
			field_count++;
		}
		if (field_count < NEURON_FIELDS)
		{
			TRACE("Number of fields read < %d, "
				"ignoring line.\n", NEURON_FIELDS);
			continue;
		}

		// Use the first field (the neuron number) to figure if this
		//  is a valid formatted line or not. Since this is the
		//  second pass we know this field is valid
		ret = sscanf(neuron_fields[NEURON_ID], "%d", &neuron_id);
		if (ret <= 0)
		{
			// Couldn't parse the neuron id field
			TRACE("Header detected, skipping.\n");
			continue;
		}

		// Now parse all the outgoing synaptic connections for this
		//  neuron
		src = neuron_ptrs[neuron_id];
		c = &(cores[src->core_id]);

		for (int i = 0; i < field_count; i++)
		{
			TRACE("nid:%ld Parsed field: %s\n", neuron_id,
							neuron_fields[i]);
		}

		for (int i = 0; i < tech->fan_out; i++)
		{
			const int curr_synapse_field = NEURON_FIELDS +
							(i*SYNAPSE_FIELDS);
			// The list of synapse input strings is null terminted
			if (neuron_fields[curr_synapse_field][0] == '\0')
			{
				TRACE("nid:%d Added %d synapses.\n",
					src->id, src->post_connection_count);
				break;
			}

			// Parse a single synapse from the csv
			ret = sscanf(neuron_fields[curr_synapse_field +
							SYNAPSE_DEST_NID],
							    "%d", &dest_id);
			if (ret <= 0)
			{
				INFO("Error: Couldn't parse synapse \"%s\".\n",
					neuron_fields[curr_synapse_field +
							SYNAPSE_DEST_NID]);
				exit(1);
			}
			else if ((dest_id < 0) || (dest_id >= neuron_count))
			{
				INFO("Error: synapse dest neuron (%d) out of "
					"range [0 <= #neuron < %d].\n", dest_id,
								neuron_count);
				exit(1);
			}

			sscanf(neuron_fields[curr_synapse_field +
								SYNAPSE_WEIGHT],
								"%f", &weight);
			dest = neuron_ptrs[dest_id];
			// Create the new synapse and add it to the end out
			//  the fan-out list
			//  core
			s = &(c->synapses[src->compartment]
						[src->post_connection_count]);
			s->pre_neuron = src;
			s->post_neuron = dest;
			s->weight = weight;
			// Now log the extra connection in the presynaptic
			//  neuron
			src->post_connection_count++;
			TRACE("Created synapse %d->%d (w:%f)\n",
				s->pre_neuron->id, s->post_neuron->id,
				s->weight);
		}
	}

	for (int i = 0; i < max_fields; i++)
	{
		free(neuron_fields[i]);
	}
	free(neuron_fields);
	free(line);
	return;
}

void network_init(struct core *cores, const int max_cores,
						const struct technology *tech)
{
	// Initialize state of neuromorphic cores and all their compartments
	INFO("Initializing %d cores.\n", max_cores);
	for (int i = 0; i < max_cores; i++)
	{
		struct core *c = &(cores[i]);

		c->id = i;
		c->spike_count = 0;
		c->compartments = 0;

		// Loihi is organised in a 2D mesh of 32 tiles (8x4)
		//  Groups of 4 cores share a router, forming a tile
		//  Routers are then connected in the mesh, allowing multi-hop
		c->x = (i / tech->cores_per_tile) % 8;
		c->y = (i / tech->cores_per_tile) / 8;

		for (int j = 0; j < tech->max_compartments; j++)
		{
			struct neuron *n;
			int neuron_id;

			n = &(c->neurons[j]);
			n->core_id = c->id;

			// Simply assign neurons in ascending order from core 0
			neuron_id = (i * tech->max_compartments) + j;
			n->id = neuron_id;
			n->compartment = j;

			n->post_connection_count = 0;
                        n->input_rate = 0;
                        n->log_spikes = 0;
                        n->log_voltage = 0;
		}
	}
}

void network_create_empty(struct core *cores, const struct technology *tech)
{
	// Initialize an empty network, where cores have no connections between
	//  neurons
	INFO("Creating empty network with %d cores.\n", tech->max_cores);

	// Initialise network with no connections
	for (int i = 0; i < tech->max_cores; i++)
	{
		struct core *c = &(cores[i]);
		for (int j = 0; j < c->compartments; j++)
		{
			struct neuron *n = &(c->neurons[j]);

			n->fired = 0;
			n->post_connection_count = 0;
			n->reset = 0.0;
			n->potential = n->reset;
			n->current = 0.0;

			// No connections so these parameters don't matter
			n->current_time_const = 1.0e-3;
			n->potential_time_const = 2.0e-3;
			n->current_decay =
				-(exp(-dt / n->current_time_const) - 1.0);
			n->potential_decay =
				-(exp(-dt / n->potential_time_const) - 1.0);

			n->threshold = -1.0; // Negative threshold always fires
			n->bias = 0.0;
		}
	}
}

