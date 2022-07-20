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

void network_read_csv(FILE *fp, const struct technology *tech,
			struct architecture *arch, struct compartment **compartment_ptrs)
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
	struct input *input_ptr;
	struct synapse *s;
	char *line;
        int neuron_count, input_count, compartment_id, ret;
	int neuron_id, dest_id, curr_input, is_input;

	// TODO: better define the max fields
	const int max_fields = NEURON_FIELDS + (4096 * SYNAPSE_FIELDS);
	const int max_compartments = arch->max_compartments;

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
		neuron_fields[i] =
				(char *) malloc(MAX_FIELD_LEN * sizeof(char));
		if (neuron_fields[i] == NULL)
		{
			INFO("Error: Failed to allocate memory for text.\n");
			exit(1);
		}
	}

	network_init(tech, arch);
	line = (char *) malloc(sizeof(char) * MAX_CSV_LINE);
	if (line == NULL)
	{
		INFO("Error: Couldn't allocate memory for text input.\n");
		exit(1);
	}

	neuron_count = 0;
	input_count = 0;

	// Step 1 - Create all neurons and map these to the hardware. Initialize
	//  the neuron and count the number of synapses to allocate.
	while (fgets(line, MAX_CSV_LINE, fp))
	{
		struct compartment *c;
		char *token;
		int synapse_count, field_count;

		if (neuron_count > max_compartments)
		{
			INFO("Error: inputting too many neurons, max is %d.\n",
							max_compartments);
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
		while (token != NULL)
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
		synapse_count = (field_count - NEURON_FIELDS) / SYNAPSE_FIELDS;
		assert(synapse_count >= 0);

		for (int i = 0; i < field_count; i++)
		{
			TRACE("nid:%d Parsed field: %s\n", neuron_count,
							neuron_fields[i]);
		}

		is_input = (neuron_fields[NEURON_ID][0] == 'i');
		if (is_input)
		{
			struct input *input_ptr =
					&(arch->external_inputs[input_count]);

			TRACE("Creating network input %d.\n", input_count);
			assert(input_count < arch->max_external_inputs);

			if (synapse_count)
			{
				input_ptr->synapses = (struct synapse *)
					malloc(synapse_count *
							sizeof(struct synapse));
				if (input_ptr->synapses == NULL)
				{
					INFO("Error: Couldn't allocate"
								"synapses.\n");
					exit(1);
				}
			}
			input_count++;
			continue;
		}
		else // This line is defining a neuron
		{
			ret = sscanf(neuron_fields[NEURON_ID], "%d",
								&neuron_id);
		}

		// Parse the neuron id
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
		}
		else if (neuron_id != neuron_count)
		{
			INFO("Error: #line (%d) != #neurons (%d).\n",
					neuron_id, neuron_count);
			exit(1);
		}

		sscanf(neuron_fields[COMPARTMENT_ID], "%d", &compartment_id);
		// Add neuron to specified neuromorphic core
		c = &(arch->compartments[compartment_id]);
		assert(c->compartment_used == 0);

		// Parse related neuron parameters
		sscanf(neuron_fields[THRESHOLD_VOLTAGE], "%lf",
							&(c->threshold));
		sscanf(neuron_fields[RESET_VOLTAGE], "%lf", &(c->reset));
		sscanf(neuron_fields[RECORD_SPIKES], "%d",
							&(c->log_spikes));
		sscanf(neuron_fields[RECORD_VOLTAGE], "%d",
							&(c->log_voltage));
		sscanf(neuron_fields[FORCE_UPDATE], "%d",
							&(c->force_update));
		c->fired = 0;
		c->potential = c->reset;
		// TODO: force update is really referring to whether there's
		//  a bias or not, I figured out
		c->update_needed = c->force_update;
		c->compartment_used = 1;

		// TODO: Hard coded LIF / CUBA time constants for now
		c->current_time_const = 1.0e-3;
		c->potential_time_const = 2.0e-3;
		c->current_decay =
			-(exp(-dt / c->current_time_const) - 1.0);
		c->potential_decay =
			-(exp(-dt / c->potential_time_const) - 1.0);

		// Allocate synaptic memory
		// TODO: use the memory block to model the space available for
		//  synaptic memory
		if (synapse_count)
		{
			c->synapses = (struct synapse *) malloc(synapse_count *
							sizeof(struct synapse));
			if (c->synapses == NULL)
			{
				INFO("Error: couldn't allocate synapse mem "
							"for nid:%u.\n", c->id);
				exit(1);
			}
		}

		// Keep track of which CSV line corresponds to which physical
		//  neuron in the hardware.  This info will be important for
		//  linking the synapse data to neuron compartments
		compartment_ptrs[neuron_count] = c;

		TRACE("Added nid:%d vt:%lf r%lf log_s:%d "
                        "log_v:%d\n", c->id, c->threshold, c->reset,
			c->log_spikes, c->log_voltage);
		neuron_count++;
	}
	INFO("Created %d inputs.\n", input_count);
	INFO("Created %d neurons.\n", neuron_count);

	curr_input = 0;
	// Next parse the whole file again, but this time read the synapse data
	fseek(fp, 0, SEEK_SET);
	while (fgets(line, MAX_CSV_LINE, fp))
	{
		struct compartment *src, *dest;
		char *token;
		int field_count, synapse_count;

		for (int i = 0; i < max_fields; i++)
		{
			neuron_fields[i][0] = '\0';
		}

		// Read all csv fields into a buffer
		field_count = 0;
		token = strtok(line, ",");
		while (token != NULL)
		{
			// This time read all the fields in the line, but we're
			//  only interested in the synapse data
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
		synapse_count = (field_count - NEURON_FIELDS) / SYNAPSE_FIELDS;

		// Use the first field (the neuron number) to figure if this
		//  is a valid formatted line or not. Since this is the
		//  second pass we know this field is valid
		is_input = (neuron_fields[NEURON_ID][0] == 'i');
		input_ptr = NULL;
		if (is_input)
		{
			// An input is a virtual connnection - it isn't
			//  associated with a neuron on the chip, but is
			//  connected to other neurons on the chip
			src = NULL;
			input_ptr = &(arch->external_inputs[curr_input]);
			TRACE("Parsing network input: %d.\n", curr_input);
			curr_input++;
		}
		else
		{
			ret = sscanf(neuron_fields[NEURON_ID], "%d",
								&neuron_id);
			if (ret <= 0)
			{
				// Couldn't parse the neuron id field
				TRACE("Header detected, skipping.\n");
				continue;
			}
			// Now parse all the outgoing synaptic connections for
			//  this neuron
			src = compartment_ptrs[neuron_id];
			for (int i = 0; i < field_count; i++)
			{
				TRACE("nid:%d Parsed field: %s\n", neuron_id,
							neuron_fields[i]);
			}
		}

		for (int i = 0; i < synapse_count; i++)
		{
			float weight;
			const int curr_synapse_field = NEURON_FIELDS +
							(i*SYNAPSE_FIELDS);

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
			dest = compartment_ptrs[dest_id];
			// Create the new synapse and add it to the end out
			//  the fan-out list core
			if (is_input)
			{
				s = &(input_ptr->synapses[
					input_ptr->post_connection_count]);
				s->pre_neuron = NULL;
				s->post_neuron = dest;
				s->weight = weight;
				input_ptr->post_connection_count++;
				TRACE("Created input synapse i->%d (w:%f)\n",
					s->post_neuron->id, s->weight);
			}
			else
			{
				s = &(src->synapses[
						src->post_connection_count]);
				s->pre_neuron = src;
				s->post_neuron = dest;
				s->weight = weight;
				src->post_connection_count++;
				TRACE("Created synapse %d->%d (w:%f)\n",
					s->pre_neuron->id, s->post_neuron->id,
					s->weight);
			}
			s->energy = 0.0;
		}
		if (is_input)
		{
			TRACE("nid:i Added %d synapses.\n",
			input_ptr->post_connection_count);
		}
		else
		{
			TRACE("nid:%d Added %d synapses.\n",
			src->id,
			src->post_connection_count);
		}
	}
	INFO("Initialized neurons and inputs.\n");

	for (int i = 0; i < max_fields; i++)
	{
		free(neuron_fields[i]);
	}
	free(neuron_fields);
	free(line);
	return;
}

void network_init(const struct technology *tech, struct architecture *arch)
{
	for (int i = 0; i < arch->max_compartments; i++)
	{
		struct compartment *c = &(arch->compartments[i]);

		c->id = i;
		c->post_connection_count = 0;
		c->log_spikes = 0;
		c->log_voltage = 0;
		c->spike_count = 0;

		c->fired = 0;
		c->potential = 0.0;
		c->current = 0.0;
		c->bias = 0.0;
		c->threshold = 0.0;
		c->reset = 0.0;
		c->update_needed = 0;
		c->compartment_used = 0;
	}

	for (int i = 0; i < arch->max_routers; i++)
	{
		struct router *r = &(arch->routers[i]);

		// TODO:
		r->east = NULL;
		r->west = NULL;
		r->north = NULL;
		r->south = NULL;
		r->energy = 0.0;
	}
}

void network_create_empty(const struct technology *tech,
						struct architecture *arch)
{
	// Initialize an empty network, where cores have no connections between
	//  neurons
	for (int i = 0; i < arch->max_compartments; i++)
	{
		struct compartment *c = &(arch->compartments[i]);

		c->fired = 0;
		c->post_connection_count = 0;
		c->reset = 0.0;
		c->potential = c->reset;
		c->current = 0.0;

		c->threshold = -1.0; // Negative threshold always fires
		c->bias = 0.0;
	}
}

