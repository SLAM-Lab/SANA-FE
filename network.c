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
			struct architecture *arch, struct neuron **neuron_ptrs)
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

	const int max_fields = NEURON_FIELDS + (4096 * SYNAPSE_FIELDS);
	const int max_neurons = arch->max_neurons;

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
		struct neuron *n;
		char *token;
		int synapse_count, field_count;

		if (neuron_count >= max_neurons)
		{
			INFO("Error: inputting too many neurons, max is %d.\n",
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
		synapse_count = field_count - NEURON_FIELDS;
		assert(synapse_count >= 0);

		for (int i = 0; i < field_count; i++)
		{
			TRACE("nid:%d Parsed field: %s\n", neuron_count,
							neuron_fields[i]);
		}

		is_input = (neuron_fields[NEURON_ID][0] == 'i');
		if (is_input)
		{
			TRACE("Creating network input %d.\n", input_count);
			input_count++;
			assert(input_count <= arch->max_external_inputs);
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
		n = &(arch->neurons[compartment_id]);

		// Parse related neuron parameters
		sscanf(neuron_fields[THRESHOLD_VOLTAGE], "%lf",
							&(n->threshold));
		sscanf(neuron_fields[RESET_VOLTAGE], "%lf", &(n->reset));
		sscanf(neuron_fields[RECORD_SPIKES], "%d", &(n->log_spikes));
		sscanf(neuron_fields[RECORD_VOLTAGE], "%d", &(n->log_voltage));
		n->fired = 0;
		n->potential = n->reset;
		n->active = 1;

		// Hard coded LIF / CUBA time constants for now
		// TODO: parameterize based on network description (csv)
		n->current_time_const = 1.0e-3;
		n->potential_time_const = 2.0e-3;
		n->current_decay = -(exp(-dt / n->current_time_const) - 1.0);
		n->potential_decay =
				-(exp(-dt / n->potential_time_const) - 1.0);

		// Allocate synaptic memory
		// TODO: use the memory block to model the space available for
		//  synaptic memory
		//assert(n->synapses == NULL);
		INFO("synapse_count: %d\n", synapse_count);
		if (synapse_count)
		{
			n->synapses = (struct synapse *) malloc(synapse_count *
							sizeof(struct synapse));
			if (n->synapses == NULL)
			{
				INFO("Error: couldn't allocate synapse mem "
							"for nid:%u.\n", n->id);
				exit(1);
			}
		}

		// Keep track of which CSV line corresponds to which physical
		//  neuron in the hardware.  This info will be important for
		//  linking the synapse data to neuron compartments
		neuron_ptrs[neuron_count] = n;

		TRACE("Added nid:%d vt:%lf r%lf log_s:%d "
                        "log_v:%d\n", n->id, n->threshold, n->reset,
			n->log_spikes, n->log_voltage);
		neuron_count++;
	}
	INFO("Created %d neurons.\n", neuron_count);

	curr_input = 0;
	// Next parse the whole file again, but this time read the synapse data
	fseek(fp, 0, SEEK_SET);
	while (fgets(line, MAX_CSV_LINE, fp))
	{
		struct neuron *src, *dest;
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
		synapse_count = field_count - NEURON_FIELDS;

		// Use the first field (the neuron number) to figure if this
		//  is a valid formatted line or not. Since this is the
		//  second pass we know this field is valid
		is_input = (neuron_fields[NEURON_ID][0] == 'i');
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
			src = neuron_ptrs[neuron_id];
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
			//  the fan-out list core
			if (is_input)
			{
				s = &(input_ptr->synapses[
					input_ptr->post_connection_count]);
				s->pre_neuron = NULL;
				s->post_neuron = dest;
				input_ptr->post_connection_count++;
			}
			else
			{
				s = &(src->synapses[
						src->post_connection_count]);
				s->pre_neuron = src;
				s->post_neuron = dest;
				src->post_connection_count++;
			}
			s->weight = weight;
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

void network_init(const struct technology *tech, struct architecture *arch)
{
	for (int i = 0; i < arch->max_neurons; i++)
	{
		struct neuron *n = &(arch->neurons[i]);

		n->id = i;
		n->post_connection_count = 0;
		n->log_spikes = 0;
		n->log_voltage = 0;

		n->fired = 0;
		n->potential = 0.0;
		n->current = 0.0;
		n->bias = 0.0;
		n->threshold = 0.0;
		n->reset = 0.0;
		n->potential_decay = 0.0;
		n->current_decay = 0.0;
		n->potential_time_const = 0.0;
		n->current_time_const = 0.0;
	}

	for (int i = 0; i < arch->max_routers; i++)
	{
		struct router *r = &(arch->routers[i]);

		// Loihi is organised in a 2D mesh of 32 tiles (8x4)
		//  Groups of 4 cores share a router, forming a tile
		// TODO: generalize (i.e. where does 8 come from)
		r->x = i % 8;
		r->y = i / 8;
	}
}

void network_create_empty(const struct technology *tech,
						struct architecture *arch)
{
	// Initialize an empty network, where cores have no connections between
	//  neurons
	for (int i = 0; i < arch->max_neurons; i++)
	{
		struct neuron *n = &(arch->neurons[i]);

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

