#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "arch.h"
#include "network.h"
#include "command.h"
#include "sim.h"

// Each hardware timestep corresponds to a simulation of the spiking network for
//  dt seconds. This relates to the LIF time constant.
// TODO: This doesn't belong in here, it belongs with the snn stuff
const double dt = 1.0e-3; // Seconds


void command_parse_command(char fields[][MAX_FIELD_LEN], const int field_count,
				struct network *net, struct architecture *arch)
{
	char command_type;

	command_type = fields[0][0];
	// Sanity check input
	if ((command_type == '\0') || (command_type == '\n') ||
							(command_type == '#'))
	{
		INFO("Warning: No command, skipping.\n");
		return;
	}

	// Process the command based on the type
	TRACE("Parsing command type: %c\n", command_type);
	switch (command_type)
	{
	case '\n':
		// Line is a comment
		break;
	case 'g': // Add neuron group
		command_parse_neuron_group(net, fields, field_count);
		break;
	case 'n': // Add neuron
		command_parse_neuron(net, fields, field_count);
		break;
	case '@': // Add network-on-chip (noc)
		command_parse_noc(arch, fields, field_count);
		break;
	case 't':
		command_parse_tile(arch, fields, field_count);
		break;
	case 'c':
		command_parse_core(arch, fields, field_count);
		break;
	/*
	case 's':
		arch_parse_synapse(arch, fields);
		break;
	case 'i':
		arch_parse_axon_input(arch, fields);
		break;
	case '+':
		arch_parse_soma(arch, fields);
		break;
	case 'o':
		arch_parse_axon_output(arch, fields);
		break;
	case 'e':
		//arch_parse_external_input(arch, id_list, line);
		break;
	*/
	case '>':
		// TODO: step simulation
		break;
	default:
		TRACE("Warning: unrecognized unit (%c) - skipping.\n",
							command_type);
		break;
	}
}

void command_parse_line(char *line, char fields[][MAX_FIELD_LEN], struct network *net, struct architecture *arch)
{
	char *token;
	int field_count;

	assert(net != NULL);
	assert(arch != NULL);
	// Zero initialize all fields, each entry is variable length
	//  The neuron fields are fixed, but this is follow by a
	//  variable number of connections (0-MAX_FANOUT).  The synaptic
	//  fields are null delimited
	for (int i = 0; i < MAX_FIELDS; i++)
	{
		fields[i][0] = '\0';
	}


	field_count = 0;
	token = strtok(line, TOKEN_SEPERATORS);
	while (token != NULL)
	{
		// Only parse the neuron fields, the rest of the fields
		//  contain the synaptic data which we'll read next pass
		strncpy(fields[field_count], token, MAX_FIELD_LEN);
		token = strtok(NULL, TOKEN_SEPERATORS);
		field_count++;
	}

	// TODO: max_fields is specific to the thing we're decoding
	/*
	if (field_count < MAX_FIELDS)
	{
		TRACE("nid:%d Number of fields read < %d, "
			"ignoring line.\n", neuron_count,
							NEURON_FIELDS);
		continue;
	}
	
	connection_count = (field_count - NEURON_FIELDS) / CONNECTION_FIELDS;
	assert(connection_count >= 0);
	*/

	command_parse_command(fields, field_count, net, arch);
}

void command_read_file(FILE *fp, struct network *net,
						struct architecture *arch)
{
	char *line;
	char (*fields)[MAX_FIELD_LEN];

	assert(net != NULL);
	assert(arch != NULL);
	
	fields = (char (*)[MAX_FIELD_LEN])
			malloc(sizeof(char[MAX_FIELD_LEN]) * MAX_FIELDS);
	line = (char *) malloc(sizeof(char) * MAX_CSV_LINE);
	if ((line == NULL) || (fields == NULL))
	{
		INFO("Error: Couldn't allocate memory for text input.\n");
		exit(1);
	}
	// Allocate memory to hold the string fields for each field
	// TODO: two rounds of malloc are excessive, we can do this with one
	/*
	for (int i = 0; i < MAX_FIELDS; i++)
	{
		fields[i] = (char *) malloc(MAX_CSV_LINE * sizeof(char));
		if (fields[i] == NULL)
		{
			INFO("Error: Failed to allocate memory for text.\n");
			exit(1);
		}
	}
	*/
		
	while (fgets(line, MAX_CSV_LINE, fp))
	{
		// Convert line to fields
		command_parse_line(line, fields, net, arch);
	}
}

int command_parse_tile(struct architecture *const arch,
			char fields[][MAX_FIELD_LEN], const int field_count)
{
	TRACE("Parsing tile.\n");
	return arch_create_tile(arch);
}

int command_parse_noc(struct architecture *arch, char fields[][MAX_FIELD_LEN],
							const int field_count)
{
	// Here we can define the links between tiles, all the stuff like whether
	//  it's a mesh, tree can be done externally. Here we just apply raw
	//  links. Maybe we can even remove NoC links dynamically. So we parse
	//  whether to create (1) or remove (0) a link between two tiles, or
	//  something like this
	/*
	struct tile *t;
	int id, ret;

	ret = sscanf(fields[1], "%u", &id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output id (%s).\n", fields[1]);
		exit(1);
	}
	assert(id < arch->tile_count);
	t = &(arch->tiles[id]);
	t->id = id;

	// TODO: change to single scanf (everywhere, not just in this function)
	ret = sscanf(fields[2], "%u", &(t->x));
	ret += sscanf(fields[3], "%u", &(t->y));
	if (ret < 2)
	{
		INFO("Error: couldn't parse router r:%d.\n", id);
		exit(1);
	}
	TRACE("Router parsed r:%u x(%d) y(%d).\n", t->id, t->x, t->y);
	*/
	return 0;
}

int command_parse_core(struct architecture *arch, char fields[][MAX_FIELD_LEN],
							const int field_count)
{
	unsigned int tile_id;
	int ret;

	TRACE("Parsing core.\n");
	if (field_count < 2)
	{
		INFO("Error: Invalid <add core> command; not enough fields.\n");
		exit(1);
	}
	ret = sscanf(fields[1], "%u", &tile_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse tile id (%s).\n", fields[1]);
		exit(1);
	}
	return arch_create_core(arch, tile_id);
}

int command_parse_neuron_group(struct network *const net,
						char fields[][MAX_FIELD_LEN],
							const int field_count)
{
	double threshold, reset;
	unsigned int neuron_count;
	int ret;

	if (field_count < 4)
	{
		INFO("Error: Invalid <add neuron group> command; "
							"not enough fields.\n");
		exit(1);
	}
	ret = sscanf(fields[1], "%u", &neuron_count);
	ret += sscanf(fields[1], "%lf", &threshold);
	ret += sscanf(fields[1], "%lf", &reset);
	
	if (ret < 3)
	{
		INFO("Error: Couldn't parse command\n");
		exit(1);
	}

	return network_create_neuron_group(net, neuron_count, threshold, reset);
}

int command_parse_neuron(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count)
{
	return 0;
}



/*
void arch_parse_axon_input(struct architecture *arch,
					char fields[][ARCH_MAX_FIELD_LEN])
{
	struct axon_input *axon_in;
	int ret, id, router_id;
	char type;

	assert(fields && fields[0]);
	type = fields[0][0];
	assert(type == 'i');

	ret = sscanf(fields[1], "%u", &id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon input id (%s).\n", fields[1]);
		exit(1);
	}
	assert(id < arch->max_axon_inputs);
	axon_in = &(arch->axon_inputs[id]);
	axon_in->id = id;

	ret = sscanf(fields[2], "%u", &router_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon router (%s).\n", fields[2]);
		exit(1);
	}
	assert(router_id < arch->max_routers);
	axon_in->r = &(arch->routers[router_id]);

	ret = sscanf(fields[3], "%d", &axon_in->fan_in);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon input fan in (%s).\n",
								fields[3]);
		exit(1);
	}

	TRACE("Axon input parsed i:%d r(%d) fan_in(%d)\n", axon_in->id,
					axon_in->r->id, axon_in->fan_in);
}

void arch_parse_axon_output(struct architecture *arch,
					char fields[][ARCH_MAX_FIELD_LEN])
{
	struct axon_output *axon_out;
	int ret, id, router_id;
	char type;

	assert(fields && fields[0]);
	type = fields[0][0];
	assert(type == 'o');

	ret = sscanf(fields[1], "%d", &id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output id (%s).\n", fields[1]);
		exit(1);
	}
	assert(id < arch->max_axon_outputs);
	axon_out = &(arch->axon_outputs[id]);
	axon_out->id = id;

	ret = sscanf(fields[2], "%d", &router_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon router (%s).\n", fields[2]);
		exit(1);
	}
	assert(router_id < arch->max_routers);
	axon_out->r = &(arch->routers[router_id]);

	ret = sscanf(fields[3], "%d", &axon_out->fan_out);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output fan out (%s).\n",
								fields[3]);
		exit(1);
	}

	TRACE("Axon output parsed o:%d r(%d) fan_out(%d)\n", axon_out->id,
					axon_out->r->id, axon_out->fan_out);
}

*/


// TODO: this is probably going to be deleted soon
/*
void command_read_csv(FILE *fp, struct network *net, struct architecture *arch)
{
	// Build arbitrary spiking network from a csv file
	//
	// SNN is defined in a single csv file, with one row per neuron.
	//  To build the network, make two passes over the file.
	//  1) Create all the neurons in the network, map these to hardware
	//  2) Add the connections and connect neurons to form a network
	//  We can't connect the neurons until we have made an initial pass
	//  to define them all.  This is why a two pass approach was needed.
	//
	// The intention is probably to have the csv machine generated.
	//  There are a number of neuron fields, followed by up to FAN_OUT
	//  connections each with a smaller number of fields.
	//
	// See network.h to see all the fields and what they mean
	struct neuron_group *group;
	struct input *input_ptr;
	struct connection *con;
	char *line;
	int neuron_count, input_count, neuron_id, ret;
	int compartment_id, dest_id, curr_input, is_input;

	// TODO: need to allocate memory for all the neurons
	// I complicated it a bit by having groups. We need to know how many
	//  neurons are in each group before allocating.

	network_init(net);


	neuron_count = 0;
	input_count = 0;

	// Step 1 - Create all neurons and map these to the hardware. Initialize
	//  the neuron and count the number of connections to allocate.
	while (fgets(line, MAX_CSV_LINE, fp))
	{
		struct neuron *c;
		char *token;
		int connection_count, field_count;

		if (neuron_count > max_neurons)
		{
			INFO("Error: inputting too many neurons, max is %d.\n",
							max_neurons);
			exit(1);
		}

		// Zero initialize all fields, each entry is variable length
		//  The neuron fields are fixed, but this is follow by a
		//  variable number of connections (0-MAX_FANOUT).  The synaptic
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
							NETWORK_MAX_FIELD_LEN);
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
		connection_count = (field_count - NEURON_FIELDS) / CONNECTION_FIELDS;
		assert(connection_count >= 0);

		for (int i = 0; i < field_count; i++)
		{
			TRACE("nid:%d Parsed field: %s\n", neuron_count,
							neuron_fields[i]);
		}

		is_input = (neuron_fields[NEURON_ID][0] == 'i');
		if (is_input)
		{
			struct input *input_ptr =
					&(net->external_inputs[input_count]);

			TRACE("Creating network input %d.\n", input_count);
			assert(input_count < net->external_input_count);

			if (connection_count)
			{
				input_ptr->connections = (struct connection *)
					malloc(connection_count *
						sizeof(struct connection));
				if (input_ptr->connections == NULL)
				{
					INFO("Error: Couldn't allocate"
								"connections.\n");
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
		c = &(arch->neurons[compartment_id]);
		assert(c->neuron_used == 0);

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
		c->neuron_used = 1;

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
		if (connection_count)
		{
			c->connections = (struct connection *) malloc(connection_count *
							sizeof(struct connection));
			if (c->connections == NULL)
			{
				INFO("Error: couldn't allocate connection mem "
							"for nid:%u.\n", c->id);
				exit(1);
			}
		}

		// Keep track of which CSV line corresponds to which physical
		//  neuron in the hardware.  This info will be important for
		//  linking the connection data to neuron neurons
		neuron_ptrs[neuron_count] = c;

		TRACE("Added nid:%d vt:%lf r%lf log_s:%d "
			"log_v:%d\n", c->id, c->threshold, c->reset,
			c->log_spikes, c->log_voltage);
		neuron_count++;
	}
	INFO("Created %d inputs.\n", input_count);
	INFO("Created %d neurons.\n", neuron_count);

	curr_input = 0;
	// Next parse the whole file again, but this time read the connection data
	fseek(fp, 0, SEEK_SET);
	while (fgets(line, NETWORK_MAX_CSV_LINE, fp))
	{
		struct neuron *src, *dest;
		char *token;
		int field_count, connection_count;

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
			//  only interested in the connection data
			strncpy(neuron_fields[field_count], token,
							NETWORK_MAX_FIELD_LEN);
			token = strtok(NULL, ",");
			field_count++;
		}
		if (field_count < NEURON_FIELDS)
		{
			TRACE("Number of fields read < %d, "
				"ignoring line.\n", NEURON_FIELDS);
			continue;
		}
		connection_count =
			(field_count - NEURON_FIELDS) / CONNECTION_FIELDS;

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
			input_ptr = &(net->external_inputs[curr_input]);
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

		for (int i = 0; i < connection_count; i++)
		{
			float weight;
			const int curr_connection_field = NEURON_FIELDS +
							(i*CONNECTION_FIELDS);

			// Parse a single connection from the csv
			ret = sscanf(neuron_fields[curr_connection_field +
							CONNECTION_DEST_NID],
							    "%d", &dest_id);
			if (ret <= 0)
			{
				INFO("Error: Couldn't parse connection \"%s\".\n",
					neuron_fields[curr_connection_field +
							CONNECTION_DEST_NID]);
				exit(1);
			}
			else if ((dest_id < 0) || (dest_id >= neuron_count))
			{
				INFO("Error: connection dest neuron (%d) out of "
					"range [0 <= #neuron < %d].\n", dest_id,
								neuron_count);
				exit(1);
			}

			sscanf(neuron_fields[curr_connection_field +
								CONNECTION_WEIGHT],
								"%f", &weight);
			dest = neuron_ptrs[dest_id];
			// Create the new connection and add it to the end out
			//  the fan-out list core
			if (is_input)
			{
				con = &(input_ptr->connections[
					input_ptr->post_connection_count]);
				con->pre_neuron = NULL;
				con->post_neuron = dest;
				con->weight = weight;
				input_ptr->post_connection_count++;
				TRACE("Created input connection i->%d (w:%f)\n",
					con->post_neuron->id, con->weight);
			}
			else
			{
				con = &(src->connections[
						src->post_connection_count]);
				con->pre_neuron = src;
				con->post_neuron = dest;
				con->weight = weight;
				src->post_connection_count++;
				TRACE("Created connection %d->%d (w:%f)\n",
					con->pre_neuron->id, con->post_neuron->id,
					con->weight);
			}
		}
		if (is_input)
		{
			TRACE("nid:i Added %d connections.\n",
			input_ptr->post_connection_count);
		}
		else
		{
			TRACE("nid:%d Added %d connections.\n",
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
*/

/*
struct architecture command_arch_init(FILE *fp)
{
	// Parse the architecture file, counting the number of units needed
	//  and then allocating the memory
	struct architecture arch;
	char line[ARCH_LINE_LEN];
	char *field;

	arch.tiles = NULL;
	//arch.max_neurons = 0;
	//arch.max_mem_blocks = 0;
	//arch.max_routers = 0;
	//arch.max_timers = 0;
	//arch.max_axon_inputs = 0;
	//arch.max_axon_outputs = 0;
	//arch.max_external_inputs = 131072; // HACK

	// TODO: refactor - arch_count_units(fp) ?
	while (fgets(line, ARCH_LINE_LEN, fp))
	{
		unsigned int unit_count;
		char block_type = line[0];

		field = strtok(line, " ");
		field = strtok(NULL, " ");
		unit_count = arch_get_count(field);
		TRACE("unit count:%d.\n", unit_count);

		// TODO: can refactor so the switch sets an enum?
		//  That doesn't really add anything though since we'd need
		//  another switch
		switch (block_type)
		{
		case '\0':
		case '#':
		case '\n':
			// Line is a comment
			continue;
		case 't':
			arch.tile_count += unit_count;
			break;
		case 'i':
			arch.max_axon_inputs += unit_count;
			break;
		case 'o':
			arch.max_axon_outputs += unit_count;
			break;
		case 'e':
			arch.max_external_inputs += unit_count;
			break;
		case 'd':
			break;
		case 's':
			break;
		default:
			TRACE("Warning: unrecognized unit (%c) - skipping.\n",
								block_type);
			break;
		}
	}

	INFO("Parsed %d neurons.\n", arch.max_neurons);

	// Reset file pointer to start of file
	fseek(fp, 0, SEEK_SET);

	// Based on the number of different units, allocate enough memory to
	//  simulate this design
	INFO("Allocating memory for %d neurons.\n", arch.max_neurons);
	*/
	/*
	arch.neurons = (struct neuron *)
		malloc(arch.max_neurons * sizeof(struct neuron));
	arch.mem_blocks = (struct mem *)
		malloc(arch.max_mem_blocks * sizeof(struct mem));
	arch.routers = (struct router *)
		malloc(arch.max_routers * sizeof(struct router));
	arch.axon_inputs = (struct axon_input *)
		malloc(arch.max_axon_inputs * sizeof(struct axon_input));
	arch.axon_outputs = (struct axon_output *)
		malloc(arch.max_axon_outputs * sizeof(struct axon_output));
	arch.timers = (double *) malloc(arch.max_timers * sizeof(double));

	if ((arch.neurons == NULL) || (arch.mem_blocks == NULL) ||
		(arch.routers == NULL) || (arch.timers == NULL) ||
		(arch.axon_inputs == NULL) || (arch.axon_outputs == NULL))
	{
		INFO("Error: Failed to allocate neuron.\n");
		exit(1);
	}
	*/

	/*
	for (int i = 0; i < arch.max_neurons; i++)
	{
		struct neuron *n = &(arch.neurons[i]);

		n->id = i;
		n->synapses = NULL;
		n->axon_in = NULL;
		n->axon_out = NULL;

		n->fired = 0;
		n->update_needed = 0;
		n->neuron_used = 0;
	}

	for (int i = 0; i < arch.max_mem_blocks; i++)
	{
		struct mem *mem = &(arch.mem_blocks[i]);

		mem->id = i;
		TRACE("Allocating synapse memory for block %d.\n", mem->id);
	}
	*/
/*
	for (int i = 0; i < arch.max_axon_outputs; i++)
	{
		struct axon_output *axon_out = &(arch.axon_outputs[i]);

		axon_out->packets_sent = (unsigned int *)
			malloc(arch.max_axon_inputs * sizeof(unsigned int));
		if (axon_out->packets_sent == NULL)
		{
			INFO("Error: Failed to allocate axon memory.\n");
			exit(1);
		}
	}

	INFO("Allocating memory for %d external inputs.\n",
						arch.max_external_inputs);
	arch.external_inputs = (struct input *)
		malloc(arch.max_external_inputs * sizeof(struct input));
	if (arch.external_inputs == NULL)
	{
		INFO("Error: Failed to allocate input memory.\n");
		exit(1);
	}

	// Zero initialize the input nodes
	for (int i = 0; i < arch.max_external_inputs; i++)
	{
		struct input *in = &(arch.external_inputs[i]);

		in->synapses = NULL;
		in->post_connection_count = 0;
		in->send_spike = 0;
	}

	for (int i = 0; i < arch.max_timers; i++)
	{
		arch.timers[i] = 0.0;
	}

	arch.initialized = 1;

	return arch;
}
*/

/* Probably don't need this code, let's see
// TODO: I don't know if this is a wasted effort. We could get away with only
//  supporting ranges, and the only place ranges are really useful is for
//  neuron neurons (of which there are many). This is a lot of complicated
//  code and structures for debatable gain. Maybe get rid of this
static int arch_parse_list(const char *field, struct range *list)
{
	char value[ARCH_MAX_VALUE_DIGITS];
	int is_range, list_entries, v, curr_char;

	list_entries = 0;
	curr_char = 0;
	is_range = 0;
	v = 0;
	for (curr_char = 0; curr_char < ARCH_MAX_FIELD_LEN; curr_char++)
	{
		// Parse a field, which is a comma-separated list containing
		// 1) values e.g. (3)
		// 2) value ranges e.g. (3..7)
		if ((field[curr_char] == ',') || (field[curr_char] == '\0'))
		{
			struct range *curr_range = &(list[list_entries]);

			value[v] = '\0';
			// Parse the previous number
			if (is_range)
			{
				int ret = sscanf(value, "%u..%u",
					&(curr_range->min), &(curr_range->max));
				if (ret < 2)
				{
					INFO("Error: Parsing range (%s).\n",
									field);
					exit(1);
				}
				TRACE("Parsed range of values (%u to %u)\n",
					curr_range->min, curr_range->max);
			}
			else
			{
				int ret = sscanf(value, "%u",
							&curr_range->min);
				if (ret < 1)
				{
					INFO("Error: Parsing value (%s).\n",
									field);
					exit(1);
				}
				// Just one value is read i.e. min == max
				curr_range->max = curr_range->min;
				TRACE("Parsed value (%u)\n", curr_range->min);
			}

			v = 0;
			is_range = 0;
			list_entries++;
			if (field[curr_char] == '\0')
			{
				// End of string, stop processing
				break;
			}
		}
		else // copy string into temporary buffer for parsing
		{
			value[v] = field[curr_char];
			if ((field[curr_char] == '.') &&
						(field[curr_char+1] == '.'))
			{
				// This entry is a range of values e.g. 1..3
				is_range = 1;
			}
			v++;
		}
	}

	return list_entries;
}

static unsigned int arch_get_count(const char *field)
{
	const char *curr = field;
	unsigned int values;
	// Count how many values the field has

	// A field can either be a number, a comma separated list,
	//  or a range of values. For simplicity, this format doesn't allow you
	//  to combine lists and ranges
	values = 0;
	while (values < ARCH_MAX_VALUES)
	{
		// Either count the number of commas for lists, but if a range
		//  parse the left and right hand values. Otherwise it's just a
		//  single value.
		if (curr == NULL)
		{
			break;
		}
		if ((curr[0] == '\n') || (curr[0] == '\0'))
		{
			values++;
			break;
		}
		else if (isdigit(curr[0]))
		{
			curr++;
		}
		else if (curr[0] == ',')
		{
			values++;
			curr++;
		}
		else if ((curr[0] == '.') && (curr[1] == '.'))
		{
			int ret;
			unsigned int lval, rval;

			ret = sscanf(field, "%d..%d", &lval, &rval);
			if (ret)
			{
				// Calculate the range (inclusive), e.g.
				// 1..3 -> 1,2,3 which is 3 values
				values = (rval - lval) + 1;
				TRACE("Parsed range (%d..%d).\n", lval, rval);
			}
			else
			{
				values = 0;
				TRACE("Warning: invalid range (%s).\n", curr);
			}
			break;
		}
		else
		{
			TRACE("Warning: invalid field (%s).\n", field);
			break;
		}
	}

	TRACE("Parsed %d values.\n", values);
	return values;
}
*/