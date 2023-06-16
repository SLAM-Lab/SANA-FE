// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// command.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "print.h"
#include "arch.h"
#include "network.h"
#include "command.h"
#include "sim.h"

int command_parse_command(char fields[][MAX_FIELD_LEN],
				const int field_count, struct network *net,
				struct architecture *arch,
				FILE *probe_spikes_fp, FILE *potential_trace_fp,
				FILE *perf_fp)
{
	int ret;
	char command_type;

	command_type = fields[0][0];
	// Sanity check input
	if ((command_type == '\0') || (command_type == '\n') ||
							(command_type == '#'))
	{
		INFO("Warning: No command, skipping.\n");
		return COMMAND_OK;
	}

	// Process the command based on the type
	switch (command_type)
	{
	case 'g': // Add neuron group
		ret = command_parse_neuron_group(net, fields, field_count);
		break;
	case 'n': // Add neuron
		ret = command_parse_neuron(net, arch, fields, field_count);
		break;
	case '@': // Add network-on-chip connections (noc)
		ret = command_parse_noc(arch, fields, field_count);
		break;
	case 't': // Add tile
		ret = command_parse_tile(arch, fields, field_count);
		break;
	case 'c': // Add core
		ret = command_parse_core(arch, fields, field_count);
		break;
	case '&': // Map neuron to hardware
		ret = command_map_hardware(net, arch, fields, field_count);
		break;
	case 'i': // Add axon input processor
		ret = command_parse_axon_input(arch, fields, field_count);
		break;
	case 's': // Add synapse processor
		ret = command_parse_synapse(arch, fields, field_count);
		break;
	case '+': // Add soma processor
		ret = command_parse_soma(arch, fields, field_count);
		break;
	case 'o': // Add axon output processor
		ret = command_parse_axon_output(arch, fields, field_count);
		break;
	case 'e': // Add group of external inputs
		ret = command_parse_extern_input_group(net, fields, field_count);
		break;
	case '<': // Add single input node
		ret = command_parse_extern_input_node(net, fields, field_count);
		break;
	case '$': // Input vector (either rates or individual spikes)
		ret = command_parse_input_spikes(net, fields, field_count);
		break;
	case '*': // Step simulator
		ret = command_parse_step_sim(net, arch, fields, field_count,
						probe_spikes_fp,
						potential_trace_fp, perf_fp);
		break;
	case 'l':  // Load commands from file
		ret = command_parse_load_commands(net, arch, fields,
							field_count,
							probe_spikes_fp,
							potential_trace_fp,
							perf_fp);
		break;
	default:
		TRACE("Warning: unrecognized unit (%c) - skipping.\n",
							command_type);
		ret = COMMAND_OK;
		break;
	}

	return ret;
}

int command_parse_line(char *line, char fields[][MAX_FIELD_LEN],
			struct network *net, struct architecture *arch,
			FILE *probe_spikes_fp, FILE *potential_trace_fp,
			FILE *perf_fp)
{
	char *token;
	int field_count, last_char_idx;

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
	last_char_idx = strlen(line)-1;
	assert(line[last_char_idx] == '\n');
	line[last_char_idx] = '\0';
	token = strtok(line, TOKEN_SEPERATORS);
	while (token != NULL)
	{
		// Only parse the neuron fields, the rest of the fields
		//  contain the synaptic data which we'll read next pass
		strncpy(fields[field_count], token, MAX_FIELD_LEN);
		token = strtok(NULL, TOKEN_SEPERATORS);
		field_count++;
	}

	return command_parse_command(fields, field_count, net, arch,
					probe_spikes_fp, potential_trace_fp,
								perf_fp);
}

int command_parse_file(FILE *fp, struct network *net,
			struct architecture *arch, FILE *probe_spikes_fp,
			FILE *potential_trace_fp, FILE *perf_fp)
{
	char *line;
	char (*fields)[MAX_FIELD_LEN];
	int ret;

	assert(net != NULL);
	assert(arch != NULL);

	fields = (char (*)[MAX_FIELD_LEN])
			malloc(sizeof(char[MAX_FIELD_LEN]) * MAX_FIELDS);
	line = (char *) malloc(sizeof(char) * MAX_CSV_LINE);
	ret = COMMAND_OK;

	if ((line == NULL) || (fields == NULL))
	{
		INFO("Error: Couldn't allocate memory for text input.\n");
		ret = COMMAND_FAIL;
	}
	else
	{
		while (fgets(line, MAX_CSV_LINE, fp) && (ret != COMMAND_FAIL))
		{
			ret = command_parse_line(line, fields, net, arch,
							probe_spikes_fp,
							potential_trace_fp,
							perf_fp);
		}
	}

	free(line);
	free(fields);

	return ret;
}

int command_parse_noc(struct architecture *arch, char fields[][MAX_FIELD_LEN],
							const int field_count)
{
	int dimensions, width, height, tile_count, ret;

	// Here we can define the links between tiles, all the stuff like whether
	//  it's a mesh, tree can be done externally. Here we just apply raw
	//  links. Maybe we can even remove NoC links dynamically. So we parse
	//  whether to create (1) or remove (0) a link between two tiles, or
	//  something like this
	INFO("Parsing NoC.\n");
	if (strncmp("mesh", fields[1], strlen("mesh")) != 0)
	{
		INFO("Error: Only 'mesh' NoC type supported (%s).\n",
								fields[1]);
		return COMMAND_FAIL;
	}

	ret = sscanf(fields[2], "%d", &dimensions);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse axon output id (%s).\n", fields[1]);
		exit(1);
	}

	if (dimensions != 2)
	{
		INFO("Error: Only 2 dimensions supported (%d).\n", dimensions);
	}

	ret = sscanf(fields[3], "%d", &width);
	ret += sscanf(fields[4], "%d", &height);
	if (ret < 2)
	{
		INFO("Error: Couldn't parse dimensions (%s).\n", fields[1]);
		exit(1);
	}

	tile_count = width * height;
	if (tile_count != arch->tile_count)
	{
		INFO("Error: width*height (%d) != tile count (%d)\n",
						tile_count, arch->tile_count);
	}

	return arch_create_noc(arch, width, height);
}

int command_parse_tile(struct architecture *const arch,
			char fields[][MAX_FIELD_LEN], const int field_count)
{
	double energy_east_west_hop, energy_north_south_hop;
	double time_east_west_hop, time_north_south_hop;
	int ret;

	TRACE("Parsing tile.\n");

	ret = sscanf(fields[3], "%lf", &energy_east_west_hop);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse energy (%s).\n", fields[3]);
		return COMMAND_FAIL;
	}

	ret = sscanf(fields[4], "%lf", &time_north_south_hop);

	if (ret < 1)
	{
		INFO("Error: Couldn't parse time (%s).\n", fields[4]);
		return COMMAND_FAIL;
	}

	ret = sscanf(fields[5], "%lf", &energy_north_south_hop);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse energy (%s).\n", fields[5]);
		return COMMAND_FAIL;
	}

	ret = sscanf(fields[6], "%lf", &time_east_west_hop);

	if (ret < 1)
	{
		INFO("Error: Couldn't parse energy (%s).\n", fields[6]);
		return COMMAND_FAIL;
	}

	return arch_create_tile(arch,
				energy_east_west_hop, energy_north_south_hop,
				time_east_west_hop, time_north_south_hop);
}

int command_parse_core(struct architecture *arch,
			char fields[][MAX_FIELD_LEN],
			const int field_count)
{
	struct tile *t;
	int tile_id;
	int ret;

	TRACE("Parsing core.\n");
	if (field_count < 2)
	{
		INFO("Error: Invalid <add core> command; not enough fields.\n");
		return COMMAND_FAIL;
	}
	ret = sscanf(fields[1], "%d", &tile_id);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse tile id (%s).\n", fields[1]);
		return COMMAND_FAIL;
	}

	if (tile_id >= arch->tile_count)
	{
		INFO("Error: Invalid tile id: %d.\n", tile_id);
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);

	return arch_create_core(arch, t);
}

int command_parse_reset_mode(const char *str)
{
	int reset_mode = -1;

	if (strcmp(str, "none") == 0)
	{
		reset_mode = NEURON_NO_RESET;
	}
	else if (strcmp(str, "soft") == 0)
	{
		reset_mode = NEURON_RESET_SOFT;
	}
	else if (strcmp(str, "hard") == 0)
	{
		reset_mode = NEURON_RESET_HARD;
	}
	else if (strcmp(str, "saturate") == 0)
	{
		reset_mode = NEURON_RESET_SATURATE;
	}
	else
	{
		INFO("Error: reset mode not recognized.");
	}

	return reset_mode;
}

int command_parse_neuron_group(struct network *const net,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	double threshold, reset, reverse_threshold, reverse_reset, leak;
	int ret, neuron_count, reset_mode, reverse_reset_mode;

	if (field_count < 8)
	{
		INFO("Error: Invalid <add neuron group> command; "
							"not enough fields.\n");
		exit(1);
	}
	ret = sscanf(fields[1], "%d", &neuron_count);
	ret += sscanf(fields[2], "%lf", &threshold);
	ret += sscanf(fields[3], "%lf", &reset);
	ret += sscanf(fields[4], "%lf", &reverse_threshold);
	ret += sscanf(fields[5], "%lf", &reverse_reset);
	ret += sscanf(fields[6], "%lf", &leak);
	reset_mode = command_parse_reset_mode(fields[7]);
	reverse_reset_mode = command_parse_reset_mode(fields[8]);

	if (ret < 6)
	{
		INFO("Error: Couldn't parse command.\n");
		return COMMAND_FAIL;
	}
	else
	{
		TRACE("Creating neuron group with %d neurons.\n", neuron_count);
		return network_create_neuron_group(net, neuron_count, threshold,
			reset, reverse_threshold, reverse_reset, leak,
						reset_mode, reverse_reset_mode);
	}
}

int command_parse_neuron(struct network *const net, struct architecture *arch,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	struct neuron_group *group;
	struct neuron *n;
	double bias;
	int neuron_group_id, neuron_id, connection_count;
	int log_spikes, log_potential, force_update, ret;

	TRACE("Parsing neuron.\n");
	if (field_count < NEURON_FIELDS)
	{
		INFO("Error: Invalid <neuron> command; (%d/%d) fields.\n",
						field_count, NEURON_FIELDS);
		return COMMAND_FAIL;
	}
	ret = sscanf(fields[NEURON_GROUP_ID], "%d", &neuron_group_id);
	ret += sscanf(fields[NEURON_ID], "%d", &neuron_id);
	ret += sscanf(fields[NEURON_BIAS], "%lf", &bias);
	ret += sscanf(fields[NEURON_RECORD_SPIKES], "%d", &log_spikes);
	ret += sscanf(fields[NEURON_RECORD_VOLTAGE], "%d", &log_potential);
	ret += sscanf(fields[NEURON_FORCE_UPDATE], "%d", &force_update);

	if (ret < (NEURON_FIELDS-1))
	{
		INFO("Error: Couldn't parse neuron command.\n");
		return COMMAND_FAIL;
	}

	if (neuron_group_id >= net->neuron_group_count)
	{
		INFO("Error: Group (%d) >= group count (%d).\n",
				neuron_group_id, net->neuron_group_count);
		return COMMAND_FAIL;
	}
	group = &(net->groups[neuron_group_id]);

	connection_count = (field_count - NEURON_FIELDS) / CONNECTION_FIELDS;
	TRACE("Parsed neuron gid:%d nid:%d log s:%d log v:%d force:%d "
		"connections:%d\n", neuron_group_id, neuron_id, log_spikes,
				log_potential, force_update, connection_count);



	if (neuron_id >= group->neuron_count)
	{
		INFO("Error: Neuron (%d) >= group neurons (%d).\n",
					neuron_id, group->neuron_count);
		return COMMAND_FAIL;
	}
	n = &(group->neurons[neuron_id]);
	ret = network_create_neuron(n, bias, log_spikes, log_potential,
				force_update, connection_count);
	if (ret == NETWORK_INVALID_NID)
	{
		return COMMAND_FAIL;
	}

	TRACE("nid:%d creating %d connections.\n", n->id, connection_count);
	// Parse each synapse between neurons, those neurons must already exit
	for (int i = 0; i < connection_count; i++)
	{
		// Parse each connection
		struct neuron_group *dest_group;
		struct neuron *src, *dest;
		struct connection *con;
		double weight;
		int dest_gid, dest_nid;
		const int curr_field = NEURON_FIELDS + (i*CONNECTION_FIELDS);
		assert(curr_field < field_count);
		con = &(n->connections_out[i]);

		ret = sscanf(fields[curr_field+CONNECTION_DEST_GID],
							"%d", &dest_gid);
		ret += sscanf(fields[curr_field+CONNECTION_DEST_NID],
							"%d", &dest_nid);
		TRACE("group id:%s\n", fields[curr_field+CONNECTION_DEST_GID]);
		TRACE("neuron id:%s\n", fields[curr_field+CONNECTION_DEST_NID]);
		TRACE("weight field:%s\n", fields[curr_field+CONNECTION_WEIGHT]);
		ret += sscanf(fields[curr_field+CONNECTION_WEIGHT],
							"%lf", &weight);
		if (ret < 3)
		{
			INFO("group id:%s\n",
				fields[curr_field+CONNECTION_DEST_GID]);
			INFO("neuron id:%s\n",
				fields[curr_field+CONNECTION_DEST_NID]);
			INFO("weight field:%s\n",
				fields[curr_field+CONNECTION_WEIGHT]);

			INFO("Error: Couldn't parse synapse.\n");
			return COMMAND_FAIL;
		}

		src = n;
		assert(dest_gid < net->neuron_group_count);
		dest_group = &(net->groups[dest_gid]);
		assert(dest_nid < dest_group->neuron_count);
		dest = &(dest_group->neurons[dest_nid]);

		con->pre_neuron = src;
		con->post_neuron = dest;
		con->weight = weight;
		TRACE("\tAdded con %d.%d->%d.%d (w:%lf)\n",
			con->pre_neuron->group->id, con->pre_neuron->id,
			con->post_neuron->group->id, con->post_neuron->id,
			con->weight);
	}

	return COMMAND_OK;
}

int command_map_hardware(struct network *const net, struct architecture *arch,
				char fields[][MAX_FIELD_LEN],
				const int field_count)
{
	struct neuron_group *group;
	struct neuron *n;
	struct tile *t;
	struct core *c;
	struct hardware_mapping map;
	int neuron_group_id, tile_id, core_id, neuron_id, ret;

	TRACE("Parsing mapping.\n");
	if (field_count < 5)
	{
		INFO("Error: Invalid mapping command; %d fields.\n",
								field_count);
		return COMMAND_FAIL;
	}
	ret = sscanf(fields[1], "%d", &neuron_group_id);
	ret += sscanf(fields[2], "%d", &neuron_id);
	ret += sscanf(fields[3], "%d", &tile_id);
	ret += sscanf(fields[4], "%d", &core_id);
	if (ret < 4)
	{
		INFO("Error: Couldn't parse command.\n");
		return COMMAND_FAIL;
	}
	group = &(net->groups[neuron_group_id]);
	n = &(group->neurons[neuron_id]);

	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);
	map.core = c;
	// TODO: support mapping to different hardware blocks
	map.axon_in = &(c->axon_in);
	map.synapse_hw = &(c->synapse);
	map.dendrite_hw = &(c->dendrite);
	map.soma_hw = &(c->soma);
	map.axon_out = &(c->axon_out);

	return arch_map_neuron(n, map);
}

int command_parse_axon_input(struct architecture *const arch,
						char fields[][MAX_FIELD_LEN],
							const int field_count)
{
	struct tile *t;
	struct core *c;
	int tile_id, core_id, ret;

	ret = sscanf(fields[1], "%d", &tile_id);
	ret += sscanf(fields[2], "%d", &core_id);
	if (ret < 2)
	{
		INFO("Error: Couldn't parse axon input.\n");
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);

	arch_create_axon_in(arch, c);
	return COMMAND_OK;
}

int command_parse_synapse(struct architecture *const arch,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	struct tile *t;
	struct core *c;
	double spike_energy, spike_time, memory_energy, memory_time;
	int tile_id, core_id, weight_bits, word_bits, ret;

	ret = sscanf(fields[1], "%d", &tile_id);
	ret += sscanf(fields[2], "%d", &core_id);
	ret += sscanf(fields[3], "%d", &weight_bits);
	ret += sscanf(fields[4], "%d", &word_bits);
	ret += sscanf(fields[5], "%lf", &spike_energy);
	ret += sscanf(fields[6], "%lf", &spike_time);
	ret += sscanf(fields[7], "%lf", &memory_energy);
	ret += sscanf(fields[8], "%lf", &memory_time);
	if (ret < 8)
	{
		INFO("Error: Couldn't parse synapse processor.\n");
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);

	arch_create_synapse(arch, c, weight_bits, word_bits,
			spike_energy, spike_time, memory_energy, memory_time);
	return COMMAND_OK;
}

int command_parse_soma(struct architecture *const arch,
			char fields[][MAX_FIELD_LEN], const int field_count)
{
	struct tile *t;
	struct core *c;
	double active_energy, active_time, inactive_energy, inactive_time;
	double spiking_energy, spiking_time;
	int tile_id, core_id, model, ret;

	ret = sscanf(fields[1], "%d", &tile_id);
	ret += sscanf(fields[2], "%d", &core_id);

	if (strcmp(fields[3], "integrate_fire") == 0)
	{
		model = NEURON_IF;
	}
	else if (strcmp(fields[3], "leaky_integrate_fire") == 0)
	{
		model = NEURON_LIF;
	}
	else if (strcmp(fields[3], "truenorth") == 0)
	{
		model = NEURON_TRUENORTH;
	}
	else
	{
		// TODO: make invalid, must specify neuron model
		model = NEURON_IF;
	}


	ret += sscanf(fields[4], "%lf", &active_energy);
	ret += sscanf(fields[5], "%lf", &active_time);
	ret += sscanf(fields[6], "%lf", &inactive_energy);
	ret += sscanf(fields[7], "%lf", &inactive_time);
	ret += sscanf(fields[8], "%lf", &spiking_energy);
	ret += sscanf(fields[9], "%lf", &spiking_time);
	if (ret < 8)
	{
		INFO("Error: Couldn't parse soma processor.\n");
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);

	arch_create_soma(arch, c, model, active_energy, active_time,
				inactive_energy, inactive_time,
				spiking_energy, spiking_time);
	return COMMAND_OK;
}

int command_parse_axon_output(struct architecture *arch,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	struct tile *t;
	struct core *c;
	double spike_energy, spike_time;
	int tile_id, core_id, ret;

	ret = sscanf(fields[1], "%d", &tile_id);
	ret += sscanf(fields[2], "%d", &core_id);
	ret += sscanf(fields[3], "%lf", &spike_energy);
	ret += sscanf(fields[4], "%lf", &spike_time);
	if (ret < 4)
	{
		INFO("Error: Couldn't parse axon input.\n");
		return COMMAND_FAIL;
	}
	t = &(arch->tiles[tile_id]);
	c = &(t->cores[core_id]);

	TRACE("Parsed axon output: spike energy:%e time:%e\n",
						spike_energy, spike_time);

	arch_create_axon_out(arch, c, spike_energy, spike_time);
	return COMMAND_OK;
}

int command_parse_extern_input_group(struct network *const net,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	int input_type, input_count, ret;

	TRACE("Parsing group of inputs.\n");
	if (net->external_inputs != NULL)
	{
		INFO("Error: Inputs already created.\n");
		return COMMAND_FAIL;
	}

	if (field_count < 3)
	{
		INFO("Error: Invalid <input group> command; (%d/3) fields.\n",
								field_count);
		return COMMAND_FAIL;
	}

	ret = sscanf(fields[1], "%d", &input_count);
	if (ret < 1)
	{
		INFO("Error: Couldn't parse input count.\n");
		return COMMAND_FAIL;
	}

	// Parse the type of inputs to use. For now just allow one type of
	//  input. If, in future, we wanted to mix and match, I could move this
	//  stuff to the input node
	if (strstr("event", fields[2]))
	{
		input_type = INPUT_EVENT;
	}
	else if (strstr("poisson", fields[2]) == 0)
	{
		input_type = INPUT_RATE;
	}
	else if (strstr("rate", fields[2]) == 0)
	{
		input_type = INPUT_POISSON;
	}
	else
	{
		INFO("Error: Invalid neuron type (%s).\n", fields[2]);
		return COMMAND_FAIL;
	}

 	return net_create_inputs(net, input_count, input_type);
}

int command_parse_extern_input_node(struct network *const net,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	struct input *in;
	int input_id, connection_count, ret;

	TRACE("Parsing input node.\n");
	if (field_count < 2)
	{
		INFO("Error: Invalid <input> command; (%d/2) fields.\n",
							field_count);
		return COMMAND_FAIL;
	}
	ret = sscanf(fields[1], "%d", &input_id);
	if ((ret < 1) || (input_id < 0) ||
					(input_id >= net->external_input_count))
	{
		INFO("Error: Invalid input ID (%s).\n", fields[1]);
		return COMMAND_FAIL;
	}
	in = &(net->external_inputs[input_id]);

	connection_count = (field_count - 2) / CONNECTION_FIELDS;
 	ret = net_create_input_node(in, connection_count);
	if (ret == NETWORK_INVALID_NID)
	{
		return COMMAND_FAIL;
	}
	TRACE("Parsed input iid:%d connections:%d\n",
					input_id, connection_count);
	if (connection_count <= 0)
	{
		INFO("Warning: input %d with no outgoing connections.\n",
								input_id);
		return COMMAND_OK;
	}

	TRACE("iid:%d creating %d connections.\n",
					in->id, in->post_connection_count);
	// Parse each synapse between neurons, those neurons must already exit
	for (int i = 0; i < connection_count; i++)
	{
		// Parse each connection
		struct neuron_group *dest_group;
		struct neuron *src, *dest;
		struct connection *con;
		double weight;
		int dest_gid, dest_nid;
		const int curr_field = 2 + (i*CONNECTION_FIELDS);
		assert(curr_field < field_count);
		con = &(in->connections[i]);

		ret = sscanf(fields[curr_field+CONNECTION_DEST_GID],
							"%d", &dest_gid);
		ret += sscanf(fields[curr_field+CONNECTION_DEST_NID],
							"%d", &dest_nid);
		TRACE("weight field:%s\n", fields[curr_field+CONNECTION_WEIGHT]);
		ret += sscanf(fields[curr_field+CONNECTION_WEIGHT],
							"%lf", &weight);
		if (ret < 3)
		{
			INFO("Error: Couldn't parse synapse.\n");
			return COMMAND_FAIL;
		}

		src = NULL; // Is not valid when using an input node to feed
		dest_group = &(net->groups[dest_gid]);
		dest = &(dest_group->neurons[dest_nid]);

		con->pre_neuron = src;
		con->post_neuron = dest;

		// TODO: clean up hacks, quantization effects on Loihi
		con->weight = (int) (weight * 64*127.0);
		con->weight = con->weight / (64*127.0);
		TRACE("\tAdded con in[%d]->%d.%d (w:%lf)\n", in->id,
			con->post_neuron->group->id, con->post_neuron->id,
			con->weight);
	}

	return COMMAND_OK;
}

int command_parse_input_spikes(struct network *const net,
						char fields[][MAX_FIELD_LEN],
						const int field_count)
{
	double val;
	int ret, input_count;

	input_count = field_count - 1;
	if (input_count > net->external_input_count)
	{
		INFO("Error: Too many inputs given.\n");
		return COMMAND_FAIL;
	}

	for (int i = 0; i < input_count; i++)
	{
		struct input *in = &(net->external_inputs[i]);
		ret = sscanf(fields[i+1], "%lf", &val);
		if (ret < 1)
		{
			INFO("Error: Couldn't parse input value (%s)\n",
								fields[i+1]);
			return COMMAND_FAIL;
		}
		else if (val < 0)
		{
			INFO("Warning: id:%d input rate < 0 (%lf)\n",
							in->id, in->spike_val);
		}
		//else if (val > 1.0)
		//{
		//	INFO("Warning: id:%d input rate > 1 (%lf)\n",
		//					in->id, in->spike_val);
		//}
		net_set_input(net, i, val);
		INFO("Parsed input %d=%lf\n", in->id, in->spike_val);
	}

	return COMMAND_OK;
}

int command_parse_step_sim(struct network *const net,
						struct architecture *const arch,
						char fields[][MAX_FIELD_LEN],
						const int field_count,
						FILE *probe_spikes_fp,
						FILE *potential_trace_fp,
						FILE *perf_fp)
{
	// TODO: should there be a way of getting stats out of this
	//  the sim_stats is normally accumulated by the simulation
	//  Probably I need to create a "simulation" struct with all of this
	//  in, including the file pointers that are copied everywhere
	// TODO: we need to know which timestep this is for the rate based stuff
	//  figure this out before I can reenable the function
	//sim_timestep(net, arch, timestep, probe_spikes_fp, potential_trace_fp,
	//							perf_fp);
	// TODO: disabled for now
	return COMMAND_OK;
}

int command_parse_load_commands(struct network *const net,
						struct architecture *const arch,
						char fields[][MAX_FIELD_LEN],
						const int field_count,
						FILE *probe_spikes_fp,
						FILE *potential_trace_fp,
						FILE *probe_perf)
{
	FILE *fp;
	int ret;

	fp = fopen(fields[1], "r");
	if (fp == NULL)
	{
		return COMMAND_FAIL;
	}
	ret = command_parse_file(fp, net, arch, probe_spikes_fp,
						potential_trace_fp, probe_perf);
	fclose(fp);

	return ret;
}
