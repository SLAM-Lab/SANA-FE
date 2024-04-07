// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.c
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <sstream>
#include <list>

#include "arch.hpp"
#include "print.hpp"
#include "network.hpp"
#include "models.hpp"

int total_connection_count = 0;

int network_create_neuron_group(struct network &net, const int neuron_count,
	const std::list<attribute> &attr)
{
	struct neuron_group *group;
	int id, ret;

	id = net.neuron_group_count;
	assert(id < NETWORK_MAX_NEURON_GROUPS);
	net.neuron_group_count++;

	INFO("Creating neuron group: %d with %d neurons\n", id, neuron_count);
	group = &(net.groups[id]);
	group->neurons.reserve(neuron_count);
	group->neuron_count = neuron_count;

	group->default_soma_hw_name = "";
	group->default_synapse_hw_name = "";
	group->default_max_connections_out = 0;
	group->default_log_potential = 0; // Disabled by default
	group->default_log_spikes = 0;
	group->default_threshold = 1.0;
	group->default_reverse_threshold = -1.0;
	group->default_reset = 0.0;
	group->default_reverse_reset = 0.0;
	group->reset_mode = NEURON_RESET_HARD;
	group->reverse_reset_mode = NEURON_NO_RESET; // Disabled by default
	group->default_force_update = 0;
	// Default is no leak (potential decay), i.e., the potential for the
	//  next timestep is 100% of the previous timestep's
	group->default_leak_decay = 1.0;
	group->default_leak_bias = 0.0;
	// TODO: change how group attributes are handled. Maybe don't bother
	//  parsing them. Just parse them over and over for each new neuron?
	//  Or create a neuron object that is the template to copy the rest
	//  from
	for (auto a: attr)
	{
		ret = 1;
		std::istringstream ss(a.value_str);
		if (a.key == "soma_hw_name")
		{
			group->default_soma_hw_name = a.value_str;
			ret = 1;
		}
		else if (a.key == "synapse_hw_name")
		{
			group->default_synapse_hw_name = a.value_str;
			ret = 1;
		}
		else if (a.key == "threshold")
		{
			ss >> group->default_threshold;
		}
		else if (a.key == "reverse_threshold")
		{
			ss >> group->default_reverse_threshold;
		}
		else if (a.key == "reset")
		{
			ss >> group->default_reset;
		}
		else if (a.key == "reverse_reset")
		{
			ss >> group->default_reverse_reset;
		}
		else if (a.key == "reset_mode")
		{
			group->reset_mode =
				network_parse_reset_mode(a.value_str);
			// Was parsed successfully if we got here
			ret = 1;
		}
		else if (a.key =="reverse_reset_mode")
		{
			group->reverse_reset_mode =
				network_parse_reset_mode(a.value_str);
			// Was parsed successfully if we got here
			ret = 1;
		}
		else if (a.key =="leak_decay")
		{
			ss >> group->default_leak_decay;
		}
		else if (a.key =="leak_bias")
		{
			ss >> group->default_leak_bias;
		}
		else if (a.key =="log_v")
		{
			ss >> group->default_log_potential;
		}
		else if (a.key =="log_spikes")
		{
			ss >> group->default_log_spikes;
		}
		else if (a.key =="connections_out")
		{
			ss >> group->default_max_connections_out;
		}
		else if (a.key =="force_update")
		{
			ss >> group->default_force_update;
		}

		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a.key.c_str(),
				a.value_str.c_str());
			exit(1);
		}
	}

	// Initialize all neurons in this group
	for (int i = 0; i < group->neuron_count; i++)
	{
		neuron n;
		n.id = i;
		n.group = group;
		n.connection_out_count = 0;
		INFO("Default soma name:%s\n", group->default_soma_hw_name.c_str());
		INFO("n id:%d\n", n.id);
		n.soma_hw_name = "";
		n.soma_hw_name = group->default_soma_hw_name,

		// Initialize neuron using group attributes
		n.log_spikes = group->default_log_spikes;
		n.log_potential = group->default_log_potential;
		n.force_update = group->default_force_update;
		n.max_connections_out = group->default_max_connections_out;

		n.fired = 0;
		// By default, dendrite current resets every timestep

		n.update_needed = 0;
		n.neuron_status = IDLE;
		n.spike_count = 0;

		// Initially the neuron is not mapped to anything
		n.core = NULL;
		n.soma_hw = NULL;

		n.maps_in = NULL;
		n.maps_out = NULL;
		n.maps_in_count = 0;
		n.maps_out_count = 0;

		n.is_init = 0;

		// Create Soma Class Instance
		// TODO: problem with this is we now make the soma name
		//  optional! Fix this
		//n.soma_model = plugin_get_soma(n.soma_hw_name);
		n.soma_model = new Loihi_Lif_Model();
		n.soma_model->set_attributes(attr);

		group->neurons.push_back(n);
	}

	INFO("Created neuron group gid:%d count:%d "
	     "threshold:%lf neg threshold:%lf "
	     "pos reset:%lf reset mode:%d "
	     "neg reset:%lf neg reset mode:%d\n",
		group->id, group->neuron_count, group->default_threshold,
		group->default_reverse_threshold, group->default_reset,
		group->reset_mode, group->default_reverse_reset,
		group->reverse_reset_mode);

	return id;
}

int network_create_neuron(
	struct neuron *const n, const std::list<attribute> &attr)
{
	// Each hardware timestep corresponds to a simulation of the spiking
	//  network for dt seconds. This relates to the LIF time constant.
	assert(n != NULL);
	if (n->is_init)
	{
		INFO("Error: Trying to redefine neuron %d.\n", n->id);
		return NETWORK_INVALID_NID;
	}

	/*** Set attributes ***/
	for (auto a: attr)
	{
		int ret = 1;

		std::istringstream ss(a.value_str);
		if (a.key == "hw_name")
		{
			n->soma_hw_name = a.value_str;
		}
		else if (a.key == "connections_out")
		{
			ss >> n->max_connections_out;
		}
		else if (a.key == "log_spikes")
		{
			ss >> n->log_spikes;
		}
		else if (a.key == "log_v")
		{
			ss >> n->log_potential;
		}
		else
		{
		 	TRACE1("Attribute %s not supported.\n", a->key);
		}

		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a.key.c_str(),
				a.value_str.c_str());
			exit(1);
		}
	}

	// Set the initial update state, no spikes can arrive before the first
	//  time-step but we can force the neuron to update, or bias it
	n->update_needed = n->force_update; // || (fabs(n->bias) > 0.0));
	n->neuron_status = IDLE;

	n->soma_last_updated = 0;
	n->dendrite_last_updated = 0;

	n->core = NULL;
	assert(n->connections_out.size() = 0);
	n->connections_out.reserve(n->max_connections_out);

	// Check if need to create Soma Class instance
	if (n->soma_model == nullptr)
	{
		// TODO: remove hack, make this user input
		INFO("Soma hw name:%s", n->soma_hw_name.c_str());
		//n->soma_model = plugin_get_soma(n->soma_hw_name);
		n->soma_model = new Loihi_Lif_Model();
		INFO("Creating new neuron %d\n", n->id);
	}
	n->soma_model->set_attributes(attr);

	TRACE1("Created neuron: gid:%d nid:%d force:%d soma:%s\n",
		n->group->id, n->id, n->force_update, n->soma_hw_name.c_str());
	n->is_init = 1;
	return n->id;
}

int network_connect_neurons(struct connection &con,
	struct neuron *const src, struct neuron *const dest,
	const std::list<attribute> &attr)
{
	INFO("dest id:%d\n", dest->id);
	INFO("dest group:%d\n", dest->group->id);
	INFO("default synapse name len:%ld\n", dest->group->default_synapse_hw_name.length());
	INFO("con hw:%s\n", con.synapse_hw_name.c_str());
	con.synapse_hw_name = dest->group->default_synapse_hw_name;
	con.pre_neuron = src;
	con.post_neuron = dest;
	con.weight = 1.0;

	for (auto a: attr)
	{
		std::istringstream ss(a.value_str);
		if ((a.key[0] == 'w') || (a.key == "weight"))
		{
			ss >> con.weight;
		}
		else if (a.key == "hw_name")
		{
			con.synapse_hw_name = a.value_str;
		}
	}

	TRACE1("\tAdded con %d.%d->%d.%d (w:%lf)\n", con.pre_neuron->group->id,
		con.pre_neuron->id, con.post_neuron->group->id,
		con.post_neuron->id, con.weight);
	return RET_OK;
}

/*
int network_create_inputs(struct network *const net, const int input_count,
							const int input_type)
{
	assert(input_count > 0);
	assert(net->external_inputs == NULL);

	net->external_inputs =
		(struct input *) malloc(sizeof(struct input) * input_count);
	if (net->external_inputs == NULL)
	{
		INFO("Error: Couldn't allocate memory for network inputs.\n");
		return -1;
	}
	net->external_input_count = input_count;

	for (int i = 0; i < input_count; i++)
	{
		struct input *in = &(net->external_inputs[i]);

		in->id = i;
		in->connections = NULL;
		in->spike_val = 0.0;
		in->post_connection_count = 0;
		in->type = input_type;
	}

	return 0;
}

int network_create_input_node(struct input *const in, const int connection_count)
{
	TRACE("Creating input node with %d connections.\n", connection_count);

	assert(in->connections == NULL);
	in->connections = (struct connection *)
		malloc(sizeof(struct connection) * connection_count);
	if (in->connections == NULL)
	{
		INFO("Error: Couldn't allocate connection memory.\n");
		return NETWORK_INVALID_NID;
	}

	// Zero initialize all connections
	in->post_connection_count = connection_count;
	for (int i = 0; i < in->post_connection_count; i++)
	{
		struct connection *con = &(in->connections[i]);
		con->pre_neuron = NULL;
		con->post_neuron = NULL;
		con->weight = 0.0;
		con->current = 0.0;
		con->synaptic_current_decay = 0.0;
		con->delay = 0.0;
	}

	return in->id;
}
*/

void network_init(struct network *const net)
{
	net->neuron_group_count = 0;
	net->external_input_count = 0;
	net->external_inputs = NULL;

	for (int i = 0; i < NETWORK_MAX_NEURON_GROUPS; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		group->id = i;
		group->neuron_count = 0;
	}
}

void network_free(struct network *const net)
{
	for (int i = 0; i < net->external_input_count; i++)
	{
		struct input *in = &(net->external_inputs[i]);
		free(in->connections);
	}
	free(net->external_inputs);

	return;
}

struct neuron *network_id_to_neuron_ptr(
	struct network *const net, const struct neuron_id id)
{
	struct neuron_group *group;
	struct neuron *neuron;

	if (id.group < net->neuron_group_count)
	{
		INFO("ERROR: Group %d > max %d.\n", id.group,
			net->neuron_group_count);
		exit(1);
	}

	group = &(net->groups[id.group]);
	if (id.neuron < group->neuron_count)
	{
		INFO("ERROR: Neuron %d > max %d.\n", id.neuron,
			group->neuron_count);
		exit(1);
	}
	neuron = &(group->neurons[id.neuron]);

	return neuron;
}

void network_set_input(
	struct network *const net, const int input_id, const double rate)
{
	struct input *in = &(net->external_inputs[input_id]);
	in->rate = rate;
}

void network_check_mapped(struct network *const net)
{
	// Check that all network neurons are mapped to a physical core
	//  If a neuron is not, print an error message and stop the simulation
	for (int i = 0; i < net->neuron_group_count; i++)
	{
		struct neuron_group *group = &(net->groups[i]);

		for (int j = 0; j < group->neuron_count; j++)
		{
			struct neuron *n = &(group->neurons[j]);

			if (n->core == NULL)
			{
				INFO("Error: Neuron %d.%d not mapped to H/W.\n",
					group->id, n->id);
				exit(1);
			}
		}
	}
}
