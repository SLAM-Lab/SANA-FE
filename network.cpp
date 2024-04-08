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

int network_create_neuron_group(Network &net, const int neuron_count,
	const std::list<Attribute> &attr)
{
	int id, ret;
	id = net.groups.size();
	INFO("Creating neuron group: %d with %d neurons\n", id, neuron_count);
	NeuronGroup &group = net.groups[id];
	group.neurons.reserve(neuron_count);

	group.default_soma_hw_name = "";
	group.default_synapse_hw_name = "";

	group.default_log_potential = 0; // Disabled by default
	group.default_log_spikes = 0;
	group.default_force_update = 0;

	for (auto a: attr)
	{
		ret = 1;
		std::istringstream ss(a.value_str);
		if (a.key == "soma_hw_name")
		{
			group.default_soma_hw_name = a.value_str;
			ret = 1;
		}
		else if (a.key == "synapse_hw_name")
		{
			group.default_synapse_hw_name = a.value_str;
			ret = 1;
		}
		else if (a.key =="log_v")
		{
			ss >> group.default_log_potential;
		}
		else if (a.key =="log_spikes")
		{
			ss >> group.default_log_spikes;
		}
		else if (a.key =="force_update")
		{
			ss >> group.default_force_update;
		}

		if (ret < 1)
		{
			INFO("Invalid attribute (%s:%s)\n", a.key.c_str(),
				a.value_str.c_str());
			exit(1);
		}
	}

	// Initialize all neurons in this group
	for (int i = 0; i < neuron_count; i++)
	{
		Neuron n;
		n.id = i;
		n.parent_group_id = group.id;
		n.connection_out_count = 0;
		INFO("Default soma name:%s\n", group.default_soma_hw_name.c_str());
		INFO("n id:%d\n", n.id);
		n.soma_hw_name = group.default_soma_hw_name,

		// Initialize neuron using group attributes
		n.log_spikes = group.default_log_spikes;
		n.log_potential = group.default_log_potential;
		n.force_update = group.default_force_update;
		n.max_connections_out = group.default_max_connections_out;

		n.fired = 0;
		// By default, dendrite current resets every timestep

		n.update_needed = 0;
		n.neuron_status = IDLE;
		n.spike_count = 0;

		// Initially the neuron is not mapped to anything
		n.core = nullptr;
		n.soma_hw = nullptr;
		n.axon_out_hw = nullptr;
		n.is_init = 0;

		// Create Soma Class Instance
		// TODO: problem with this is we now make the soma name
		//  optional! Fix this
		//n.soma_model = plugin_get_soma(n.soma_hw_name);
		n.model = new LoihiLifModel();
		n.model->set_attributes(attr);

		group.neurons.push_back(n);
	}

	INFO("Created neuron group gid:%d", group.id);

	return id;
}

int network_create_neuron(Neuron &n, const std::list<Attribute> &attr)
{
	// Each hardware timestep corresponds to a simulation of the spiking
	//  network for dt seconds. This relates to the LIF time constant.
	if (n.is_init)
	{
		INFO("Error: Trying to redefine neuron %d.\n", n.id);
		return NETWORK_INVALID_NID;
	}

	/*** Set attributes ***/
	for (auto a: attr)
	{
		int ret = 1;

		std::istringstream ss(a.value_str);
		if (a.key == "hw_name")
		{
			n.soma_hw_name = a.value_str;
		}
		else if (a.key == "connections_out")
		{
			ss >> n.max_connections_out;
		}
		else if (a.key == "log_spikes")
		{
			ss >> n.log_spikes;
		}
		else if (a.key == "log_v")
		{
			ss >> n.log_potential;
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
	n.update_needed = n.force_update; // || (fabs(n.bias) > 0.0));
	n.neuron_status = IDLE;

	n.soma_last_updated = 0;
	n.dendrite_last_updated = 0;

	n.core = NULL;
	assert(n.connections_out.size() = 0);
	n.connections_out.reserve(n.max_connections_out);

	// Check if need to create Soma Class instance
	if (n.model == nullptr)
	{
		// TODO: remove hack, make this user input
		INFO("Soma hw name:%s", n.soma_hw_name.c_str());
		//n.soma_model = plugin_get_soma(n.soma_hw_name);
		n.model = new LoihiLifModel();
		INFO("Creating new neuron %d\n", n.id);
	}
	n.model->set_attributes(attr);

	TRACE1("Created neuron: gid:%d nid:%d force:%d soma:%s\n",
		n.parent_group_id, n.id, n.force_update,
		n.soma_hw_name.c_str());
	n.is_init = 1;
	return n.id;
}

int network_connect_neurons(Connection &con,
	Neuron &src, Neuron &dest,
	const std::list<Attribute> &attr)
{
	INFO("dest id:%d.%d\n", dest.id, dest.parent_group_id);
	// TODO: set based on group defaults
	con.pre_neuron = &src;
	con.post_neuron = &dest;
	con.weight = 1.0;
	con.last_updated = 0;

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

	INFO("\tAdded con %d.%d->%d.%d (w:%lf)\n",
		con.pre_neuron->parent_group_id,
		con.pre_neuron->id,
		con.post_neuron->parent_group_id,
		con.post_neuron->id, con.weight);
	return RET_OK;
}

void network_check_mapped(Network &net)
{
	// Check that all network neurons are mapped to a physical core
	//  If a neuron is not, print an error message and stop the simulation
	for (auto &group: net.groups)
	{
		for (auto &n: group.neurons)
		{
			if (n.core == nullptr)
			{
				INFO("Error: Neuron %d.%d not mapped to H/W.\n",
					group.id, n.id);
				exit(1);
			}
		}
	}
}

/*
int network_create_inputs(Network *const net, const int input_count,
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
	in->connections = (Connection *)
		malloc(sizeof(Connection) * connection_count);
	if (in->connections == NULL)
	{
		INFO("Error: Couldn't allocate connection memory.\n");
		return NETWORK_INVALID_NID;
	}

	// Zero initialize all connections
	in->post_connection_count = connection_count;
	for (int i = 0; i < in->post_connection_count; i++)
	{
		Connection *con = &(in->connections[i]);
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

/*
Neuron *network_id_to_neuron_ptr(
	Network *const net, const struct neuron_id id)
{
	struct neuron_group *group;
	Neuron *neuron;

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
*/