// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.c
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <sstream>
#include <list>
#include <memory>
#include <unordered_map>

#include "arch.hpp"
#include "print.hpp"
#include "network.hpp"
#include "models.hpp"

using namespace sanafe;

Connection::Connection(const int connection_id)
{
	id = connection_id;
	pre_neuron = nullptr;
	post_neuron = nullptr;
	synapse_hw = nullptr;
	last_updated = 0;

	current = 0.0;
	weight = 0.0;
	delay = 0.0;
	synaptic_current_decay = 0.0;
}

sanafe::NeuronGroup::NeuronGroup(const size_t group_id, const int neuron_count)
{
	INFO("Creating neuron group: %lu with %d neurons\n",
		group_id, neuron_count);
	id = group_id;
	default_max_connections_out = 0;
	default_soma_hw_name = "";
	default_synapse_hw_name = "";
	default_log_potential = false;
	default_log_spikes = false;
	default_force_update = false;

	// Reserve space for the neurons to go
	neurons.reserve(neuron_count);
}

sanafe::Neuron::Neuron(const size_t neuron_id)
{
	id = id;
	log_spikes = false;
	log_potential = false;
	neuron_status = sanafe::IDLE;
	spike_count = 0;
	max_connections_out = 0;

	// Initially the neuron is not mapped to anything
	core = nullptr;
	soma_hw = nullptr;
	axon_out_hw = nullptr;

	soma_last_updated = 0;
	dendrite_last_updated = 0;
	core = nullptr;
}

NeuronGroup &sanafe::Network::create_neuron_group(const int neuron_count,
	const std::unordered_map<std::string, std::string> &attr)
{
	const auto id = groups.size();

	groups.push_back(NeuronGroup(id, neuron_count));
	NeuronGroup &group = groups.back();

	for (auto a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		std::istringstream ss(value_str);
		if (key == "soma_hw_name")
		{
			group.default_soma_hw_name = value_str;
		}
		else if (key == "synapse_hw_name")
		{
			group.default_synapse_hw_name = value_str;
		}
		else if (key =="log_v")
		{
			ss >> group.default_log_potential;
		}
		else if (key =="log_spikes")
		{
			ss >> group.default_log_spikes;
		}
		else if (key =="force_update")
		{
			ss >> group.default_force_update;
		}
	}

	// Initialize all neurons in this group
	group.default_attributes = attr;
	for (int i = 0; i < neuron_count; i++)
	{
		Neuron n(i);
		n.id = i;
		n.parent_group_id = group.id;
		TRACE1("Default soma name:%s\n",
			group.default_soma_hw_name.c_str());
		TRACE1("nid:%d\n", n.id);
		n.soma_hw_name = group.default_soma_hw_name,

		// Initialize neuron using group attributes
		n.log_spikes = group.default_log_spikes;
		n.log_potential = group.default_log_potential;
		n.max_connections_out = group.default_max_connections_out;
		n.parent_net = this;
		group.neurons.push_back(n);
	}

	INFO("Created neuron group gid:%d\n", group.id);

	return group;
}

void sanafe::Network::load_net_file(
	const std::string &filename, Architecture &arch)
{
	std::ifstream network_fp;
	network_fp.open(filename);
	if (network_fp.fail())
	{
		throw std::runtime_error("Error: Network file failed to open.");
	}
	INFO("Reading network from file.\n");
	int ret = description_parse_net_file(network_fp, *this, arch);
	network_fp.close();
	if (ret == RET_FAIL)
	{
		throw std::invalid_argument("Error: Invalid network file.");
	}
	network_check_mapped(*this);
	arch_create_axons(arch);
}

void sanafe::Neuron::set_attributes(
	const std::unordered_map<std::string, std::string> &attr)
{
	// Each hardware timestep corresponds to a simulation of the spiking
	//  network for dt seconds. This relates to the LIF time constant.
	/*** Set attributes ***/
	for (auto a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		std::istringstream ss(value_str);
		if (key == "hw_name")
		{
			soma_hw_name = value_str;
		}
		else if (key == "connections_out")
		{
			ss >> max_connections_out;
		}
		else if (key == "log_spikes")
		{
			ss >> log_spikes;
		}
		else if (key == "log_v")
		{
			ss >> log_potential;
		}
		else
		{
		 	TRACE1("Attribute %s not supported.\n", key.str());
		}
	}

	assert(connections_out.size() == 0);
	connections_out.reserve(max_connections_out);

	// Check if need to create Soma Class instance
	if (core != nullptr)
	{
		// TODO: remove hack, make this user input
		TRACE1("Soma hw name: %s", n.soma_hw_name.c_str());
		//n.soma_model = plugin_get_soma(n.soma_hw_name);
		model = std::shared_ptr<SomaModel>(
			new LoihiLifModel(parent_group_id, id));
		TRACE1("Creating new neuron %d\n", n.id);
	}
	if (model != nullptr)
	{
		model->set_attributes(attr);
	}
	attributes.insert(attr.begin(), attr.end());

	TRACE1("Created neuron: gid:%d nid:%d soma:%s\n",
		parent_group_id, id, soma_hw_name.c_str());
}

void sanafe::NeuronGroup::set_attribute_multiple(
	const std::string &attr, const std::vector<std::string> &values)
{
	if (values.size() != neurons.size())
	{
		INFO("Error: Attribute values must be defined for all neurons "
			"(%lu!=%lu).\n", attr.size(), neurons.size());
	}
	size_t nid = 0;
	std::unordered_map<std::string, std::string> map;
	for (const std::string &v: values)
	{
		auto &n = neurons[nid];
		map[attr] = v;
		n.set_attributes(map);
		nid++;
	}
}

void sanafe::NeuronGroup::connect_neurons(
	const NeuronGroup &dest_group,
	const std::vector<size_t> dest_neurons,
	const std::unordered_map<std::string, std::string> &attr)
{
	return;
}


void sanafe::Neuron::connect_to_neuron(
	Neuron &dest,
	const std::unordered_map<std::string, std::string> &attr)
{
	connections_out.push_back(Connection(connections_out.size()));
	Connection &con = connections_out.back();
	con.pre_neuron = this;
	con.post_neuron = &dest;

	for (auto a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		std::istringstream ss(value_str);
		if ((key[0] == 'w') || (key == "weight"))
		{
			ss >> con.weight;
		}
		else if (key == "hw_name")
		{
			con.synapse_hw_name = value_str;
		}
	}

	TRACE1("\tAdded con %d.%d->%d.%d (w:%lf)\n",
		con.pre_neuron->parent_group_id,
		con.pre_neuron->id,
		con.post_neuron->parent_group_id,
		con.post_neuron->id, con.weight);
	return;
}

void sanafe::network_check_mapped(Network &net)
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