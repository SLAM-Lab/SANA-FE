// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.cpp
#include <cctype>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <list>
#include <set>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <functional> // For std::reference_wrapper

#include "print.hpp"
#include "arch.hpp"
#include "network.hpp"

using namespace sanafe;

sanafe::Architecture::Architecture()
{
	noc_buffer_size = 0;
	noc_init = false;

	return;
}

int sanafe::Architecture::get_core_count()
{
	int core_count = 0;
	for (Tile &tile: tiles)
	{
		core_count += tile.cores.size();
	}

	return core_count;
}

int sanafe::Architecture::set_noc_attributes(
	const std::unordered_map<std::string, std::string> &attr)
{
	const size_t tile_count = tiles.size();
	if (tile_count == 0)
	{
		// The NoC interconnect is defined after tiles are all defined
		//  This is because we link the tiles together in the NoC mesh
		INFO("Error: NoC must be defined after tiles.\n");
		throw std::runtime_error(
			"Error: NoC must be defined after tiles");
	}

	// Default values
	noc_width = 1;
	noc_height = 1;
	noc_buffer_size = 0;

	for (auto a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		std::istringstream ss(value_str);
		if (key == "width")
		{
			ss >> noc_width;
		}
		else if (key == "height")
		{
			ss >> noc_height;
		}
		else if (key == "link_buffer_size")
		{
			ss >> noc_buffer_size;
		}
	}
	assert((noc_height * noc_width) <= ARCH_MAX_TILES);
	noc_init = true;
	TRACE1("NoC created, mesh, width:%d height:%d.\n", noc_width,
		noc_height);
	return 0;
}

sanafe::Message::Message()
{
	// Initialize message variables. Mark most fields as invalid either
	//  using NaN or -Inf values where possible.
	src_neuron = nullptr;

	dummy_message = false;
	generation_delay = 0.0;
	network_delay = NAN;
	receive_delay = NAN;
	blocked_latency = 0.0;
	hops = -1;
	spikes = -1;
	sent_timestamp = -INFINITY;
	received_timestamp = -INFINITY;
	processed_timestamp = -INFINITY;
	timestep = -1;
	next = nullptr;

	src_x = -1;
	src_y = -1;
	dest_x = -1;
	dest_y = -1;
	dest_core_offset = -1;
	dest_tile_id = -1;
	dest_axon_id = -1;
}

sanafe::Tile::Tile(const int tile_id)
{
	id = tile_id;
	energy = 0.0;
	x = 0;
	y = 0;
	hops = 0L;
	messages_received = 0L;
	total_neurons_fired = 0L;
	east_hops = 0L;
	west_hops = 0L;
	north_hops = 0L;
	south_hops = 0L;
}

sanafe::Core::Core(const int core_id, const int tile_id, const int core_offset)
{
	id = core_id;
	parent_tile_id = tile_id;
	offset = core_offset;
	energy = 0.0;
	latency_after_last_message = 0.0;
}

sanafe::Tile &sanafe::Architecture::create_tile(
	const std::unordered_map<std::string, std::string> &attr)

{
	tiles.push_back(Tile(tiles.size()));
	Tile &tile = tiles.back();
	tiles_vec.push_back(tile);

	// Set attributes
	tile.energy_east_hop = 0.0;
	tile.latency_east_hop = 0.0;
	tile.energy_north_hop = 0.0;
	tile.latency_north_hop = 0.0;
	tile.energy_west_hop = 0.0;
	tile.latency_west_hop = 0.0;
	tile.energy_south_hop = 0.0;
	tile.latency_south_hop = 0.0;
	for (auto a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		std::istringstream ss(value_str);
		if (key == "energy_east")
		{
			ss >> tile.energy_east_hop;
		}
		else if (key == "latency_east")
		{
			ss >> tile.latency_east_hop;
		}
		else if (key == "energy_west")
		{
			ss >> tile.energy_west_hop;
		}
		else if (key == "latency_west")
		{
			ss >> tile.latency_west_hop;
		}
		else if (key == "energy_north")
		{
			ss >> tile.energy_north_hop;
		}
		else if (key == "latency_north")
		{
			ss >> tile.latency_north_hop;
		}
		else if (key == "energy_south")
		{
			ss >> tile.energy_south_hop;
		}
		else if (key == "latency_south")
		{
			ss >> tile.latency_south_hop;
		}
	}

	return tile;
}

void Architecture::load_arch_file(const std::string &filename)
{
	std::ifstream arch_fp(filename);
	if (arch_fp.fail())
	{
		throw std::invalid_argument(
			"Error: Architecture file failed to open.");
	}
	int ret = description_parse_arch_file(arch_fp, *this);
	arch_fp.close();
	if (ret == RET_FAIL)
	{
		throw std::invalid_argument(
			"Error: Invalid architecture file.");
	}
}

sanafe::AxonInUnit::AxonInUnit(const std::string &axon_in_name):
	name(axon_in_name)
{
	return;
}

sanafe::SynapseUnit::SynapseUnit(const std::string &synapse_name):
	name(synapse_name)
{
	energy = 0.0;
	time = 0.0;
	return;
}

sanafe::SomaUnit::SomaUnit(const std::string &soma_name):
	name(soma_name)
{
	neuron_updates = 0;
	neurons_fired = 0;
	neuron_count = 0;
	energy = 0.0;
	time = 0.0;
	return;
}

sanafe::AxonOutUnit::AxonOutUnit(const std::string &axon_out_name):
	name(axon_out_name)
{
	energy = 0.0;
	time = 0.0;
	return;
}

sanafe::Core &sanafe::Architecture::create_core(
	const size_t tile_id,
	const std::unordered_map<std::string, std::string> &attr)
{
	if (tile_id > tiles.size())
	{
		throw std::invalid_argument("Error: Tile ID > total tiles");
	}
	Tile &tile = tiles_vec[tile_id];
	const int core_offset = tile.cores.size();
	tile.cores.push_back(Core(get_core_count(), tile.id, core_offset));
	Core &c = tile.cores.back();
	cores_vec.push_back(c);
	tile.cores_vec.push_back(c);

	// *** Set attributes ***
	c.buffer_pos = BUFFER_SOMA;
	c.max_neurons = 1024;
	for (const auto &a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		if (key == "buffer_before")
		{
			if (value_str == "soma")
			{
				c.buffer_pos = BUFFER_SOMA;
			}
		}
		else if (key == "max_neurons")
		{
			std::istringstream ss(value_str);
			ss >> c.max_neurons;
		}
	}

	// Initialize core state
	c.energy = 0.0;
	c.next_message = Message();

	TRACE1("Core created id:%d.%d (tile:%d).\n", c.parent_tile_id, c.id);
	return c;
}

sanafe::AxonInUnit &sanafe::Core::create_axon_in(
	const std::string &name,
	const std::unordered_map<std::string, std::string> &attr)
{
	axon_in_hw.push_back(AxonInUnit(name));
	struct AxonInUnit &in = axon_in_hw.back();

	in.energy = 0.0;
	in.time = 0.0;
	in.parent_tile_id = parent_tile_id;

	in.energy_spike_message = 0.0;
	in.latency_spike_message = 0.0;
	for (const auto &curr: attr)
	{
		const std::string &key = curr.first;
		const std::string &value_str = curr.second;
		std::istringstream ss(value_str);
		if (key == "energy_message")
		{
			ss >> in.energy_spike_message;
		}
		else if (key == "latency_message")
		{
			ss >> in.latency_spike_message;
		}
	}

	TRACE2("Axon input created (c:%d.%d)\n", c.t->id, c.id);

	return in;
}

sanafe::SynapseUnit &sanafe::Core::create_synapse(
	const std::string &name,
	const std::unordered_map<std::string, std::string> &attr)
{
	synapse.push_back(SynapseUnit(name));
	SynapseUnit &s = synapse.back();

	/**** Set attributes ****/
	s.energy_memory_access = 0.0;
	s.latency_memory_access = 0.0;
	s.energy_spike_op = 0.0;
	s.latency_spike_op = 0.0;
	s.weight_bits = 8;
	s.name = name;
	for (const auto &curr: attr)
	{
		const std::string &key = curr.first;
		const std::string &value_str = curr.second;
		std::istringstream ss(value_str);
		if (key == "model")
		{
			s.model = arch_parse_synapse_model(value_str);
		}
		else if (key == "energy_memory")
		{
			ss >> s.energy_memory_access;
		}
		else if (key == "latency_memory")
		{
			ss >> s.latency_memory_access;
		}
		else if (key == "energy_spike")
		{
			ss >> s.energy_spike_op;
		}
		else if (key == "latency_spike")
		{
			ss >> s.latency_spike_op;
		}
	}

	TRACE1("Synapse processor created (c:%d.%d)\n", c.parent_tile_id, c.id);

	return s;
}

sanafe::SomaUnit &sanafe::Core::create_soma(
	const std::string &name,
	const std::unordered_map<std::string, std::string> &attr)
{
	INFO("cid:%d creating soma sid:%lu with %lu attributes\n",
		id, soma.size(), attr.size());
	soma.push_back(SomaUnit(name));
	SomaUnit &s = soma.back();

	/*** Set attributes ***/
	s.energy_access_neuron = 0.0;
	s.latency_access_neuron = 0.0;
	s.energy_update_neuron = 0.0;
	s.latency_update_neuron = 0.0;
	s.energy_spiking = 0.0;
	s.latency_spiking = 0.0;
	s.leak_towards_zero = 1;
	s.noise_type = NOISE_NONE;

	for (const auto &a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		std::istringstream ss(value_str);
		INFO("Soma attribute k:%s v%s\n", key.c_str(),
			value_str.c_str());
		if (key == "energy_update_neuron")
		{
			ss >> s.energy_update_neuron;
		}
		else if (key == "latency_update_neuron")
		{
			ss >> s.latency_update_neuron;
		}
		else if (key == "energy_access_neuron")
		{
			ss >> s.energy_access_neuron;
		}
		else if (key == "latency_access_neuron")
		{
			ss >> s.latency_access_neuron;
		}
		else if (key == "energy_spike_out")
		{
			ss >> s.energy_spiking;
		}
		else if (key == "latency_spike_out")
		{
			ss >> s.latency_spiking;
		}
		else if (key == "noise")
		{
			s.noise_type = NOISE_FILE_STREAM;
			s.noise_stream = fopen(value_str.c_str(), "r");
			TRACE1("Opening noise str: %s\n", value_str.c_str());
			if (s.noise_stream == NULL)
			{
				INFO("Error: Failed to open noise stream: %s.\n",
					value_str.c_str());
				exit(1);
			}
		}
	}

	TRACE1("Soma processor created (c:%d.%d)\n", c.parent_tile_id,
		c.offset);
	return s;
}

sanafe::AxonOutUnit &sanafe::Core::create_axon_out(
	const std::string &name,
	const std::unordered_map<std::string, std::string> &attr)
{
	axon_out_hw.push_back(AxonOutUnit(name));
	AxonOutUnit &out = axon_out_hw.back();

	out.packets_out = 0;
	out.energy = 0.0;
	out.time = 0.0;
	out.parent_tile_id = parent_tile_id;

	/*** Set attributes ***/
	out.energy_access = 0.0;
	out.latency_access = 0.0;
	for (const auto &curr: attr)
	{
		const std::string &key = curr.first;
		const std::string &value_str = curr.second;
		std::istringstream ss(value_str);
		if (key == "energy")
		{
			ss >> out.energy_access;
		}
		else if (key == "latency")
		{
			ss >> out.latency_access;
		}
	}

	// Track the tile the axon interfaces with
	TRACE1("Axon output created (c:%d.%d)\n", c.parent_tile_id, c.offset);

	return out;
}

void sanafe::arch_create_axons(Architecture &arch)
{
	TRACE1("Creating all connection maps.\n");
	for (Tile &tile: arch.tiles)
	{
		for (Core &c: tile.cores)
		{
			for (auto n_ptr: c.neurons)
			{
				arch_map_neuron_connections(*n_ptr);
			}
		}
	}

	TRACE1("Finished creating connection maps.\n");
	arch_print_axon_summary(arch);
}

void sanafe::arch_print_axon_summary(Architecture &arch)
{
	int in_count, out_count, core_used;
	in_count = 0;
	out_count = 0;

	INFO("** Mapping summary **\n");
	for (Tile &tile: arch.tiles)
	{
		// For debug only, print the axon maps
		for (Core &c: tile.cores)
		{
			core_used = 0;
			for (std::vector<Neuron *>::size_type k = 0;
				k < c.neurons.size(); k++)
			{
#ifdef DEBUG
				Neuron *n = c->neurons[k];
				TRACE2("\tnid:%d.%d ", n->group->id, n->id);
				TRACE2("i:%d o:%d\n", n->maps_in_count,
					n->maps_out_count);
#endif
				core_used = 1;
			}

			if (core_used)
			{
				// TODO: update these counts of axons
				//in_count += c.axon_in_hw.axons;
				//out_count += c.axon_out.map_count;
			}
		}
	}
	const int core_count = arch.get_core_count();
	INFO("Total cores: %d\n", core_count);
	INFO("Average in map count: %lf\n", (double) in_count / core_count);
	INFO("Average out map count: %lf\n", (double) out_count / core_count);

	return;
}

void sanafe::arch_map_neuron_connections(Neuron &pre_neuron)
{
	// Setup the connections between neurons and map them to hardware
	assert(pre_neuron.core != nullptr);

	// Figure out the unique set of cores that this neuron broadcasts to
	TRACE1("Counting connections for neuron nid:%d\n", pre_neuron.id);
	std::vector<bool> core_inserted();
	std::set<Core *> cores_out;
	for (Connection &curr_connection: pre_neuron.connections_out)
	{
		TRACE1("Looking at connection id: %d\n", curr_connection.id);
		Core *dest_core = curr_connection.post_neuron->core;
		cores_out.insert(dest_core);
		TRACE1("Connected to dest core: %d\n", dest_core->id);
	}

	TRACE1("Creating connections for neuron nid:%d to %lu core(s)\n",
		pre_neuron.id, cores_out.size());
	for (Core *dest_core: cores_out)
	{
		// Create the axon, and add it to both the destination and
		//  source cores
		arch_allocate_axon(pre_neuron, *dest_core);
	}
	TRACE3("Counted all maps for nid:%d count: %d\n",
		pre_neuron.id);

	for (Connection &curr_connection: pre_neuron.connections_out)
	{
		// Add every connection to the axon. Also link to the map in the
		//  post synaptic core / neuron
		Core &post_core = *(curr_connection.post_neuron->core);
		TRACE1("Adding connection:%d\n", curr_connection.id);
		arch_add_connection_to_axon(curr_connection, post_core);
	}
	TRACE1("Finished mapping connections to hardware for nid:%d.%d.\n",
		pre_neuron.parent_group_id, pre_neuron.id);

	return;
}

void sanafe::Core::map_neuron(Neuron &n)
{
	TRACE1("Mapping nid:%d to core: %d\n", n.id, id);
	// Map the neuron to hardware units
	if (n.core != nullptr)
	{
		throw std::runtime_error(
			"Error: Neuron already mapped, dynamic remapping not "
			"supported.");
	}

	n.core = this;
	TRACE1("Mapping neuron %d to core %d\n", n.id, id);
	neurons.push_back(&n);

	if (neurons.size() > max_neurons)
	{
		INFO("Error: Exceeded maximum neurons per core (%lu)",
			max_neurons);
		throw std::runtime_error(
			"Error: Exceeded maximum neurons per core.");
	}

	// Map neuron model to soma hardware unit in this core. Search through
	//  all neuron models implemented by this core and return the one that
	//  matches. If no soma hardware is specified, default to the first
	//  one defined
	if (soma.size() == 0)
	{
		INFO("Error: No soma units defined for cid:%d\n", id);
		throw std::runtime_error("Error: No soma units defined");
	}
	n.soma_hw = &(soma[0]);
	if (n.soma_hw_name.length() > 0)
	{
		bool soma_found = false;
		for (auto &soma_hw: soma)
		{
			if (n.soma_hw_name == soma_hw.name)
			{
				n.soma_hw = &soma_hw;
				soma_found = true;
			}
		}
		if (!soma_found)
		{
			INFO("Error: Could not map neuron nid:%d (hw:%s) "
				"to any soma h/w.\n", n.id,
				n.soma_hw_name.c_str());
			throw std::runtime_error(
				"Error: Could not map neuron to soma h/w");
		}
	}
	n.soma_hw->neuron_count++;

	// TODO: support multiple axon outputs
	if (axon_out_hw.size() == 0)
	{
		INFO("Error: No axon out units defined for cid:%d\n", id);
		throw std::runtime_error("Error: No axon out units defined");
	}
	n.axon_out_hw = &(axon_out_hw[0]);

	// Pass all the model specific arguments
	if (n.model == nullptr)
	{
		//n.soma_model = plugin_get_soma(n.soma_hw_name);
		n.model = std::shared_ptr<SomaModel>(
			new LoihiLifModel(n.parent_group_id, n.id));
		const NeuronGroup &group = n.parent_net->groups_vec[
			n.parent_group_id];
		// First set the group's default attribute values, and then
		//  any defined by the neuron
		n.model->set_attributes(group.default_attributes);
		n.model->set_attributes(n.attributes);
	}

	return;
}

void sanafe::arch_allocate_axon(Neuron &pre_neuron, Core &post_core)
{
	// Create a new input axon at a receiving (destination) core
	//  Then create the output axon at the sending core. Finally
	//  update the presynaptic neuron and postsynaptic neuron

	Core &pre_core = *(pre_neuron.core);

	TRACE3("Adding connection to core.\n");
	// Allocate the axon and its connections at the post-synaptic core
	AxonInModel in;
	in.active_synapses = 0;
	in.last_updated = -1;
	post_core.axons_in.push_back(in);
	const int new_axon_in_address = post_core.axons_in.size() - 1;

	// Add the axon at the sending, pre-synaptic core
	TRACE1("axon in address:%d for core:%d.%d\n",
		new_axon_in_address, post_core.parent_tile_id, post_core.id);
	AxonOutModel out;
	out.dest_axon_id = new_axon_in_address;
	out.dest_core_offset = post_core.offset;
	out.dest_tile_id = post_core.parent_tile_id;
	out.src_neuron_id = pre_neuron.id;
	pre_core.axons_out.push_back(out);
	const int new_axon_out_address = pre_core.axons_out.size() - 1;

	// Then add the output axon to the sending pre-synaptic neuron
	pre_neuron.axon_out_addresses.push_back(new_axon_out_address);
	TRACE1("nid:%d.%d cid:%d.%d added one output axon address %d.\n",
		pre_neuron.parent_group_id, pre_neuron.id,
		pre_core.parent_tile_id, pre_core.offset,
		new_axon_out_address);

	return;
}

void sanafe::arch_add_connection_to_axon(Connection &con, Core &post_core)
{
	// Add a given connection to the axon in the post-synaptic
	//  (destination) core
	TRACE3("Adding to connection to axon:%lu\n",
		post_core.axons_out.size()-1);

	post_core.synapses.push_back(&con);
	const std::vector<Connection *>::size_type synapse_address =
		post_core.synapses.size() - 1;

	// Access the most recently created axon in for the post-synaptic core
	AxonInModel &last_added_target_axon = post_core.axons_in.back();
	last_added_target_axon.synapse_addresses.push_back(synapse_address);

	// Map the connection to the synapse hardware
	const std::vector<SynapseUnit>::size_type default_hw_id = 0;
	// Default to the first defined hardware unit (there must be at least
	// one hardware unit defined)
	SynapseUnit &synapse_hw = post_core.synapse[default_hw_id];
	if (con.synapse_hw_name.length() > 0)
	{
		bool mapped = false;
		// Search for the specified synapse hardware
		for (const auto &s: post_core.synapse)
		{
			if (con.synapse_hw_name == s.name)
			{
				synapse_hw = s;
				mapped = true;
				break;
			}
		}
		if (!mapped)
		{
			INFO("Error: Could not map connection to synapse h/w.\n");
			exit(1);
		}
	}
	con.synapse_hw = &synapse_hw;

	return;
}

int sanafe::arch_parse_synapse_model(const std::string &model_str)
{
	int model;

	if (model_str == "current_based")
	{
		model = SYNAPSE_CUBA;
	}
	else
	{
		INFO("Error: No synapse model specified (%s)\n",
			model_str.c_str());
		exit(1);
	}

	return model;
}
