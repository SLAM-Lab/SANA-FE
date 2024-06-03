// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.cpp
#include <cctype>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <filesystem> // For std::filesystem::path
#include <functional> // For std::reference_wrapper
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <sstream>

#include "arch.hpp"
#include "description.hpp"
#include "network.hpp"
#include "models.hpp"
#include "plugins.hpp"
#include "print.hpp"

sanafe::Architecture::Architecture()
{
	noc_buffer_size = 0;
	noc_init = false;
	max_cores_per_tile = 0U;

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
	const std::map<std::string, std::string> &attr)
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

	for (const auto &a: attr)
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

	// Set the x and y coordinates of the NoC tiles
	for (Tile &tile: tiles)
	{
		tile.x = tile.id / noc_height;
		tile.y = tile.id % noc_height;
	}

	noc_init = true;
	TRACE1("NoC created, mesh, width:%d height:%d.\n", noc_width,
		noc_height);
	return 0;
}

std::string sanafe::Architecture::info()
{
	std::ostringstream ss;
	ss << "sanafe::Architecture(tiles=" << tiles.size();
	ss << ", cores=" << cores_vec.size() << ")";

	return ss.str();
}

std::string sanafe::Architecture::description() const
{
	// TODO: change to just a regular map - this way attributes will be
	//  sorted and printed alphabetically. This will make for more
	//  consistent and deterministic behavior
	std::map<std::string, std::string> attributes;
	attributes["link_buffer_size"] = print_int(noc_buffer_size);
	attributes["topology"] = "mesh";
	attributes["width"] = print_int(noc_width);
	attributes["height"] = print_int(noc_height);

	std::ostringstream ss;
	ss << '@' << print_format_attributes(attributes) << std::endl;
	return ss.str();
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
	sent_timestamp = - std::numeric_limits<double>::infinity();
	received_timestamp = - std::numeric_limits<double>::infinity();
	processed_timestamp = - std::numeric_limits<double>::infinity();
	timestep = -1;
	in_noc = false;

	src_x = -1;
	src_y = -1;
	dest_x = -1;
	dest_y = -1;
	dest_core_offset = -1;
	dest_core_id = -1;
	dest_tile_id = -1;
	dest_axon_id = -1;
	dest_axon_hw = -1;
}

sanafe::Tile::Tile(const std::string &name, const int tile_id) : name(name)
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

std::string sanafe::Tile::info() const
{
	std::ostringstream ss;
	ss << "sanafe::Tile(tile=" << id << " cores=";
	ss << cores_vec.size() << ")";

	return ss.str();
}

std::string sanafe::Tile::description() const
{
	std::map<std::string, std::string> attributes;
	attributes["energy_east"] = print_float(energy_east_hop);
	attributes["energy_west"] = print_float(energy_west_hop);
	attributes["energy_north"] = print_float(energy_north_hop);
	attributes["energy_south"] = print_float(energy_south_hop);
	attributes["latency_east"] = print_float(energy_east_hop);
	attributes["latency_west"] = print_float(energy_west_hop);
	attributes["latency_north"] = print_float(energy_north_hop);
	attributes["latency_south"] = print_float(energy_south_hop);

	std::ostringstream ss;
	ss << "t " << name << print_format_attributes(attributes) << std::endl;
	return ss.str();
}

sanafe::Core::Core(
	const std::string &name, const int core_id, const int tile_id,
	const int core_offset) : name(name)
{
	id = core_id;
	parent_tile_id = tile_id;
	offset = core_offset;
	energy = 0.0;
	latency_after_last_message = 0.0;
}

sanafe::Tile &sanafe::Architecture::create_tile(
	const std::string &name,
	const std::map<std::string, std::string> &attr)

{
	tiles.push_back(Tile(name, tiles.size()));
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
	for (const auto &a: attr)
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

void sanafe::Architecture::load_arch_description(
	const std::filesystem::path &filename)
{
	std::ifstream arch_fp(filename);
	if (arch_fp.fail())
	{
		throw std::invalid_argument(
			"Error: Architecture file failed to open.");
	}
	int ret = description_parse_arch_file(arch_fp, *this);
	arch_fp.close();
	if (ret == sanafe::RET_FAIL)
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

sanafe::DendriteUnit::DendriteUnit(const std::string &dendrite_name):
	name(dendrite_name)
{
	energy = 0.0;
	time = 0.0;
	energy_access = 0.0;
	latency_access = 0.0;
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
	const std::string &name,
	const size_t tile_id,
	const std::map<std::string, std::string> &attr)
{
	if (tile_id > tiles.size())
	{
		throw std::invalid_argument("Error: Tile ID > total tiles");
	}
	Tile &tile = tiles_vec[tile_id];
	const int core_offset = tile.cores.size();
	tile.cores.push_back(
		Core(name, get_core_count(), tile.id, core_offset));
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

	max_cores_per_tile = std::max<size_t>(
		max_cores_per_tile, tile.cores.size());
	TRACE1("Core created id:%d.%d (tile:%d).\n", c.parent_tile_id, c.id);
	return c;
}

std::string sanafe::Core::info() const
{
	std::map<std::string, std::string> attributes;
	std::ostringstream ss;
	ss << "sanafe::Core(name= " << name << " tile=" << parent_tile_id << ")";
	return ss.str();
}

std::string sanafe::Core::description() const
{
	std::ostringstream ss;
	ss << "c " << name << ' ' << parent_tile_id << std::endl;
	return ss.str();
}

std::string sanafe::AxonInUnit::description() const
{
	std::map<std::string, std::string> attributes;
	attributes["energy_message"] = print_float(energy_spike_message);
	attributes["latency_message"] = print_float(latency_spike_message);
	std::ostringstream ss;
	ss << "i " << name << ' ' << parent_tile_id;
	ss << ' ' << parent_core_offset;
	ss << print_format_attributes(attributes) << std::endl;
	return ss.str();
}

std::string sanafe::SynapseUnit::description() const
{
	std::map<std::string, std::string> attributes;
	attributes["energy_spike"] = print_float(energy_spike_op);
	attributes["latency_spike"] = print_float(latency_spike_op);
	attributes["model"] = model;
	std::ostringstream ss;
	ss << "s " << name << ' ' << parent_tile_id << ' ' << parent_core_offset;
	ss << print_format_attributes(attributes) << std::endl;
	return ss.str();
}

std::string sanafe::DendriteUnit::description() const
{
	std::map<std::string, std::string> attributes;
	attributes["model"] = model;
	attributes["energy"] = print_float(energy_access);
	attributes["latency"] = print_float(latency_access);
	std::ostringstream ss;
	ss << "d " << name << ' ' << parent_tile_id;
	ss << ' ' << parent_core_offset;
	ss << print_format_attributes(attributes) << std::endl;
	return ss.str();
}

std::string sanafe::SomaUnit::description() const
{
	std::map<std::string, std::string> attributes;
	attributes["model"] = model;
	attributes["energy_access_neuron"] = print_float(energy_access_neuron);
	attributes["latency_access_neuron"] =
		print_float(latency_access_neuron);
	attributes["energy_update_neuron"] = print_float(energy_update_neuron);
	attributes["latency_update_neuron"] =
		print_float(latency_update_neuron);
	attributes["energy_spike_out"] = print_float(energy_access_neuron);
	attributes["latency_spike_out"] = print_float(latency_access_neuron);
	std::ostringstream ss;
	ss << "+ " << name << ' ' << parent_tile_id;
	ss << ' ' << parent_core_offset;
	ss << print_format_attributes(attributes) << std::endl;
	return ss.str();
}

std::string sanafe::AxonOutUnit::description() const
{
	std::map<std::string, std::string> attributes;
	attributes["energy"] = print_float(energy_access);
	attributes["latency"] = print_float(latency_access);
	std::ostringstream ss;
	ss << "o " << name << ' ' << parent_tile_id;
	ss << ' ' << parent_core_offset;
	ss << print_format_attributes(attributes) << std::endl;
	return ss.str();
}

sanafe::AxonInUnit &sanafe::Core::create_axon_in(
	const std::string &name,
	const std::map<std::string, std::string> &attr)
{
	axon_in_hw.push_back(AxonInUnit(name));
	struct AxonInUnit &in = axon_in_hw.back();

	in.energy = 0.0;
	in.time = 0.0;
	in.parent_tile_id = parent_tile_id;
	in.parent_core_offset = offset;

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
	const std::map<std::string, std::string> &attr)
{
	synapse.push_back(SynapseUnit(name));
	SynapseUnit &s = synapse.back();
	s.parent_tile_id = parent_tile_id;
	s.parent_core_offset = offset;

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
			s.model = value_str;
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

	TRACE1("Synapse processor created (c:%d.%d)\n",
		c.parent_tile_id, c.offset);

	return s;
}

sanafe::DendriteUnit &sanafe::Core::create_dendrite(
	const std::string &name,
	const std::map<std::string, std::string> &attr)
{
	dendrite.push_back(DendriteUnit(name));
	DendriteUnit &d = dendrite.back();
	d.parent_tile_id = parent_tile_id;
	d.parent_core_offset = offset;

	/**** Set attributes ****/
	for (const auto &a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;
		std::istringstream ss(value_str);
		if (key == "model")
		{
			d.model = value_str;
		}
		else if (key == "energy")
		{
			ss >> d.energy_access;
		}
		else if (key == "latency")
		{
			ss >> d.latency_access;
		}
	}
	TRACE1("Dendrite processor created (c:%d.%d)\n",
		c.parent_tile_id, c.offset);

	return d;
}

sanafe::SomaUnit &sanafe::Core::create_soma(
	const std::string &name,
	const std::map<std::string, std::string> &attr)
{
	TRACE1("cid:%d creating soma sid:%lu with %lu attributes\n",
		id, soma.size(), attr.size());
	soma.push_back(SomaUnit(name));
	SomaUnit &s = soma.back();
	s.parent_tile_id = parent_tile_id;
	s.parent_core_offset = offset;

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
		TRACE2("Soma attribute k:%s v%s\n", key.c_str(),
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
			/*
			s.noise_type = NOISE_FILE_STREAM;
			s.noise_stream = fopen(value_str.c_str(), "r");
			TRACE1("Opening noise str: %s\n", value_str.c_str());
			if (s.noise_stream == NULL)
			{
				INFO("Error: Failed to open noise stream: %s.\n",
					value_str.c_str());
				exit(1);
			}
			*/
		}
		else if (key == "model")
		{
			s.model = value_str;
		}
		else if (key == "plugin_lib")
		{
			s.plugin_lib = std::filesystem::path(value_str);
		}
	}

	TRACE1("Soma processor created (c:%d.%d)\n", c.parent_tile_id,
		c.offset);
	return s;
}

sanafe::AxonOutUnit &sanafe::Core::create_axon_out(
	const std::string &name,
	const std::map<std::string, std::string> &attr)
{
	axon_out_hw.push_back(AxonOutUnit(name));
	AxonOutUnit &out = axon_out_hw.back();
	out.parent_tile_id = parent_tile_id;
	out.parent_core_offset = offset;

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
	std::vector<bool> core_inserted;
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

	// Map neuron model to dendrite and soma hardware units in this core.
	//  Search through all models implemented by this core and return the
	//  one that matches. If no dendrite / soma hardware is specified,
	//  default to the first one defined
	if (dendrite.size() == 0)
	{
		INFO("Error: No dendrite units defined for cid:%d\n", id);
		throw std::runtime_error("Error: No dendrite units defined");
	}
	n.dendrite_hw = &(dendrite[0]);
	if (n.dendrite_hw_name.length() > 0)
	{
		bool dendrite_found = false;
		for (auto &dendrite_hw: dendrite)
		{
			if (n.dendrite_hw_name == dendrite_hw.name)
			{
				n.dendrite_hw = &dendrite_hw;
				dendrite_found = true;
			}
		}
		if (!dendrite_found)
		{
			INFO("Error: Could not map neuron nid:%d (hw:%s) "
				"to any dendrite h/w.\n", n.id,
				n.dendrite_hw_name.c_str());
			throw std::runtime_error(
				"Error: Could not map neuron to dendrite h/w");
		}
	}

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
	if (n.soma_model == nullptr)
	{
		// Setup the soma model
		TRACE1("Soma hw name: %s", soma_hw_name.c_str());
		assert(n.soma_hw != nullptr);
		if (n.soma_hw->plugin_lib.empty())
		{
			// Use built in models
			INFO("Creating soma built-in model %s.\n",
				n.soma_hw->model.c_str());
			n.soma_model = sanafe::model_get_soma(
				n.soma_hw->model, n.parent_group_id, n.id);
		}
		else
		{
			INFO("Creating soma from plugin %s.\n",
				n.soma_hw->plugin_lib.c_str());
			n.soma_model = plugin_get_soma(
				n.soma_hw->model, n.parent_group_id, n.id,
				n.soma_hw->plugin_lib);
		}
		const NeuronGroup &group = n.parent_net->groups_vec[
			n.parent_group_id];
		assert(n.soma_model != nullptr);
		// First set the group's default attribute values, and then
		//  any defined by the neuron
		n.soma_model->set_attributes(group.default_attributes);
		n.soma_model->set_attributes(n.attributes);
	}

	if (n.dendrite_model == nullptr)
	{
		// Setup the soma model
		TRACE1("Dendrite hw name: %s", dendrite_hw_name.c_str());
		if (n.dendrite_hw->plugin_lib.empty())
		{
			// Use built in models
			INFO("Creating dendrite built-in model %s.\n",
				n.dendrite_hw->model.c_str());
			n.dendrite_model = sanafe::model_get_dendrite(
				n.dendrite_hw->model);
		}
		else
		{
			INFO("Creating dendrite from plugin %s.\n",
				n.dendrite_hw->plugin_lib.c_str());
			n.dendrite_model = sanafe::plugin_get_dendrite(
				n.dendrite_hw->model,
				n.dendrite_hw->plugin_lib);
		}
		const NeuronGroup &group = n.parent_net->groups_vec[
			n.parent_group_id];
		assert(n.dendrite_model != nullptr);
		// Global attributes for all compartments in the dendrite
		n.dendrite_model->set_attributes(group.default_attributes);
		n.dendrite_model->set_attributes(n.attributes);
		// Set compartment specific attributes
		for (auto &compartment: n.dendrite_compartments)
		{
			n.dendrite_model->set_attributes(
				compartment.id, compartment.attributes);
		}
		// Set branch specific attributes
		for (auto &branch: n.dendrite_branches)
		{
			n.dendrite_model->set_attributes(
				branch.src_compartment_id,
				branch.dest_compartment_id,
				branch.attributes);
		}
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
	const size_t new_axon_in_address = post_core.axons_in.size() - 1;

	// Add the axon at the sending, pre-synaptic core
	TRACE1("axon in address:%d for core:%d.%d\n",
		new_axon_in_address, post_core.parent_tile_id, post_core.id);
	AxonOutModel out;
	out.dest_axon_id = new_axon_in_address;
	out.dest_core_offset = post_core.offset;
	out.dest_tile_id = post_core.parent_tile_id;
	out.src_neuron_id = pre_neuron.id;
	pre_core.axons_out.push_back(out);
	const size_t new_axon_out_address = pre_core.axons_out.size() - 1;

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
	con.synapse_hw = &(post_core.synapse[default_hw_id]);
	if (con.synapse_hw_name.length() > 0)
	{
		bool mapped = false;
		// Search for the specified synapse hardware
		for (auto &s: post_core.synapse)
		{
			if (con.synapse_hw_name == s.name)
			{
				con.synapse_hw = &s;
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

	return;
}

void sanafe::Architecture::save_arch_description(
	const std::filesystem::path &path)
{
	std::ofstream out(path);
	if (!out.is_open())
	{
		INFO("Error: Couldn't open arch file to save to: %s\n",
			path.c_str());
		throw std::invalid_argument(
			"Error: Couldn't open arch file to save to.");
	}

	assert(tiles_vec.size() == tiles.size());
	INFO("tiles vector size:%lu\n", tiles_vec.size());
	INFO("tiles size:%lu\n", tiles.size());

	for (const Tile &tile: tiles_vec)
	{
		out << tile.description();
		for (const Core &core: tile.cores)
		{
			out << core.description();
			for (const AxonInUnit &in: core.axon_in_hw)
			{
				out << in.description();
			}
			for (const SynapseUnit &s: core.synapse)
			{
				out << s.description();
			}
			for (const DendriteUnit &d: core.dendrite)
			{
				out << d.description();
			}
			for (const SomaUnit &s: core.soma)
			{
				out << s.description();
			}
			for (const AxonOutUnit &o: core.axon_out_hw)
			{
				out << o.description();
			}
		}
	}
	out << description();
}
