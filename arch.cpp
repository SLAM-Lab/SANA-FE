// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.c
#include <cctype>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <list>
#include <set>
#include <iostream>
#include <sstream>

#include "print.hpp"
#include "arch.hpp"
#include "network.hpp"
#include "description.hpp"

Architecture::Architecture()
{
	noc_buffer_size = 0;
	noc_init = false;

	// TODO: remove these
	spike_vector_on = false;
	return;
}

int Architecture::get_core_count()
{
	int core_count = 0;
	for (auto &tile: tiles)
	{
		core_count += tile.cores.size();
	}

	return core_count;
}

Message::Message()
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
	dest_core_id = -1;
	dest_tile_id = -1;
	dest_axon_id = -1;
}

int arch_create_noc(struct Architecture &arch, const std::list<Attribute> &attr)
{
	const int tile_count = arch.tiles.size();
	if (tile_count <= 0)
	{
		// The NoC interconnect is built after tiles are all defined
		//  This is because we link the tiles together in the NoC mesh
		INFO("Error: NoC must be built after tiles defined.\n");
		exit(1);
	}

	// Default values
	arch.noc_width = 1;
	arch.noc_height = 1;
	arch.noc_buffer_size = 0;

	for (auto a: attr)
	{
		std::istringstream ss(a.value_str);
		if (a.key == "width")
		{
			ss >> arch.noc_width;
		}
		else if (a.key == "height")
		{
			ss >> arch.noc_height;
		}
		else if (a.key == "link_buffer_size")
		{
			ss >> arch.noc_buffer_size;
		}
	}
	assert((arch.noc_height * arch.noc_width) <= ARCH_MAX_TILES);
	arch.noc_init = 1;
	TRACE1("NoC created, mesh, width:%d height:%d.\n", arch.noc_width,
		arch.noc_height);
	return 0;
}

int arch_create_tile(Architecture &arch, const std::list<Attribute> &attr)

{
	Tile tile;
	const int tile_count = arch.tiles.size();

	tile.id = tile_count;
	tile.energy = 0.0;
	tile.x = 0;
	tile.y = 0;

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
		std::istringstream ss(a.value_str);
		if (a.key == "energy_east")
		{
			ss >> tile.energy_east_hop;
		}
		else if (a.key == "latency_east")
		{
			ss >> tile.latency_east_hop;
		}
		else if (a.key == "energy_west")
		{
			ss >> tile.energy_west_hop;
		}
		else if (a.key == "latency_west")
		{
			ss >> tile.latency_west_hop;
		}
		else if (a.key == "energy_north")
		{
			ss >> tile.energy_north_hop;
		}
		else if (a.key == "latency_north")
		{
			ss >> tile.latency_north_hop;
		}
		else if (a.key == "energy_south")
		{
			ss >> tile.energy_south_hop;
		}
		else if (a.key == "latency_south")
		{
			ss >> tile.latency_south_hop;
		}
	}

	arch.tiles.push_back(tile);
	return tile.id;
}

int arch_create_core(
	Architecture &arch, Tile &tile, const std::list<Attribute> &attr)
{
	const int core_offset = tile.cores.size();
	Core c;
	c.offset = core_offset;
	c.id = arch.get_core_count();
	c.parent_tile_id = tile.id;
	c.max_neurons = 1024;

	// *** Set attributes ***
	c.buffer_pos = BUFFER_SOMA;
	for (auto a: attr)
	{
		if (a.key == "buffer_before")
		{
			if (a.value_str == "soma")
			{
				c.buffer_pos = BUFFER_SOMA;
			}
		}
		else if (a.key == "max_neurons")
		{
			std::istringstream ss(a.value_str);
			ss >> c.max_neurons;
		}
	}

	// Initialize core state
	c.soma_count = 0;
	c.synapse_count = 0;
	c.energy = 0.0;

	c.next_message = Message();
	tile.cores.push_back(c);

	TRACE1("Core created id:%d.%d (tile:%d).\n", c.parent_tile_id, .id);
	return c.id;
}

void arch_create_axon_in(
	Core &c, const std::string &name, const std::list<Attribute> &attr)
{
	struct AxonInUnit in;
	in.energy = 0.0;
	in.time = 0.0;
	in.parent_tile_id = c.parent_tile_id;

	in.energy_spike_message = 0.0;
	in.latency_spike_message = 0.0;
	for (auto curr: attr)
	{
		std::istringstream ss(curr.value_str);
		if (curr.key == "name")
		{
			in.name = curr.value_str;
		}
		else if (curr.key == "energy_message")
		{
			ss >> in.energy_spike_message;
		}
		else if (curr.key == "latency_message")
		{
			ss >> in.latency_spike_message;
		}
	}

	c.axon_in_hw.push_back(in);
	TRACE2("Axon input created (c:%d.%d)\n", c->t->id, c->id);

	return;
}

void arch_create_synapse(struct Core &c, const std::string &name,
	const std::list<Attribute> &attr)
{
	SynapseUnit s;
	s.name = name;
	s.energy = 0.0;
	s.time = 0.0;

	/**** Set attributes ****/
	s.energy_memory_access = 0.0;
	s.latency_memory_access = 0.0;
	s.energy_spike_op = 0.0;
	s.latency_spike_op = 0.0;
	s.weight_bits = 8;
	for (auto curr: attr)
	{
		std::istringstream ss(curr.value_str);
		if (curr.key == "name")
		{
			s.name = curr.value_str;
		}
		else if (curr.key == "model")
		{
			s.model = arch_parse_synapse_model(curr.value_str);
		}
		else if (curr.key == "energy_memory")
		{
			ss >> s.energy_memory_access;
		}
		else if (curr.key == "latency_memory")
		{
			ss >> s.latency_memory_access;
		}
		else if (curr.key == "energy_spike")
		{
			ss >> s.energy_spike_op;
		}
		else if (curr.key == "latency_spike")
		{
			ss >> s.latency_spike_op;
		}
	}
	c.synapse.push_back(s);

	TRACE1("Synapse processor created (c:%d.%d)\n", c.parent_tile_id, c.id);

	return;
}

void arch_create_soma(struct Core &c, const std::string &name,
	const std::list<Attribute> &attr)
{
	SomaUnit s;
	s.name = name;
	s.neuron_updates = 0;
	s.neurons_fired = 0;
	s.neuron_count = 0;
	s.energy = 0.0;
	s.time = 0.0;

	/*** Set attributes ***/
	s.energy_access_neuron = 0.0;
	s.latency_access_neuron = 0.0;
	s.energy_update_neuron = 0.0;
	s.latency_update_neuron = 0.0;
	s.energy_spiking = 0.0;
	s.latency_spiking = 0.0;
	s.leak_towards_zero = 1;
	s.noise_type = NOISE_NONE;

	for (auto a: attr)
	{
		std::istringstream ss(a.value_str);
		if (a.key == "energy_update_neuron")
		{
			ss >> s.energy_update_neuron;
		}
		else if (a.key == "latency_update_neuron")
		{
			ss >> s.latency_update_neuron;
		}
		else if (a.key == "energy_access_neuron")
		{
			ss >> s.energy_access_neuron;
		}
		else if (a.key == "latency_access_neuron")
		{
			ss >> s.latency_access_neuron;
		}
		else if (a.key == "energy_spike_out")
		{
			ss >> s.energy_spiking;
		}
		else if (a.key == "latency_spike_out")
		{
			ss >> s.latency_spiking;
		}
		else if (a.key == "noise")
		{
			s.noise_type = NOISE_FILE_STREAM;
			s.noise_stream = fopen(a.value_str.c_str(), "r");
			TRACE1("Opening noise str: %s\n", a.value_str.c_str());
			if (s.noise_stream == NULL)
			{
				INFO("Error: Failed to open noise stream: %s.\n",
					a.value_str.c_str());
				exit(1);
			}
		}
	}
	c.soma.push_back(s);

	TRACE1("Soma processor created (c:%d.%d)\n", c.parent_tile_id,
		c.offset);
	return;
}

void arch_create_axon_out(struct Core &c, const std::list<Attribute> &attr)
{
	AxonOutUnit out;
	out.packets_out = 0;
	out.energy = 0.0;
	out.time = 0.0;

	/*** Set attributes ***/
	out.energy_access = 0.0;
	out.latency_access = 0.0;
	for (auto curr: attr)
	{
		std::istringstream ss(curr.value_str);
		if (curr.key == "energy")
		{
			ss >> out.energy_access;
		}
		else if (curr.key == "latency")
		{
			ss >> out.latency_access;
		}
	}

	// Track the tile the axon interfaces with
	out.parent_tile_id = c.parent_tile_id;
	c.axon_out_hw.push_back(out);
	TRACE1("Axon output created (c:%d.%d)\n", c.parent_tile_id, c.offset);

	return;
}

void arch_create_axons(struct Architecture &arch)
{
	TRACE1("Creating all connection maps.\n");
	for (auto &tile: arch.tiles)
	{
		for (auto &c: tile.cores)
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

void arch_print_axon_summary(struct Architecture &arch)
{
	int in_count, out_count, core_used;
	in_count = 0;
	out_count = 0;

	INFO("** Mapping summary **\n");
	for (auto &tile: arch.tiles)
	{
		// For debug only, print the axon maps
		for (auto &c: tile.cores)
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

void arch_map_neuron_connections(Neuron &pre_neuron)
{
	// Setup the connections between neurons and map them to hardware
	assert(pre_neuron.core != nullptr);

	// Figure out the unique set of cores that this neuron broadcasts to
	INFO("Counting connections for neuron nid:%d\n", pre_neuron.id);
	std::set<Core *> cores_out;
	for (Connection &curr_connection: pre_neuron.connections_out)
	{
		INFO("Looking at connection id: %d\n", curr_connection.id);
		Core *dest_core = curr_connection.post_neuron->core;
		cores_out.insert(dest_core);
		INFO("Connected to dest core: %d\n", dest_core->id);
	}

	INFO("Creating connections for neuron nid:%d to %lu core(s)\n",
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
		INFO("Adding connection:%d\n", curr_connection.id);
		arch_add_connection_to_axon(curr_connection, post_core);
	}
	INFO("Finished mapping connections to hardware for nid:%d.%d.\n",
		pre_neuron.parent_group_id, pre_neuron.id);

	return;
}

int arch_map_neuron(Neuron &n, Core &c)
{
	// Map the neuron to hardware units
	assert(n.core == NULL);

	n.core = &c;
	TRACE1("Mapping neuron %d to core %d\n", n.id, c->id);
	c.neurons.push_back(&n);

	if (c.neurons.size() > c.max_neurons)
	{
		INFO("Error: Exceeded maximum neurons per core (%lu)",
			c.max_neurons);
		exit(1);
	}

	// Map neuron model to soma hardware unit in this core. Search through
	//  all neuron models implemented by this core and return the one that
	//  matches. If no soma hardware is specified, default to the first
	//  one defined
	n.soma_hw = &(c.soma[0]);
	if (n.soma_hw_name.length() > 0)
	{
		bool soma_found = false;
		for (auto &soma_hw: c.soma)
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
			exit(1);
		}
	}
	n.soma_hw->neuron_count++;

	// TODO: support multiple axon outputs
	n.axon_out_hw = &(c.axon_out_hw[0]);

	return RET_OK;
}

void arch_allocate_axon(Neuron &pre_neuron, Core &post_core)
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
	INFO("axon in address:%d for core:%d.%d\n",
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
	INFO("nid:%d.%d cid:%d.%d added one output axon address %d.\n",
		pre_neuron.parent_group_id, pre_neuron.id,
		pre_core.parent_tile_id, pre_core.offset,
		new_axon_out_address);

	return;
}

void arch_add_connection_to_axon(Connection &con, Core &post_core)
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
		for (auto &synapse_hw: post_core.synapse)
		{
			if (con.synapse_hw_name == synapse_hw.name)
			{
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

int arch_parse_synapse_model(const std::string &model_str)
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
