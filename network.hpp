// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.h - (spiking) neural network functionality. Spiking neural
//  networks are represented as groups of neurons. A neuron group might have a
//  bunch of neurons all with the same properties (and common hardware).
//  Each neuron has its own state and a set of connections to other neurons.
//  These structures have links to hardware for performance simulation.
//  Here we include different neuron, synapse and dendrite models.
#ifndef NETWORK_HEADER_INCLUDED_
#define NETWORK_HEADER_INCLUDED_

#include <cstdint>
#include <list>
#include <memory>
#include <unordered_map>
#include "plugins.hpp"
#include "models.hpp"

namespace sanafe
{
enum ConnectionConfigFormat
{
	// This is the structure of the CSV format for specifying synaptic
	//  connections in a network for this simulator.  Each row represents
	//  a unique connection.
	CONNECTION_DEST_GID = 0,
	CONNECTION_DEST_NID,
	CONNECTION_WEIGHT,
	CONNECTION_FIELDS,
};

// Forward declarations
struct Core;
struct Neuron;
struct NeuronGroup;
struct SomaUnit;
struct SomaModel;
struct SynapseUnit;
struct AxonOutUnit;

struct Connection
{
	Neuron *post_neuron, *pre_neuron;
	SynapseUnit *synapse_hw;
	std::string synapse_hw_name;
	double weight, current, synaptic_current_decay;
	int id, delay, last_updated;

	Connection(const int connection_id);
};

struct Neuron
{
	std::vector<Connection> connections_out;
	std::vector<int> axon_out_addresses;
	std::unordered_map<std::string, std::string> attributes;

	// Mapped hardware
	Network *parent_net;
	Core *core, *post_synaptic_cores;
	SomaUnit *soma_hw;
	AxonOutUnit *axon_out_hw;
	std::string soma_hw_name;

	std::shared_ptr<SomaModel> model;

	// Track the timestep each hardware unit was last updated
	bool is_init, fired, force_update, log_spikes, log_potential;
	bool update_needed;
	int id, parent_group_id;
	int spike_count;
	int soma_last_updated, dendrite_last_updated;
	int max_connections_out, maps_in_count, maps_out_count;

	double dendritic_current_decay, processing_latency;
	double current, charge;
	NeuronStatus neuron_status;
	int forced_spikes;

	Neuron(const size_t neuron_id);
	void set_attributes(const std::unordered_map<std::string, std::string> &attr);
	Connection &connect_to_neuron(Neuron &dest, const std::unordered_map<std::string, std::string> &attr);
};

class NeuronGroup
{
public:
	// A neuron group is a collection of neurons that share common
	//  parameters. All neurons must be based on the same neuron model.
	std::vector<Neuron> neurons;
	std::string default_soma_hw_name;
	std::string default_synapse_hw_name;
	std::unordered_map<std::string, std::string> default_attributes;

	int id;
	int default_max_connections_out;
	bool default_log_potential, default_log_spikes, default_force_update;

	NeuronGroup(const size_t group_id, const int neuron_count);
	Neuron &define_neuron(const size_t id, const std::unordered_map<std::string, std::string> &attr);
};

class Network
{
public:
	std::vector<NeuronGroup> groups;
	Network() {};
	NeuronGroup &create_neuron_group(const int neuron_count, const std::unordered_map<std::string, std::string> &attr);
private:
	Network(const Network &copy);
};

void network_check_mapped(Network &net);
}

#endif
