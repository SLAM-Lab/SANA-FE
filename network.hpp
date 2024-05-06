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

#define NETWORK_MAX_NEURON_GROUPS 1024
#define NETWORK_INVALID_NID -1
#define MAX_FIELDs 128
#define MAX_FIELD_LEN 64

#include <cstdint>
#include <list>
#include <memory>
#include "plugins.hpp"
#include "models.hpp"

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
	std::vector<Attribute> attributes;

	// Mapped hardware
	Core *core, *post_synaptic_cores;
	SomaUnit *soma_hw;
	AxonOutUnit *axon_out_hw;
	std::string soma_hw_name;

	std::shared_ptr<SomaModel> model;

	// Track the timestep each hardware unit was last updated
	int id, parent_group_id, is_init, fired;
	int max_connections_out;
	int force_update, spike_count;
	int soma_last_updated, dendrite_last_updated;
	int maps_in_count, maps_out_count;
	bool log_spikes, log_potential, update_needed;

	double dendritic_current_decay, processing_latency;
	double current, charge;
	sanafe::NeuronStatus neuron_status;
	int forced_spikes;
};

struct NeuronGroup
{
	// A neuron group is a collection of neurons that share common
	//  parameters. All neurons must be based on the same neuron model.
	std::vector<Neuron> neurons;
	std::string default_soma_hw_name;
	std::string default_synapse_hw_name;
	std::vector<Attribute> default_attributes;

	int id;
	int default_max_connections_out;
	bool default_log_potential, default_log_spikes, default_force_update;
};

struct Network
{
	std::vector<NeuronGroup> groups;
};

struct Architecture;
struct Core;

int network_create_neuron(Network &net, Neuron &n, const std::vector<Attribute> &attr);
int network_create_neuron_group(Network &net, const int neuron_count, const std::vector<Attribute> &attr);
//Neuron *network_id_to_neuron_ptr(Network *const net, const NeuronId id);
int network_connect_neurons(Connection &con, Neuron &src, Neuron &dest, const std::vector<Attribute> &attr);
void network_check_mapped(Network &net);



#endif
