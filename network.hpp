// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.hpp - (spiking) neural network functionality. Spiking neural
//  networks are represented as groups of neurons. A neuron group might have a
//  bunch of neurons all with the same properties (and common hardware).
//  Each neuron has its own state and a set of connections to other neurons.
//  These structures have links to hardware for performance simulation.
//  Here we include different neuron, synapse and dendrite models.
#ifndef NETWORK_HEADER_INCLUDED_
#define NETWORK_HEADER_INCLUDED_

#include <cstdint>
#include <filesystem>
#include <functional> // For std::reference_wrapper
#include <list>
#include <map>
#include <memory>
#include <optional>

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
// Architecture
class Architecture;
class Core;
struct SynapseUnit;
struct DendriteUnit;
struct SomaUnit;
struct AxonOutUnit;
// Network
class NeuronGroup;
class Neuron;
struct Connection;
struct Synapse;
struct Compartment;
struct Branch;
// Models
class SynapseModel;
class SomaModel;
class DendriteModel;

enum NeuronStatus: int;

class Network
{
public:
    std::list<NeuronGroup> groups;
    std::vector<std::reference_wrapper<NeuronGroup> > groups_vec;
    Network() {};
    NeuronGroup &create_neuron_group(const int neuron_count, const std::map<std::string, std::string> &attr);
    void load_net_description(const std::string &filename, Architecture &arch);
    std::string info() const;
    void save_net_description(const std::filesystem::path &path, const bool save_mapping=true) const;
    void check_mapped() const;

private:
    // Do *NOT* allow Network objects to be copied
    //  This is because Neuron objects link back to their parent Network
    //  (and need to be aware of the parent NeuronGroup). Linking to parent
    //  objects allows us to efficiently store Attributes for neurons i.e.,
    //  by avoiding duplication of shared attributes.
    //  If the Network was moved or copied, all parent links in Neurons
    //  would be invalid.
    Network(const Network &copy);
};

class NeuronGroup
{
public:
    // A neuron group is a collection of neurons that share common
    //  parameters. All neurons must be based on the same neuron model.
    std::vector<Neuron> neurons;
    std::string default_soma_hw_name;
    std::string default_synapse_hw_name;
    std::map<std::string, std::string> default_attributes;

    std::filesystem::path default_soma_plugin;
    int id;
    int default_max_connections_out, default_max_compartments;
    bool default_log_potential, default_log_spikes, default_force_update;

    int get_id() { return id; }
    NeuronGroup(const size_t group_id, const int neuron_count);
    void set_attribute_multiple(const std::string &attr, const std::vector<std::string> &values);
    void connect_neurons(NeuronGroup &dest_group, const std::vector<std::pair<int, int> > &src_dest_id_pairs, const std::map<std::string, std::vector<std::string> > &attr_lists);
    std::string info() const;
    std::string description() const;
};

class Neuron
{
public:
    std::vector<Connection> connections_out;
    std::vector<Compartment> dendrite_compartments;
    std::vector<Branch> dendrite_branches;
    std::vector<int> axon_out_addresses;
    std::map<std::string, std::string> attributes;

    // Mapped hardware
    Network *parent_net;
    Core *core, *post_synaptic_cores;
    DendriteUnit *dendrite_hw;
    SomaUnit *soma_hw;
    AxonOutUnit *axon_out_hw;
    std::string soma_hw_name, default_synapse_hw_name, dendrite_hw_name;

    std::shared_ptr<SomaModel> soma_model;
    std::shared_ptr<DendriteModel> dendrite_model;

    bool force_update, log_spikes, log_potential;
    int id, parent_group_id, forced_spikes, spike_count;
    int soma_last_updated, dendrite_last_updated;
    int max_connections_out, max_compartments;
    int maps_in_count, maps_out_count;
    double processing_latency;
    NeuronStatus status;

    // Inputs to H/W units
    std::vector<Synapse> dendrite_input_synapses;
    double soma_input_charge;
    bool axon_out_input_fired;

    explicit Neuron(const size_t neuron_id);
    int get_id() { return id; }
    void set_attributes(const std::map<std::string, std::string> &attr);
    void connect_to_neuron(Neuron &dest, const size_t dest_compartment_id, const std::map<std::string, std::string> &attr);
    void create_compartment(const std::map<std::string, std::string> &compartment_attr);
    void create_branch(const size_t src_compartment_id, const size_t dest_compartment_id, const std::map<std::string, std::string> &branch_attr);
    std::string info() const;
    std::string description(const bool write_mapping=true) const;
};

struct Connection
{
    std::map<std::string, std::string> attributes;
    std::shared_ptr<SynapseModel> synapse_model;
    Neuron *post_neuron, *pre_neuron;
    SynapseUnit *synapse_hw;
    std::string synapse_hw_name;
    size_t dest_compartment;
    int id, delay, last_updated;

    explicit Connection(const int connection_id);
    std::string description() const;
};

struct Synapse
{
    double current;
    int dest_compartment;
};

struct Compartment
{
    std::map<std::string, std::string> attributes;
    size_t id, parent_neuron_id;
    Compartment(const size_t compartment_id) : id(compartment_id) {}
};

struct Branch
{
    size_t id, src_compartment_id, dest_compartment_id;
    std::map<std::string, std::string> attributes;
    Branch(const size_t branch_id) : id(branch_id) {}
};

} // namespace

#endif
