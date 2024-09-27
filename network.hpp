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

#include <any>
#include <cstdint>
#include <filesystem>
#include <functional> // For std::reference_wrapper
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>
#include <variant>

#include "print.hpp"

namespace sanafe
{
// Forward declarations
// Architecture
class Architecture;
class Core;
struct SynapseUnit;
struct DendriteUnit;
struct SomaUnit;
struct AxonOutUnit;
// Network
class Network;
class NeuronGroup;
class Neuron;
struct Connection;
struct Synapse;
// Models
struct ModelParam;

struct NeuronTemplate;
enum NeuronStatus : int
{
    INVALID_NEURON_STATE,
    IDLE,
    UPDATED,
    FIRED
};

struct NeuronTemplate
{
    std::map<std::string, ModelParam> soma_model_params;
    std::map<std::string, ModelParam> dendrite_model_params;
    std::string soma_hw_name{""};
    std::string default_synapse_hw_name{""};
    std::string dendrite_hw_name{""};
    bool log_spikes{false};
    bool log_potential{false};
    bool force_synapse_update{false};
    bool force_dendrite_update{false};
    bool force_soma_update{false};
};

struct Conv2DParameters
{
    int input_width{};
    int input_height{};
    int input_channels{};
    int kernel_width{};
    int kernel_height{};
    int kernel_count{1};
    int stride_width{1};
    int stride_height{1};
};

class Neuron
{
public:
    std::vector<Connection> connections_out;
    std::vector<int> axon_out_addresses;
    std::map<std::string, ModelParam> soma_model_params;
    std::map<std::string, ModelParam> dendrite_model_params;

    // Mapped hardware
    Network &parent_net;
    Core *core{nullptr};
    Core *post_synaptic_cores{nullptr};
    DendriteUnit *dendrite_hw{nullptr};
    SomaUnit *soma_hw{nullptr};
    AxonOutUnit *axon_out_hw{nullptr};
    std::string soma_hw_name;
    std::string default_synapse_hw_name;
    std::string dendrite_hw_name;

    std::string parent_group_id;
    size_t id{};
    size_t mapped_address{};
    int spike_count{0};
    int maps_in_count{0};
    int maps_out_count{0};
    NeuronStatus status{sanafe::IDLE};
    bool force_synapse_update;
    bool force_dendrite_update;
    bool force_soma_update;
    bool log_spikes;
    bool log_potential;

    // Inputs to H/W units
    std::vector<Synapse> dendrite_input_synapses;
    double soma_input_charge{0.0};
    bool axon_out_input_spike{false};

    explicit Neuron(size_t neuron_id, Network &parent_net, const std::string parent_group_id, const NeuronTemplate &config);
    [[nodiscard]] size_t get_id() const { return id; }
    void set_attributes(const NeuronTemplate &attributes);
    Connection &connect_to_neuron(Neuron &dest);
    [[nodiscard]] std::string info() const;
    //std::string description(const bool write_mapping=true) const;
};

class NeuronGroup
{
public:
    // A neuron group is a collection of neurons that share common
    //  parameters. All neurons must be based on the same neuron model.
    std::vector<Neuron> neurons;
    NeuronTemplate default_neuron_config;
    Network &parent_net;
    std::string name;
    size_t order_created{0};
    size_t position_defined{0};
    [[nodiscard]] std::string get_id() const { return name; }
    explicit NeuronGroup(const std::string group_name, Network &parent_net, size_t neuron_count, const NeuronTemplate &default_config);

    //void set_attribute_multiple(const std::string &attr, const std::vector<std::any> &values);
    void connect_neurons_dense(NeuronGroup &dest_group, const std::map<std::string, std::vector<ModelParam>> &attribute_lists);
    void connect_neurons_sparse(NeuronGroup &dest_group, const std::map<std::string, std::vector<ModelParam>> &attribute_lists, const std::vector<std::pair<size_t, size_t> > &source_dest_id_pairs);
    void connect_neurons_conv2d(NeuronGroup &dest_group, const std::map<std::string, std::vector<ModelParam>> &attribute_lists, const Conv2DParameters &convolution);
    [[nodiscard]] std::string info() const;
    //std::string description() const;
};

class Network
{
public:
    std::unordered_map<std::string, NeuronGroup> groups;
    std::string name;
    bool record_attributes{true};

    explicit Network(std::string net_name = "") : name(std::move(net_name)) {};
    ~Network() = default;
    Network(Network &&) = default;
    Network &operator=(Network &&) = default;
    // Do *NOT* allow Network objects to be copied
    //  This is because Neuron objects link back to their parent Network
    //  (and need to be aware of the parent NeuronGroup). Linking to parent
    //  objects allows us to efficiently store Attributes for neurons i.e.,
    //  by avoiding duplication of shared attributes.
    //  If the Network was moved or copied, all parent links in Neurons
    //  would be invalid.
    Network(const Network &) = delete;
    Network &operator=(const Network &) = delete;

    NeuronGroup &create_neuron_group(const std::string name, size_t neuron_count, const NeuronTemplate &default_config);
    [[nodiscard]] std::string info() const;
    //void save_net_description(const std::filesystem::path &path, const bool save_mapping=true) const;
    void check_mapped() const;
};

Network load_net(const std::filesystem::path &path, Architecture &arch, bool use_netlist_format = false);

struct Connection
{
    std::map<std::string, ModelParam> synapse_params;
    std::map<std::string, ModelParam> dendrite_params;
    Neuron *post_neuron;
    Neuron *pre_neuron;
    SynapseUnit *synapse_hw;
    std::string synapse_hw_name;
    size_t synapse_address{0UL};
    int id;

    explicit Connection(int connection_id);
    //std::string description() const;
};

struct Synapse
{
    double current;
    Connection &con;
};

struct NeuronAddress
{
    std::string group_name{};
    std::optional<size_t> neuron_id{std::nullopt};
};

// Alternative ModelParameter implementations. Keep for now.. maybe find
//  somewhere better for this
// Simpler NeuronAttribute implementation that supports one level of lists
//  This is basically just a specialization of the variant class
//using NeuronAttribute =
//        std::variant<bool, int, double, std::string, std::vector<bool>,
//                std::vector<int>, std::vector<double>>;

/*
// More complex attribute / parameter implementation. However, this relies on
//  undefined behavior as std::map is being defined with an incomplete type.
//  This seems to work for GCC but isn't portable.
template<typename T>
class IncompleteTypeWrapper
{
public:
    IncompleteTypeWrapper(const T &t)
    {
        value = std::make_unique<T>(t);
    }
    IncompleteTypeWrapper(const IncompleteTypeWrapper &other)
    {
        value = std::make_unique<T>(*(other.value.get()));
    }
    ~IncompleteTypeWrapper() = default;

    auto operator== (const IncompleteTypeWrapper &other) const
    {
        return value == other.value;
    }
    auto operator!= (const IncompleteTypeWrapper &other) const
    {
        return value != other.value;
    }
    operator T() { return value.get(); }
private:
    std::unique_ptr<T> value;
};

template <typename Attribute>
using AttributeBase =
        std::variant<bool, int, double, std::string, std::vector<Attribute>,
                IncompleteTypeWrapper<std::map<std::string, Attribute>>>;

// The "using" keyword cannot be used recursively
struct AttributePrototype
{
    using type = AttributeBase<AttributePrototype>;
};
using NeuronAttribute3 = AttributePrototype::type;
*/

} // namespace

#endif
