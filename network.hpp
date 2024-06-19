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
#include <variant>

#include "print.hpp"

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
// Models
class SynapseModel;
class SomaModel;
class DendriteModel;

struct NeuronTemplate;
enum NeuronStatus: int;

class Network
{
public:
    // Use a vector of dynamically allocated groups, so we get vectors random
    //  access, but do not reallocate objects when growing the vector
    std::vector<std::unique_ptr<NeuronGroup>> groups;
    std::string name;
    explicit Network(const std::string &net_name) : name(net_name) {};
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

    NeuronGroup &create_neuron_group(const std::string &name, const size_t neuron_count, const NeuronTemplate &default_config);
    std::string info() const;
    //void save_net_description(const std::filesystem::path &path, const bool save_mapping=true) const;
    void check_mapped() const;
};

Network load_net(const std::filesystem::path &path, Architecture &arch);

// An attribute can contain a scalar value, or either a list or named set of
//  attributes i.e., attributes are recursively defined attributes. However,
//  in C++, variants cannot be defined recursively. Use metaprogramming to
//  define the attribute type, with some C++ magic. Note that vectors *can* be
//  created with incomplete types, whereas maps *cannot*.
struct ModelParam
{
    operator bool() const { return std::get<bool>(value); }
    operator int() const { return std::get<int>(value); }
    operator double() const
    {
        if (std::holds_alternative<double>(value))
        {
            return std::get<double>(value);
        }
        else if (std::holds_alternative<int>(value))
        {
            // Assume it is safe to convert from any integer to double
            INFO("Warning: Casting integer value to double type.\n");
            return static_cast<double>(std::get<int>(value));
        }
        else
        {
            std::string error = "Error: Attribute ";
            if (name.has_value())
            {
                error += name.value();
            }
            error += " cannot be cast to a double";
            throw std::runtime_error(error);
        }
    }
    operator std::string() const { return std::get<std::string>(value); }
    template <typename T> operator std::vector<T>() const
    {
        std::vector<T> cast_vector;
        const auto &value_vector =
                std::get<std::vector<ModelParam>>(value);
        cast_vector.reserve(value_vector.size());

        for (const auto &element : value_vector)
        {
            cast_vector.push_back(std::get<T>(element.value));
        }
        return cast_vector;
    }
    template <typename T>
    operator std::map<std::string, ModelParam>() const
    {
        std::map<std::string, ModelParam> cast_map;
        const auto value_vector =
                std::get<std::vector<ModelParam>>(value);
        for (const auto &element : value_vector)
        {
            cast_map[element.name.value()] = static_cast<T>(element.value);
        }
        return cast_map;
    }
    bool operator==(const ModelParam &rhs) const
    {
        return (value == rhs.value);
    }

    std::variant<bool, int, double, std::string, std::vector<ModelParam>>
            value;
    std::optional<std::string> name;
    // In C++17, we cannot use std::map (which would be the natural choice) with
    //  incomplete types i.e., cannot use std::map in such a recursive
    //  structure. Considering this, and the fact that performance is not as
    //  important for this struct, label every attribute with a name and if the
    //  user wants to use "map" style lookups e.g., foo = attribute["key"]
    //  then support casting the struct to a std::map.
    //  There have been other discussions on this topic e.g., for implementing
    //  JSON and YAML parsers, but they end up either requiring Boost or other
    //  dependencies, and / or rely on undefined C++ behavior and generally
    //  require complex solutions.
};

struct NeuronTemplate
{
    std::map<std::string, ModelParam> soma_model_params;
    std::map<std::string, ModelParam> dendrite_model_params;
    std::string soma_hw_name, default_synapse_hw_name, dendrite_hw_name;
    bool log_spikes, log_potential, force_update;

    explicit NeuronTemplate(const std::string &soma_hw_name = "", const std::string &default_synapse_hw_name = "", const std::string &dendrite_hw_name = "", const bool log_spikes = false, const bool log_potential = false, const bool force_update = false);
};

class NeuronGroup
{
public:
    // A neuron group is a collection of neurons that share common
    //  parameters. All neurons must be based on the same neuron model.
    std::vector<Neuron> neurons;
    NeuronTemplate default_neuron_config;
    std::string name;
    int id;
    int get_id() { return id; }
    explicit NeuronGroup(const std::string &group_name, const size_t group_id, const size_t neuron_count, const NeuronTemplate &default_config);
    //void set_attribute_multiple(const std::string &attr, const std::vector<std::any> &values);
    //void connect_neurons(NeuronGroup &dest_group, const std::vector<std::pair<int, int> > &src_dest_id_pairs, const std::map<std::string, std::vector<std::any>> &attr_lists);
    std::string info() const;
    //std::string description() const;
};

class Neuron
{
public:
    std::vector<Connection> connections_out;
    std::vector<int> axon_out_addresses;
    std::map<std::string, ModelParam> soma_model_params;
    std::map<std::string, ModelParam> dendrite_model_params;

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
    int maps_in_count, maps_out_count;
    NeuronStatus status;

    // Inputs to H/W units
    std::vector<Synapse> dendrite_input_synapses;
    double soma_input_charge;
    bool axon_out_input_spike;

    explicit Neuron(const size_t neuron_id, const size_t parent_group_id, const NeuronTemplate &config);
    int get_id() { return id; }
    void set_attributes(const NeuronTemplate &attributes);
    void connect_to_neuron(Neuron &dest, const std::map<std::string, ModelParam> &synapse_params, const std::map<std::string, ModelParam> &dendrite_params,  const std::optional<std::string> &synapse_hw_name = std::nullopt);
    std::string info() const;
    //std::string description(const bool write_mapping=true) const;
};

struct Connection
{
    std::map<std::string, ModelParam> synapse_params;
    std::map<std::string, ModelParam> dendrite_params;
    std::shared_ptr<SynapseModel> synapse_model;
    Neuron *post_neuron, *pre_neuron;
    SynapseUnit *synapse_hw;
    std::string synapse_hw_name;
    int id, delay, last_updated;

    explicit Connection(const int connection_id);
    //std::string description() const;
};

struct Synapse
{
    double current;
    // TODO: prepare the dendrite parameters into an object? would be much
    //  more efficient. This can be an any pointer? Or to a base dendrite class
    std::map<std::string, ModelParam> dendrite_params;
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
