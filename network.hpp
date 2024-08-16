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
class SynapseModel;
class SomaModel;
class DendriteModel;

struct NeuronTemplate;

enum NeuronStatus: int { INVALID_NEURON_STATE, IDLE, UPDATED, FIRED };

// An attribute can contain a scalar value, or either a list or named set of
//  attributes i.e., attributes are recursively defined attributes. However,
//  in C++, variants cannot be defined recursively.
struct ModelParam
{
    operator bool() const
    {
        if (std::holds_alternative<bool>(value))
        {
            return std::get<bool>(value);
        }
        else if (std::holds_alternative<int>(value))
        {
            TRACE1("Warning: Casting integer value to bool type.\n");
            return (std::get<int>(value) != 0);
        }

        std::string error = "Error: Attribute ";
        if (name.has_value())
        {
            error += name.value();
        }
        error += " cannot be cast to a bool";
        throw std::runtime_error(error);
    }
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
            TRACE1("Warning: Casting integer value to double type.\n");
            return static_cast<double>(std::get<int>(value));
        }

        std::string error = "Error: Attribute ";
        if (name.has_value())
        {
            error += name.value();
        }
        error += " cannot be cast to a double";
        throw std::runtime_error(error);
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
            cast_vector.push_back(static_cast<T>(element));
        }
        return cast_vector;
    }
    template <typename T>
    operator std::map<std::string, ModelParam>() const
    {
        std::map<std::string, ModelParam> cast_map;
        const auto &value_vector =
                std::get<std::vector<ModelParam>>(value);
        for (const auto &element : value_vector)
        {
            cast_map[element.name.value()] = static_cast<T>(element);
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
    bool log_spikes, log_potential, force_synapse_update, force_dendrite_update;
    bool force_soma_update;

    explicit NeuronTemplate(std::string soma_hw_name = "", std::string default_synapse_hw_name = "", std::string dendrite_hw_name = "", bool log_spikes = false, bool log_potential = false, bool force_synapse_update = false, bool force_dendrite_update = false, bool force_soma_update = false);
};

class Neuron
{
public:
    std::vector<Connection> connections_out;
    std::vector<int> axon_out_addresses;
    std::map<std::string, ModelParam> soma_model_params;
    std::map<std::string, ModelParam> dendrite_model_params;

    // Mapped hardware
    Network *parent_net{nullptr};
    Core *core{nullptr};
    Core *post_synaptic_cores{nullptr};
    DendriteUnit *dendrite_hw{nullptr};
    SomaUnit *soma_hw{nullptr};
    AxonOutUnit *axon_out_hw{nullptr};
    std::string soma_hw_name;
    std::string default_synapse_hw_name;
    std::string dendrite_hw_name;

    std::shared_ptr<SomaModel> soma_model{nullptr};
    std::shared_ptr<DendriteModel> dendrite_model{nullptr};

    std::string parent_group_id;
    size_t id{};
    int forced_spikes{0};
    int spike_count{0};
    int soma_last_updated{0};
    int dendrite_last_updated{0};
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
    Network *parent_net{nullptr};
    std::string name;
    size_t order_created{0};
    size_t position_defined{0};
    [[nodiscard]] std::string get_id() const { return name; }
    explicit NeuronGroup(const std::string group_name, Network &parent_net, size_t neuron_count, const NeuronTemplate &default_config);

    //void set_attribute_multiple(const std::string &attr, const std::vector<std::any> &values);
    //void connect_neurons(NeuronGroup &dest_group, const std::vector<std::pair<int, int> > &src_dest_id_pairs, const std::map<std::string, std::vector<std::any>> &attr_lists);
    [[nodiscard]] std::string info() const;
    //std::string description() const;
};

class Network
{
public:
    // Use a vector of dynamically allocated groups, so we get vectors random
    //  access, but do not reallocate objects when growing the vector
    std::map<std::string, NeuronGroup> groups;
    std::string name;
    explicit Network(std::string net_name) : name(std::move(net_name)) {};
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

Network load_net(const std::filesystem::path &path, Architecture &arch);

struct Connection
{
    std::map<std::string, ModelParam> synapse_params;
    std::map<std::string, ModelParam> dendrite_params;
    std::shared_ptr<SynapseModel> synapse_model;
    Neuron *post_neuron;
    Neuron *pre_neuron;
    SynapseUnit *synapse_hw;
    std::string synapse_hw_name;
    size_t synapse_address{0UL};
    int id;
    int last_updated{0};

    explicit Connection(int connection_id);
    //std::string description() const;
};

struct Synapse
{
    double current;
    // TODO: prepare the dendrite parameters into an object? would be much
    //  more efficient. This can be an any pointer? Or to a base dendrite class
    //std::map<std::string, ModelParam> dendrite_params;
    Connection &con;
};

struct NeuronAddress
{
    std::string group_name{};
    size_t neuron_id{};
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
