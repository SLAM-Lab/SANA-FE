// Copyright (c) 2025 - The University of Texas at Austin
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

#include "param.hpp"
#include "print.hpp"
#include "chip.hpp"

namespace sanafe
{
// Forward declarations
// Architecture
class Architecture;
struct CoreConfiguration;
// Network
class SpikingNetwork;
class NeuronGroup;
class Neuron;
struct Connection;
struct Synapse;
struct NeuronAddress;
// Models
struct ModelAttribute;

struct NeuronConfiguration
{
    std::map<std::string, ModelAttribute> model_attributes{};
    std::optional<std::string> soma_hw_name{};
    std::optional<std::string> default_synapse_hw_name{};
    std::optional<std::string> dendrite_hw_name{};
    std::optional<bool> log_spikes{};
    std::optional<bool> log_potential{};
    std::optional<bool> force_synapse_update{};
    std::optional<bool> force_dendrite_update{};
    std::optional<bool> force_soma_update{};
};

struct NeuronAddress
{
    std::string group_name{};
    std::optional<size_t> neuron_offset{std::nullopt};
    [[nodiscard]] std::string info() const;
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
    std::vector<Connection> edges_out;
    std::map<std::string, ModelAttribute> model_attributes;
    std::string soma_hw_name{};
    std::string default_synapse_hw_name{};
    std::string dendrite_hw_name{};
    std::string parent_group_name;
    SpikingNetwork &parent_net;
    size_t offset{};
    std::optional<CoreAddress> core_address{std::nullopt};
    size_t mapping_order{};
    // Optionally set flags for updating and traces
    bool force_synapse_update{false};
    bool force_dendrite_update{false};
    bool force_soma_update{false};
    bool log_spikes{false};
    bool log_potential{false};

    explicit Neuron(size_t neuron_offset, SpikingNetwork &net, const std::string parent_group_id, const NeuronConfiguration &config);
    [[nodiscard]] size_t get_id() const { return offset; }
    size_t connect_to_neuron(Neuron &dest);
    void map_to_core(const CoreConfiguration &core);
    void set_attributes(const NeuronConfiguration &attributes);
    [[nodiscard]] std::string info() const;
};

class NeuronGroup
{
public:
    // A neuron group is a collection of neurons that share common attributes
    std::vector<Neuron> neurons;
    NeuronConfiguration default_neuron_config;
    std::string name;

    [[nodiscard]] std::string get_name() const { return name; }
    explicit NeuronGroup(const std::string group_name, SpikingNetwork &net, size_t neuron_count, const NeuronConfiguration &default_config);

    void connect_neurons_dense(NeuronGroup &dest_group, const std::map<std::string, std::vector<ModelAttribute>> &attribute_lists);
    void connect_neurons_sparse(NeuronGroup &dest_group, const std::map<std::string, std::vector<ModelAttribute>> &attribute_lists, const std::vector<std::pair<size_t, size_t> > &source_dest_id_pairs);
    void connect_neurons_conv2d(NeuronGroup &dest_group, const std::map<std::string, std::vector<ModelAttribute>> &attribute_lists, const Conv2DParameters &convolution);
    [[nodiscard]] std::string info() const;
};

class SpikingNetwork
{
public:
    std::map<std::string, NeuronGroup> groups;
    std::string name;

    explicit SpikingNetwork(std::string net_name = "") : name(std::move(net_name)) {};
    ~SpikingNetwork() = default;
    SpikingNetwork(SpikingNetwork &&) = default;
    SpikingNetwork &operator=(SpikingNetwork &&) = default;
    // Do *NOT* allow Network objects to be copied
    //  This is because Neuron objects link back to their parent Network
    //  (and need to be aware of the parent NeuronGroup). Linking to parent
    //  objects allows us to efficiently store Attributes for neurons i.e.,
    //  by avoiding duplication of shared attributes.
    //  If the Network was moved or copied, all parent links in Neurons
    //  would be invalid.
    SpikingNetwork(const SpikingNetwork &) = delete;
    SpikingNetwork &operator=(const SpikingNetwork &) = delete;

    NeuronGroup &create_neuron_group(const std::string name, size_t neuron_count, const NeuronConfiguration &default_config);
    [[nodiscard]] std::string info() const;
    void save(const std::filesystem::path &path, const bool use_netlist_format = false) const;

    size_t update_mapping_count();
private:
    void save_netlist(const std::filesystem::path &path) const;
    void save_yaml(const std::filesystem::path &path) const;

    size_t mapping_count{0};
};

struct Connection
{
    std::map<std::string, ModelAttribute> synapse_attributes;
    std::map<std::string, ModelAttribute> dendrite_attributes;
    std::string synapse_hw_name{};
    NeuronAddress pre_neuron{};
    NeuronAddress post_neuron{};
    int id;

    Connection(size_t id) : id(id) {}
    [[nodiscard]] std::string info() const;
};

SpikingNetwork load_net(const std::filesystem::path &path, Architecture &arch, bool use_netlist_format = false);

} // namespace

#endif
