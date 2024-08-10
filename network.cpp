// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.cpp
#include <any>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem> // For std::filesystem::path
#include <fstream>
#include <functional> // For std::reference_wrapper
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <sstream>

#include "arch.hpp"
#include "description.hpp"
#include "models.hpp"
#include "network.hpp"
#include "print.hpp"

sanafe::NeuronTemplate::NeuronTemplate(std::string soma_hw_name,
        std::string default_synapse_hw_name,
        std::string dendrite_hw_name, const bool log_spikes,
        const bool log_potential, const bool force_update)
        : soma_hw_name(std::move(soma_hw_name))
        , default_synapse_hw_name(std::move(default_synapse_hw_name))
        , dendrite_hw_name(std::move(dendrite_hw_name))
        , log_spikes(log_spikes)
        , log_potential(log_potential)
        , force_update(force_update)
{
}

sanafe::Connection::Connection(const int connection_id)
        : post_neuron(nullptr)
        , pre_neuron(nullptr)
        , synapse_hw(nullptr)
        , id(connection_id)
{
}

/*
std::string sanafe::Connection::description() const
{
    assert(pre_neuron != nullptr);
    assert(post_neuron != nullptr);

    std::ostringstream ss;
    ss << pre_neuron->parent_group_id << '.' << pre_neuron->id;
    ss << "->";
    ss << post_neuron->parent_group_id << '.' << post_neuron->id;
    ss << print_format_attributes(attributes) << std::endl;
    ;
    return ss.str();
}
*/

sanafe::NeuronGroup::NeuronGroup(std::string group_name,
        Network &parent_net, const NeuronTemplate &default_config)
        : default_neuron_config(default_config)
        , parent_net(&parent_net)
        , name(std::move(group_name))
{
}

sanafe::Neuron::Neuron(std::string neuron_id,
        Network &parent_net, std::string parent_group_id,
        const NeuronTemplate &config)
        : parent_net(&parent_net)
        , soma_hw_name(config.soma_hw_name)
        , default_synapse_hw_name(config.default_synapse_hw_name)
        , dendrite_hw_name(config.dendrite_hw_name)
        , id(std::move(neuron_id))
        , parent_group_id(std::move(parent_group_id))
        , force_update(config.force_update)
        , log_spikes(config.log_spikes)
        , log_potential(config.log_potential)
{
    soma_hw_name = config.soma_hw_name;
    default_synapse_hw_name = config.default_synapse_hw_name;
    dendrite_hw_name = config.dendrite_hw_name;

    soma_model_params.insert(
            config.soma_model_params.begin(), config.soma_model_params.end());
    dendrite_model_params.insert(config.dendrite_model_params.begin(),
            config.dendrite_model_params.end());
}

std::string sanafe::Neuron::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Neuron(nid=" << parent_group_id << '.' << id;
    ss << " connections_out=" << connections_out.size();
    //ss << " attributes={" << print_format_attributes(attributes) << "})";
    return ss.str();
}

/*
std::string sanafe::Neuron::description(const bool write_mapping) const
{
    std::ostringstream ss;
    ss << "n " << parent_group_id << '.' << id;
    ss << print_format_attributes(attributes);
    ss << std::endl;
    if (write_mapping && (core != nullptr))
    {
        ss << "& " << parent_group_id << '.' << id;
        ss << '@' << core->parent_tile_id << '.' << core->id;
        ss << std::endl;
    }
    return ss.str();
}
*/

sanafe::NeuronGroup &sanafe::Network::create_neuron_group(
        const std::string &name, const NeuronTemplate &default_config)
{
    groups.emplace(name, NeuronGroup(name, *this, default_config));
    INFO("Created neuron group gid:%s\n", name.c_str());

    return groups.at(name);
}

/*
std::string sanafe::NeuronGroup::description() const
{
    std::ostringstream ss;
    ss << "g " << neurons.size();
    ss << print_format_attributes(default_attributes) << std::endl;
    return ss.str();
}
*/

std::string sanafe::NeuronGroup::info() const
{
    std::ostringstream ss;
    ss << "sanafe::NeuronGroup(gid=" << name;
    ss << " neurons=" << neurons.size();
    //ss << " attributes={" << print_format_attributes(default_attributes);
    ss << "})";
    return ss.str();
}

sanafe::Neuron &sanafe::NeuronGroup::create_neuron(
        const std::string &id, const NeuronTemplate &config)
{
    neurons.emplace(id, Neuron(id, *parent_net, name, config));
    //INFO("Created neuron nid:%s.%s\n", name.c_str(), id.c_str());

    return neurons.at(id);
}

void sanafe::Neuron::connect_to_neuron(Neuron &dest,
        const std::map<std::string, ModelParam> &synapse_params,
        const std::map<std::string, ModelParam> &dendrite_params,
        const std::optional<std::string> &synapse_hw_name)
{
    connections_out.emplace_back(Connection(connections_out.size()));
    Connection &con = connections_out.back();
    con.pre_neuron = this;
    con.post_neuron = &dest;
    if (synapse_hw_name.has_value())
    {
        con.synapse_hw_name = synapse_hw_name.value();
    }
    else
    {
        con.synapse_hw_name = default_synapse_hw_name;
    }
    con.synapse_params.insert(synapse_params.begin(), synapse_params.end());
    con.dendrite_params.insert(dendrite_params.begin(), dendrite_params.end());

    //INFO("\tAdded con %s.%s->%s.%s (w:%lf)\n", con.pre_neuron->parent_group_id.c_str(),
    //        con.pre_neuron->id.c_str(), con.post_neuron->parent_group_id.c_str(),
    //        con.post_neuron->id.c_str(), static_cast<double>(con.synapse_params["w"]));
}

sanafe::Network sanafe::load_net(
        const std::filesystem::path &path, Architecture &arch)
{
    std::ifstream network_fp;

    network_fp.open(path);
    if (network_fp.fail())
    {
        const std::string error = "Error: Network file: failed to open (" +
                std::string(path) + ").";
        throw std::invalid_argument(error);
    }
    INFO("Loading network from file: %s\n", path.c_str());
    Network net = description_parse_network_file_yaml(network_fp, arch);
    network_fp.close();

    return net;
}

std::string sanafe::Network::info() const
{
    return "sanafe::Network(groups=" + std::to_string(groups.size()) + ")";
}

void sanafe::Network::check_mapped() const
{
    // Check that all network neurons are mapped to a physical core
    //  If a neuron is not, print an error message and stop the simulation
    for (const auto &group_entry : groups)
    {
        const NeuronGroup &group = group_entry.second;
        for (const auto &neuron_entry : group.neurons)
        {
            const Neuron &neuron = neuron_entry.second;
            if (neuron.core == nullptr)
            {
                INFO("Error: Neuron %s.%s not mapped to H/W.\n",
                        group.name.c_str(), neuron.id.c_str());
                throw std::runtime_error("Error: Neuron isn't mapped");
            }
        }
    }
}

/*
void sanafe::Network::save_net_description(
        const std::filesystem::path &path, const bool save_mapping) const
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        INFO("Error: Couldn't open net file to save to: %s\n", path.c_str());
        throw std::invalid_argument(
                "Error: Couldn't open net file to save to.");
    }

    // Save all groups first
    for (const auto &group : groups)
    {
        out << group->description();
    }

    // Now save all neurons and connections
    for (const auto &group : groups)
    {
        for (const Neuron &n : group->neurons)
        {
            // Save neuron description
            out << n.description(save_mapping);
            // Save all edges for this neuron
            for (const Connection &con : n.connections_out)
            {
                out << con.description();
            }
        }
    }

    return;
}
*/

// Reimplement setting functions for multiple neurons and edges
/*
void sanafe::NeuronGroup::set_attribute_multiple(
        const std::string &attr, const std::vector<std::string> &values)
{
    if (values.size() != neurons.size())
    {
        INFO("Error: Attribute values must be defined for all neurons "
             "(%lu!=%lu).\n",
                attr.size(), neurons.size());
    }
    size_t nid = 0;
    std::map<std::string, std::string> map;
    for (const std::string &v : values)
    {
        auto &n = neurons[nid];
        map[attr] = v;
        n.set_attributes(map);
        nid++;
    }
}
*/

/*
void sanafe::NeuronGroup::connect_neurons(NeuronGroup &dest_group,
        const std::vector<std::pair<int, int>> &src_dest_id_pairs,
        const std::map<std::string, std::vector<std::string>> &attr_lists)
{
    int edge_id = 0;
    for (auto [src_id, dest_id] : src_dest_id_pairs)
    {
        TRACE2("Connecting neurons, neurons.size=%lu\n", neurons.size());
        if ((src_id < 0) || (static_cast<size_t>(src_id) >= neurons.size()))
        {
            INFO("src_id:%d out of range (0 <= id <= %lu).\n", src_id,
                    neurons.size());
            throw std::invalid_argument("Error: src id is out of range.");
        }
        if ((dest_id < 0) ||
                (static_cast<size_t>(dest_id) >= dest_group.neurons.size()))
        {
            INFO("dest_id:%d out of range (0 <= id <= %lu).\n", dest_id,
                    dest_group.neurons.size());
            throw std::invalid_argument("Error: dest nid is out of range.");
        }

        // Create attributes map for this neuron
        std::map<std::string, std::string> attr;
        for (auto &[key, value_list] : attr_lists)
        {
            if (value_list.size() != src_dest_id_pairs.size())
            {
                INFO("Error: Length of attribute list != "
                     "number of defined edges. (%lu!=%lu).\n",
                        value_list.size(), src_dest_id_pairs.size());
                throw std::invalid_argument(
                        "Error: Length of attribute list != "
                        "number of defined edges.");
            }
            attr[key] = value_list[src_id];
        }

        Neuron &src = neurons[src_id];
        Neuron &dest = dest_group.neurons[dest_id];
        src.connect_to_neuron(dest, attr);
        edge_id++;
    }
}
*/
