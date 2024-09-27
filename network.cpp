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
//#include "models.hpp"
#include "network.hpp"
#include "print.hpp"

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

sanafe::NeuronGroup::NeuronGroup(const std::string group_name,
        Network &parent_net, const size_t neuron_count,
        const NeuronTemplate &default_config)
        : default_neuron_config(default_config)
        , parent_net(&parent_net)
        , name(std::move(group_name))
{
    neurons.reserve(neuron_count);
    for (size_t nid = 0; nid < neuron_count; nid++)
    {
        neurons.emplace_back(
                Neuron(nid, parent_net, group_name, default_config));
    }
}

sanafe::Neuron::Neuron(const size_t neuron_id, Network &parent_net,
        const std::string parent_group_id, const NeuronTemplate &config)
        : parent_net(&parent_net)
        , soma_hw_name(config.soma_hw_name)
        , default_synapse_hw_name(config.default_synapse_hw_name)
        , dendrite_hw_name(config.dendrite_hw_name)
        , parent_group_id(std::move(parent_group_id))
        , id(neuron_id)
        , force_synapse_update(config.force_synapse_update)
        , force_dendrite_update(config.force_dendrite_update)
        , force_soma_update(config.force_soma_update)
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

void sanafe::Neuron::set_attributes(const NeuronTemplate &attributes)
{
    // TODO: make NeuronTemplate into NeuronParameters and make each field optional
    log_spikes = attributes.log_spikes;
    log_potential = attributes.log_potential;
    force_dendrite_update = attributes.force_dendrite_update;
    force_soma_update = attributes.force_soma_update;
    force_synapse_update = attributes.force_synapse_update;

    // TODO
    for (auto &[key, param] : attributes.soma_model_params)
    {
        if (param.forward_to_dendrite)
        {
            if (dendrite_hw != nullptr)
            {
                dendrite_hw->set_attribute(mapped_address, key, param);
            }
            soma_model_params.insert(attributes.soma_model_params.begin(),
                    attributes.soma_model_params.end());
        }
        if (param.forward_to_soma)
        {
            if (soma_hw != nullptr)
            {
                soma_hw->set_attribute(mapped_address, key, param);
            }
            // TODO: enable and disable recording attributes
            //  This can be automatically disabled when running in kernel mode
            //   and loading from description file.
            //  This can automatically be enabled when operating via the Python
            //   interface, since in that mode it might be useful to make the
            //   network saveable
            soma_model_params.insert(attributes.soma_model_params.begin(),
                    attributes.soma_model_params.end());
        }
    }

    // TODO: remove the separate soma and dendrite parameters, now its supported
    //  by just one
    for (auto &[key, param] : attributes.dendrite_model_params)
    {
        if (param.forward_to_soma)
        {
            if (param.forward_to_dendrite)
            {
                if (dendrite_hw != nullptr)
                {
                    dendrite_hw->set_attribute(mapped_address, key, param);
                }
                soma_model_params.insert(attributes.soma_model_params.begin(),
                        attributes.soma_model_params.end());
            }
            if (soma_hw != nullptr)
            {
                soma_hw->set_attribute(mapped_address, key, param);
            }
            // TODO: enable and disable recording attributes
            //  This can be automatically disabled when running in kernel mode
            //   and loading from description file.
            //  This can automatically be enabled when operating via the Python
            //   interface, since in that mode it might be useful to make the
            //   network saveable
            soma_model_params.insert(attributes.soma_model_params.begin(),
                    attributes.soma_model_params.end());
        }

    }
}

std::string sanafe::Neuron::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Neuron(nid=" << parent_group_id << '.' << id;
    ss << " connections_out=" << connections_out.size() << ")";
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
        const std::string name, const size_t neuron_count,
        const NeuronTemplate &default_config)
{
    groups.emplace(
            name, NeuronGroup(name, *this, neuron_count, default_config));
    INFO("Created neuron group gid:%s with %zu neurons\n", name.c_str(),
            neuron_count);

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

sanafe::Connection &sanafe::Neuron::connect_to_neuron(Neuron &dest)
{
    connections_out.emplace_back(Connection(connections_out.size()));
    Connection &con = connections_out.back();
    con.pre_neuron = this;
    con.post_neuron = &dest;
    //if (synapse_hw_name.has_value())
    //{
    //    con.synapse_hw_name = synapse_hw_name.value();
    //}
    //else
    {
        con.synapse_hw_name = default_synapse_hw_name;
    }

    TRACE1("\tAdded con %s.%s->%s.%s (w:%lf)\n",
            con.pre_neuron->parent_group_id.c_str(), con.pre_neuron->id.c_str(),
            con.post_neuron->parent_group_id.c_str(),
            con.post_neuron->id.c_str(),
            static_cast<double>(con.synapse_params["w"]));

    return con;
}

sanafe::Network sanafe::load_net(const std::filesystem::path &path,
        Architecture &arch, const bool use_netlist_format)
{
    std::ifstream network_fp;

    network_fp.open(path);
    if (network_fp.fail())
    {
        const std::string error = "Error: Network file: failed to open (" +
                std::string(path) + ").";
        throw std::invalid_argument(error);
    }

    Network net;
    if (use_netlist_format)
    {
        INFO("Loading network from netlist file (legacy): %s\n", path.c_str());
        // Fall back to the original netlist based format used by SANA-FE v1.
        //  This is supported mainly for back-compatibility
        net = description_parse_network_file_netlist(network_fp, arch);
    }
    else
    {
        INFO("Loading network from YAML file: %s\n", path.c_str());
        net = description_parse_network_file_yaml(network_fp, arch);
    }
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
        for (const auto &neuron : group.neurons)
        {
            if (neuron.core == nullptr)
            {
                INFO("Error: Neuron %s.%zu not mapped to H/W.\n",
                        group.name.c_str(), neuron.id);
                throw std::runtime_error("Error: Neuron isn't mapped");
            }
        }
    }
}

void sanafe::NeuronGroup::connect_neurons_sparse(NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelParam>> &attribute_lists,
        const std::vector<std::pair<size_t, size_t>> &source_dest_id_pairs)
{
    for (auto [source_id, dest_id] : source_dest_id_pairs)
    {
        TRACE2("Connecting neurons, neurons.size=%lu\n", neurons.size());
        if (source_id >= neurons.size())
        {
            INFO("source_id:%zu out of range (0 <= id <= %zu).\n", source_id,
                    neurons.size());
            throw std::invalid_argument("Error: src id is out of range.");
        }
        if (dest_id >= dest_group.neurons.size())
        {
            INFO("dest_id:%zu out of range (0 <= id <= %zu).\n", dest_id,
                    dest_group.neurons.size());
            throw std::invalid_argument("Error: dest nid is out of range.");
        }

        Neuron &source = neurons[source_id];
        Neuron &dest = dest_group.neurons[dest_id];
        Connection &con = source.connect_to_neuron(dest);

        // Create attributes map for this neuron
        std::map<std::string, ModelParam> attributes;
        for (auto &[key, value_list] : attribute_lists)
        {
            if (value_list.size() != source_dest_id_pairs.size())
            {
                INFO("Error: Length of attribute list != "
                     "number of defined edges. (%zu!=%zu).\n",
                        value_list.size(), source_dest_id_pairs.size());
                throw std::invalid_argument(
                        "Error: Length of attribute list != "
                        "number of defined edges.");
            }
            attributes[key] = value_list[source_id];
        }

        // TODO: create a set attributes function for the connection that
        //  filters when forwarding to hardware unit
        con.dendrite_params = attributes;
        con.synapse_params = attributes;
    }
}

void sanafe::NeuronGroup::connect_neurons_conv2d(NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelParam>> &attribute_lists,
        const Conv2DParameters &convolution)
{
    // Only support channel last storage
    // TODO:There are two things I need to consider -
    //  1) is the order that we map to the neurons i.e. how do we flatten
    //  the two layers from 3d/4d tensors down to 1d. This affects how to calculate the src and dest index
    //  2) is how we represent the filter/kernel attributes in the description file. Maybe its easier
    //    not to store it flattened for the user. In sana-fe we can assume a
    //    internal representation like here.
    //
    // TODO: This code is based on the SNNToolbox. Confusingly,
    //  by default, the SNNToolbox uses channels_last (tensorflow default)
    //  when saving the filter but uses channels first when actually
    //  mapping layer inputs/outputs to the layer of neurons. We should support
    //  filter_channels_last and input_channels_last and output_channels_last
    //  separately, to be as flexible as possible. For now hard-code.
    const int pad_width = 0;
    const int pad_height = 0;

    int output_width = ((convolution.input_width + (2 * pad_width) -
                                convolution.kernel_width) /
                               convolution.stride_width) +
            1;
    int output_height = ((convolution.input_height + (2 * pad_height) -
                                 convolution.kernel_height) /
                                convolution.stride_height) +
            1;
    int output_channels = convolution.kernel_count;

    size_t expected_input_size = convolution.input_channels *
            convolution.input_width * convolution.input_height;
    size_t expected_output_size =
            output_channels * output_width * output_height;

    if (expected_input_size != neurons.size())
    {
        const std::string error = "Expected " +
                std::to_string(expected_output_size) +
                " neurons in source group for convolution but there are " +
                std::to_string(dest_group.neurons.size()) + " neurons.\n";
        INFO("Error: %s", error.c_str());
        throw std::invalid_argument(error);
    }
    if (expected_output_size != dest_group.neurons.size())
    {
        const std::string error = "Expected " +
                std::to_string(expected_output_size) +
                " neurons in dest group for convolution but there are " +
                std::to_string(dest_group.neurons.size()) + " neurons.\n";
        INFO("Error: %s", error.c_str());
        throw std::invalid_argument(error);
    }

    // Create the convolutional connections
    for (int c_out = 0; c_out < output_channels; ++c_out)
    {
        for (int y_out = 0; y_out < output_height; ++y_out)
        {
            for (int x_out = 0; x_out < output_width; ++x_out)
            {
                int dest_idx = c_out * output_width * output_height;
                dest_idx += y_out * output_width;
                dest_idx += x_out;

                Neuron &dest = dest_group.neurons[dest_idx];
                for (int c_in = 0; c_in < convolution.input_channels; ++c_in)
                {
                    for (int y_filter = 0; y_filter < convolution.kernel_height;
                            ++y_filter)
                    {
                        const int y_position =
                                (y_out * convolution.stride_height) + y_filter;
                        if ((y_position < 0) ||
                                (y_position >= convolution.input_height))
                        {
                            continue;
                        }
                        for (int x_filter = 0;
                                x_filter < convolution.kernel_width; ++x_filter)
                        {
                            const int x_position =
                                    (x_out * convolution.stride_width) +
                                    x_filter;
                            if (x_position < 0 ||
                                    x_position >= convolution.input_width)
                            {
                                continue;
                            }

                            // Calculate the index of the src neuron
                            int source_idx = c_in * convolution.input_width *
                                    convolution.input_height;
                            source_idx += ((y_out * convolution.stride_height) +
                                                  y_filter) *
                                    convolution.input_width;
                            source_idx += ((x_out * convolution.stride_width) +
                                    x_filter);

                            Neuron &source = neurons[source_idx];

                            // weight = kernels[y_kernel, x_kernel, c_in, c_out]
                            int filter_idx = y_filter *
                                    convolution.kernel_width *
                                    convolution.input_channels *
                                    output_channels;
                            filter_idx += x_filter *
                                    convolution.input_channels *
                                    output_channels;
                            filter_idx += c_in * output_channels;
                            filter_idx += c_out;

                            // Create the connection
                            Connection &con = source.connect_to_neuron(dest);

                            // Set the attributes for this connection, using
                            //  the list of attributes
                            for (auto &[key, attribute_list] : attribute_lists)
                            {
                                if (attribute_list.size() <=
                                        static_cast<size_t>(filter_idx))
                                {
                                    INFO("Error: Not enough entries defined "
                                         "for attribute: %s\n",
                                            key.c_str());
                                    throw std::invalid_argument(
                                            "Not enough entries defined "
                                            "for attribute");
                                }
                                const ModelParam &attribute =
                                        attribute_list[filter_idx];
                                if (attribute.forward_to_synapse)
                                {
                                    con.synapse_params[key] = attribute;
                                }
                                if (attribute.forward_to_dendrite)
                                {
                                    con.dendrite_params[key] = attribute;
                                }
                            }
                            TRACE1("%s.%zu->%s.%zu w=%d\n",
                                    source.parent_group_id.c_str(), source.id,
                                    dest.parent_group_id.c_str(), dest.id,
                                    (int) con.synapse_params["weight"]);
                        }
                    }
                }
            }
        }
    }

    // TODO: hack, remove
    /*
    for (auto &n : neurons)
    {
        for (auto &con : n.connections_out)
        {
            printf("e %s.%zu->%s.%zu w=%0.1lf\n",
                    con.pre_neuron->parent_group_id.c_str(), con.pre_neuron->id,
                    con.post_neuron->parent_group_id.c_str(),
                    con.post_neuron->id, (double) con.synapse_params["weight"]);
        }
    }
    */
}

void sanafe::NeuronGroup::connect_neurons_dense(NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelParam>> &attribute_lists)
{
    for (size_t source_index = 0; source_index < neurons.size(); ++source_index)
    {
        Neuron &source = neurons[source_index];

        for (size_t dest_index = 0; dest_index < dest_group.neurons.size();
                ++dest_index)
        {
            Neuron &dest = dest_group.neurons[dest_index];
            const size_t list_index =
                    (source_index * dest_group.neurons.size()) + dest_index;
            Connection &con = source.connect_to_neuron(dest);
            for (auto &[key, attribute_list] : attribute_lists)
            {
                if (attribute_list.size() <= list_index)
                {
                    INFO("Error: Not enough entries defined "
                         "for attribute: %s\n",
                            key.c_str());
                    throw std::invalid_argument("Not enough entries defined "
                                                "for attribute");
                }
                const ModelParam &attribute = attribute_list[list_index];
                if (attribute.forward_to_synapse)
                {
                    con.synapse_params[key] = attribute;
                }
                if (attribute.forward_to_dendrite)
                {
                    con.dendrite_params[key] = attribute;
                }
            }
        }
    }
}

/*
void sanafe::Network::save_netlist(
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