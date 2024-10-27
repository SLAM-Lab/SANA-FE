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
#include "network.hpp"
#include "print.hpp"

sanafe::NeuronGroup::NeuronGroup(const std::string group_name,
        SpikingNetwork &net,
        const size_t neuron_count,
        const NeuronConfiguration &default_config)
        : default_neuron_config(default_config)
        , name(std::move(group_name))
{
    neurons.reserve(neuron_count);
    for (size_t nid = 0; nid < neuron_count; nid++)
    {
        neurons.emplace_back(
                Neuron(nid, net, group_name, default_config));
    }
}

sanafe::Neuron::Neuron(const size_t neuron_id,
        SpikingNetwork &net,
        const std::string parent_group_id, const NeuronConfiguration &config)
        : parent_group_id(std::move(parent_group_id))
        , parent_net(net)
        , id(neuron_id)
{
    configure(config);
}

void sanafe::Neuron::map_to_core(const CoreConfiguration &core)
{
    this->core_id = core.address.id;
    mapping_order = parent_net.update_mapping_count();
    TRACE1(NET, "Mapping order for nid:%s.%zu = %zu\n", parent_group_id.c_str(), id,
            mapping_order);

    return;
}

void sanafe::Neuron::configure(const NeuronConfiguration &config)
{
    if (config.default_synapse_hw_name.has_value())
    {
        default_synapse_hw_name = config.default_synapse_hw_name.value();
    }
    if (config.dendrite_hw_name.has_value())
    {
        dendrite_hw_name = config.dendrite_hw_name.value();
    }
    if (config.soma_hw_name.has_value())
    {
        soma_hw_name = config.soma_hw_name.value();
    }
    if (config.log_spikes.has_value())
    {
        log_spikes = config.log_spikes.value();
    }
    if (config.log_potential.has_value())
    {
        log_potential = config.log_potential.value();
    }
    if (config.force_dendrite_update)
    {
        force_dendrite_update = config.force_dendrite_update.value();
    }
    if (config.force_soma_update.has_value())
    {
        force_soma_update = config.force_soma_update.value();
    }
    if (config.force_synapse_update)
    {
        force_synapse_update = config.force_synapse_update.value();
    }

    model_parameters.insert(config.model_parameters.begin(),
            config.model_parameters.end());
}

std::string sanafe::Neuron::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Neuron(nid=" << parent_group_id << '.' << id;
    ss << " edges_out=" << edges_out.size() << ")";
    //ss << " attributes={" << print_format_attributes(attributes) << "})";
    return ss.str();
}

sanafe::NeuronGroup &sanafe::SpikingNetwork::create_neuron_group(
        const std::string name, const size_t neuron_count,
        const NeuronConfiguration &default_config)
{
    groups.emplace(
            name, NeuronGroup(name, *this, neuron_count, default_config));
    TRACE1(NET, "Created neuron group gid:%s with %zu neurons\n", name.c_str(),
            neuron_count);

    return groups.at(name);
}

size_t sanafe::SpikingNetwork::update_mapping_count()
{
    ++mapping_count;
    return mapping_count;
}

std::string sanafe::NeuronGroup::info() const
{
    std::ostringstream ss;
    ss << "sanafe::NeuronGroup(gid=" << name;
    ss << " neurons=" << neurons.size();
    //ss << " attributes={" << print_format_attributes(default_attributes);
    ss << "})";
    return ss.str();
}

size_t sanafe::Neuron::connect_to_neuron(Neuron &dest)
{
    edges_out.emplace_back(edges_out.size());
    Connection &edge = edges_out.back();
    edge.pre_neuron.neuron_id = this->id;
    edge.pre_neuron.group_name = this->parent_group_id;
    edge.post_neuron.neuron_id = dest.id;
    edge.post_neuron.group_name = dest.parent_group_id;
    //if (synapse_hw_name.has_value())
    //{
    //    con.synapse_hw_name = synapse_hw_name.value();
    //}
    //else
    {
        edge.synapse_hw_name = dest.default_synapse_hw_name;
    }

    TRACE1(NET, "\tAdded con %s.%zu->%s.%zu\n",
            edge.pre_neuron.group_name.c_str(),
            edge.pre_neuron.neuron_id.value(),
            edge.post_neuron.group_name.c_str(),
            edge.post_neuron.neuron_id.value());

    return edge.id;
}

sanafe::SpikingNetwork sanafe::load_net(const std::filesystem::path &path,
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

    SpikingNetwork net;
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

std::string sanafe::SpikingNetwork::info() const
{
    return "sanafe::Network(groups=" + std::to_string(groups.size()) + ")";
}

// TODO: change this function to map network to hardware
//  First, map every group and neuron onto the hardware. Any unmapped neurons can raise an error
//  Second, map all edges as connections on the hardware
/*
void sanafe::SpikingNetwork::check_mapped() const
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
                INFO("neuron in sanafe: %p\n", (void *) &neuron);
                std::string error = "Error: Neuron " + group.name + "." +
                        std::to_string(neuron.id) + " not mapped to H/W.\n";
                throw std::runtime_error(error);
            }
        }
    }
}
*/

void sanafe::NeuronGroup::connect_neurons_sparse(NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelParam>> &attribute_lists,
        const std::vector<std::pair<size_t, size_t>> &source_dest_id_pairs)
{
    for (auto [source_id, dest_id] : source_dest_id_pairs)
    {
        TRACE2(NET, "Connecting neurons, neurons.size=%lu\n", neurons.size());
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
        size_t connection_idx = source.connect_to_neuron(dest);
        Connection &con = source.edges_out[connection_idx];

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

        con.synapse_params = attributes;
        con.dendrite_params = attributes;
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
                            const size_t con_idx =
                                    source.connect_to_neuron(dest);
                            Connection &con = source.edges_out[con_idx];

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
                                if (attribute.forward_to_dendrite)
                                {
                                    con.dendrite_params[key] = attribute;
                                }
                                if (attribute.forward_to_synapse)
                                {
                                    con.synapse_params[key] = attribute;
                                }
                            }
                            TRACE1(NET, "%s.%zu->%s.%zu w=%d\n",
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
            const size_t con_idx = source.connect_to_neuron(dest);
            Connection &con = source.edges_out[con_idx];
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
