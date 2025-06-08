// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// network.cpp
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem> // For std::filesystem::path
#include <fstream>
#include <functional> // For std::reference_wrapper
#include <map>
#include <optional>
#include <sstream>

#include "arch.hpp"
#include "attribute.hpp"
#include "netlist.hpp"
#include "network.hpp"
#include "print.hpp"
#include "yaml_snn.hpp"
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

std::string sanafe::NeuronAddress::info() const
{
    std::ostringstream ss;
    // Output the address as a string rather than an object, making it easier
    //  to integrate to other formats later e.g. YAML
    assert(!group_name.empty());
    ss << group_name;
    if (neuron_offset.has_value())
    {
        ss << '.' << neuron_offset.value();
    }

    return ss.str();
}

std::string sanafe::Connection::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Connection(pre_neuron=" << pre_neuron.group_name;
    if (pre_neuron.neuron_offset.has_value())
    {
        ss << "." << pre_neuron.neuron_offset.value();
    }

    ss << " post_neuron=" << post_neuron.group_name;
    if (post_neuron.neuron_offset.has_value())
    {
        ss << "." << post_neuron.neuron_offset.value();
    }
    ss << ")";

    return ss.str();
}

sanafe::NeuronGroup::NeuronGroup(const std::string group_name,
        SpikingNetwork &net, const size_t neuron_count,
        const NeuronConfiguration &default_config)
        : default_neuron_config(default_config)
        , name(group_name)
{
    neurons.reserve(neuron_count);
    for (size_t neuron_offset = 0; neuron_offset < neuron_count;
            ++neuron_offset)
    {
        neurons.emplace_back(neuron_offset, net, group_name, default_config);
    }
}

sanafe::Neuron::Neuron(const size_t neuron_offset, SpikingNetwork &net,
        std::string parent_group_name, const NeuronConfiguration &config)
        : parent_group_name(std::move(parent_group_name))
        , parent_net(net)
        , offset(neuron_offset)
{
    set_attributes(config);
}

void sanafe::Neuron::map_to_core(const CoreConfiguration &core)
{
    this->core_address = core.address;
    SpikingNetwork &net = parent_net;
    mapping_order = net.update_mapping_count();
    TRACE1(NET, "Mapping order for nid:%s.%zu = %zu\n",
            parent_group_name.c_str(), offset, mapping_order);
}

void sanafe::Neuron::set_attributes(const NeuronConfiguration &config)
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

    model_attributes.insert(
            config.model_attributes.begin(), config.model_attributes.end());
}

std::string sanafe::Neuron::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Neuron(nid=" << parent_group_name << '.' << offset;
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
    edge.pre_neuron.neuron_offset = this->offset;
    edge.pre_neuron.group_name = this->parent_group_name;
    edge.post_neuron.neuron_offset = dest.offset;
    edge.post_neuron.group_name = dest.parent_group_name;
    edge.synapse_hw_name = dest.default_synapse_hw_name;

    TRACE1(NET, "\tAdded con %s.%zu->%s.%zu\n",
            edge.pre_neuron.group_name.c_str(),
            edge.pre_neuron.neuron_offset.value(),
            edge.post_neuron.group_name.c_str(),
            edge.post_neuron.neuron_offset.value());

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
        net = netlist_parse_file(network_fp, arch);
    }
    else
    {
        INFO("Loading network from YAML file: %s\n", path.c_str());
        net = yaml_parse_network_file(network_fp, arch);
    }
    network_fp.close();

    return net;
}

std::string sanafe::SpikingNetwork::info() const
{
    return "sanafe::Network(groups=" + std::to_string(groups.size()) + ")";
}

void sanafe::NeuronGroup::connect_neurons_sparse(NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelAttribute>>
                &attribute_lists,
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
        const size_t connection_idx = source.connect_to_neuron(dest);
        Connection &con = source.edges_out[connection_idx];

        // Create attributes map for this neuron
        std::map<std::string, ModelAttribute> attributes;
        for (const auto &[key, value_list] : attribute_lists)
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

        con.synapse_attributes = attributes;
        con.dendrite_attributes = attributes;
    }
}

// 2D Convolution Neural Network Connection Algorithm
//

void sanafe::NeuronGroup::connect_neurons_conv2d(NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelAttribute>>
                &attribute_lists,
        const Conv2DParameters &convolution)
{
    // This algorithm creates synaptic connections between two neuron groups to
    //  implement a 2D convolutional layer, similar to CNNs in deep learning
    //  frameworks.
    //
    // Algorithm:
    // ==============
    // 1. Setup phase:
    //  - Calculate output dimensions based on input size, kernel size, and
    //     stride
    //  - Validate that neuron group sizes match expected dimensions
    //
    // 2. Connection phase (nested loops):
    //  For each output position (channel, y, x) i.e., neuron:
    //    For each input channel:
    //      For each kernel position (ky, kx):
    //        - Calculate corresponding input neuron position
    //        - Skip if position is outside input bounds (boundary handling)
    //        - Create synaptic connection from input neuron to output neuron
    //        - Set connection weights/attributes from the kernel filter values
    // Calculate output dimensions and validate inputs
    // 1. Setup phase
    const auto dims = conv2d_calculate_dimensions(convolution);
    conv2d_validate_neuron_counts(dest_group, dims);

    // 2. Connection phase - Create the convolutional connections for each
    //  output neuron
    for (int c_out = 0; c_out < dims.output_channels; ++c_out)
    {
        for (int y_out = 0; y_out < dims.output_height; ++y_out)
        {
            for (int x_out = 0; x_out < dims.output_width; ++x_out)
            {
                const Conv2DCoordinate out{c_out, y_out, x_out};
                conv2d_create_output_neuron_connections(
                        dest_group, attribute_lists, convolution, dims, out);
            }
        }
    }
}

void sanafe::NeuronGroup::conv2d_create_output_neuron_connections(
        NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelAttribute>>
                &attribute_lists,
        const Conv2DParameters &convolution, const Conv2DOutputDimensions &dims,
        const Conv2DCoordinate &output_coordinate)
{
    int dest_idx =
            output_coordinate.channel * dims.output_width * dims.output_height;
    dest_idx += output_coordinate.y * dims.output_width;
    dest_idx += output_coordinate.x;
    Neuron &dest = dest_group.neurons[dest_idx];

    // Connect to source/input neurons through kernel
    for (int c_in = 0; c_in < convolution.input_channels; ++c_in)
    {
        conv2d_create_kernel_connections(dest, attribute_lists, convolution,
                dims, output_coordinate, c_in);
    }
}

void sanafe::NeuronGroup::conv2d_create_kernel_connections(Neuron &dest,
        const std::map<std::string, std::vector<ModelAttribute>>
                &attribute_lists,
        const Conv2DParameters &convolution, const Conv2DOutputDimensions &dims,
        const Conv2DCoordinate &out, int c_in)
{
    for (int y_filter = 0; y_filter < convolution.kernel_height; ++y_filter)
    {
        const int y_position = (out.y * convolution.stride_height) + y_filter;
        if (!conv2d_is_position_valid(y_position, convolution.input_height))
        {
            continue;
        }

        for (int x_filter = 0; x_filter < convolution.kernel_width; ++x_filter)
        {
            const int x_position =
                    (out.x * convolution.stride_width) + x_filter;
            if (!conv2d_is_position_valid(x_position, convolution.input_width))
            {
                continue;
            }

            const Conv2DPosition position = {out, c_in, y_filter, x_filter};
            const auto indices =
                    conv2d_calculate_indices(convolution, dims, position);

            Neuron &source = neurons[indices.source_idx];
            conv2d_create_and_configure_connection(
                    source, dest, attribute_lists, indices.filter_idx);
        }
    }
}

sanafe::Conv2DOutputDimensions sanafe::NeuronGroup::conv2d_calculate_dimensions(
        const Conv2DParameters &convolution)
{
    // Only support channels-last storage for now
    // Inputs must be flattened (C-style) to 1D arrays
    //
    // TODO: This code is based on the SNNToolbox. Confusingly,
    //  by default, the SNNToolbox uses channels_last (tensorflow default)
    //  when saving the filter but uses channels first when actually
    //  mapping layer inputs/outputs to the layer of neurons. We should support
    //  filter_channels_last and input_channels_last and output_channels_last
    //  separately, to be as flexible as possible. For now hard-code.
    const int pad_width = 0;
    const int pad_height = 0;

    Conv2DOutputDimensions dims{};
    // Standard convolution output size formula: (W - K + 2P) / S + 1"
    dims.output_width = ((convolution.input_width + (2 * pad_width) -
                                 convolution.kernel_width) /
                                convolution.stride_width) +
            1;
    dims.output_height = ((convolution.input_height + (2 * pad_height) -
                                  convolution.kernel_height) /
                                 convolution.stride_height) +
            1;
    dims.output_channels = convolution.kernel_count;

    dims.expected_input_size = static_cast<size_t>(convolution.input_channels) *
            static_cast<size_t>(convolution.input_width) *
            static_cast<size_t>(convolution.input_height);
    dims.expected_output_size = static_cast<size_t>(convolution.kernel_count) *
            static_cast<size_t>(dims.output_width) *
            static_cast<size_t>(dims.output_height);

    return dims;
}

void sanafe::NeuronGroup::conv2d_validate_neuron_counts(
        const NeuronGroup &dest_group, const Conv2DOutputDimensions &dims) const
{
    // Check the total neuron counts for the source and destination neuron
    //  groups against the expected number of neurons, given the convolution
    //  dimensions.
    if (dims.expected_input_size != neurons.size())
    {
        const std::string error = "Expected " +
                std::to_string(dims.expected_input_size) +
                " neurons in source group for convolution but there are " +
                std::to_string(neurons.size()) + " neurons.\n";
        INFO("Error: %s", error.c_str());
        throw std::invalid_argument(error);
    }

    if (dims.expected_output_size != dest_group.neurons.size())
    {
        const std::string error = "Expected " +
                std::to_string(dims.expected_output_size) +
                " neurons in dest group for convolution but there are " +
                std::to_string(dest_group.neurons.size()) + " neurons.\n";
        INFO("Error: %s", error.c_str());
        throw std::invalid_argument(error);
    }
}

sanafe::Conv2DIndices sanafe::NeuronGroup::conv2d_calculate_indices(
        const Conv2DParameters &convolution, const Conv2DOutputDimensions &dims,
        const Conv2DPosition &pos)
{
    // Calculate indices for the neurons within the neuron groups, and for the
    //  correct index into the flattened filter array
    Conv2DIndices indices{};
    indices.dest_idx = pos.output_coordinate.channel * dims.output_width *
            dims.output_height;
    indices.dest_idx += pos.output_coordinate.y * dims.output_width;
    indices.dest_idx += pos.output_coordinate.x;

    // Calculate source index
    indices.source_idx =
            pos.c_in * convolution.input_width * convolution.input_height;
    indices.source_idx +=
            ((pos.output_coordinate.y * convolution.stride_height) +
                    pos.y_filter) *
            convolution.input_width;
    indices.source_idx +=
            ((pos.output_coordinate.x * convolution.stride_width) +
                    pos.x_filter);

    // Calculate filter index
    indices.filter_idx = pos.y_filter * convolution.kernel_width *
            convolution.input_channels * convolution.kernel_count;
    indices.filter_idx += pos.x_filter * convolution.input_channels *
            convolution.kernel_count;
    indices.filter_idx += pos.c_in * convolution.kernel_count;
    indices.filter_idx += pos.output_coordinate.channel;
    // Filter laid out as [y][x][input_channel][output_channel]

    return indices;
}

bool sanafe::NeuronGroup::conv2d_is_position_valid(
        int position, int max_size) noexcept
{
    return (position >= 0) && (position < max_size);
}

void sanafe::NeuronGroup::conv2d_create_and_configure_connection(Neuron &source,
        Neuron &dest,
        const std::map<std::string, std::vector<ModelAttribute>>
                &attribute_lists,
        int filter_idx)
{
    // Create the connection
    const size_t con_idx = source.connect_to_neuron(dest);
    Connection &con = source.edges_out[con_idx];

    // Set the attributes for this connection
    for (const auto &[key, attribute_list] : attribute_lists)
    {
        if (attribute_list.size() <= static_cast<size_t>(filter_idx))
        {
            INFO("Error: Not enough entries defined for attribute (%zu): %s\n",
                    attribute_list.size(), key.c_str());
            throw std::invalid_argument(
                    "Not enough entries defined for attribute");
        }

        const ModelAttribute &attribute = attribute_list[filter_idx];
        if (attribute.forward_to_dendrite)
        {
            con.dendrite_attributes[key] = attribute;
        }
        if (attribute.forward_to_synapse)
        {
            con.synapse_attributes[key] = attribute;
        }
    }

    TRACE1(NET, "%s.%zu->%s.%zu w=%d\n", source.parent_group_name.c_str(),
            source.offset, dest.parent_group_name.c_str(), dest.offset,
            (int) con.synapse_attributes["weight"]);
}

void sanafe::NeuronGroup::connect_neurons_dense(NeuronGroup &dest_group,
        const std::map<std::string, std::vector<ModelAttribute>>
                &attribute_lists)
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
            for (const auto &[key, attribute_list] : attribute_lists)
            {
                if (attribute_list.size() <= list_index)
                {
                    INFO("Error: Not enough entries defined "
                         "for attribute: %s\n",
                            key.c_str());
                    throw std::invalid_argument("Not enough entries defined "
                                                "for attribute");
                }
                const ModelAttribute &attribute = attribute_list[list_index];
                if (attribute.forward_to_synapse)
                {
                    con.synapse_attributes[key] = attribute;
                }
                if (attribute.forward_to_dendrite)
                {
                    con.dendrite_attributes[key] = attribute;
                }
            }
        }
    }
}

void sanafe::SpikingNetwork::save_netlist(
        const std::filesystem::path &path) const
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        INFO("Error: Couldn't open net file to save to: %s\n", path.c_str());
        throw std::invalid_argument(
                "Error: Couldn't open net file to save to.");
    }

    const auto group_name_to_id = create_group_name_to_id_mapping();

    save_groups_to_netlist(out);
    save_neurons_to_netlist(out, group_name_to_id);
    save_mappings_to_netlist(out, group_name_to_id);
}

std::map<std::string, size_t>
sanafe::SpikingNetwork::create_group_name_to_id_mapping() const
{
    std::map<std::string, size_t> group_name_to_id;
    size_t id = 0UL;

    for (const auto &[group_name, group] : groups)
    {
        group_name_to_id[group_name] = id;
        ++id;
    }

    return group_name_to_id;
}

void sanafe::SpikingNetwork::save_groups_to_netlist(std::ofstream &out) const
{
    TRACE1(NET, "Saving groups\n");

    for (const auto &[group_name, group] : groups)
    {
        TRACE1(NET, "Saving group:%s\n", group_name.c_str());
        out << netlist_group_to_netlist(group) << "\n";
    }
}

void sanafe::SpikingNetwork::save_neurons_to_netlist(std::ofstream &out,
        const std::map<std::string, size_t> &group_name_to_id) const
{
    TRACE1(NET, "Saving neurons\n");

    for (const auto &[group_name, group] : groups)
    {
        for (const Neuron &neuron : group.neurons)
        {
            // Save neuron description
            out << netlist_neuron_to_netlist(neuron, *this, group_name_to_id)
                << "\n";

            // Save all connections for this neuron
            for (const Connection &connection : neuron.edges_out)
            {
                out << netlist_connection_to_netlist(
                               connection, group_name_to_id)
                    << "\n";
            }
        }
    }
}

void sanafe::SpikingNetwork::save_mappings_to_netlist(std::ofstream &out,
        const std::map<std::string, size_t> &group_name_to_id) const
{
    TRACE1(NET, "Saving mappings\n");

    // Collect all neurons from all groups
    std::vector<std::reference_wrapper<const Neuron>> all_neurons;
    for (const auto &[group_name, group] : groups)
    {
        const auto &neurons = group.neurons;
        all_neurons.insert(all_neurons.end(), neurons.begin(), neurons.end());
    }

    // Sort by mapping_order
    std::sort(all_neurons.begin(), all_neurons.end(),
            [](const Neuron &a, const Neuron &b) {
                return a.mapping_order < b.mapping_order;
            });

    // Save sorted mappings
    for (const Neuron &neuron : all_neurons)
    {
        out << netlist_mapping_to_netlist(neuron, group_name_to_id) << "\n";
    }
}

void sanafe::SpikingNetwork::save_yaml(const std::filesystem::path &path) const
{
    yaml_write_network(path, *this);
    yaml_write_mappings_file(path, *this);
}

void sanafe::SpikingNetwork::save(
        const std::filesystem::path &path, const bool use_netlist_format) const
{
    if (use_netlist_format)
    {
        save_netlist(path);
    }
    else
    {
        save_yaml(path);
    }
}
