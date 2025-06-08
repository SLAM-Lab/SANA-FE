// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef YAML_SNN_HEADER_INCLUDED_
#define YAML_SNN_HEADER_INCLUDED_

#include "c4/yml/fwd.hpp"
#include "c4/yml/node.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <list>
#include <map>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "fwd.hpp"
#include "network.hpp"

namespace sanafe
{
// YAML network description
SpikingNetwork yaml_parse_network_file(std::ifstream &fp, Architecture &arch);
SpikingNetwork yaml_parse_network_section(const ryml::Parser &parser, ryml::ConstNodeRef net_node);
void yaml_parse_neuron_group_section(const ryml::Parser &parser, ryml::ConstNodeRef groups_node, SpikingNetwork &net);
size_t description_count_neurons(const ryml::Parser &parser, ryml::ConstNodeRef neuron_node);
void yaml_parse_edges_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef edges_node, SpikingNetwork &net);
void yaml_parse_group(const ryml::Parser &parser, ryml::ConstNodeRef neuron_group_node, SpikingNetwork &net);
void yaml_parse_neuron_section(const ryml::Parser &parser, ryml::ConstNodeRef neuron_node, NeuronGroup &neuron_group);
// Neuron to Hardware Mapping
void description_parse_mapping_section_yaml(const ryml::Parser &parser,
        ryml::ConstNodeRef mappings_node, Architecture &arch,
        SpikingNetwork &net);
void description_parse_mapping(const ryml::Parser &parser,
        ryml::ConstNodeRef mapping_info, Architecture &arch,
        SpikingNetwork &net);
void description_parse_mapping_info(const ryml::Parser &parser,
        ryml::ConstNodeRef info, Neuron &n, std::string &core_name);
void description_parse_neuron(const std::string &id, const ryml::Parser &parser,
        ryml::ConstNodeRef attributes, NeuronGroup &neuron_group);
NeuronConfiguration yaml_parse_neuron_attributes(const ryml::Parser &parser,
        ryml::ConstNodeRef attributes,
        const NeuronConfiguration &default_template = NeuronConfiguration());
void description_parse_edge(const std::string &description,
        const ryml::Parser &parser, ryml::ConstNodeRef attributes_node,
        SpikingNetwork &network);
void description_parse_neuron_connection(const NeuronAddress &source_address,
        const NeuronAddress &target_address, const ryml::Parser &parser,
        ryml::ConstNodeRef attributes_node, SpikingNetwork &net);

void description_parse_hyperedge(const NeuronAddress &source_address,
        const NeuronAddress &target_address, const ryml::Parser &parser,
        ryml::ConstNodeRef hyperedge_node, SpikingNetwork &net);
void yaml_parse_conv2d(NeuronGroup &source_group, const ryml::Parser &parser,
        ryml::ConstNodeRef hyperedge_node, NeuronGroup &target_group);
void yaml_parse_sparse(NeuronGroup &source_group, const ryml::Parser &parser,
        ryml::ConstNodeRef hyperedge_node, NeuronGroup &target_group);
void yaml_parse_dense(NeuronGroup &source_group, const ryml::Parser &parser,
        ryml::ConstNodeRef hyperedge_node, NeuronGroup &target_group);
bool yaml_parse_conv2d_attributes(ryml::ConstNodeRef attribute, Conv2DParameters &convolution);
void yaml_parse_unit_specific_attributes(const ryml::Parser &parser, ryml::ConstNodeRef parent_node, std::map<std::string, std::vector<ModelAttribute>> &attribute_lists);

void description_parse_edge_attributes(Connection &edge,
        const ryml::Parser &parser, ryml::ConstNodeRef attributes_node);
std::tuple<NeuronAddress, NeuronAddress> description_parse_edge_description(const std::string_view &description);

// Functions for writing YAML
void yaml_write_network(
        std::filesystem::path path, const sanafe::SpikingNetwork &network);
void yaml_write_mappings_file(
        std::filesystem::path path, const SpikingNetwork &network);
void yaml_create_mappings(ryml::NodeRef &node, std::vector<std::reference_wrapper<const Neuron>> &all_neurons, std::list<std::string> &strings);
c4::yml::NodeRef yaml_serialize_network(c4::yml::NodeRef root, const sanafe::SpikingNetwork &network, std::list<std::string> &strings);
c4::yml::NodeRef yaml_serialize_neuron_group(c4::yml::NodeRef parent, const NeuronGroup &group, std::list<std::string> &strings);
ryml::NodeRef yaml_serialize_neuron_run(ryml::NodeRef neurons_node, const std::tuple<int, int> &neuron_run,const NeuronGroup &group, std::list<std::string> &strings);
c4::yml::NodeRef yaml_serialize_model_attributes(const std::map<std::string, ModelAttribute> &default_values, c4::yml::NodeRef parent, const std::map<std::string, ModelAttribute> &attributes);

std::string write_edge_format(const Connection &connection);

}

#endif