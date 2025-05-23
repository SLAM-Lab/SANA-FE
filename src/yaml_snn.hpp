// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef YAML_SNN_HEADER_INCLUDED_
#define YAML_SNN_HEADER_INCLUDED_

#include <cstdio>
#include <fstream>
#include <iostream>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "fwd.hpp"
#include "network.hpp"

namespace sanafe
{
// YAML network description
SpikingNetwork description_parse_network_file_yaml(std::ifstream &fp, Architecture &arch);
SpikingNetwork description_parse_network_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef net_node);
void description_parse_neuron_group_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef groups_node, SpikingNetwork &net);
size_t description_count_neurons(const ryml::Parser &parser, ryml::ConstNodeRef neuron_node);
void description_parse_edges_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef edges_node, SpikingNetwork &net);
void description_parse_group(const ryml::Parser &parser, ryml::ConstNodeRef neuron_group_node, SpikingNetwork &net);
void description_parse_neuron_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef neuron_node, NeuronGroup &neuron_group);
// Neuron to Hardware Mapping
void description_parse_mapping_section_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef mappings_node, Architecture &arch, SpikingNetwork &net);
void description_parse_mapping(const ryml::Parser &parser, const ryml::ConstNodeRef mapping_info, Architecture &arch, SpikingNetwork &net);
void description_parse_mapping_info(const ryml::Parser &parser, const ryml::ConstNodeRef info, Neuron &n, std::string &core_name);
void description_parse_neuron(const std::string &id, const ryml::Parser &parser, const ryml::ConstNodeRef attributes, NeuronGroup &neuron_group);
NeuronConfiguration description_parse_neuron_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes, const NeuronConfiguration &default_template = NeuronConfiguration());
void description_parse_edge(const std::string &description, const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node, SpikingNetwork &network);
void description_parse_neuron_connection(const NeuronAddress &source_address, const NeuronAddress &target_address, const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node, SpikingNetwork &net);
void description_parse_hyperedge(const NeuronAddress &source_address, const NeuronAddress &target_address, const ryml::Parser &parser, const ryml::ConstNodeRef hyperedge_node, SpikingNetwork &net);
void description_parse_edge_attributes(Connection &edge, const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node);
std::tuple<NeuronAddress, NeuronAddress> description_parse_edge_description(const std::string_view &description);

// Functions for writing YAML
void description_write_network_yaml(const std::filesystem::path path, const sanafe::SpikingNetwork &network);
void description_write_mappings_yaml(const std::filesystem::path path, const SpikingNetwork &network);
c4::yml::NodeRef description_serialize_network_to_yaml(c4::yml::NodeRef root, const sanafe::SpikingNetwork &network, std::list<std::string> &strings);
c4::yml::NodeRef description_serialize_neuron_group_to_yaml(c4::yml::NodeRef parent, const NeuronGroup &group, std::list<std::string> &strings);
c4::yml::NodeRef description_serialize_neuron_to_yaml(c4::yml::NodeRef parent, const Neuron &neuron, std::list<std::string> &strings);
c4::yml::NodeRef description_serialize_connection_to_yaml(c4::yml::NodeRef parent, const Connection &connection, std::list<std::string> &strings);
c4::yml::NodeRef description_serialize_model_attributes_to_yaml(c4::yml::NodeRef parent, const std::map<std::string, ModelAttribute> &attributes, const std::map<std::string, ModelAttribute> &default_values);
c4::yml::NodeRef description_serialize_model_attribute_to_yaml(c4::yml::NodeRef parent, const ModelAttribute &attribute);
}

#endif