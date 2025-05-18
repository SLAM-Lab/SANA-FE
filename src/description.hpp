// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef DESCRIPTION_HEADER_INCLUDED_
#define DESCRIPTION_HEADER_INCLUDED_

#include <cstdio>
#include <fstream>
#include <iostream>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

#include <ryml.hpp>

#include "arch.hpp"
#include "network.hpp"

namespace sanafe
{

class DescriptionParsingError : public std::invalid_argument
{
    // Custom exception class for architecture and network description parsing
public:
    DescriptionParsingError(const std::string &error, const ryml::Parser &parser, const ryml::ConstNodeRef &node);
    [[nodiscard]] const char *what() const noexcept override;
private:
    std::string message;
};

template <typename T>
T description_required_field(const ryml::Parser &parser, ryml::ConstNodeRef node, const std::string &key);

std::string description_get_type_string(const std::type_info &type);

// Architecture description
Architecture description_parse_arch_file_yaml(std::ifstream &fp);
Architecture description_parse_arch_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef arch_node);
void description_parse_tile_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef tile_node, Architecture &arch);
void description_parse_core_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef core_node, size_t parent_tile_id, Architecture &arch);
void description_parse_axon_in_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef axon_in_node, CoreConfiguration &parent_core);
AxonInPowerMetrics description_parse_axon_in_attributes_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
void description_parse_synapse_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef synapse_node, CoreConfiguration &parent_core);
ModelInfo description_parse_synapse_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes);
void description_parse_dendrite_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef dendrite_node, CoreConfiguration &parent_core);
void description_parse_soma_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef soma_node, CoreConfiguration &parent_core);
void description_parse_axon_out_section(const ryml::Parser &parser, ryml::ConstNodeRef axon_out_node, CoreConfiguration &parent_core);
std::pair<size_t, size_t> description_parse_range_yaml(const std::string &range_str);
CorePipelineConfiguration description_parse_core_pipeline_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
TilePowerMetrics description_parse_tile_metrics_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
NetworkOnChipConfiguration description_parse_noc_configuration_yaml(const ryml::Parser &parser, ryml::ConstNodeRef noc_attributes);

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
void description_parse_edge(const std::string &description, const ryml::Parser &parser, const ryml::ConstNodeRef edge_node, SpikingNetwork &network);
void description_parse_neuron_connection(const NeuronAddress &source_address, const NeuronAddress &dest_address, const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node, SpikingNetwork &net);
void description_parse_hyperedge(const NeuronAddress &source_address, const NeuronAddress &dest_address, const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node, SpikingNetwork &net);
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

// Helper for YAML serialization of variant types

// Helper functions
std::map<std::string, ModelAttribute> description_parse_model_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node);
ModelAttribute description_parse_attribute_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
void check_key(const ryml::Parser &parser, const ryml::ConstNodeRef node, const std::string &key);
ryml::NodeRef description_serialize_variant_value_to_yaml(ryml::NodeRef node, const std::variant<bool, int, double, std::string, std::vector<sanafe::ModelAttribute>> &value);

// ********** Netlist (SANA-FE v1) network description format  **********
SpikingNetwork netlist_parse_file(std::ifstream &fp, Architecture &arch);
void netlist_get_fields(std::vector<std::string_view> &fields, const std::string &line);
size_t field_to_int(const std::string_view &field);
std::pair<std::string, size_t> netlist_parse_neuron_field(const std::string_view &neuron_field);
std::pair<size_t, size_t> netlist_parse_core_field(const std::string_view &core_field);
std::tuple<std::string, size_t, std::string, size_t> netlist_parse_edge_field(const std::string_view &edge_field);
std::tuple<std::string, size_t, size_t, size_t> netlist_parse_mapping_field(const std::string_view &mapping_field);
void netlist_read_network_entry(const std::vector<std::string_view> &fields, Architecture &arch, SpikingNetwork &net, const int line_number);

std::string netlist_group_to_netlist(const NeuronGroup &group);
std::string netlist_neuron_to_netlist(const Neuron &neuron, const SpikingNetwork &net, const std::map<std::string, size_t> &group_name_to_id);
std::string netlist_mapping_to_netlist(const Neuron &neuron, const std::map<std::string, size_t> &group_name_to_id);
std::string netlist_connection_to_netlist(const Connection &con, const std::map<std::string, size_t> &group_name_to_id);
std::string netlist_attributes_to_netlist(const std::map<std::string, ModelAttribute> &model_attributes, const std::map<std::string, ModelAttribute> &default_attributes);

void netlist_read_group(const std::vector<std::string_view> &fields, SpikingNetwork &net, const int line_number);
void netlist_read_neuron(const std::vector<std::string_view> &fields, SpikingNetwork &net, const int line_number);
void netlist_read_edge(const std::vector<std::string_view> &fields, SpikingNetwork &net, const int line_number);
void netlist_read_mapping(const std::vector<std::string_view> &fields, Architecture &arch, SpikingNetwork &net, const int line_number);

std::map<std::string, ModelAttribute> netlist_parse_attributes(const std::vector<std::string_view> &attribute_fields, const int line_number);
std::map<std::string, ModelAttribute> netlist_parse_embedded_json(const std::vector<std::string_view> &attribute_fields, const int line_number);
}

#endif
