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
void description_parse_axon_in_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef axon_in_node, Core &parent_core);
AxonInPowerMetrics description_parse_axon_in_attributes_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
void description_parse_synapse_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef synapse_node, Core &parent_core);
std::pair<SynapsePowerMetrics, ModelInfo> description_parse_synapse_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes);
void description_parse_dendrite_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef dendrite_node, Core &parent_core);
void description_parse_soma_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef soma_node, Core &parent_core);
void description_parse_axon_out_section(const ryml::Parser &parser, ryml::ConstNodeRef axon_out_node, Core &parent_core);
std::pair<size_t, size_t> description_parse_range_yaml(const std::string &range_str);
CorePipelineConfiguration description_parse_core_pipeline_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
TilePowerMetrics description_parse_tile_metrics_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
NetworkOnChipConfiguration description_parse_noc_configuration_yaml(const ryml::Parser &parser, ryml::ConstNodeRef noc_attributes);

// Network description
Network description_parse_network_file_yaml(std::ifstream &fp, Architecture &arch);
Network description_parse_network_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef net_node);
void description_parse_neuron_group_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef groups_node, Network &net);
void description_parse_edges_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef edges_node, Network &net);
void description_parse_group(const ryml::Parser &parser, ryml::ConstNodeRef neuron_group_node, Network &net);
void description_parse_neuron_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef neuron_node, NeuronGroup &neuron_group);
// Neuron to Hardware Mapping
void description_parse_mapping_section_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef mappings_node, Architecture &arch, Network &net);
void description_parse_mapping(const ryml::Parser &parser, Neuron &neuron, const ryml::ConstNodeRef mapping_info, Architecture &arch);

void description_parse_neuron(const std::string &id, const ryml::Parser &parser, const ryml::ConstNodeRef attributes, NeuronGroup &neuron_group);
NeuronTemplate description_parse_neuron_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes, const NeuronTemplate &default_template = NeuronTemplate());

void description_parse_edge(const std::string &description, const ryml::Parser &parser, const ryml::ConstNodeRef edge_node, Network &network);
void description_parse_edge_attributes(Connection &con, const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node);
std::tuple<NeuronAddress, NeuronAddress> description_parse_edge_description(const std::string_view &description);

// Helper functions
std::map<std::string, ModelParam> description_parse_model_parameters_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef parameters_node);
ModelParam description_parse_parameter_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
void check_key(const ryml::Parser &parser, const ryml::ConstNodeRef node, const std::string &key);
}

#endif
