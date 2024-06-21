#ifndef DESCRIPTION_HEADER_INCLUDED_
#define DESCRIPTION_HEADER_INCLUDED_

#include <cstdio>
#include <fstream>
#include <iostream>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "arch.hpp"
#include "network.hpp"

namespace sanafe
{

class DescriptionParsingError : public std::invalid_argument
{
    // Custom exception class for architecture and network description parsing
public:
    DescriptionParsingError(const std::string &error, const YAML::Mark &pos);
    [[nodiscard]] const char *what() const noexcept override;
private:
    std::string message;
};

template <typename T>
T description_required_field(const YAML::Node &node, const std::string &key);

template <typename T>
std::string description_get_type_string(const T &value);

Architecture description_parse_arch_file_yaml(std::ifstream &fp);
Architecture description_parse_arch_section_yaml(const YAML::Node &arch_node);
void description_parse_tile_section_yaml(const YAML::Node &tile_node, Architecture &arch);
void description_parse_core_section_yaml(const YAML::Node &core_node, size_t parent_tile_id, Architecture &arch);
void description_parse_axon_in_section_yaml(const YAML::Node &axon_in_node, Core &parent_core);
AxonInPowerMetrics description_parse_axon_in_attributes_yaml(const YAML::Node &attributes, const Core &parent_core);
void description_parse_synapse_section_yaml(const YAML::Node &synapse_node, Core &parent_core);
std::pair<SynapsePowerMetrics, ModelInfo> description_parse_synapse_attributes_yaml(const YAML::Node &synapse_node, Core &parent_core);
void description_parse_dendrite_section_yaml(const YAML::Node &dendrite_node, Core &parent_core);
void description_parse_soma_section_yaml(const YAML::Node &soma_node, Core &parent_core);
void description_parse_axon_out_section(const YAML::Node &axon_out_node, Core &parent_core);
std::pair<size_t, size_t> description_parse_range_yaml(const std::string &range_str);
CorePipelineConfiguration description_parse_core_pipeline_yaml(const YAML::Node &attributes);
TilePowerMetrics description_parse_tile_metrics_yaml(const YAML::Node &attributes);
NetworkOnChipConfiguration description_parse_noc_configuration_yaml(const YAML::Node &noc_attributes);

Network description_parse_net_file_yaml(std::ifstream &fp, Architecture &arch);
Network description_parse_net_section_yaml(const YAML::Node &net_node, Architecture &arch);
void description_parse_neuron_group_section_yaml(const YAML::Node &neuron_group_node, Architecture &arch, Network &net);
void description_parse_edge_section_yaml(const YAML::Node &edge_node, Network &net);
void description_parse_neuron_section_yaml(const YAML::Node &neuron_node, Architecture &arch, NeuronGroup &neuron_group);
NeuronTemplate description_parse_neuron_attributes_yaml(const YAML::Node &attributes, const NeuronTemplate &default_template = NeuronTemplate());
void description_parse_neuron_mapping_subsection_yaml(const YAML::Node &mapping_node, Architecture &arch, Neuron &neuron);
void description_parse_edge_mapping_subsection_yaml(const YAML::Node &mapping_node, Connection &edge);

std::map<std::string, ModelParam> description_parse_model_parameters_yaml(const YAML::Node &parameters_node);
ModelParam description_parse_parameter_yaml(const YAML::Node &attribute_node);
size_t description_count_neurons_yaml(const YAML::Node &neurons_node);
std::string description_yaml_parsing_error(const std::string &detail, const YAML::Mark &location);

}

#endif
