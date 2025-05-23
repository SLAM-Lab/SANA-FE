// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef YAML_ARCH_HEADER_INCLUDED_
#define YAML_ARCH_HEADER_INCLUDED_

#include <cstdio>
#include <fstream>
#include <iostream>
#include <optional>
#include <string_view>
#include <variant>
#include <vector>

#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "arch.hpp"
#include "fwd.hpp"

namespace sanafe
{
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
CorePipelineConfiguration description_parse_core_pipeline_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
TilePowerMetrics description_parse_tile_metrics_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
NetworkOnChipConfiguration description_parse_noc_configuration_yaml(const ryml::Parser &parser, ryml::ConstNodeRef noc_attributes);
}

#endif