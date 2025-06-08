// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef YAML_ARCH_HEADER_INCLUDED_
#define YAML_ARCH_HEADER_INCLUDED_


#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <string_view>

#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)
#include <c4/yml/fwd.hpp>
#include <c4/yml/node.hpp>

#include "arch.hpp"
#include "fwd.hpp"

namespace sanafe
{

struct PipelineUnitSectionInfo
{
    std::string_view name;
    std::function<void(const ryml::Parser &, const ryml::ConstNodeRef &, CoreConfiguration &, const std::string_view &, const std::string &)> parsing_function;
};

// Architecture description
Architecture description_parse_arch_file_yaml(std::ifstream &fp);
Architecture description_parse_arch_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef arch_node);

void description_parse_tile_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef tile_node, Architecture &arch);
void description_parse_core_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef core_node, size_t parent_tile_id, Architecture &arch);
void description_parse_core_yaml(const ryml::Parser &parser, ryml::ConstNodeRef core_node, size_t parent_tile_id, Architecture &arch, const std::string_view &name);

void yaml_parse_axon_in(const ryml::Parser &parser, ryml::ConstNodeRef axon_in_node, CoreConfiguration &parent_core,  const std::string_view & /*type*/, const std::string &name);
void yaml_parse_axon_out(const ryml::Parser &parser, ryml::ConstNodeRef axon_out_node, CoreConfiguration &parent_core, const std::string_view & /*type*/, const std::string &name);
void yaml_parse_processing_unit(const ryml::Parser &parser, ryml::ConstNodeRef synapse_node, CoreConfiguration &parent_core, const std::string_view &type, const std::string &name);

AxonInPowerMetrics yaml_parse_axon_in_attributes(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
AxonOutPowerMetrics yaml_parse_axon_out_attributes(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
ModelInfo yaml_parse_processing_unit_attributes(const ryml::Parser &parser, const ryml::ConstNodeRef &attributes);

void yaml_merge_or_create_hardware_unit(CoreConfiguration &parent_core, const std::string &name, ModelInfo &model_details, const std::string_view &section);
void yaml_set_implements_flag(PipelineUnitConfiguration &hw, const std::string_view &section);

template<typename ParseFunc>
void yaml_parse_pipeline_entry(const ryml::Parser &parser, const ryml::ConstNodeRef &unit_node, CoreConfiguration &parent_core, const std::string_view &type, ParseFunc parsing_function);

//void description_parse_synapse_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef synapse_node, CoreConfiguration &parent_core);
//ModelInfo description_parse_synapse_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes);
//void description_parse_dendrite_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef dendrite_node, CoreConfiguration &parent_core);
//void description_parse_dendrite_yaml(const ryml::Parser &parser, ryml::ConstNodeRef dendrite_node, CoreConfiguration &parent_core, std::string name);
//sanafe::ModelInfo description_parse_dendrite_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes);
//void description_parse_soma_section_yaml(const ryml::Parser &parser, ryml::ConstNodeRef soma_node, CoreConfiguration &parent_core);

CorePipelineConfiguration description_parse_core_pipeline_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
TilePowerMetrics description_parse_tile_metrics_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes);
NetworkOnChipConfiguration description_parse_noc_configuration_yaml(const ryml::Parser &parser, ryml::ConstNodeRef noc_attributes);
LookupTable<double> yaml_parse_sync_delay_table(const ryml::Parser &parser, const ryml::ConstNodeRef &noc_attributes, const std::string &model_type);
}

#endif