// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef YAML_COMMON_HEADER_INCLUDED_
#define YAML_COMMON_HEADER_INCLUDED_

#include <map>

#include <ryml.hpp>

#include "attribute.hpp"

namespace sanafe
{
class YamlDescriptionParsingError : public std::invalid_argument
{
    // Custom exception class for architecture and network description parsing
public:
    YamlDescriptionParsingError(const std::string &error, const ryml::Parser &parser, const ryml::ConstNodeRef &node);
    [[nodiscard]] const char *what() const noexcept override;
private:
    std::string message;
};

template <typename T>
T yaml_required_field(const ryml::Parser &parser, ryml::ConstNodeRef node, const std::string &key);
std::string yaml_get_type_string(const std::type_info &type);

// Helper functions
std::map<std::string, ModelAttribute> description_parse_model_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node);
ModelAttribute yaml_parse_attribute(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
std::vector<ModelAttribute> yaml_parse_attribute_list(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
std::vector<ModelAttribute> yaml_parse_attribute_map(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
std::variant<bool, int, double, std::string, std::vector<ModelAttribute>> yaml_parse_attribute_scalar(const ryml::ConstNodeRef attribute_node);
void check_key(const ryml::Parser &parser, const ryml::ConstNodeRef node, const std::string &key);
ryml::NodeRef description_serialize_variant_value_to_yaml(ryml::NodeRef node, const std::variant<bool, int, double, std::string, std::vector<sanafe::ModelAttribute>> &value);
std::pair<size_t, size_t> yaml_parse_range(const std::string &range_str);
}

#endif
