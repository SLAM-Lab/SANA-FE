// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef DESCRIPTION_HEADER_INCLUDED_
#define DESCRIPTION_HEADER_INCLUDED_

#include <map>

#include <ryml.hpp>

#include "attribute.hpp"

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

// Helper functions
std::map<std::string, ModelAttribute> description_parse_model_attributes_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node);
ModelAttribute description_parse_attribute_yaml(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
std::vector<ModelAttribute> description_parse_attribute_list(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
std::vector<ModelAttribute> description_parse_attribute_map(const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node);
std::variant<bool, int, double, std::string, std::vector<ModelAttribute>> description_parse_attribute_scalar(const ryml::ConstNodeRef attribute_node);
void check_key(const ryml::Parser &parser, const ryml::ConstNodeRef node, const std::string &key);
ryml::NodeRef description_serialize_variant_value_to_yaml(ryml::NodeRef node, const std::variant<bool, int, double, std::string, std::vector<sanafe::ModelAttribute>> &value);
std::pair<size_t, size_t> description_parse_range_yaml(const std::string &range_str);
}

#endif
