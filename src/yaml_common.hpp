// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef YAML_COMMON_HEADER_INCLUDED_
#define YAML_COMMON_HEADER_INCLUDED_

#include <cstddef>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <variant>
#include <vector>

#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <c4/yml/fwd.hpp>
#include <c4/yml/node.hpp>

#include "fwd.hpp"

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

std::string yaml_get_type_string(const std::type_info &type);

template <typename T>
T yaml_required_field(const ryml::Parser &parser,
        const ryml::ConstNodeRef node, const std::string &key)
{
    // Wrapper around YAML library for field=map[key], adding more error prints
    if (node.invalid())
    {
        const std::string message = "Invalid node when looking up key: " + key;
        throw std::runtime_error(message);
    }
    const ryml::ConstNodeRef field_node = node.find_child(key.c_str());
    if (field_node.invalid())
    {
        const std::string message = "Key '" + key + "' does not exist";
        throw YamlDescriptionParsingError(message, parser, node);
    }
    if (!field_node.has_val())
    {
        const std::string message = "'" + key + "' value should be a scalar";
        throw YamlDescriptionParsingError(message, parser, field_node);
    }

    T field{};
    // Efficiently convert to type T by trying the RapidYAML reader.
    //  If read() fails, it returns false and execution falls through
    if (c4::yml::read(field_node, &field))
    {
        return field; // type T
    }

    const std::string message = "Could not cast field '" + key +
            "' to type: " + yaml_get_type_string(typeid(field));
    throw YamlDescriptionParsingError(message, parser, field_node);
}

template <typename T>
std::optional<T> yaml_optional_field(
        ryml::ConstNodeRef node, const std::string &key)
{
    // A terser way of getting optionally set values from YAML
    std::optional<T> optional_value;
    if (!(node.find_child(key.c_str()).invalid()))
    {
        // A few redundant steps needed so that RapidYAML doesn't complain
        T value;
        node[key.c_str()] >> value;
        optional_value = value;
    }

    return optional_value;
}

// Helper functions
std::map<std::string, ModelAttribute> description_parse_model_attributes_yaml(const ryml::Parser &parser, ryml::ConstNodeRef attributes_node);
ModelAttribute yaml_parse_attribute(const ryml::Parser &parser, ryml::ConstNodeRef attribute_node);
std::vector<ModelAttribute> yaml_parse_attribute_list(const ryml::Parser &parser, ryml::ConstNodeRef attribute_node);
std::vector<ModelAttribute> yaml_parse_attribute_map(const ryml::Parser &parser, ryml::ConstNodeRef attribute_node);
std::variant<bool, int, double, std::string, std::vector<ModelAttribute>>
yaml_parse_attribute_scalar(ryml::ConstNodeRef attribute_node);
void check_key(const ryml::Parser &parser, ryml::ConstNodeRef node, const std::string &key);
ryml::NodeRef description_serialize_variant_value_to_yaml(ryml::NodeRef node, const std::variant<bool, int, double, std::string, std::vector<sanafe::ModelAttribute>> &value);
std::pair<size_t, size_t> yaml_parse_range(const std::string &range_str);
}

#endif
