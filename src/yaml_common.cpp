// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>
#include <utility>
#include <variant>
#include <vector>

// False positives: clang-tidy is unable to see RapidYAML headers
#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "attribute.hpp"
#include "print.hpp"
#include "yaml_common.hpp"

// NOLINTBEGIN(misc-include-cleaner)
sanafe::YamlDescriptionParsingError::YamlDescriptionParsingError(
        const std::string &error, const ryml::Parser &parser,
        const ryml::ConstNodeRef &node)
        : std::invalid_argument(error)
{
    const ryml::Location pos = parser.location(node);
    message = "Error: " + error + " (Line " + std::to_string(pos.line + 1) +
            ':' + std::to_string(pos.col + 1) + ").\n";
}
// NOLINTEND(misc-include-cleaner)

const char *sanafe::YamlDescriptionParsingError::what() const noexcept
{
    return message.c_str();
}

void sanafe::check_key(const ryml::Parser &parser,
        const ryml::ConstNodeRef node, const std::string &key)
{
    if (!node.is_map())
    {
        throw YamlDescriptionParsingError(
                "Node should be a mapping\n. For more info on YAML mappings "
                "refer to the YAML 1.2 specification, 7.4.2 'Flow Mappings' "
                "and 8.2.2 'Block Mappings'",
                parser, node);
    }
    const ryml::ConstNodeRef child = node.find_child(key.c_str());
    if (child.invalid())
    {
        const std::string message = "Value for key '" + key + "' not defined";
        throw YamlDescriptionParsingError(message, parser, node);
    }
}

std::string sanafe::yaml_get_type_string(const std::type_info &type)
{
    if (type == typeid(bool))
    {
        return "bool";
    }
    if (type == typeid(int))
    {
        return "int";
    }
    if (type == typeid(size_t))
    {
        return "size_t";
    }
    if (type == typeid(double))
    {
        return "double";
    }
    if (type == typeid(std::string))
    {
        return "string";
    }

    // Not a scalar type; fall back to default name which may be mangled
    return type.name();
}

std::map<std::string, sanafe::ModelAttribute>
sanafe::description_parse_model_attributes_yaml( // NOLINT(misc-no-recursion)
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node)
{
    std::map<std::string, ModelAttribute> model_attributes;

    if (attributes_node.is_seq())
    {
        for (const auto &mapping_list_node : attributes_node)
        {
            // Recursive call to flatten list of mappings
            auto new_attributes = description_parse_model_attributes_yaml(
                    parser, mapping_list_node);
            model_attributes.insert(
                    new_attributes.begin(), new_attributes.end());
        }
    }
    else if (attributes_node.is_map())
    {
        for (const auto &node : attributes_node)
        {
            std::string key_str;
            node >> ryml::key(key_str); // NOLINT(misc-include-cleaner)
            std::set<std::string> unit_specific_keys = {
                    "synapse", "dendrite", "soma"};
            if (unit_specific_keys.find(key_str) == unit_specific_keys.end())
            {
                //INFO("Parsing attribute: %s\n", key.c_str());
                model_attributes[key_str] =
                        yaml_parse_attribute(parser, node);
            }
        }
    }
    else
    {
        throw std::invalid_argument(
                "Error: Model attributes must be an ordered map or mapping of "
                "named attributes.\n");
    }

    return model_attributes;
}

sanafe::ModelAttribute sanafe::yaml_parse_attribute(
        const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node)
{
    ModelAttribute attribute;

    if (attribute_node.is_seq())
    {
        attribute.value =
                yaml_parse_attribute_list(parser, attribute_node);
    }
    else if (attribute_node.is_map())
    {
        attribute.value =
                yaml_parse_attribute_map(parser, attribute_node);
    }
    else
    {
        attribute.value = yaml_parse_attribute_scalar(attribute_node);
    }

    return attribute;
}

std::vector<sanafe::ModelAttribute> sanafe::yaml_parse_attribute_list(
        const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node)
{
    // Create an list of unnamed attributes
    std::vector<ModelAttribute> attribute_list;
    for (const auto &node : attribute_node)
    {
        TRACE2(DESCRIPTION, "Parsing sub-attribute in list.\n");
        ModelAttribute curr = yaml_parse_attribute(parser, node);
        attribute_list.push_back(std::move(curr));
    }
    TRACE2(DESCRIPTION,
            "Setting attribute to an list of %zu unnamed attributes\n",
            attribute_list.size());

    return attribute_list;
}

std::vector<sanafe::ModelAttribute> sanafe::yaml_parse_attribute_map(
        const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node)
{
    std::vector<ModelAttribute> attribute_map;
    for (const auto &node : attribute_node)
    {
        TRACE2(DESCRIPTION, "Parsing mapping of attributes.\n");
        std::string key;
        node >> ryml::key(key);

        // Recursively parse YAML attribute NOLINTNEXTLINE(misc-no-recursion)
        ModelAttribute curr = yaml_parse_attribute(parser, node);
        curr.name = key;
        TRACE2(DESCRIPTION, "Saving to key: %s\n", key.c_str());
        attribute_map.push_back(std::move(curr));
    }
    TRACE1(DESCRIPTION, "Setting attribute to a list of %zu named attributes\n",
            attribute_map.size());

    return attribute_map;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
// False positive: Sequential type-checking pattern is simple despite multiple
//  return points
sanafe::AttributeVariant sanafe::yaml_parse_attribute_scalar(
        const ryml::ConstNodeRef attribute_node)
{
    if (attribute_node.invalid())
    {
        throw std::invalid_argument("Invalid node.\n");
    }
    //  Parse YAML scalar by attempting type conversions in priority order.
    //   YAML scalars are ambiguous - the same text can often be parsed a
    //   multiple types (e.g., "1" -> bool, int, double, or string;
    //   "true" -> bool or string).
    //
    //   Parsing order reflects type specificity:
    //   1. int - most restrictive numeric type
    //   2. double - accepts ints but preserves decimal values
    //   3. bool - handles true/false/yes/no values
    //   4. string - accepts anything as fallback
    //   This ensures we capture the most specific type possible while
    //   maintaining compatibility with YAML's flexible scalar representation.
    //
    // NOLINTBEGIN(misc-include-cleaner)
    // False positive: clang-tidy is unable to see the RapidYAML headers, so
    //  suppress warnings for c4::yml::read
    int decoded_int = 0;
    if (c4::yml::read(attribute_node, &decoded_int))
    {
        TRACE1(DESCRIPTION, "Parsed as int: %d.\n", decoded_int);
        return decoded_int;
    }

    double decoded_double = std::numeric_limits<double>::quiet_NaN();
    if (c4::yml::read(attribute_node, &decoded_double))
    {
        TRACE1(DESCRIPTION, "Parsed as float: %lf.\n", decoded_double);
        return decoded_double;
    }

    bool decoded_bool = false;
    if (c4::yml::read(attribute_node, &decoded_bool))
    {
        TRACE2(DESCRIPTION, "Parsed as bool: %d.\n", decoded_bool);
        return decoded_bool;
    }

    std::string decoded_str;
    if (c4::yml::read(attribute_node, &decoded_str))
    {
        TRACE2(DESCRIPTION, "Parsed as string: %s.\n", decoded_str.c_str());
        return std::move(decoded_str);
    }
    // NOLINTEND(misc-include-cleaner)

    // All parsing attempts failed - this indicates a fundamental issue
    //  with the YAML node structure or RapidYAML library state
    const std::string error_msg =
            "Failed to parse YAML scalar node as any supported type";
    INFO("Error: %s.\n", error_msg.c_str());
    throw std::runtime_error(error_msg);
}
// NOLINTEND(readability-function-cognitive-complexity)

std::pair<size_t, size_t> sanafe::yaml_parse_range(
        const std::string &range_str)
{
    constexpr std::string_view range_delimiter = "..";
    const size_t delimiter_pos = range_str.find(range_delimiter);
    if (delimiter_pos == std::string::npos)
    {
        throw std::runtime_error("Range delimiter '..' not found");
    }

    // Locate the mandatory ".." delimiter
    const auto bracket_start = range_str.find('[');
    const auto bracket_end = range_str.find(']');

    // Determine parsing boundaries (brackets are optional)
    const size_t start_pos =
            (bracket_start != std::string::npos) ? (bracket_start + 1UL) : 0UL;
    const size_t end_pos = (bracket_end != std::string::npos) ?
            bracket_end :
            range_str.length();

    // Validate that we have a reasonable range to parse
    if ((end_pos <= start_pos) || (delimiter_pos <= start_pos) ||
            (delimiter_pos >= end_pos))
    {
        throw std::runtime_error("Invalid range format");
    }

    // Extract the numeric parts before and after the delimiter
    const std::string first_str =
            range_str.substr(start_pos, delimiter_pos - start_pos);
    const std::string last_str =
            range_str.substr((delimiter_pos + range_delimiter.length()),
                    (end_pos - delimiter_pos) - range_delimiter.length());

    // Finally, convert strings to numbers
    size_t range_first = 0UL;
    size_t range_last = 0UL;
    try
    {
        range_first = std::stoull(first_str);
        range_last = std::stoull(last_str);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Invalid range string, failed to convert");
    }
    TRACE1(DESCRIPTION, "Range: %zu to %zu\n", range_first, range_last);
    if (range_first > range_last)
    {
        throw std::runtime_error("Invalid range; first > last");
    }

    return std::make_pair(range_first, range_last);
}
