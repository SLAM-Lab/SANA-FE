// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <algorithm>
#include <cassert>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <variant>
#include <vector>

#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "attribute.hpp"
#include "arch.hpp"
#include "chip.hpp"
#include "description.hpp"
#include "network.hpp"
#include "print.hpp"

sanafe::DescriptionParsingError::DescriptionParsingError(
        const std::string &error, const ryml::Parser &parser,
        const ryml::ConstNodeRef &node)
        : std::invalid_argument(error)
{
    const ryml::Location pos = parser.location(node);
    message = "Error: " + error + " (Line " + std::to_string(pos.line + 1) +
            ':' + std::to_string(pos.col + 1) + ").\n";
}

const char *sanafe::DescriptionParsingError::what() const noexcept
{
    return message.c_str();
}

void sanafe::check_key(const ryml::Parser &parser,
        const ryml::ConstNodeRef node, const std::string &key)
{
    if (!node.is_map())
    {
        throw DescriptionParsingError(
                "Node should be a mapping\n. For more info on YAML mappings "
                "refer to the YAML 1.2 specification, 7.4.2 'Flow Mappings' "
                "and 8.2.2 'Block Mappings'",
                parser, node);
    }
    const ryml::ConstNodeRef child = node.find_child(key.c_str());
    if (child.invalid())
    {
        const std::string message = "Value for key '" + key + "' not defined";
        throw DescriptionParsingError(message, parser, node);
    }
}

template <typename T>
T sanafe::description_required_field(const ryml::Parser &parser,
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
        throw DescriptionParsingError(message, parser, node);
    }
    if (!field_node.has_val())
    {
        const std::string message = "'" + key + "' value should be a scalar";
        throw DescriptionParsingError(message, parser, field_node);
    }

    T field;
    // Efficiently convert to type T by trying the RapidYAML reader.
    //  If read() fails, it returns false and execution falls through
    if (c4::yml::read(field_node, &field))
    {
        return field; // type T
    }

    const std::string message = "Could not cast field '" + key +
            "' to type: " + description_get_type_string(typeid(field));
    throw DescriptionParsingError(message, parser, field_node);
}

// Template instantiations. The alternative is to have the implemented function
//  in the header file, but I prefer it this way, as the instantiations should
//  be limited to a few basic types
template bool sanafe::description_required_field<bool>(
        const ryml::Parser &parser, ryml::ConstNodeRef node,
        const std::string &key);
template int sanafe::description_required_field<int>(const ryml::Parser &parser,
        ryml::ConstNodeRef node, const std::string &key);
template double sanafe::description_required_field<double>(
        const ryml::Parser &parser, ryml::ConstNodeRef node,
        const std::string &key);
template std::string sanafe::description_required_field<std::string>(
        const ryml::Parser &parser, ryml::ConstNodeRef node,
        const std::string &key);

std::string sanafe::description_get_type_string(const std::type_info &type)
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
sanafe::description_parse_model_attributes_yaml(
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
            node >> ryml::key(key_str);
            std::set<std::string> unit_specific_keys = {
                    "synapse", "dendrite", "soma"};
            if (unit_specific_keys.find(key_str) == unit_specific_keys.end())
            {
                //INFO("Parsing attribute: %s\n", key.c_str());
                model_attributes[key_str] =
                        description_parse_attribute_yaml(parser, node);
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


sanafe::ModelAttribute sanafe::description_parse_attribute_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attribute_node)
{
    ModelAttribute attribute;

    if (attribute_node.is_seq())
    {
        // Create an list of unnamed attributes
        std::vector<ModelAttribute> attribute_list;
        for (const auto &node : attribute_node)
        {
            TRACE2(DESCRIPTION, "Parsing sub-attribute in list.\n");
            ModelAttribute curr =
                    description_parse_attribute_yaml(parser, node);
            attribute_list.push_back(std::move(curr));
        }
        TRACE2(DESCRIPTION,
                "Setting attribute to an list of %zu unnamed attributes\n",
                attribute_list.size());
        attribute.value = attribute_list;
    }
    else if (attribute_node.is_map())
    {
        // Create a list of named attributes
        std::vector<ModelAttribute> attribute_list;
        for (const auto &node : attribute_node)
        {
            TRACE2(DESCRIPTION, "Parsing mapping of attributes.\n");
            ModelAttribute curr =
                    description_parse_attribute_yaml(parser, node);
            std::string key;
            node >> ryml::key(key);
            curr.name = key;
            TRACE2(DESCRIPTION, "Saving to key: %s\n", key.c_str());
            attribute_list.push_back(std::move(curr));
        }
        TRACE1(DESCRIPTION,
                "Setting attribute to a list of %zu named attributes\n",
                attribute_list.size());
        attribute.value = attribute_list;
    }
    else
    {
        if (attribute_node.invalid())
        {
            throw std::invalid_argument("Invalid node.\n");
        }
        int decoded_int = 0;
        double decoded_double = std::numeric_limits<double>::quiet_NaN();
        bool decoded_bool = false;
        std::string decoded_str;
        // BEGINNOLINT(misc-include-cleaner)
        if (c4::yml::read(attribute_node, &decoded_int))
        {
            TRACE1(DESCRIPTION, "Parsed int: %d.\n", decoded_int);
            attribute.value = decoded_int;
        }
        else if (c4::yml::read(attribute_node, &decoded_double))
        {
            TRACE1(DESCRIPTION, "Parsed float: %lf.\n", decoded_double);
            attribute.value = decoded_double;
        }
        else if (c4::yml::read(attribute_node, &decoded_bool))
        {
            TRACE2(DESCRIPTION, "Parsed bool: %d.\n", decoded_bool);
            attribute.value = decoded_bool;
        }
        else if (c4::yml::read(attribute_node, &decoded_str))
        {
            TRACE2(DESCRIPTION, "Parsed string: %s.\n", decoded_str.c_str());
            attribute.value = std::move(decoded_str);
        }
        // ENDNOLINT(misc-include-cleaner)
    }

    return attribute;
}

std::pair<size_t, size_t> sanafe::description_parse_range_yaml(
        const std::string &range_str)
{
    const size_t delimiter_pos = range_str.find("..");

    size_t start_pos = 0;
    size_t end_pos = 0;
    if (range_str.find('[') == std::string::npos)
    {
        start_pos = 0;
    }
    else
    {
        start_pos = range_str.find('[') + 1;
    }

    if (range_str.find(']') == std::string::npos)
    {
        end_pos = range_str.length();
    }
    else
    {
        end_pos = range_str.find(']');
    }

    if ((end_pos <= start_pos) || (delimiter_pos == std::string::npos) ||
            (delimiter_pos <= start_pos) || (delimiter_pos >= end_pos))
    {
        throw std::runtime_error("Invalid range format");
    }

    const std::string first_str =
            range_str.substr(start_pos, delimiter_pos - start_pos);
    const std::string last_str =
            range_str.substr(delimiter_pos + 2, (end_pos - delimiter_pos) - 2);

    size_t first = 0;
    size_t last = 0;
    try
    {
        first = std::stoull(first_str);
        last = std::stoull(last_str);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Invalid range string");
    }
    TRACE1(DESCRIPTION, "Range: %zu to %zu\n", first, last);
    if (first > last)
    {
        throw std::runtime_error("Invalid range; first > last");
    }

    return std::make_pair(first, last);
}
