// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef PARAM_HEADER_INCLUDED_
#define PARAM_HEADER_INCLUDED_

#include <optional>
#include <variant>

#include "print.hpp"

namespace sanafe
{
// An attribute can contain a scalar value, or either a list or named set of
//  attributes i.e., attributes can be recursively defined. However,
//  in C++, variants cannot be defined recursively, so create this new class.
struct ModelParam
{
    operator bool() const
    {
        if (std::holds_alternative<bool>(value))
        {
            return std::get<bool>(value);
        }
        if (std::holds_alternative<int>(value))
        {
            TRACE1(ARCH, "Warning: Casting integer value to bool type.\n");
            return (std::get<int>(value) != 0);
        }

        std::string error = "Error: Attribute ";
        if (name.has_value())
        {
            error += name.value();
        }
        error += " cannot be cast to a bool";
        throw std::runtime_error(error);
    }
    operator int() const
    {
        return std::get<int>(value);
    }
    operator double() const
    {
        if (std::holds_alternative<double>(value))
        {
            return std::get<double>(value);
        }
        if (std::holds_alternative<int>(value))
        {
            // Assume it is safe to convert from any integer to double
            TRACE1(ARCH, "Warning: Casting integer value to double type.\n");
            return static_cast<double>(std::get<int>(value));
        }

        std::string error = "Error: Attribute ";
        if (name.has_value())
        {
            error += name.value();
        }
        error += " cannot be cast to a double";
        throw std::runtime_error(error);
    }

    operator std::string() const
    {
        return std::get<std::string>(value);
    }
    template <typename T> operator std::vector<T>() const
    {
        std::vector<T> cast_vector;
        const auto &value_vector = std::get<std::vector<ModelParam>>(value);
        cast_vector.reserve(value_vector.size());

        for (const auto &element : value_vector)
        {
            cast_vector.push_back(static_cast<T>(element));
        }
        return cast_vector;
    }
    template <typename T> operator std::map<std::string, ModelParam>() const
    {
        std::map<std::string, ModelParam> cast_map;
        const auto &value_vector = std::get<std::vector<ModelParam>>(value);
        for (const auto &element : value_vector)
        {
            cast_map[element.name.value()] = static_cast<T>(element);
        }
        return cast_map;
    }
    bool operator==(const ModelParam &rhs) const
    {
        return (value == rhs.value);
    }

    // In C++17, we cannot use std::map (which would be the natural choice) with
    //  incomplete types i.e., cannot use std::map in such a recursive
    //  structure. Considering this, and the fact that performance is not as
    //  important for this struct, label every attribute with a name and if the
    //  user wants to use "map" style lookups e.g., foo = attribute["key"]
    //  then support casting the struct to a std::map.
    //  There have been other discussions on this topic e.g., for implementing
    //  JSON and YAML parsers, but they end up either requiring Boost or other
    //  dependencies, and / or rely on undefined C++ behavior and generally
    //  require complex solutions.
    std::variant<bool, int, double, std::string, std::vector<ModelParam>> value;
    std::optional<std::string> name;

    // Filters control which hardware units can receive this parameter
    bool forward_to_synapse{true};
    bool forward_to_dendrite{true};
    bool forward_to_soma{true};
};
}
#endif
