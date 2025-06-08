// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// print.cpp
#include <any>
#include <map>
#include <string>

#include "print.hpp"

std::string sanafe::print_format_attributes(
        const std::map<std::string, std::any> &attr)
{
    std::string attr_str;

    for (const auto &a : attr)
    {
        const std::string &key = a.first;
        const auto value_str = std::any_cast<std::string>(a.second);

        attr_str += " ";
        attr_str += key;
        attr_str += "=";
        attr_str += value_str;
    }
    return attr_str;
}
