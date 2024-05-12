// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <sstream>
#include <string>
#include <unordered_map>

#include "print.hpp"

std::string sanafe::print_format_attributes(
	const std::unordered_map<std::string, std::string> &attr)
{
	std::string attr_str;

	for (const auto &a: attr)
	{
		const std::string &key = a.first;
		const std::string &value_str = a.second;

		attr_str += ' ' + key + '=' + value_str;
	}
	return attr_str;
}

std::string sanafe::print_float(const double val)
{
	std::ostringstream ss;
	ss << val;
	return ss.str();
}

std::string sanafe::print_int(const long int val)
{
	std::ostringstream ss;
	ss << val;
	return ss.str();
}
