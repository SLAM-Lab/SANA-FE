// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// utils.hpp
#ifndef UTILS_HEADER_INCLUDED_
#define UTILS_HEADER_INCLUDED_

#include <chrono>
#include <cstddef>
#include <map>
#include <stdexcept>

namespace sanafe
{
double calculate_elapsed_time(const std::chrono::time_point<std::chrono::high_resolution_clock> &ts_start, const std::chrono::time_point<std::chrono::high_resolution_clock> &ts_end);
size_t abs_diff(size_t a, size_t b);

template <typename T> struct LookupTable
{
    std::map<size_t, T> values;

    T get(const int x) const
    {
        return get(static_cast<size_t>(x));
    }

    T get(const size_t x) const
    {
        if (values.empty())
        {
            throw std::runtime_error("Table is empty");
        }

        // Find the largest key that is <= x
        auto it = values.upper_bound(x);
        if (it == values.begin())
        {
            // x is smaller than all keys, return first value
            return values.begin()->second;
        }
        --it; // Move to the largest key <= x
        return it->second;
    }
};
}

#endif