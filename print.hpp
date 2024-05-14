// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// print.hpp
#ifndef PRINT_HEADER_INCLUDED_
#define PRINT_HEADER_INCLUDED_

#include <map>
#include <string>
#include <cstdio>

namespace sanafe
{
std::string print_format_attributes(const std::map<std::string, std::string> &attr);
std::string print_float(const double val);
std::string print_int(const long int val);
}

// DEBUG and TRACE print with source annotations, TRACE is only enabled for
//  verbose debug printing.
#define INFO(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
#ifdef DEBUG
#define TRACE1_PYBIND(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
#define SIM_TRACE1(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
#define TRACE1(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
#define TRACE2(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
//#define TRACE3(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
#define TRACE3(...) do {} while (0)
#else
#define TRACE1_PYBIND(...) do {} while (0)
#define SIM_TRACE1(...) do {} while (0)
#define TRACE1(...) do {} while (0)
#define TRACE2(...) do {} while (0)
#define TRACE3(...) do {} while (0)
#endif

#endif
