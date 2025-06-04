// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// print.hpp
#ifndef PRINT_HEADER_INCLUDED_
#define PRINT_HEADER_INCLUDED_

#include <any>
#include <cstdio>
#include <map>
#include <string>

namespace sanafe
{
std::string print_format_attributes(const std::map<std::string, std::any> &attr);
}

// DEBUG and TRACE<n> print with source annotations, TRACE<n> levels offer more
//  detailed and verbose output for debugging. If you want to add a new category
//  for tracing you must follow the following steps:
//  1) Below, add your new category macro e.g., #define DEBUG_CATEGORY_FOO <N+1>
//  2) In CMakeLists.txt, set the default debug level for the new category
//  3) In CMakeLists.txt, add the new category to the validation loop, which
//      adds a compiler flag for each category i.e., add your category inside
//      foreach(category <append to this list>)

// Define categories
// NOLINTBEGIN
#define DEBUG_CATEGORY_CHIP 0
#define DEBUG_CATEGORY_ARCH 1
#define DEBUG_CATEGORY_NET 2
#define DEBUG_CATEGORY_PYMODULE 3
#define DEBUG_CATEGORY_DESCRIPTION 4
#define DEBUG_CATEGORY_MODELS 5
#define DEBUG_CATEGORY_SCHEDULER 6
#define DEBUG_CATEGORY_PLUGINS 7

// Default debug levels; CMake can override these defaults later
#ifndef DEBUG_LEVEL_ARCH
#define DEBUG_LEVEL_ARCH 0
#endif
#ifndef DEBUG_LEVEL_NET
#define DEBUG_LEVEL_NET 0
#endif
#ifndef DEBUG_LEVEL_PYMODULE
#define DEBUG_LEVEL_PYMODULE 0
#endif
#ifndef DEBUG_LEVEL_DESCRIPTION
#define DEBUG_LEVEL_DESCRIPTION 0
#endif
#ifndef DEBUG_LEVEL_MODELS
#define DEBUG_LEVEL_MODELS 0
#endif
#ifndef DEBUG_LEVEL_PLUGINS
#define DEBUG_LEVEL_PLUGINS 0
#endif
#ifndef DEBUG_LEVEL_SCHEDULER
#define DEBUG_LEVEL_SCHEDULER 0
#endif
#ifndef DEBUG_LEVEL_CHIP
#define DEBUG_LEVEL_CHIP 0
#endif
// NOLINTEND

// Conditional source info based on build configuration
//NOLINTBEGIN
#ifdef ENABLE_SOURCE_INFO
    #define SOURCE_INFO() fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__)
#else
    #define SOURCE_INFO() do {} while (0)
#endif
// NOLINTEND

//NOLINTBEGIN
#ifdef ENABLE_DEBUG_PRINTS
    #define INFO(...) do { \
        SOURCE_INFO(); \
        fprintf(stdout, __VA_ARGS__); \
    } while (0)

    #define TRACE1(category, ...) do { \
        if (DEBUG_LEVEL_##category >= 1) { \
            SOURCE_INFO(); \
            fprintf(stdout, __VA_ARGS__); \
        } \
    } while (0)

    #define TRACE2(category, ...) do { \
        if (DEBUG_LEVEL_##category >= 2) { \
            SOURCE_INFO(); \
            fprintf(stdout, __VA_ARGS__); \
        } \
    } while (0)

    #define TRACE3(category, ...) do { \
        if (DEBUG_LEVEL_##category >= 3) { \
            SOURCE_INFO(); \
            fprintf(stdout, __VA_ARGS__); \
        } \
    } while (0)
#else
    // No-op versions for stripped builds
    #define INFO(...) do {} while (0)
    #define TRACE1(category, ...) do {} while (0)
    #define TRACE2(category, ...) do {} while (0)
    #define TRACE3(category, ...) do {} while (0)
#endif
// NOLINTEND

#endif
