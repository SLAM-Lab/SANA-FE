#ifndef PRINT_HEADER_INCLUDED_
#define PRINT_HEADER_INCLUDED_

#include "stdio.h"

// DEBUG and TRACE print with source annotations, TRACE is only enabled for
//  verbose debug printing
#define INFO(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
#ifdef DEBUG
#define TRACE(...) do { fprintf(stdout, "[%s:%d:%s()] ", __FILE__, __LINE__, __func__); fprintf(stdout, __VA_ARGS__); } while (0)
#else
#define TRACE(...) do {} while (0)
#endif

#endif