// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// plugins.hpp
#ifndef PLUGINS_HEADER_INCLUDED_
#define PLUGINS_HEADER_INCLUDED_
#include <string>

namespace sanafe
{
class SomaModel;

void plugin_init_soma(char* name);
SomaModel *plugin_get_soma(const std::string &name);
}

#endif