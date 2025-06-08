// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// plugins.hpp
#ifndef PLUGINS_HEADER_INCLUDED_
#define PLUGINS_HEADER_INCLUDED_
#include <filesystem>
#include <memory>
#include <string>

#include "fwd.hpp"

namespace sanafe
{
struct DlHandleDeleter
{
    void operator()(void *handle) const;
};
void plugin_init_hw(const std::string &model_name, const std::filesystem::path &plugin_path);
std::shared_ptr<PipelineUnit> plugin_get_hw(const std::string &model_name, const std::filesystem::path &plugin_path);
}

#endif
