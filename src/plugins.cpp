// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  plugins.cpp
#include <filesystem>
#include <iostream>
#include <map>
#include <string>

#include <dlfcn.h>

#include "models.hpp"
#include "plugins.hpp"
#include "print.hpp"

using _create_hw = sanafe::PipelineUnit *();

std::map<std::string, _create_hw *> plugin_create_hw;

void sanafe::plugin_init_hw(
        const std::string &model_name, const std::filesystem::path &plugin_path)
{
    const std::string create = "create_" + model_name;

    // Load the soma library
    INFO("Loading synapse plugin:%s\n", plugin_path.c_str());
    void *hw = dlopen(plugin_path.c_str(), RTLD_LAZY | RTLD_GLOBAL);
    if (hw == nullptr)
    {
        INFO("Error: Couldn't load library %s\n", plugin_path.c_str());
        throw std::runtime_error("Error: Could not load library.\n");
    }

    // Reset DLL errors
    dlerror();

    // Function to create an instance of the Soma class
    INFO("Loading function: %s\n", create.c_str());
    auto *create_func = (_create_hw *) dlsym(hw, create.c_str());
    plugin_create_hw[model_name] = create_func;

    const char *dlsym_error = dlerror();
    if (dlsym_error != nullptr)
    {
        INFO("Error: Couldn't load symbol %s: %s\n", create.c_str(),
                dlsym_error);
        throw std::runtime_error("Error: Could not load symbol.\n");
    }
    INFO("Loaded plugin symbols for %s.\n", model_name.c_str());
}

std::shared_ptr<sanafe::PipelineUnit> sanafe::plugin_get_hw(
        const std::string &model_name, const std::filesystem::path &plugin_path)
{
    if (plugin_path.empty())
    {
        throw std::runtime_error("No plugin path given.");
    }

    TRACE1(PLUGINS, "Getting synapse:%s\n", model_name.c_str());
    if (plugin_create_hw.count(model_name) == 0)
    {
        plugin_init_hw(model_name, plugin_path);
    }

    return std::shared_ptr<PipelineUnit>(plugin_create_hw[model_name]());
}
