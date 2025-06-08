// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  plugins.cpp
#include <filesystem>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <dlfcn.h>

#include "pipeline.hpp"
#include "plugins.hpp"
#include "print.hpp"

using create_hw = sanafe::PipelineUnit *();

namespace // anonymous
{
// Manage the different plugins and their corresponding factory routines. For
//  now, use a couple of global maps (ignoring any clang lint warnings).
//  Probably not the cleanest or most modern, but it works and should be self-
//  contained in this file.
std::map<std::string, create_hw *>
        plugin_create_hw; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

// Use a unique_ptr with the custom deleter to automatically manage the library
//  handle
using DlHandlePtr = std::unique_ptr<void, sanafe::DlHandleDeleter>;
std::unordered_map<std::string, DlHandlePtr>
        plugin_handles; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)
}

void sanafe::DlHandleDeleter::operator()(void *handle) const
{
    if (handle != nullptr)
    {
        dlclose(handle);
    }
}

void sanafe::plugin_init_hw(
        const std::string &model_name, const std::filesystem::path &plugin_path)
{
    const std::string create = "create_" + model_name;

    // Load the soma library
    INFO("Loading plugin:%s\n", plugin_path.c_str());
    void *hw = dlopen(plugin_path.c_str(), RTLD_LAZY | RTLD_GLOBAL);
    plugin_handles[model_name] = DlHandlePtr(hw);
    if (hw == nullptr)
    {
        INFO("Error: Couldn't load library %s\n", plugin_path.c_str());
        throw std::runtime_error("Error: Could not load library.\n");
    }

    // Reset DLL errors
    dlerror();

    // Function to create an instance of the Soma class
    INFO("Loading function: %s\n", create.c_str());
    // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
    auto *create_func =
            reinterpret_cast<create_hw *>(dlsym(hw, create.c_str()));
    plugin_create_hw[model_name] = create_func;
    // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

    const char *dlsym_error = dlerror();
    if (dlsym_error != nullptr)
    {
        INFO("Error: Couldn't load symbol %s: %s\n", create.c_str(),
                dlsym_error);
        // This will also automatically close the library through its unique_ptr
        plugin_handles.erase(model_name);
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

    TRACE1(PLUGINS, "Getting model:%s\n", model_name.c_str());
    if (plugin_create_hw.count(model_name) == 0)
    {
        plugin_init_hw(model_name, plugin_path);
    }

    return std::shared_ptr<PipelineUnit>(plugin_create_hw[model_name]());
}
