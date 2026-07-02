// Copyright (c) 2026 - The University of Texas at Austin
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

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include "pipeline.hpp"
#include "plugins.hpp"
#include "print.hpp"

using create_hw = sanafe::PipelineUnit *();

namespace // anonymous
{
#ifdef _WIN32
inline void *plugin_dlopen(const char *path)
{
    return reinterpret_cast<void *>(LoadLibraryA(path));
}
inline void *plugin_dlsym(void *handle, const char *symbol)
{
    return reinterpret_cast<void *>(
            GetProcAddress(reinterpret_cast<HMODULE>(handle), symbol));
}
inline void plugin_dlclose(void *handle)
{
    FreeLibrary(reinterpret_cast<HMODULE>(handle));
}
inline std::string plugin_dlerror()
{
    const DWORD err = GetLastError();
    if (err == 0)
    {
        return {};
    }
    LPSTR buf = nullptr;
    const size_t len = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER |
                    FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
            nullptr, err, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
            reinterpret_cast<LPSTR>(&buf), 0, nullptr);
    std::string msg(buf, len);
    LocalFree(buf);
    return msg;
}
#else
inline void *plugin_dlopen(const char *path)
{
    return dlopen(path, RTLD_LAZY | RTLD_LOCAL);
}
inline void *plugin_dlsym(void *handle, const char *symbol)
{
    return dlsym(handle, symbol);
}
inline void plugin_dlclose(void *handle)
{
    dlclose(handle);
}
[[maybe_unused]] inline std::string plugin_dlerror()
{
    const char *err = dlerror();
    return (err != nullptr) ? std::string(err) : std::string();
}
#endif

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
} // anonymous namespace

void sanafe::DlHandleDeleter::operator()(void *handle) const
{
    if (handle != nullptr)
    {
        plugin_dlclose(handle);
    }
}

void sanafe::plugin_init_hw(
        const std::string &model_name, const std::filesystem::path &plugin_path)
{
    const std::string create = "create_" + model_name;

    // Load the soma library
    INFO("Loading plugin:%s\n", plugin_path.string().c_str());
    void *hw = plugin_dlopen(plugin_path.string().c_str());
    plugin_handles[model_name] = DlHandlePtr(hw);
    if (hw == nullptr)
    {
        INFO("Error: Couldn't load library %s\n", plugin_path.string().c_str());
        throw std::runtime_error("Error: Could not load library.\n");
    }

    // Function to create an instance of the Soma class
    INFO("Loading function: %s\n", create.c_str());
    // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
    auto *create_func =
            reinterpret_cast<create_hw *>(plugin_dlsym(hw, create.c_str()));
    plugin_create_hw[model_name] = create_func;
    // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

    // cppcheck-suppress knownConditionTrueFalse
    if (hw == nullptr)
    {
        INFO("Error: Couldn't load library %s: %s\n",
                plugin_path.string().c_str(), plugin_dlerror().c_str());
        throw std::runtime_error("Error: Could not load library.\n");
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
