// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// plugins.hpp
#ifndef PLUGINS_HEADER_INCLUDED_
#define PLUGINS_HEADER_INCLUDED_
#include <filesystem>
#include <string>

namespace sanafe
{
class SynapseModel;
class DendriteModel;
class SomaModel;

void plugin_init_synapse(const std::string &model_name, const std::filesystem::path &plugin_path);
void plugin_init_dendrite(const std::string &model_name, const std::filesystem::path &plugin_path);
void plugin_init_soma(const std::string &model_name, const std::filesystem::path &plugin_path);
std::shared_ptr<SynapseModel> plugin_get_synapse(const std::string &model_name, std::filesystem::path plugin_path="");
std::shared_ptr<DendriteModel> plugin_get_dendrite(const std::string &model_name, std::filesystem::path plugin_path="");
std::shared_ptr<SomaModel> plugin_get_soma(const std::string &model_name, const std::string &gid, const std::string &nid, std::filesystem::path plugin_path="");
}

#endif
