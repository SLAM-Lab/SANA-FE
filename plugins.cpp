// Copyright (c) 2024 - The University of Texas at Austin
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

typedef sanafe::SomaModel *_create_soma(const int gid, const int nid);
//typedef void _destroy_soma(std::shared_ptr<sanafe::SomaModel> model);
std::map<std::string, _create_soma *> plugin_create_soma;
//std::map<std::string, _destroy_soma *> destroy_soma;

void sanafe::plugin_init_soma(
	const std::string &model_name, const std::filesystem::path &plugin_path)
{
	const std::string create = "create_" + model_name;
	//const std::string destroy = "destroy_" + model_name;

	// Load the soma library
	INFO("Loading soma plugin:%s\n", plugin_path.c_str());
	void *soma = dlopen(plugin_path.c_str(), RTLD_LAZY | RTLD_GLOBAL);
	if (soma == nullptr)
	{
		INFO("Error: Couldn't load library %s\n", plugin_path.c_str());
		exit(1);
	}

	// Reset DLL errors
	dlerror();

	// Function to create an instance of the Soma class
	INFO("Loading function: %s\n", create.c_str());
	_create_soma *create_func = (_create_soma *) dlsym(soma, create.c_str());
	plugin_create_soma[model_name] = create_func;
	const char *dlsym_error = dlerror();
	if (dlsym_error)
	{
		INFO("Error: Couldn't load symbol %s: %s\n", create.c_str(),
			dlsym_error);
		exit(1);
	}

	// Function to destroy an instance of the Soma class
	//destroy_soma[model_name] =
	//	(_destroy_soma *) dlsym(soma, destroy.c_str());

	//dlsym_error = dlerror();
	//if (dlsym_error)
	//{
	//	INFO("Error: Couldn't load symbol %s: %s\n", destroy.c_str(),
	//		dlsym_error);
	//	exit(1);
	//}
	INFO("Loaded plugin symbols for %s.\n", model_name.c_str());
}

std::shared_ptr<sanafe::SomaModel> sanafe::plugin_get_soma(
	const std::string &model_name, const int gid, const int nid,
	std::filesystem::path plugin_path)
{
	if (plugin_path.empty())
	{
		// Use default path, which is in the plugin directory and
		const std::string model_filename = model_name + ".so";
		plugin_path =
			std::filesystem::path(".") /
			std::filesystem::path("plugins") /
			model_filename;
		INFO("No path given; setting path to default: %s\n",
			plugin_path.c_str());
	}

	INFO("Getting soma:%s\n", model_name.c_str());
	if (plugin_create_soma.count(model_name) == 0)
	{
		plugin_init_soma(model_name, plugin_path);
	}

	return std::shared_ptr<SomaModel>(
		plugin_create_soma[model_name](gid, nid));
}
