// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  plugins.cpp
#include <string>
#include <iostream>
#include <map>
#include <dlfcn.h>

#include "plugins.hpp"
#include "models.hpp"
#include "print.hpp"

using namespace sanafe;

typedef SomaModel *_create_soma(void);
typedef void _destroy_soma(SomaModel *model);
std::map<std::string, _create_soma *> create_soma;
std::map<std::string, _destroy_soma *> destroy_soma;

void plugin_init_soma(const std::string &name)
{
	const std::string lib_path = "./plugins/" + name + ".so";
	const std::string create = "create_" + name;
	const std::string destroy = "destroy_" + name;

	// Load the soma library
	void *soma = dlopen(lib_path.c_str(), RTLD_LAZY | RTLD_GLOBAL);
	if (!soma)
	{
		INFO("Error: Couldn't load library %s\n", lib_path.c_str());
		exit(1);
	}

	// Reset DLL errors
	dlerror();

	// Function to create an instance of the Soma class
	_create_soma *create_func = (_create_soma *) dlsym(soma,create.c_str());
	create_soma[name] = create_func;
	const char *dlsym_error = dlerror();
	if (dlsym_error)
	{
		INFO("Error: Couldn't load symbol %s: %s\n", create.c_str(),
			dlsym_error);
		exit(1);
	}

	// Function to destroy an instance of the Soma class
	destroy_soma[name] =
		(_destroy_soma *) dlsym(soma, destroy.c_str());

	dlsym_error = dlerror();
	if (dlsym_error)
	{
		INFO("Error: Couldn't load symbol %s: %s\n", destroy.c_str(),
			dlsym_error);
		exit(1);
	}
	INFO("Loaded plugin symbols for %s.\n", name.c_str());
}

SomaModel *plugin_get_soma(const std::string &name)
{
	INFO("Getting soma:%s\n", name.c_str());

	if (!create_soma.count(name))
	{
		plugin_init_soma(name);
	}

	return create_soma[name]();
}
