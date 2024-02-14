// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  plugins.cpp
#include <dlfcn.h>
#include <string.h>
#include <iostream>
#include <map>
#include "plugins.hpp"
#include "print.hpp"

std::map<std::string, _create_soma *> create_soma;
std::map<std::string, _destroy_soma *> destroy_soma;

void init_soma(char* name){
    char lib_path[11 + MAX_SOMA_LEN] = "./plugins/";
    strcat(lib_path, name);
    strcat(lib_path, ".so");
    char create[8 + MAX_SOMA_LEN] = "create_";
    strcat(create, name);
    char destroy[9 + MAX_SOMA_LEN] = "destroy_";
    strcat(destroy, name);

    // load the soma library
    void* soma = dlopen(lib_path, RTLD_LAZY);
    if (!soma) {
        INFO("Error: Couldn't load library %s\n", lib_path);
        exit(1);
    }

    // reset errors
    dlerror();
    
    // Function to create an instance of the Soma class
    _create_soma* create_func = (_create_soma*) dlsym(soma, create);
    create_soma[std::string(name)] = create_func;
    const char* dlsym_error = dlerror();
    if (dlsym_error) {
        INFO("Error: Couldn't load symbol %s: %s\n", create, dlsym_error);
        exit(1);
    }
    
    // Function to destroy an instance of the Soma class
    destroy_soma[std::string(name)] = (_destroy_soma*) dlsym(soma, destroy);
    dlsym_error = dlerror();
    if (dlsym_error) {
        INFO("Error: Couldn't load symbol %s: %s\n", destroy, dlsym_error);
        exit(1);
    }
    INFO("Loaded plugin symbols for %s.\n", name);
}

Base_Soma* get_soma(char* name){
    std::string name_s = std::string(name);
    if(!create_soma.count(name_s)){
        init_soma(name);
    }
    return create_soma[name_s]();
}