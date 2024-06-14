// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// pymodule.cpp
#include <any>
#include <map>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "arch.hpp"
#include "network.hpp"
#include "print.hpp"
#include "sim.hpp"

#define PYBIND11_DETAILED_ERROR_MESSAGES

std::string pyobject_to_str(const pybind11::object &value);
std::vector<std::string> pylist_to_str_vec(const pybind11::iterable &value);
std::map<std::string, std::string> dict_to_str_map(const pybind11::dict &vals);
pybind11::dict run_data_to_dict(const sanafe::RunData &run);
void py_connect_neurons(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group,
        const std::vector<std::pair<int, int>> &src_dest_id_pairs,
        const pybind11::dict &vals);

std::vector<std::string> pylist_to_str_vec(const pybind11::iterable &value)
{
    // Convert any iterable python object e.g., a list or tuple, to a
    //  C++ vector of strings. This can then be used by the kernel.
    std::vector<std::string> value_vec_str;

    for (auto iter = pybind11::iter(value);
            iter != pybind11::iterator::sentinel(); ++iter)
    {
        value_vec_str.push_back(
                pyobject_to_str(pybind11::cast<pybind11::object>(*iter)));
    }

    return value_vec_str;
}

std::map<std::string, std::vector<std::string>> pylist_or_object_to_str_vec(
        const size_t count, const pybind11::dict &vals)
{
    std::map<std::string, std::vector<std::string>> attr_lists;
    for (auto &v : vals)
    {
        const std::string key = pybind11::cast<std::string>(v.first);
        if (!pybind11::isinstance<pybind11::iterable>(v.second))
        {
            // Only a single value was defined for this attribute,
            //  this should be duplicated for all connections.
            //  First get the number of new connections to create
            std::vector<std::string> vec(count,
                    pyobject_to_str(
                            pybind11::cast<pybind11::object>(v.second)));
            attr_lists[key] = vec;
        }
        else
        {
            attr_lists[key] = pylist_to_str_vec(
                    pybind11::cast<pybind11::iterable>(v.second));
        }
    }

    return attr_lists;
}

std::map<std::string, std::string> dict_to_str_map(const pybind11::dict &vals)
{
    // Convert a Python dict to a set of C++ attribute strings,
    //  ready to be used by the kernel
    std::map<std::string, std::string> map;
    for (const auto &v : vals)
    {
        TRACE1_PYBIND("Adding dict val: dict['%s']\n",
                pybind11::cast<std::string>(v.first).c_str());
        const std::string key = pybind11::cast<std::string>(v.first);
        const std::string value_str =
                pyobject_to_str(pybind11::cast<pybind11::object>(v.second));
        map[key] = value_str;
        TRACE1_PYBIND("Set map[%s]=%s\n", key.c_str(), map[key].c_str());
    }
    TRACE1_PYBIND("Converted map.size()=%lu\n", map.size());

    return map;
}

std::string pyobject_to_str(const pybind11::object &value)
{
    // Convert a Python object (string, int or float) to a C++ string.
    //  This is required by the kernel. Although this isn't the
    //  most efficient way to share data between Python and the kernel, it
    //  does enforce a representation that will always be possible to save
    //  and load (and be reproducible)
    std::string value_str;
    if (pybind11::isinstance<pybind11::str>(value))
    {
        value_str = pybind11::cast<std::string>(value);
        TRACE1_PYBIND("Converted str to str:%s\n", value_str.c_str());
    }
    else if (pybind11::isinstance<pybind11::int_>(value))
    {
        value_str = std::to_string(pybind11::cast<int>(value));
        TRACE1_PYBIND("Converted int to str:%s\n", value_str.c_str());
    }
    else if (pybind11::isinstance<pybind11::float_>(value))
    {
        std::ostringstream os;
        os << pybind11::cast<float>(value);
        value_str = os.str();
        TRACE1_PYBIND("Converted float to str:%s\n", value_str.c_str());
    }
    else
    {
        throw std::invalid_argument("Error: dict has unsupported type");
    }

    return value_str;
}

pybind11::dict run_data_to_dict(const sanafe::RunData &run)
{
    pybind11::dict run_data_dict;
    run_data_dict["timestep_start"] = run.timestep_start;
    run_data_dict["timesteps_executed"] = run.timesteps_executed;
    run_data_dict["energy"] = run.energy;
    run_data_dict["sim_time"] = run.sim_time;
    run_data_dict["spikes"] = run.spikes;
    run_data_dict["packets_sent"] = run.packets_sent;
    run_data_dict["neurons_fired"] = run.neurons_fired;

    return run_data_dict;
}

PYBIND11_MODULE(sanafecpp, m)
{
    m.doc() = R"pbdoc(
    SANA-FE CPP Kernel Module
    --------------------------------

    .. currentmodule:: sanafecpp

    .. autosummary::
       :toctree: _generate

           SANA_FE
    )pbdoc";

    m.def("load_arch", &sanafe::load_arch);
    m.def("load_net", &sanafe::load_net);
    pybind11::register_exception<std::runtime_error>(m, "KernelRuntimeError");
    pybind11::register_exception<std::invalid_argument>(
            m, "KernelInvalidArgument");

    pybind11::class_<sanafe::Network>(m, "Network")
            .def(pybind11::init<>())
            .def("__repr__", &sanafe::Network::info)
            .def("save_net_description", &sanafe::Network::save_net_description,
                    pybind11::arg("path"),
                    pybind11::arg("write_mappings") = true)
            .def(
                    "create_neuron_group",
                    [](sanafe::Network *self, const int neuron_count,
                            pybind11::dict &attr) -> sanafe::NeuronGroup & {
                        return self->create_neuron_group(
                                neuron_count, dict_to_str_map(attr));
                    },
                    pybind11::return_value_policy::reference_internal)
            .def_readwrite("groups", &sanafe::Network::groups);

    pybind11::class_<sanafe::NeuronGroup>(m, "NeuronGroup")
            .def(pybind11::init<int, int>(), pybind11::arg("id"),
                    pybind11::arg("neuron_count"))
            .def("__repr__", &sanafe::NeuronGroup::info)
            .def("set_attribute_multiple",
                    [](sanafe::NeuronGroup *self, const std::string &attr_name,
                            pybind11::iterable &values) {
                        return self->set_attribute_multiple(
                                attr_name, pylist_to_str_vec(values));
                    })
            .def("connect_neurons",
                    [](sanafe::NeuronGroup *self,
                            sanafe::NeuronGroup &dest_group,
                            const std::vector<std::pair<int, int>>
                                    &src_dest_id_pairs,
                            pybind11::dict &attr) {
                        return self->connect_neurons(dest_group,
                                src_dest_id_pairs,
                                pylist_or_object_to_str_vec(
                                        src_dest_id_pairs.size(), attr));
                    })
            .def("get_id", &sanafe::NeuronGroup::get_id)
            .def_readwrite("neurons", &sanafe::NeuronGroup::neurons);

    pybind11::class_<sanafe::Neuron>(m, "Neuron")
            .def(pybind11::init<int>(), pybind11::arg("id"))
            .def("__repr__", &sanafe::Neuron::info)
            .def("set_attributes", &sanafe::Neuron::set_attributes)
            .def("connect_to_neuron",
                    [](sanafe::Neuron *self, sanafe::Neuron &dest,
                            std::optional<pybind11::dict> attr =
                                    pybind11::none()) {
                        if (!attr.has_value())
                        {
                            attr = pybind11::dict();
                        }
                        return self->connect_to_neuron(dest,
                                dict_to_str_map(attr.value()));
                    })
            .def("get_id", &sanafe::Neuron::get_id);

    pybind11::class_<sanafe::NetworkOnChipConfiguration>(
            m, "NetworkOnChipConfiguration")
            .def(pybind11::init<int, int, int>(),
                    pybind11::arg("width_in_tiles") = 1,
                    pybind11::arg("height_in_tiles") = 1,
                    pybind11::arg("link_buffer_size") = 0);

    pybind11::class_<sanafe::TilePowerMetrics>(m, "TilePowerMetrics")
            .def(pybind11::init<double, double, double, double, double, double,
                         double, double>(),
                    pybind11::arg("energy_north") = 0.0,
                    pybind11::arg("latency_north") = 0.0,
                    pybind11::arg("energy_east") = 0.0,
                    pybind11::arg("latency_east") = 0.0,
                    pybind11::arg("energy_south") = 0.0,
                    pybind11::arg("latency_south") = 0.0,
                    pybind11::arg("energy_west") = 0.0,
                    pybind11::arg("latency_west") = 0.0);

    pybind11::class_<sanafe::CorePipelineConfiguration>(
            m, "CorePipelineConfiguration")
            .def(pybind11::init<const std::string, const size_t>(),
                    pybind11::arg("buffer_pos") = "soma",
                    pybind11::arg("max_neurons_supported") = 1024);

    pybind11::class_<sanafe::CoreAddress>(m, "CoreAddress")
            .def(pybind11::init<size_t, size_t, size_t>(),
                    pybind11::arg("id"),
                    pybind11::arg("parent_tile_id"),
                    pybind11::arg("offset_within_tile"));

    pybind11::class_<sanafe::Architecture>(m, "Architecture")
            .def(pybind11::init<std::string,
                    sanafe::NetworkOnChipConfiguration>())
            .def("__repr__", &sanafe::Architecture::info)
            .def("save_arch_description",
                    &sanafe::Architecture::save_arch_description)
            .def("create_tile", &sanafe::Architecture::create_tile,
                    pybind11::return_value_policy::reference_internal)
            .def("create_core", &sanafe::Architecture::create_core,
                    pybind11::return_value_policy::reference_internal)
            .def_readwrite("tiles", &sanafe::Architecture::tiles);

    pybind11::class_<sanafe::Tile>(m, "Tile")
            .def(pybind11::init<const std::string, int,
                    sanafe::TilePowerMetrics>())
            .def("get_id", &sanafe::Tile::get_id)
            .def_readwrite("cores", &sanafe::Tile::cores);

    pybind11::class_<sanafe::Core, std::shared_ptr<sanafe::Core>>(m, "Core")
            .def(pybind11::init<const std::string, const sanafe::CoreAddress &,
                         const sanafe::CorePipelineConfiguration>(),
                    pybind11::arg("name"), pybind11::arg("address"),
                    pybind11::arg("pipeline"))
            .def("__repr__", &sanafe::Core::info)
            // TODO: re-enable this at some point
            /*
            .def(
                    "create_axon_in",
                    [](sanafe::Core *self, const std::string &name,
                            pybind11::dict &attr) {
                        return self->create_axon_in(
                                name, dict_to_str_map(attr));
                    },
                    pybind11::return_value_policy::reference_internal)
            .def(
                    "create_dendrite",
                    [](sanafe::Core *self, const std::string &name,
                            pybind11::dict &attr) {
                        return self->create_dendrite(
                                name, dict_to_str_map(attr));
                    },
                    pybind11::return_value_policy::reference_internal)
            .def(
                    "create_synapse",
                    [](sanafe::Core *self, const std::string &name,
                            pybind11::dict &attr) {
                        return self->create_synapse(
                                name, dict_to_str_map(attr));
                    },
                    pybind11::return_value_policy::reference_internal)
            .def(
                    "create_soma",
                    [](sanafe::Core *self, const std::string &name,
                            pybind11::dict &attr) {
                        return self->create_soma(name, dict_to_str_map(attr));
                    },
                    pybind11::return_value_policy::reference_internal)
            .def(
                    "create_axon_out",
                    [](sanafe::Core *self, const std::string &name,
                            pybind11::dict &attr) {
                        return self->create_axon_out(
                                name, dict_to_str_map(attr));
                    },
                    pybind11::return_value_policy::reference_internal)
            */
            .def("map_neuron", &sanafe::Core::map_neuron)
            .def("get_id", &sanafe::Core::get_id)
            .def("get_offset", &sanafe::Core::get_offset);
    /*
            .def_readwrite("axon_in_hw", &sanafe::Core::axon_in_hw)
            .def_readwrite("synapse", &sanafe::Core::synapse)
            .def_readwrite("dendrite", &sanafe::Core::dendrite)
            .def_readwrite("soma", &sanafe::Core::soma)
            .def_readwrite("axon_out_hw", &sanafe::Core::axon_out_hw);
            */

    //pybind11::class_<sanafe::AxonInUnit>(m, "AxonInUnit")
    //        .def(pybind11::init<std::string>());

    //pybind11::class_<sanafe::DendriteUnit>(m, "DendriteUnit")
    //        .def(pybind11::init<std::string>());

    //pybind11::class_<sanafe::SynapseUnit>(m, "SynapseUnit")
    //        .def(pybind11::init<std::string>());

    //pybind11::class_<sanafe::SomaUnit>(m, "SomaUnit")
    //        .def(pybind11::init<std::string>());

    //pybind11::class_<sanafe::AxonOutUnit>(m, "AxonOutUnit")
    //        .def(pybind11::init<std::string>());

    pybind11::class_<sanafe::Simulation>(m, "Simulation")
            .def(pybind11::init<sanafe::Architecture &, sanafe::Network &,
                         std::string, bool, bool, bool, bool>(),
                    pybind11::arg("arch"), pybind11::arg("net"),
                    pybind11::arg("output_dir") = ".",
                    pybind11::arg("record_spikes") = false,
                    pybind11::arg("record_potentials") = false,
                    pybind11::arg("record_perf") = false,
                    pybind11::arg("record_messages") = false)
            .def(
                    "run",
                    [](sanafe::Simulation *self, const long int timesteps,
                            const long int heartbeat) {
                        return run_data_to_dict(
                                self->run(timesteps, heartbeat));
                    },
                    pybind11::arg("timesteps") = 1,
                    pybind11::arg("heartbeat") = 100)
            .def("get_power", &sanafe::Simulation::get_power)
            .def("get_run_summary", [](sanafe::Simulation *self) {
                return run_data_to_dict(self->get_run_summary());
            });
}
