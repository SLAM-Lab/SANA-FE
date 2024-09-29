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

// TODO: support conv, dense and sparse hyperedges here
// TODO: support set_attributes for neurons and for groups

std::map<std::string, sanafe::ModelParam> pydict_to_model_parameters(
        const pybind11::dict &dictionary, const bool forward_to_synapse = true,
        const bool forward_to_dendrite = true,
        const bool forward_to_soma = false);
sanafe::ModelParam pyobject_to_model_parameter(const pybind11::object &value);
pybind11::dict run_data_to_dict(const sanafe::RunData &run);

std::map<std::string, sanafe::ModelParam> pydict_to_model_parameters(
        const pybind11::dict &dictionary, const bool forward_to_synapse,
        const bool forward_to_dendrite, const bool forward_to_soma)
{
    // Convert a Python dict to a set of C++ attribute strings,
    //  ready to be used by the kernel
    std::map<std::string, sanafe::ModelParam> map;
    for (const auto &key_value_pair : dictionary)
    {
        TRACE1_PYBIND("Adding dict val: dict['%s']\n",
                pybind11::cast<std::string>(v.first).c_str());
        const std::string key =
                pybind11::cast<std::string>(key_value_pair.first);

        pybind11::object value =
                pybind11::cast<pybind11::object>(key_value_pair.second);
        sanafe::ModelParam parameter = pyobject_to_model_parameter(value);

        parameter.forward_to_synapse = forward_to_synapse;
        parameter.forward_to_dendrite = forward_to_dendrite;
        parameter.forward_to_soma = forward_to_soma;
        parameter.name = key;
        map[key] = parameter;
        TRACE1_PYBIND("Set map[%s]=%s\n", key.c_str(), map[key].c_str());
        INFO("Set map[%s]\n", key.c_str());
    }
    TRACE1_PYBIND("Converted map.size()=%zu\n", map.size());

    return map;
}

sanafe::ModelParam pyobject_to_model_parameter(const pybind11::object &value)
{
    // Convert a Python object (string, int or float) to a C++ string.
    //  This is required by the kernel. Although this isn't the
    //  most efficient way to share data between Python and the kernel, it
    //  does enforce a representation that will always be possible to save
    //  and load (and be reproducible)
    sanafe::ModelParam parameter;
    if (pybind11::isinstance<pybind11::str>(value))
    {
        parameter.value = pybind11::cast<std::string>(value);
        TRACE1_PYBIND("Converted str to str:%s\n", value_str.c_str());
    }
    else if (pybind11::isinstance<pybind11::int_>(value))
    {
        parameter.value = pybind11::cast<int>(value);
        TRACE1_PYBIND("Converted int to str:%s\n", value_str.c_str());
    }
    else if (pybind11::isinstance<pybind11::float_>(value))
    {
        parameter.value = pybind11::cast<float>(value);
        TRACE1_PYBIND("Converted float to str:%s\n", value_str.c_str());
    }
    else if (pybind11::isinstance<pybind11::dict>(value))
    {
        std::vector<sanafe::ModelParam> model_parameters;
        for (auto item : pybind11::cast<pybind11::dict>(value))
        {
            sanafe::ModelParam param = pyobject_to_model_parameter(
                    pybind11::cast<pybind11::object>(item.second));
            param.name = pybind11::cast<std::string>(item.first);
        }
        parameter.value = model_parameters;
    }
    else if (pybind11::isinstance<pybind11::iterable>(value))
    {
        std::vector<sanafe::ModelParam> list;
        for (auto iter = pybind11::iter(value);
                iter != pybind11::iterator::sentinel(); ++iter)
        {
            list.push_back(pyobject_to_model_parameter(
                    pybind11::cast<pybind11::object>(*iter)));
        }
        parameter.value = list;
    }
    else
    {
        throw std::invalid_argument("Error: dict has unsupported type");
    }

    return parameter;
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

sanafe::NeuronGroup pycreate_neuron_group(sanafe::Network *self,
        const std::string &group_name, const int neuron_count,
        pybind11::dict &model_dict, std::string default_synapse_hw_name,
        std::string dendrite_hw_name, bool force_dendrite_update,
        bool force_soma_update, bool force_synapse_update, bool log_potential,
        bool log_spikes, std::string soma_hw_name)
{
    sanafe::NeuronTemplate default_neuron_config;
    default_neuron_config.default_synapse_hw_name = default_synapse_hw_name;
    default_neuron_config.dendrite_hw_name = dendrite_hw_name;
    default_neuron_config.force_dendrite_update = force_dendrite_update;
    default_neuron_config.force_synapse_update = force_synapse_update;
    default_neuron_config.force_soma_update = force_soma_update;
    default_neuron_config.log_potential = log_potential;
    default_neuron_config.log_spikes = log_spikes;
    default_neuron_config.soma_hw_name = soma_hw_name;

    std::map<std::string, sanafe::ModelParam> model_parameters =
            pydict_to_model_parameters(model_dict);
    default_neuron_config.model_parameters = model_parameters;

    return self->create_neuron_group(
            group_name, neuron_count, default_neuron_config);
}

std::unique_ptr<sanafe::Neuron> pycreate_neuron(size_t neuron_id,
        sanafe::Network &parent_net, std::string parent_group_id,
        std::string soma_hw_name, std::string default_synapse_hw_name,
        std::string dendrite_hw_name, bool log_spikes, bool log_potential,
        bool force_synapse_update, bool force_dendrite_update,
        bool force_soma_update)
{
    sanafe::NeuronTemplate neuron_config;
    neuron_config.soma_hw_name = soma_hw_name;
    neuron_config.default_synapse_hw_name = default_synapse_hw_name;
    neuron_config.dendrite_hw_name = dendrite_hw_name;
    neuron_config.log_spikes = log_spikes;
    neuron_config.log_potential = log_potential;
    neuron_config.force_synapse_update = force_synapse_update;
    neuron_config.force_dendrite_update = force_dendrite_update;
    neuron_config.force_soma_update = force_soma_update;

    return std::make_unique<sanafe::Neuron>(
            neuron_id, parent_net, parent_group_id, neuron_config);
}

void pyset_attributes(sanafe::Neuron *self,
        std::optional<std::string> soma_hw_name,
        std::optional<std::string> default_synapse_hw_name,
        std::optional<std::string> dendrite_hw_name,
        std::optional<bool> log_spikes, std::optional<bool> log_potential,
        std::optional<bool> force_synapse_update,
        std::optional<bool> force_dendrite_update,
        std::optional<bool> force_soma_update, pybind11::dict model_parameters,
        pybind11::dict dendrite_specific_parameters,
        pybind11::dict soma_specific_parameters)
{
    sanafe::NeuronTemplate neuron_template{};

    neuron_template.soma_hw_name = soma_hw_name;
    neuron_template.default_synapse_hw_name = default_synapse_hw_name;
    neuron_template.dendrite_hw_name = dendrite_hw_name;
    neuron_template.log_spikes = log_spikes;
    neuron_template.log_potential = log_potential;
    neuron_template.force_synapse_update = force_synapse_update;
    neuron_template.force_dendrite_update = force_dendrite_update;
    neuron_template.force_soma_update = force_soma_update;

    auto parsed = pydict_to_model_parameters(model_parameters);
    neuron_template.model_parameters.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_parameters(
            dendrite_specific_parameters, false, true, false);
    neuron_template.model_parameters.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_parameters(
            soma_specific_parameters, false, false, true);
    neuron_template.model_parameters.insert(parsed.begin(), parsed.end());

    self->set_attributes(neuron_template);

    return;
}

PYBIND11_MODULE(sanafe, m)
{
    m.doc() = R"pbdoc(
    SANA-FE Kernel Module
    --------------------------------

    .. currentmodule:: sanafecpp

    .. autosummary::
       :toctree: _generate

           SANA_FE
    )pbdoc";

    m.def("load_arch", &sanafe::load_arch);
    m.def("load_net", &sanafe::load_net, pybind11::arg("path"),
            pybind11::arg("arch"), pybind11::arg("use_netlist_format") = false);
    pybind11::register_exception<std::runtime_error>(m, "KernelRuntimeError");
    pybind11::register_exception<std::invalid_argument>(
            m, "KernelInvalidArgument");

    pybind11::class_<sanafe::Network>(m, "Network")
            .def(pybind11::init<>())
            .def("__repr__", &sanafe::Network::info)
            .def("create_neuron_group", &pycreate_neuron_group,
                    pybind11::return_value_policy::reference_internal,
                    pybind11::arg("group_name"), pybind11::arg("neuron_count"),
                    pybind11::arg("model_parameters") = pybind11::dict(),
                    pybind11::arg("default_synapse_hw_name") = "",
                    pybind11::arg("default_dendrite_hw_name") = "",
                    pybind11::arg("force_dendrite_update") = false,
                    pybind11::arg("force_synapse_update") = false,
                    pybind11::arg("force_soma_update") = false,
                    pybind11::arg("log_potential") = false,
                    pybind11::arg("log_spikes") = false,
                    pybind11::arg("soma_hw_name") = "")
            .def_readwrite("groups", &sanafe::Network::groups);

    pybind11::class_<sanafe::NeuronGroup>(m, "NeuronGroup")
            .def("__repr__", &sanafe::NeuronGroup::info)
            .def("get_id", &sanafe::NeuronGroup::get_id)
            .def_property(
                    "neurons",
                    [](sanafe::NeuronGroup &self)
                            -> std::vector<sanafe::Neuron> & {
                        return self.neurons;
                    },
                    nullptr); // Read-only property

    pybind11::class_<sanafe::Neuron>(m, "Neuron")
            .def(pybind11::init(&pycreate_neuron), pybind11::arg("neuron_id"),
                    pybind11::arg("parent_net"),
                    pybind11::arg("parent_group_id"),
                    pybind11::arg("default_synapse_hw_name") = "",
                    pybind11::arg("dendrite_hw_name") = "",
                    pybind11::arg("force_dendrite_update") = false,
                    pybind11::arg("force_synapse_update") = false,
                    pybind11::arg("force_soma_update") = false,
                    pybind11::arg("log_potential") = false,
                    pybind11::arg("log_spikes") = false,
                    pybind11::arg("soma_hw_name") = "")
            .def("__repr__", &sanafe::Neuron::info)
            .def("set_attributes", &pyset_attributes,
                    pybind11::arg("soma_hw_name") = pybind11::none(),
                    pybind11::arg("default_synapse_hw_name") = pybind11::none(),
                    pybind11::arg("dendrite_hw_name") = pybind11::none(),
                    pybind11::arg("log_spikes") = pybind11::none(),
                    pybind11::arg("log_potential") = pybind11::none(),
                    pybind11::arg("force_synapse_update") = pybind11::none(),
                    pybind11::arg("force_dendrite_update") = pybind11::none(),
                    pybind11::arg("force_soma_update") = pybind11::none(),
                    pybind11::arg("model_parameters") = pybind11::dict(),
                    pybind11::arg("soma_parameters") = pybind11::dict(),
                    pybind11::arg("dendrite_parameters") = pybind11::dict())
            .def("connect_to_neuron",
                    [](sanafe::Neuron *self, sanafe::Neuron &dest,
                            std::optional<pybind11::dict> attr =
                                    pybind11::none()) {
                        if (!attr.has_value())
                        {
                            attr = pybind11::dict();
                        }
                        sanafe::Connection &con = self->connect_to_neuron(dest);
                        const auto parameters =
                                pydict_to_model_parameters(attr.value());
                        for (auto &[key, parameter] : parameters)
                        {
                            if (parameter.forward_to_synapse)
                            {
                                con.synapse_params[key] = parameter;
                            }
                            if (parameter.forward_to_dendrite)
                            {
                                con.dendrite_params[key] = parameter;
                            }
                        }

                        return con;
                    })
            .def("get_id", &sanafe::Neuron::get_id);

    pybind11::class_<sanafe::CoreAddress>(m, "CoreAddress")
            .def(pybind11::init<size_t, size_t, size_t>(), pybind11::arg("id"),
                    pybind11::arg("parent_tile_id"),
                    pybind11::arg("offset_within_tile"));

    pybind11::class_<sanafe::Architecture>(m, "Architecture")
            .def(pybind11::init<std::string,
                    sanafe::NetworkOnChipConfiguration>())
            .def("__repr__", &sanafe::Architecture::info)
            // TODO: support programming new architectures from python again
            //.def("create_tile", &sanafe::Architecture::create_tile,
            //        pybind11::return_value_policy::reference_internal)
            //.def("create_core", &sanafe::Architecture::create_core,
            //        pybind11::return_value_policy::reference_internal)
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
            // TODO: adding hardware units
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
