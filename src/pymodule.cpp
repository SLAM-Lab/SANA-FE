// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// pymodule.cpp
#include <any>
#include <map>
#include <sstream>
#include <unordered_map>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "arch.hpp"
#include "hardware.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"
#include "sim.hpp"

#define PYBIND11_DETAILED_ERROR_MESSAGES

// TODO: support conv, dense and sparse hyperedges here

// Forward declarations
std::map<std::string, sanafe::ModelParam> pydict_to_model_parameters(
        const pybind11::dict &dictionary, const bool forward_to_synapse = true,
        const bool forward_to_dendrite = true,
        const bool forward_to_soma = true);
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
        const std::string key =
                pybind11::cast<std::string>(key_value_pair.first);
        TRACE1(PYMODULE, "Adding dict val: dict['%s']\n", key.c_str());

        pybind11::object value =
                pybind11::cast<pybind11::object>(key_value_pair.second);
        sanafe::ModelParam parameter = pyobject_to_model_parameter(value);

        parameter.forward_to_synapse = forward_to_synapse;
        parameter.forward_to_dendrite = forward_to_dendrite;
        parameter.forward_to_soma = forward_to_soma;
        parameter.name = key;
        map[key] = parameter;
    }
    TRACE1(PYMODULE, "Converted map.size()=%zu\n", map.size());

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
        TRACE1(PYMODULE, "Converted str to parameter\n");
    }
    else if (pybind11::isinstance<pybind11::int_>(value))
    {
        parameter.value = pybind11::cast<int>(value);
        TRACE1(PYMODULE, "Converted int to parameter\n");
    }
    else if (pybind11::detail::make_caster<float>().load(value, true))
    {
        // Value can be successfully cast to a float e.g., a float or a NumPy
        //  float32 scalar. This was chosen to be more flexible than checking
        //  against the built-in float_ type
        parameter.value = pybind11::cast<float>(value);
        TRACE1(PYMODULE, "Converted float to parameter\n");
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

sanafe::NeuronGroup &pycreate_neuron_group(sanafe::SpikingNetwork *self,
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
        sanafe::SpikingNetwork &parent_net, std::string parent_group_id,
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

sanafe::TileConfiguration &pycreate_tile(sanafe::Architecture *self,
        std::string name, double energy_north_hop, double latency_north_hop,
        double energy_east_hop, double latency_east_hop,
        double energy_south_hop, double latency_south_hop,
        double energy_west_hop, double latency_west_hop)
{
    sanafe::TilePowerMetrics tile_power_metrics{};
    tile_power_metrics.energy_north_hop = energy_north_hop;
    tile_power_metrics.latency_north_hop = latency_north_hop;
    tile_power_metrics.energy_east_hop = energy_east_hop;
    tile_power_metrics.latency_east_hop = latency_east_hop;
    tile_power_metrics.energy_south_hop = energy_south_hop;
    tile_power_metrics.latency_south_hop = latency_south_hop;
    tile_power_metrics.energy_west_hop = energy_west_hop;
    tile_power_metrics.latency_west_hop = latency_west_hop;

    return self->create_tile(name, tile_power_metrics);
}

std::unique_ptr<sanafe::TileConfiguration> pyconstruct_tile(std::string name,
        size_t tile_id, double energy_north_hop, double latency_north_hop,
        double energy_east_hop, double latency_east_hop,
        double energy_south_hop, double latency_south_hop,
        double energy_west_hop, double latency_west_hop)
{
    sanafe::TilePowerMetrics tile_power_metrics{};
    tile_power_metrics.energy_north_hop = energy_north_hop;
    tile_power_metrics.latency_north_hop = latency_north_hop;
    tile_power_metrics.energy_east_hop = energy_east_hop;
    tile_power_metrics.latency_east_hop = latency_east_hop;
    tile_power_metrics.energy_south_hop = energy_south_hop;
    tile_power_metrics.latency_south_hop = latency_south_hop;
    tile_power_metrics.energy_west_hop = energy_west_hop;
    tile_power_metrics.latency_west_hop = latency_west_hop;

    return std::make_unique<sanafe::TileConfiguration>(
            name, tile_id, tile_power_metrics);
}

std::unique_ptr<sanafe::CoreConfiguration> pyconstruct_core(std::string name,
        size_t parent_tile_id, size_t offset_within_tile, size_t core_id,
        std::string buffer_position, size_t max_neurons_supported)
{
    sanafe::CorePipelineConfiguration pipeline_config{};

    pipeline_config.buffer_position =
            sanafe::pipeline_parse_buffer_pos_str(buffer_position);
    pipeline_config.max_neurons_supported = max_neurons_supported;

    sanafe::CoreAddress core_address{};
    core_address.parent_tile_id = parent_tile_id;
    core_address.offset_within_tile = offset_within_tile;
    core_address.id = core_id;

    return std::make_unique<sanafe::CoreConfiguration>(
            name, core_address, pipeline_config);
}

sanafe::CoreConfiguration &pycreate_core(sanafe::Architecture *self,
        std::string name, size_t parent_tile_id, std::string buffer_position,
        size_t max_neurons_supported)
{
    sanafe::CorePipelineConfiguration pipeline_config{};

    pipeline_config.buffer_position =
            sanafe::pipeline_parse_buffer_pos_str(buffer_position);
    pipeline_config.max_neurons_supported = max_neurons_supported;

    return self->create_core(name, parent_tile_id, pipeline_config);
}

// For now support setting attributes of both Neuron and MappedNeuron
// TODO: in the future, how I handle these may change (might simplify the Neuron
//  struct)
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

    TRACE1(PYMODULE, "Setting neuron attributes\n");
    if (DEBUG_LEVEL_PYMODULE > 0)
    {
        for (auto &[key, value] : neuron_template.model_parameters)
        {
                TRACE1(PYMODULE, "\tkey: %s\n", key.c_str());
        }
    }
    self->set_attributes(neuron_template);

    return;
}

void pyset_attributes_mapped(sanafe::MappedNeuron *self,
        pybind11::dict model_parameters,
        pybind11::dict dendrite_specific_parameters,
        pybind11::dict soma_specific_parameters)
{
    sanafe::NeuronTemplate neuron_template{};

    // Most neuron attributes cannot be set dynamically. Only forward model
    //  parameters
    // TODO: reorg how we deal with NeuronTemplates, Neuron, and MappedNeuron
    //  to hopefully simplify this and setting neuron simulator attributes and
    //  model parameters

    auto parsed = pydict_to_model_parameters(model_parameters);
    neuron_template.model_parameters.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_parameters(
            dendrite_specific_parameters, false, true, false);
    neuron_template.model_parameters.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_parameters(
            soma_specific_parameters, false, false, true);
    neuron_template.model_parameters.insert(parsed.begin(), parsed.end());

    TRACE1(PYMODULE, "Setting neuron attributes\n");
    if (DEBUG_LEVEL_PYMODULE > 0)
    {
        for (auto &[key, value] : neuron_template.model_parameters)
        {
                TRACE1(PYMODULE, "\tkey: %s\n", key.c_str());
        }
    }
    self->set_attributes(neuron_template);

    return;
}

PYBIND11_MODULE(sanafecpp, m)
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

    pybind11::enum_<sanafe::BufferPosition>(m, "BufferPosition")
            .value("BUFFER_BEFORE_DENDRITE_UNIT",
                    sanafe::BUFFER_BEFORE_DENDRITE_UNIT)
            .value("BUFFER_BEFORE_SOMA_UNIT", sanafe::BUFFER_BEFORE_SOMA_UNIT)
            .value("BUFFER_BEFORE_AXON_OUT_UNIT",
                    sanafe::BUFFER_BEFORE_AXON_OUT_UNIT)
            .value("BUFFER_POSITIONS", sanafe::BUFFER_POSITIONS);

    pybind11::class_<sanafe::SpikingNetwork>(m, "Network")
            .def(pybind11::init<>())
            .def("__repr__", &sanafe::SpikingNetwork::info)
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
            //.def_readwrite("groups", &sanafe::Network::groups);
            .def_property(
                    "groups",
                    [](sanafe::SpikingNetwork &self)
                            -> std::map<std::string, sanafe::NeuronGroup> & {
                        return self.groups;
                    },
                    nullptr, pybind11::return_value_policy::reference_internal)
            .def(
                    "__getitem__",
                    [](sanafe::SpikingNetwork &self,
                            std::string g) -> sanafe::NeuronGroup & {
                        if (self.groups.find(g) == self.groups.end())
                        {
                            throw pybind11::index_error();
                        }
                        return self.groups.at(g);
                    },
                    pybind11::return_value_policy::reference_internal);

    pybind11::class_<sanafe::NeuronGroup>(m, "NeuronGroup")
            .def("__repr__", &sanafe::NeuronGroup::info)
            .def("get_id", &sanafe::NeuronGroup::get_id)
            //.def_readwrite("neurons", &sanafe::NeuronGroup::neurons);
            .def_property("neurons",
                    pybind11::cpp_function(
                            [](sanafe::NeuronGroup &self)
                                    -> std::vector<sanafe::Neuron> & {
                                return self.neurons;
                            },
                            pybind11::return_value_policy::reference_internal),
                    nullptr) // Read-only property
            .def(
                    "__getitem__",
                    [](sanafe::NeuronGroup &self,
                            size_t i) -> sanafe::Neuron & {
                        if (i >= self.neurons.size())
                        {
                            throw pybind11::index_error();
                        }
                        return self.neurons[i];
                    },
                    pybind11::return_value_policy::reference_internal)
            .def(
                    "__iter__",
                    [](sanafe::NeuronGroup &self) {
                        return pybind11::make_iterator(self.neurons.begin(),
                                self.neurons.end(),
                                pybind11::return_value_policy::
                                        reference_internal);
                    },
                    pybind11::keep_alive<0, 1>());
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
            .def(
                    "connect_to_neuron",
                    [](sanafe::Neuron *self, sanafe::Neuron &dest,
                            std::optional<pybind11::dict> attr =
                                    pybind11::none()) {
                        if (!attr.has_value())
                        {
                            attr = pybind11::dict();
                        }
                        const size_t con_idx = self->connect_to_neuron(dest);
                        sanafe::Connection &con =
                                self->edges_out[con_idx];
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

                        return con_idx;
                    },
                    pybind11::return_value_policy::reference_internal)
            .def("map_to_core", &sanafe::Neuron::map_to_core)
            .def("get_id", &sanafe::Neuron::get_id,
                    pybind11::return_value_policy::reference_internal);

    pybind11::class_<sanafe::Architecture>(m, "Architecture")
            .def(pybind11::init<std::string,
                    sanafe::NetworkOnChipConfiguration>())
            .def("__repr__", &sanafe::Architecture::info)
            .def("create_tile", &pycreate_tile,
                    pybind11::return_value_policy::reference_internal,
                    pybind11::arg("name"),
                    pybind11::arg("energy_north_hop") = 0.0,
                    pybind11::arg("latency_north_hop") = 0.0,
                    pybind11::arg("energy_east_hop") = 0.0,
                    pybind11::arg("latency_east_hop") = 0.0,
                    pybind11::arg("energy_south_hop") = 0.0,
                    pybind11::arg("latency_south_hop") = 0.0,
                    pybind11::arg("energy_west_hop") = 0.0,
                    pybind11::arg("latency_west_hop") = 0.0)
            .def("create_core", &pycreate_core, pybind11::arg("name"),
                    pybind11::arg("parent_tile_id"),
                    pybind11::arg("buffer_position") =
                            sanafe::BUFFER_BEFORE_SOMA_UNIT,
                    pybind11::arg("max_neurons_supported") =
                            sanafe::default_max_neurons,
                    pybind11::return_value_policy::reference_internal)
            .def_readwrite("tiles", &sanafe::Architecture::tiles);

    pybind11::class_<sanafe::TileConfiguration>(m, "Tile")
            .def(pybind11::init(&pyconstruct_tile), pybind11::arg("name"),
                    pybind11::arg("id"),
                    pybind11::arg("energy_north_hop") = 0.0,
                    pybind11::arg("latency_north_hop") = 0.0,
                    pybind11::arg("energy_east_hop") = 0.0,
                    pybind11::arg("latency_east_hop") = 0.0,
                    pybind11::arg("energy_south_hop") = 0.0,
                    pybind11::arg("latency_south_hop") = 0.0,
                    pybind11::arg("energy_west_hop") = 0.0,
                    pybind11::arg("latency_west_hop") = 0.0)
            .def_readwrite("cores", &sanafe::TileConfiguration::cores);

    pybind11::class_<sanafe::CoreConfiguration,
            std::shared_ptr<sanafe::CoreConfiguration>>(m, "Core")
            .def(pybind11::init(&pyconstruct_core), pybind11::arg("name"),
                    pybind11::arg("parent_tile_id"),
                    pybind11::arg("offset_within_tile"),
                    pybind11::arg("core_id"),
                    pybind11::arg("buffer_position") =
                            sanafe::BUFFER_BEFORE_SOMA_UNIT,
                    pybind11::arg("max_neurons_supported") =
                            sanafe::default_max_neurons);
            //.def("__repr__", &sanafe::CoreConfiguration::info)
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
    pybind11::class_<sanafe::SpikingHardware>(m, "SpikingHardware")
            .def_property(
                    "mapped_neuron_groups",
                    [](sanafe::SpikingHardware &self)
                            -> std::map<std::string,
                                    std::vector<sanafe::MappedNeuron *>> & {
                        return self.mapped_neuron_groups;
                    },
                    nullptr, pybind11::return_value_policy::reference_internal)
            .def(pybind11::init<sanafe::Architecture &, std::string, bool, bool,
                         bool, bool>(),
                    pybind11::arg("arch"), pybind11::arg("output_dir") = ".",
                    pybind11::arg("record_spikes") = false,
                    pybind11::arg("record_potentials") = false,
                    pybind11::arg("record_perf") = false,
                    pybind11::arg("record_messages") = false)
            .def("load", &sanafe::SpikingHardware::load)
            .def(
                    "sim",
                    [](sanafe::SpikingHardware *self, const long int timesteps,
                            const long int heartbeat) {
                        return run_data_to_dict(
                                self->sim(timesteps, heartbeat));
                    },
                    pybind11::arg("timesteps") = 1,
                    pybind11::arg("heartbeat") = 100)
            .def("get_power", &sanafe::SpikingHardware::get_power)
            .def("get_run_summary", [](sanafe::SpikingHardware *self) {
                return run_data_to_dict(self->get_run_summary());
            })
            .def("reset", &sanafe::SpikingHardware::reset);
    pybind11::class_<sanafe::MappedNeuron>(m, "MappedNeuron")
            .def("set_attributes", &pyset_attributes_mapped,
                    pybind11::arg("model_parameters") = pybind11::dict(),
                    pybind11::arg("soma_parameters") = pybind11::dict(),
                    pybind11::arg("dendrite_parameters") = pybind11::dict());
}
