// Copyright (c) 2025 - The University of Texas at Austin
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

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "arch.hpp"
#include "network.hpp"
#include "print.hpp"
#include "chip.hpp"

#define PYBIND11_DETAILED_ERROR_MESSAGES

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

std::map<std::string, std::vector<sanafe::ModelParam>>
pydict_to_attribute_lists(pybind11::dict attributes)
{
    std::map<std::string, std::vector<sanafe::ModelParam>> map{};

    for (const auto &key_value_pair : attributes)
    {
        const std::string attribute_name =
                pybind11::cast<std::string>(key_value_pair.first);
        if (!pybind11::isinstance<pybind11::iterable>(key_value_pair.second))
        {
             std::string error("Error: Each attribute must be provided as a "
                    "1D list/array of values. Multi-dimensional arrays must be "
                    "flattened in C-order and storing channels as the last dim");
            INFO("%s\n", error.c_str());
            throw std::invalid_argument(error);
        }

        map[attribute_name] = std::vector<sanafe::ModelParam>();
        for (auto iter = pybind11::iter(key_value_pair.second);
            iter != pybind11::iterator::sentinel(); ++iter)
        {
            map[attribute_name].push_back(pyobject_to_model_parameter(
                    pybind11::cast<pybind11::object>(*iter)));
        }
    }
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
    bool is_int{false};
    bool is_float{false};
    bool is_string = pybind11::isinstance<pybind11::str>(value);

    bool is_numpy = pybind11::hasattr(value, "dtype");
    if (is_numpy)
    {
        auto dtype = value.attr("dtype");
        std::string code = dtype.attr("kind").cast<std::string>();
        is_int = (code == "i" || code == "u");
        is_float = (code == "f");
    }
    else
    {
        is_int = pybind11::isinstance<pybind11::int_>(value);
        is_float = pybind11::isinstance<pybind11::float_>(value);
        //is_float = pybind11::detail::make_caster<float>().load(value, true)
    }

    // Now check against each type and make the appropriate cast
    if (is_string)
    {
        parameter.value = pybind11::cast<std::string>(value);
        TRACE1(PYMODULE, "Converted str to parameter\n");
    }
    else if (is_int)
    {
        parameter.value = pybind11::cast<int>(value);
        TRACE1(PYMODULE, "Converted int to parameter\n");
    }
    else if (is_float)
    {
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
void pyconnect_neurons_sparse(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group,
        const pybind11::dict attributes, const pybind11::list &src_dest_id_pairs)
{
    const std::map<std::string, std::vector<sanafe::ModelParam>>
            attribute_lists = pydict_to_attribute_lists(attributes);

    // Convert python list of tuples to SANA-FE format
    if (!pybind11::isinstance<pybind11::iterable>(src_dest_id_pairs))
    {
        std::string error("Error: must provide connectivity as a list of "
            "source/destination pairs, providing the offsets within the "
            "groups.");
        INFO("%s\n", error.c_str());
        throw std::invalid_argument(error);
    }

    std::vector<std::pair<size_t, size_t>> src_dest_id_pairs_vec{};
    for (auto iter = pybind11::iter(src_dest_id_pairs);
        iter != pybind11::iterator::sentinel(); ++iter)
    {
        if (!pybind11::isinstance<pybind11::iterable>(*iter))
        {
            std::string error("Error: each entry in the src/dest id list must be a 2-tuple, "
                "representing the neuron's offsets within the src and dest neuron groups");
            INFO("%s\n", error.c_str());
            throw pybind11::value_error(error);
        }

        pybind11::sequence seq = pybind11::cast<pybind11::sequence>(*iter);
        if (seq.size() != 2)
        {
            throw pybind11::value_error("Expected a 2-tuple");
        }
        const size_t src = pybind11::cast<size_t>(seq[0]);
        const size_t dest = pybind11::cast<size_t>(seq[1]);
        src_dest_id_pairs_vec.push_back(std::make_pair(src, dest));
    }

    self->connect_neurons_sparse(dest_group, attribute_lists, src_dest_id_pairs_vec);
    return;
}

void pyconnect_neurons_conv2d(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group,
        const pybind11::dict attributes,
        int input_width, int input_height, int input_channels, int kernel_width,
        int kernel_height, int kernel_count, int stride_width,
        int stride_height)
{
    const std::map<std::string, std::vector<sanafe::ModelParam>>
            attribute_lists = pydict_to_attribute_lists(attributes);
    sanafe::Conv2DParameters config{};
    config.input_width = input_width;
    config.input_height = input_height;
    config.input_channels = input_channels;

    config.kernel_width = kernel_width;
    config.kernel_height = kernel_height;
    config.kernel_count = kernel_count;

    self->connect_neurons_conv2d(dest_group, attribute_lists, config);
    return;
}

void pyconnect_neurons_dense(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group,
        const pybind11::dict attributes)
{
    const std::map<std::string, std::vector<sanafe::ModelParam>>
            attribute_lists = pydict_to_attribute_lists(attributes);

    self->connect_neurons_dense(dest_group, attribute_lists);
    return;
}

pybind11::dict run_data_to_dict(const sanafe::RunData &run)
{
    pybind11::dict energy_dict;
    energy_dict["total"] = run.total_energy;
    energy_dict["synapse"] = run.synapse_energy;
    energy_dict["dendrite"] = run.dendrite_energy;
    energy_dict["soma"] = run.soma_energy;
    energy_dict["network"] = run.network_energy;

    pybind11::dict run_data_dict;
    run_data_dict["timestep_start"] = run.timestep_start;
    run_data_dict["timesteps_executed"] = run.timesteps_executed;
    run_data_dict["energy"] = energy_dict;
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
    sanafe::NeuronConfiguration default_neuron_config;
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
    sanafe::NeuronConfiguration neuron_config;
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
        double energy_west_hop, double latency_west_hop, bool log_energy,
        bool log_latency)
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
    tile_power_metrics.log_energy = log_energy;
    tile_power_metrics.log_latency = log_latency;

    return self->create_tile(name, tile_power_metrics);
}

std::unique_ptr<sanafe::TileConfiguration> pyconstruct_tile(std::string name,
        size_t tile_id, double energy_north_hop, double latency_north_hop,
        double energy_east_hop, double latency_east_hop,
        double energy_south_hop, double latency_south_hop,
        double energy_west_hop, double latency_west_hop, bool log_energy,
        bool log_latency)
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
    tile_power_metrics.log_energy = log_energy;
    tile_power_metrics.log_latency = log_latency;

    return std::make_unique<sanafe::TileConfiguration>(
            name, tile_id, tile_power_metrics);
}

std::unique_ptr<sanafe::CoreConfiguration> pyconstruct_core(std::string name,
        size_t parent_tile_id, size_t offset_within_tile, size_t core_id,
        std::string buffer_position, bool buffer_inside_unit,
        size_t max_neurons_supported, bool log_energy, bool log_latency)
{
    sanafe::CorePipelineConfiguration pipeline_config{};

    pipeline_config.buffer_position =
            sanafe::pipeline_parse_buffer_pos_str(buffer_position, buffer_inside_unit);
    pipeline_config.max_neurons_supported = max_neurons_supported;
    pipeline_config.log_energy = log_energy;
    pipeline_config.log_latency = log_latency;

    sanafe::CoreAddress core_address{};
    core_address.parent_tile_id = parent_tile_id;
    core_address.offset_within_tile = offset_within_tile;
    core_address.id = core_id;

    return std::make_unique<sanafe::CoreConfiguration>(
            name, core_address, pipeline_config);
}

sanafe::CoreConfiguration &pycreate_core(sanafe::Architecture *self,
        std::string name, size_t parent_tile_id, std::string buffer_position,
        bool buffer_inside_unit, size_t max_neurons_supported, bool log_energy,
        bool log_latency)
{
    sanafe::CorePipelineConfiguration pipeline_config{};

    pipeline_config.buffer_position = sanafe::pipeline_parse_buffer_pos_str(
            buffer_position, buffer_inside_unit);
    pipeline_config.max_neurons_supported = max_neurons_supported;
    pipeline_config.log_energy = log_energy;
    pipeline_config.log_latency = log_latency;

    return self->create_core(name, parent_tile_id, pipeline_config);
}

void pyconfigure(sanafe::Neuron *self,
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
    sanafe::NeuronConfiguration neuron_template{};

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
    self->configure(neuron_template);

    return;
}

void pyset_model_attributes(sanafe::MappedNeuron *self,
        pybind11::dict model_parameters,
        pybind11::dict dendrite_specific_parameters,
        pybind11::dict soma_specific_parameters)
{
    std::map<std::string, sanafe::ModelParam> converted_parameters;

    auto parsed = pydict_to_model_parameters(model_parameters);
    converted_parameters.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_parameters(
            dendrite_specific_parameters, false, true, false);
    converted_parameters.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_parameters(
            soma_specific_parameters, false, false, true);
    converted_parameters.insert(parsed.begin(), parsed.end());

    TRACE1(PYMODULE, "Setting neuron attributes\n");
    if (DEBUG_LEVEL_PYMODULE > 0)
    {
        for (auto &[key, value] : converted_parameters)
        {
                TRACE1(PYMODULE, "\tkey: %s\n", key.c_str());
        }
    }
    self->set_model_attributes(converted_parameters);

    return;
}

pybind11::dict pysim(sanafe::SpikingChip *self, const long int timesteps,
        const long int heartbeat, std::string timing_model_str, int nthreads)
{
    sanafe::TimingModel timing_model;
    if (timing_model_str == "simple")
    {
        timing_model = sanafe::TIMING_MODEL_SIMPLE;
    }
    else if (timing_model_str == "detailed")
    {
        timing_model = sanafe::TIMING_MODEL_DETAILED;
    }
    else if (timing_model_str == "cycle")
    {
        timing_model = sanafe::TIMING_MODEL_CYCLE_ACCURATE;
    }
    else
    {
        throw std::invalid_argument("Error: unsupoorted timing model");
    }

#ifdef HAVE_OPENMP
    // If nthreads <= 0, just leave it at the system default
    if (nthreads > 0)
    {
        omp_set_num_threads(nthreads);
    }
#endif

    return run_data_to_dict(self->sim(timesteps, heartbeat, timing_model));
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
            .def("save", &sanafe::SpikingNetwork::save, pybind11::arg("path"),
                    pybind11::arg("use_netlist_format") = false)
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
            .def("get_name", &sanafe::NeuronGroup::get_name)
            .def("connect_neurons_dense", &pyconnect_neurons_dense)
            .def("connect_neurons_sparse", &pyconnect_neurons_sparse)
            .def("connect_neurons_conv2d", &pyconnect_neurons_conv2d,
                pybind11::arg("dest_group"),
                pybind11::arg("attributes"),
                pybind11::arg("input_width"),
                pybind11::arg("input_height"),
                pybind11::arg("input_channels"),
                pybind11::arg("kernel_width"),
                pybind11::arg("kernel_height"),
                pybind11::arg("kernel_count") = 1,
                pybind11::arg("stride_width") = 1,
                pybind11::arg("stride_height") = 1)
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
            .def("configure", &pyconfigure,
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
                    pybind11::arg("latency_west_hop") = 0.0,
                    pybind11::arg("log_energy") = false,
                    pybind11::arg("log_latency") = false)
            .def("create_core", &pycreate_core, pybind11::arg("name"),
                    pybind11::arg("parent_tile_id"),
                    pybind11::arg("buffer_position") =
                            sanafe::BUFFER_BEFORE_SOMA_UNIT,
                    pybind11::arg("buffer_inside_unit") = false,
                    pybind11::arg("max_neurons_supported") =
                            sanafe::default_max_neurons,
                    pybind11::arg("log_energy") = false,
                    pybind11::arg("log_latency") = false,
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
                    pybind11::arg("latency_west_hop") = 0.0,
                    pybind11::arg("log_energy") = false,
                    pybind11::arg("log_latency") = false)
            .def_readwrite("cores", &sanafe::TileConfiguration::cores);

    pybind11::class_<sanafe::CoreConfiguration,
            std::shared_ptr<sanafe::CoreConfiguration>>(m, "Core")
            .def(pybind11::init(&pyconstruct_core), pybind11::arg("name"),
                    pybind11::arg("parent_tile_id"),
                    pybind11::arg("offset_within_tile"),
                    pybind11::arg("core_id"),
                    pybind11::arg("buffer_position") =
                            sanafe::BUFFER_BEFORE_SOMA_UNIT,
                    pybind11::arg("buffer_inside_unit") = false,
                    pybind11::arg("max_neurons_supported") =
                            sanafe::default_max_neurons,
                    pybind11::arg("log_energy") = false,
                    pybind11::arg("log_latency") = false);
    pybind11::class_<sanafe::SpikingChip>(m, "SpikingChip")
            .def_property(
                    "mapped_neuron_groups",
                    [](sanafe::SpikingChip &self)
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
            .def("load", &sanafe::SpikingChip::load)
            .def("sim", &pysim,
                    pybind11::arg("timesteps") = 1,
                    pybind11::arg("heartbeat") = 100,
                    pybind11::arg("timing_model") = "detailed",
                    pybind11::arg("nthreads") = 0)
            .def("get_power", &sanafe::SpikingChip::get_power)
            .def("get_run_summary", [](sanafe::SpikingChip *self) {
                return run_data_to_dict(self->get_run_summary());
            })
            .def("reset", &sanafe::SpikingChip::reset);
    pybind11::class_<sanafe::MappedNeuron>(m, "MappedNeuron")
            .def("set_model_attributes", &pyset_model_attributes,
                    pybind11::arg("model_parameters") = pybind11::dict(),
                    pybind11::arg("soma_parameters") = pybind11::dict(),
                    pybind11::arg("dendrite_parameters") = pybind11::dict());
}
