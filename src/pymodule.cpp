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
#include "chip.hpp"
#include "network.hpp"
#include "print.hpp"

#define PYBIND11_DETAILED_ERROR_MESSAGES

// Forward declarations
std::map<std::string, sanafe::ModelAttribute> pydict_to_model_attributes(
        const pybind11::dict &dictionary, const bool forward_to_synapse = true,
        const bool forward_to_dendrite = true,
        const bool forward_to_soma = true);
sanafe::ModelAttribute pyobject_to_model_attribute(
        const pybind11::object &value);
pybind11::dict run_data_to_dict(const sanafe::RunData &run);

std::map<std::string, sanafe::ModelAttribute> pydict_to_model_attributes(
        const pybind11::dict &dictionary, const bool forward_to_synapse,
        const bool forward_to_dendrite, const bool forward_to_soma)
{
    // Convert a Python dict to a set of C++ attribute strings,
    //  ready to be used by the kernel
    std::map<std::string, sanafe::ModelAttribute> map;
    for (const auto &key_value_pair : dictionary)
    {
        const std::string key =
                pybind11::cast<std::string>(key_value_pair.first);
        TRACE1(PYMODULE, "Adding dict val: dict['%s']\n", key.c_str());

        pybind11::object value =
                pybind11::cast<pybind11::object>(key_value_pair.second);
        sanafe::ModelAttribute attribute = pyobject_to_model_attribute(value);

        attribute.forward_to_synapse = forward_to_synapse;
        attribute.forward_to_dendrite = forward_to_dendrite;
        attribute.forward_to_soma = forward_to_soma;
        attribute.name = key;
        map[key] = attribute;
    }
    TRACE1(PYMODULE, "Converted map.size()=%zu\n", map.size());

    return map;
}

std::map<std::string, std::vector<sanafe::ModelAttribute>>
pydict_to_attribute_lists(pybind11::dict attributes)
{
    std::map<std::string, std::vector<sanafe::ModelAttribute>> map{};

    for (const auto &key_value_pair : attributes)
    {
        const std::string attribute_name =
                pybind11::cast<std::string>(key_value_pair.first);
        if (!pybind11::isinstance<pybind11::iterable>(key_value_pair.second))
        {
            std::string error(
                    "Error: Each attribute must be provided as a "
                    "1D list/array of values. Multi-dimensional arrays must be "
                    "flattened in C-order and storing channels as the last dim");
            INFO("%s\n", error.c_str());
            throw std::invalid_argument(error);
        }

        map[attribute_name] = std::vector<sanafe::ModelAttribute>();
        for (auto iter = pybind11::iter(key_value_pair.second);
                iter != pybind11::iterator::sentinel(); ++iter)
        {
            map[attribute_name].push_back(pyobject_to_model_attribute(
                    pybind11::cast<pybind11::object>(*iter)));
        }
    }
    return map;
}

sanafe::ModelAttribute pyobject_to_model_attribute(
        const pybind11::object &value)
{
    // Convert a Python object (string, int or float) to a C++ string.
    //  This is required by the kernel. Although this isn't the
    //  most efficient way to share data between Python and the kernel, it
    //  does enforce a representation that will always be possible to save
    //  and load (and be reproducible)
    sanafe::ModelAttribute attribute;
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
        attribute.value = pybind11::cast<std::string>(value);
        TRACE1(PYMODULE, "Converted str to attribute\n");
    }
    else if (is_int)
    {
        attribute.value = pybind11::cast<int>(value);
        TRACE1(PYMODULE, "Converted int to attribute\n");
    }
    else if (is_float)
    {
        attribute.value = pybind11::cast<float>(value);
        TRACE1(PYMODULE, "Converted float to attribute\n");
    }
    else if (pybind11::isinstance<pybind11::dict>(value))
    {
        std::vector<sanafe::ModelAttribute> model_attributes;
        for (auto item : pybind11::cast<pybind11::dict>(value))
        {
            sanafe::ModelAttribute attribute = pyobject_to_model_attribute(
                    pybind11::cast<pybind11::object>(item.second));
            attribute.name = pybind11::cast<std::string>(item.first);
        }
        attribute.value = model_attributes;
    }
    else if (pybind11::isinstance<pybind11::iterable>(value))
    {
        std::vector<sanafe::ModelAttribute> list;
        for (auto iter = pybind11::iter(value);
                iter != pybind11::iterator::sentinel(); ++iter)
        {
            list.push_back(pyobject_to_model_attribute(
                    pybind11::cast<pybind11::object>(*iter)));
        }
        attribute.value = list;
    }
    else
    {
        throw std::invalid_argument("Error: dict has unsupported type");
    }

    return attribute;
}
void pyconnect_neurons_sparse(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group, const pybind11::dict attributes,
        const pybind11::list &src_dest_id_pairs)
{
    const std::map<std::string, std::vector<sanafe::ModelAttribute>>
            attribute_lists = pydict_to_attribute_lists(attributes);

    // Convert python list of tuples to SANA-FE format
    if (!pybind11::isinstance<pybind11::iterable>(src_dest_id_pairs))
    {
        std::string error(
                "Error: must provide connectivity as a list of "
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
            std::string error(
                    "Error: each entry in the src/dest id list must be a 2-tuple, "
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

    self->connect_neurons_sparse(
            dest_group, attribute_lists, src_dest_id_pairs_vec);
    return;
}

void pyconnect_neurons_conv2d(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group, const pybind11::dict attributes,
        int input_width, int input_height, int input_channels, int kernel_width,
        int kernel_height, int kernel_count, int stride_width,
        int stride_height)
{
    const std::map<std::string, std::vector<sanafe::ModelAttribute>>
            attribute_lists = pydict_to_attribute_lists(attributes);
    sanafe::Conv2DParameters config{};
    config.input_width = input_width;
    config.input_height = input_height;
    config.input_channels = input_channels;

    config.kernel_width = kernel_width;
    config.kernel_height = kernel_height;
    config.kernel_count = kernel_count;

    config.stride_width = stride_width;
    config.stride_height = stride_height;

    self->connect_neurons_conv2d(dest_group, attribute_lists, config);
    return;
}

void pyconnect_neurons_dense(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group, const pybind11::dict attributes)
{
    const std::map<std::string, std::vector<sanafe::ModelAttribute>>
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

    std::map<std::string, sanafe::ModelAttribute> model_attributes =
            pydict_to_model_attributes(model_dict);
    default_neuron_config.model_attributes = model_attributes;

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

    pipeline_config.buffer_position = sanafe::pipeline_parse_buffer_pos_str(
            buffer_position, buffer_inside_unit);
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

void pyset_attributes(sanafe::Neuron *self,
        std::optional<std::string> soma_hw_name,
        std::optional<std::string> default_synapse_hw_name,
        std::optional<std::string> dendrite_hw_name,
        std::optional<bool> log_spikes, std::optional<bool> log_potential,
        std::optional<bool> force_synapse_update,
        std::optional<bool> force_dendrite_update,
        std::optional<bool> force_soma_update, pybind11::dict model_attributes,
        pybind11::dict dendrite_specific_attributes,
        pybind11::dict soma_specific_attributes)
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

    auto parsed = pydict_to_model_attributes(model_attributes);
    neuron_template.model_attributes.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_attributes(
            dendrite_specific_attributes, false, true, false);
    neuron_template.model_attributes.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_attributes(
            soma_specific_attributes, false, false, true);
    neuron_template.model_attributes.insert(parsed.begin(), parsed.end());

    TRACE1(PYMODULE, "Setting neuron attributes\n");
    if (DEBUG_LEVEL_PYMODULE > 0)
    {
        for (auto &[key, value] : neuron_template.model_attributes)
        {
            TRACE1(PYMODULE, "\tkey: %s\n", key.c_str());
        }
    }
    self->set_attributes(neuron_template);

    return;
}

void pyset_model_attributes(sanafe::MappedNeuron *self,
        pybind11::dict model_attributes,
        pybind11::dict dendrite_specific_attributes,
        pybind11::dict soma_specific_attributes)
{
    std::map<std::string, sanafe::ModelAttribute> converted_attributes;

    auto parsed = pydict_to_model_attributes(model_attributes);
    converted_attributes.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_attributes(
            dendrite_specific_attributes, false, true, false);
    converted_attributes.insert(parsed.begin(), parsed.end());
    parsed = pydict_to_model_attributes(
            soma_specific_attributes, false, false, true);
    converted_attributes.insert(parsed.begin(), parsed.end());

    TRACE1(PYMODULE, "Setting neuron attributes\n");
    if (DEBUG_LEVEL_PYMODULE > 0)
    {
        for (auto &[key, value] : converted_attributes)
        {
            TRACE1(PYMODULE, "\tkey: %s\n", key.c_str());
        }
    }
    self->set_model_attributes(converted_attributes);

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

// Custom Python wrapper class for the Neuron class
//  While most objects are managed fine, pybind struggles to manage object
//  lifetime for Neurons, which can be accessed as slices and returned inside
//  containers
class PyNeuronRef
{
private:
    pybind11::object group_ref_; // Keeps the Python object for the group alive
    sanafe::Neuron *neuron_; // Points to the neuron in the group

public:
    PyNeuronRef(pybind11::object group_ref, sanafe::Neuron *neuron)
            : group_ref_(std::move(group_ref))
            , neuron_(neuron)
    {
    }

    sanafe::Neuron *get() const
    {
        return neuron_;
    }

    const std::vector<sanafe::Connection> &edges_out() const
    {
        return neuron_->edges_out;
    }
};

// Custom iterator class for iterating over NeuronRef objects efficiently
class NeuronGroupIterator
{
private:
    pybind11::object group_ref_;
    sanafe::NeuronGroup *group_;
    size_t current_;
    size_t size_;

public:
    NeuronGroupIterator(pybind11::object group_ref, sanafe::NeuronGroup *group)
            : group_ref_(std::move(group_ref))
            , group_(group)
            , current_(0)
            , size_(group->neurons.size())
    {
    }

    PyNeuronRef next()
    {
        if (current_ >= size_)
        {
            throw pybind11::stop_iteration();
        }
        return PyNeuronRef(group_ref_, &group_->neurons[current_++]);
    }
};

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
                    pybind11::arg("model_attributes") = pybind11::dict(),
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
                    // Keep the "self" object alive while the returned
                    //  reference is in use
                    pybind11::keep_alive<0, 1>(),
                    pybind11::return_value_policy::reference_internal);

    pybind11::class_<sanafe::NeuronGroup>(m, "NeuronGroup")
            .def("__repr__", &sanafe::NeuronGroup::info)
            .def("get_name", &sanafe::NeuronGroup::get_name)
            .def("connect_neurons_dense", &pyconnect_neurons_dense)
            .def("connect_neurons_sparse", &pyconnect_neurons_sparse)
            .def("connect_neurons_conv2d", &pyconnect_neurons_conv2d,
                    pybind11::arg("dest_group"), pybind11::arg("attributes"),
                    pybind11::arg("input_width"), pybind11::arg("input_height"),
                    pybind11::arg("input_channels"),
                    pybind11::arg("kernel_width"),
                    pybind11::arg("kernel_height"),
                    pybind11::arg("kernel_count") = 1,
                    pybind11::arg("stride_width") = 1,
                    pybind11::arg("stride_height") = 1)
            .def("__getitem__",
                    [](pybind11::object self_obj,
                            pybind11::object index) -> pybind11::object {
                        sanafe::NeuronGroup &self =
                                pybind11::cast<sanafe::NeuronGroup &>(self_obj);

                        if (pybind11::isinstance<pybind11::int_>(index))
                        {
                            // Integer indexing
                            size_t i = pybind11::cast<size_t>(index);
                            if (i >= self.neurons.size())
                            {
                                throw pybind11::index_error();
                            }
                            // Return a NeuronRef that keeps the group reference alive
                            return pybind11::cast(
                                    PyNeuronRef(self_obj, &self.neurons[i]));
                        }
                        else if (pybind11::isinstance<pybind11::slice>(index))
                        {
                            // Handle slice indexing
                            pybind11::slice slice =
                                    pybind11::cast<pybind11::slice>(index);
                            size_t start, stop, step, slice_length;
                            if (!slice.compute(self.neurons.size(), &start,
                                        &stop, &step, &slice_length))
                            {
                                throw pybind11::error_already_set();
                            }

                            // Create a Python list of NeuronRef objects
                            pybind11::list result;
                            for (size_t i = 0; i < slice_length; ++i)
                            {
                                size_t idx = start + i * step;
                                result.append(pybind11::cast(PyNeuronRef(
                                        self_obj, &self.neurons[idx])));
                            }
                            return result;
                        }
                        throw pybind11::type_error(
                                "Index must be int or slice");
                    })
            .def("__len__",
                    [](const sanafe::NeuronGroup &self) {
                        return self.neurons.size();
                    })
            .def("__iter__", [](pybind11::object self_obj) {
                sanafe::NeuronGroup &self =
                        pybind11::cast<sanafe::NeuronGroup &>(self_obj);
                return NeuronGroupIterator(self_obj, &self);
            });
    // Neuron bindings are a little more complicated, because we use a wrapper
    //  reference class to help with managing object lifetimes. This is partly
    //  due to the fact PyBind11 struggles to keep the parent alive when
    //  returning lists of them e.g., when we want to return a slice of Neuron objects
    pybind11::class_<PyNeuronRef>(m, "Neuron")
            .def("__repr__",
                    [](const PyNeuronRef &ref) {
                        return ref.get()
                                ->info(); // Delegate to neuron's info method
                    })
            // Forward all methods from Neuron to PyNeuronRef
            .def("get_id",
                    [](const PyNeuronRef &ref) { return ref.get()->get_id(); })
            .def("map_to_core",
                    [](const PyNeuronRef &ref,
                            const sanafe::CoreConfiguration &core_configuration) {
                        return ref.get()->map_to_core(core_configuration);
                    })
            .def(
                    "set_attributes",
                    [](const PyNeuronRef &ref,
                            std::optional<std::string> soma_hw_name,
                            std::optional<std::string> default_synapse_hw_name,
                            std::optional<std::string> dendrite_hw_name,
                            std::optional<bool> log_spikes,
                            std::optional<bool> log_potential,
                            std::optional<bool> force_synapse_update,
                            std::optional<bool> force_dendrite_update,
                            std::optional<bool> force_soma_update,
                            pybind11::dict model_attributes,
                            pybind11::dict dendrite_specific_attributes,
                            pybind11::dict soma_specific_attributes) {
                        pyset_attributes(ref.get(), soma_hw_name,
                                default_synapse_hw_name, dendrite_hw_name,
                                log_spikes, log_potential, force_synapse_update,
                                force_dendrite_update, force_soma_update,
                                model_attributes, dendrite_specific_attributes,
                                soma_specific_attributes);
                    },
                    pybind11::arg("soma_hw_name") = pybind11::none(),
                    pybind11::arg("default_synapse_hw_name") = pybind11::none(),
                    pybind11::arg("dendrite_hw_name") = pybind11::none(),
                    pybind11::arg("log_spikes") = pybind11::none(),
                    pybind11::arg("log_potential") = pybind11::none(),
                    pybind11::arg("force_synapse_update") = pybind11::none(),
                    pybind11::arg("force_dendrite_update") = pybind11::none(),
                    pybind11::arg("force_soma_update") = pybind11::none(),
                    pybind11::arg("model_attributes") = pybind11::dict(),
                    pybind11::arg("soma_attributes") = pybind11::dict(),
                    pybind11::arg("dendrite_attributes") = pybind11::dict())
            .def("connect_to_neuron",
                    [](const PyNeuronRef &ref, const PyNeuronRef &dest_ref,
                            std::optional<pybind11::dict> attr =
                                    pybind11::none()) -> size_t {
                        // Extract the actual neuron from the destination ref
                        sanafe::Neuron &dest = *(dest_ref.get());

                        if (!attr.has_value())
                        {
                            attr = pybind11::dict();
                        }
                        const size_t con_idx =
                                ref.get()->connect_to_neuron(dest);
                        sanafe::Connection &con = ref.get()->edges_out[con_idx];
                        const auto attributes =
                                pydict_to_model_attributes(attr.value());
                        for (auto &[key, attribute] : attributes)
                        {
                            if (attribute.forward_to_synapse)
                            {
                                con.synapse_attributes[key] = attribute;
                            }
                            if (attribute.forward_to_dendrite)
                            {
                                con.dendrite_attributes[key] = attribute;
                            }
                        }
                        return con_idx;
                    })
            // Expose edges_out as a property
            .def_property_readonly("edges_out",
                    [](const PyNeuronRef &ref) { return ref.edges_out(); });
    pybind11::class_<NeuronGroupIterator>(m, "NeuronGroupIterator")
            .def("__iter__",
                    [](NeuronGroupIterator &it) -> NeuronGroupIterator & {
                        return it;
                    })
            .def("__next__", &NeuronGroupIterator::next);
    pybind11::class_<sanafe::Connection>(m, "Connection")
            .def_readonly("pre_neuron", &sanafe::Connection::pre_neuron)
            .def_readonly("post_neuron", &sanafe::Connection::post_neuron)
            .def_readonly(
                    "synapse_hw_name", &sanafe::Connection::synapse_hw_name)
            .def("__repr__", &sanafe::Connection::info);

    pybind11::class_<sanafe::NeuronAddress>(m, "NeuronAddress")
            .def_readonly("group_name", &sanafe::NeuronAddress::group_name)
            .def_readonly(
                    "neuron_offset", &sanafe::NeuronAddress::neuron_offset)
            .def("__repr__", &sanafe::NeuronAddress::info);

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
            .def_readwrite("tiles", &sanafe::Architecture::tiles)
            .def("cores", &sanafe::Architecture::cores);

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
            .def("sim", &pysim, pybind11::arg("timesteps") = 1,
                    pybind11::arg("heartbeat") = 100,
                    pybind11::arg("timing_model") = "detailed",
                    pybind11::arg("nthreads") = 0)
            .def("get_power", &sanafe::SpikingChip::get_power)
            .def("reset", &sanafe::SpikingChip::reset);
    pybind11::class_<sanafe::MappedNeuron>(m, "MappedNeuron")
            .def("set_model_attributes", &pyset_model_attributes,
                    pybind11::arg("model_attributes") = pybind11::dict(),
                    pybind11::arg("soma_attributes") = pybind11::dict(),
                    pybind11::arg("dendrite_attributes") = pybind11::dict());
}
