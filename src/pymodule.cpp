// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// pymodule.cpp
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ios>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/gil.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pyerrors.h>
// clang-tidy tries to remove these... not sure why, but it'll break if it does
#include <pybind11/stl.h> // NOLINT(misc-include-cleaner)
#include <pybind11/stl/filesystem.h> // NOLINT(misc-include-cleaner)

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "arch.hpp"
#include "attribute.hpp"
#include "chip.hpp"
#include "docstrings.hpp"
#include "mapped.hpp"
#include "message.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"
#include "schedule.hpp"
#include "timestep.hpp"

#define PYBIND11_DETAILED_ERROR_MESSAGES

// Normally we don't suppress these warnings, but in this file, Python allows for
//  named arguments - removing the risk of accidentally swapping args and
//  making it easier to specify a large number of args
//  PyBind11 also relies on a super long function macro, so suppress this too
// NOLINTBEGIN(bugprone-easily-swappable-parameters,readability-function-size)

// Forward declarations
namespace // anonymous
{
std::map<std::string, sanafe::ModelAttribute> pydict_to_model_attributes(
        const pybind11::dict &dictionary, bool forward_to_synapse = true,
        bool forward_to_dendrite = true, bool forward_to_soma = true);
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
        const auto key = pybind11::cast<std::string>(key_value_pair.first);
        TRACE1(PYMODULE, "Adding dict val: dict['%s']\n", key.c_str());

        const auto value =
                pybind11::cast<pybind11::object>(key_value_pair.second);
        sanafe::ModelAttribute attribute = pyobject_to_model_attribute(value);

        attribute.forward_to_synapse = forward_to_synapse;
        attribute.forward_to_dendrite = forward_to_dendrite;
        attribute.forward_to_soma = forward_to_soma;
        attribute.name = key;
        map[key] = std::move(attribute);
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
        const auto attribute_name =
                pybind11::cast<std::string>(key_value_pair.first);
        if (!pybind11::isinstance<pybind11::iterable>(key_value_pair.second))
        {
            const std::string error(
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity,readability-function-size)
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
    const bool is_string = pybind11::isinstance<pybind11::str>(value);

    const bool is_numpy = pybind11::hasattr(value, "dtype");
    if (is_numpy)
    {
        auto dtype = value.attr("dtype");
        const auto code = dtype.attr("kind").cast<std::string>();
        is_int = (code == "i" || code == "u");
        is_float = (code == "f");
    }
    else
    {
        is_int = pybind11::isinstance<pybind11::int_>(value);
        is_float = pybind11::isinstance<pybind11::float_>(value);
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

pybind11::object pymodel_attribute_to_pyobject(
        const sanafe::ModelAttribute &attribute)
{
    if (std::holds_alternative<bool>(attribute.value))
    {
        auto bool_value = static_cast<bool>(attribute);
        return pybind11::cast(bool_value);
    }
    if (std::holds_alternative<int>(attribute.value))
    {
        auto int_value = static_cast<int>(attribute);
        return pybind11::cast(int_value);
    }
    if (std::holds_alternative<double>(attribute.value))
    {
        auto float_val = static_cast<double>(attribute);
        return pybind11::cast(float_val);
    }
    if (std::holds_alternative<std::string>(attribute.value))
    {
        auto str_value = static_cast<std::string>(attribute);
        return pybind11::cast(str_value);
    }
    if (std::holds_alternative<std::vector<sanafe::ModelAttribute>>(
                attribute.value))
    {
        auto attribute_vec =
                static_cast<std::vector<sanafe::ModelAttribute>>(attribute);
        // If the vector of sub-attributes is either empty, or seems to be
        //  named attributes (i.e., a mapping)
        if (!attribute_vec.empty() && attribute_vec[0].name.has_value())
        {
            pybind11::dict pydict{};
            for (const auto &sub_attribute : attribute_vec)
            {
                if (!sub_attribute.name.has_value())
                {
                    INFO("Error: Sub-attribute is unnamed in mapping.\n");
                    throw std::runtime_error(
                            "Error: Sub-attribute is unnamed.\n");
                }
                pydict[pybind11::str(sub_attribute.name.value())] =
                        pymodel_attribute_to_pyobject(sub_attribute);
            }
            return pydict;
        }

        pybind11::list pylist{};
        for (const auto &sub_attribute : attribute_vec)
        {
            pylist.append(pymodel_attribute_to_pyobject(sub_attribute));
        }
        return pylist;
    }

    throw std::runtime_error("Unrecognized model attribute type\n");
}

pybind11::dict pymodel_attributes_to_pydict(
        const std::map<std::string, sanafe::ModelAttribute> &model_attributes)
{
    pybind11::dict attribute_dict{};

    for (const auto &[key, attribute] : model_attributes)
    {
        attribute_dict[pybind11::str(key)] =
                pymodel_attribute_to_pyobject(attribute);
    }

    return attribute_dict;
}

using PerfStatistic = std::variant<size_t, long int, double>;

std::map<std::string, PerfStatistic> timestep_data_to_map(
        const sanafe::Timestep &ts)
{
    std::map<std::string, PerfStatistic> ts_map;

    ts_map["timestep"] = ts.timestep;
    ts_map["fired"] = ts.neurons_fired;
    ts_map["updated"] = ts.neurons_updated;
    ts_map["hops"] = ts.total_hops;
    ts_map["spikes"] = ts.spike_count;
    ts_map["sim_time"] = ts.sim_time;
    ts_map["total_energy"] = ts.total_energy;

    return ts_map;
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
    run_data_dict["neurons_updated"] = run.neurons_updated;
    run_data_dict["neurons_fired"] = run.neurons_fired;

    return run_data_dict;
}

pybind11::dict message_to_dict(const sanafe::Message &m)
{
    pybind11::dict message_dict;

    message_dict["generation_delay"] = m.generation_delay;
    message_dict["network_delay"] = m.network_delay;
    message_dict["receive_delay"] = m.receive_delay;
    message_dict["blocked_delay"] = m.blocked_delay;
    message_dict["sent_timestamp"] = m.sent_timestamp;
    message_dict["received_timestamp"] = m.received_timestamp;
    message_dict["processed_timestamp"] = m.processed_timestamp;
    message_dict["timestep"] = m.timestep;
    message_dict["mid"] = m.mid;
    message_dict["spikes"] = m.spikes;
    message_dict["hops"] = m.hops;
    message_dict["src_neuron_offset"] = m.src_neuron_offset;
    message_dict["src_neuron_group_id"] = m.src_neuron_group_id;
    message_dict["src_x"] = m.src_x;
    message_dict["dest_x"] = m.dest_x;
    message_dict["src_y"] = m.src_y;
    message_dict["dest_y"] = m.dest_y;
    message_dict["src_tile_id"] = m.src_tile_id;
    message_dict["src_core_id"] = m.src_core_id;
    message_dict["src_core_offset"] = m.src_core_offset;
    message_dict["dest_tile_id"] = m.dest_tile_id;
    message_dict["dest_core_id"] = m.dest_core_id;
    message_dict["dest_core_offset"] = m.dest_core_offset;
    message_dict["dest_axon_hw"] = m.dest_axon_hw;
    message_dict["dest_axon_id"] = m.dest_axon_id;
    message_dict["placeholder"] = m.placeholder;

    return message_dict;
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
        const std::string error(
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
            const std::string error(
                    "Error: each entry in the src/dest id list must be a 2-tuple, "
                    "representing the neuron's offsets within the src and dest neuron groups");
            INFO("%s\n", error.c_str());
            throw pybind11::value_error(error);
        }

        const auto seq = pybind11::cast<pybind11::sequence>(*iter);
        if (seq.size() != 2)
        {
            throw pybind11::value_error("Expected a 2-tuple");
        }
        const size_t src = pybind11::cast<size_t>(seq[0]);
        const size_t dest = pybind11::cast<size_t>(seq[1]);
        src_dest_id_pairs_vec.emplace_back(src, dest);
    }

    self->connect_neurons_sparse(
            dest_group, attribute_lists, src_dest_id_pairs_vec);
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
}

void pyconnect_neurons_dense(sanafe::NeuronGroup *self,
        sanafe::NeuronGroup &dest_group, const pybind11::dict attributes)
{
    const std::map<std::string, std::vector<sanafe::ModelAttribute>>
            attribute_lists = pydict_to_attribute_lists(attributes);

    self->connect_neurons_dense(dest_group, attribute_lists);
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
    default_neuron_config.soma_hw_name = std::move(soma_hw_name);

    std::map<std::string, sanafe::ModelAttribute> model_attributes =
            pydict_to_model_attributes(model_dict);
    default_neuron_config.model_attributes = std::move(model_attributes);

    return self->create_neuron_group(
            group_name, neuron_count, default_neuron_config);
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

    return self->create_tile(std::move(name), tile_power_metrics);
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

    return self->create_core(std::move(name), parent_tile_id, pipeline_config);
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
#if (DEBUG_LEVEL_PYMODULE > 0)
    for (auto &[key, value] : neuron_template.model_attributes)
    {
        TRACE1(PYMODULE, "\tkey: %s\n", key.c_str());
    }
#endif
    self->set_attributes(neuron_template);
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
#if (DEBUG_LEVEL_PYMODULE > 0)
    for (auto &[key, value] : converted_attributes)
    {
        TRACE1(PYMODULE, "\tkey: %s\n", key.c_str());
        INFO("\tkey: %s\n", key.c_str());
    }
#endif
    self->set_model_attributes(converted_attributes);
}

// Determine spike trace mode
enum TraceMode : uint8_t
{
    trace_none = 0UL, // No tracing
    trace_file = 1UL, // Write to file
    trace_memory = 2UL, // Store in memory for return
};

class PyTrace
{
public:
    PyTrace(sanafe::SpikingChip *chip, pybind11::object &trace_obj,
            const bool overwrite_trace)
            : parent_chip(chip)
    {
        open(trace_obj, overwrite_trace);
    }
    virtual ~PyTrace()
    {
        if (file.is_open())
        {
            file.close();
        }
    }
    PyTrace(const PyTrace &copy) = delete;
    PyTrace(const PyTrace &&move) = delete;
    PyTrace &operator=(const PyTrace &copy) = delete;
    PyTrace &operator=(PyTrace &&move) = delete;

    void write_header();
    void update(const sanafe::Timestep &ts);
    virtual pybind11::object get_python_object() const = 0;

protected:
    sanafe::SpikingChip *parent_chip;
    pybind11::object obj{pybind11::none()};
    TraceMode mode{TraceMode::trace_none};
    std::ofstream file;

private:
    void open(pybind11::object trace_obj, bool overwrite_trace);

    virtual void get_trace_string(
            std::ostringstream &ss, const sanafe::Timestep &ts) = 0;
    virtual void get_header_string(std::ostringstream &ss) = 0;
    virtual void append_trace_data(const sanafe::Timestep &ts) = 0;
};

void PyTrace::open(pybind11::object trace_obj, const bool overwrite_trace)
{
    mode = TraceMode::trace_none;
    // The trace object can be None, a boolean (True/False) enabling memory
    //  based traces, a string specifying a trace filename or a Python file
    //  object
    obj = trace_obj;

    if (!trace_obj.is_none())
    {
        // Check if it's a boolean (True)
        if (pybind11::isinstance<pybind11::bool_>(trace_obj))
        {
            const bool trace_enabled = trace_obj.cast<bool>();
            if (trace_enabled)
            {
                mode = TraceMode::trace_memory;
            }
        }
        // Check if it's a file-like object (i.e., has a 'write' method)
        else if (pybind11::hasattr(trace_obj, "write"))
        {
            mode = TraceMode::trace_file;
            // File handle will be used directly via Python calls

            if (overwrite_trace)
            {
                obj.attr("seek")(0);
            }
        }
        // Check if it's a string (filename)
        else if (pybind11::isinstance<pybind11::str>(trace_obj))
        {
            mode = TraceMode::trace_file;
            const auto filename = trace_obj.cast<std::string>();

            const std::ios::openmode flags = overwrite_trace ?
                    (std::ios::out | std::ios::trunc) :
                    (std::ios::out | std::ios::app);
            file.open(filename, flags);
            if (!file.is_open())
            {
                throw std::runtime_error(
                        "Failed to open trace file: " + filename);
            }
        }
        else
        {
            throw std::invalid_argument(
                    "trace_obj must be None, True, a filename string, "
                    "or a file-like object");
        }
    }
}

void PyTrace::update(const sanafe::Timestep &ts)
{
    if (mode == TraceMode::trace_memory)
    {
        // Note we must always acquire the GIL before modifying Python objects
        const pybind11::gil_scoped_acquire acquire;
        // Make call to appropriate overriden virtual routine
        append_trace_data(ts);
    }
    else if (mode == TraceMode::trace_file)
    {
        std::ostringstream trace_ss;
        // Make call to appropriate overriden virtual routine
        get_trace_string(trace_ss, ts);
        if (file.is_open())
        {
            // Write to C++ file stream
            file << trace_ss.str();
            file.flush(); // Ensure data is written
        }
        else if (!obj.is_none())
        {
            // Write to Python file handle
            const pybind11::gil_scoped_acquire acquire_for_write;
            obj.attr("write")(trace_ss.str());
            if (pybind11::hasattr(obj, "flush"))
            {
                obj.attr("flush")();
            }
        }
    }
    // else do nothing
}

void PyTrace::write_header()
{
    if (mode == TraceMode::trace_file)
    {
        std::ostringstream trace_ss;
        get_header_string(trace_ss);
        if (file.is_open())
        {
            // Write to C++ file stream
            file << trace_ss.str();
            file.flush(); // Ensure data is written
        }
        else if (!obj.is_none())
        {
            // Write to Python file handle
            const pybind11::gil_scoped_acquire acquire_for_write;
            obj.attr("write")(trace_ss.str());
            if (pybind11::hasattr(obj, "flush"))
            {
                obj.attr("flush")();
            }
        }
    }
}

class PySpikeTrace : public PyTrace
{
public:
    PySpikeTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyTrace(chip, trace_obj, overwrite_trace)
    {
    }
    ~PySpikeTrace() override = default;
    void get_header_string(std::ostringstream &ss) override
    {
        sanafe::SpikingChip::sim_trace_write_spike_header(ss);
    }
    void get_trace_string(
            std::ostringstream &ss, const sanafe::Timestep &ts) override
    {
        parent_chip->sim_trace_record_spikes(ss, ts.timestep);
    }
    void append_trace_data(const sanafe::Timestep & /*ts*/) override
    {
        auto spikes = parent_chip->get_spikes();
        data.emplace_back(std::move(spikes));
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            pybind11::list data_list;
            for (const auto &spikes : data)
            {
                pybind11::list spike_list;
                for (const sanafe::NeuronAddress &spike : spikes)
                {
                    spike_list.append(pybind11::cast(spike));
                }
                data_list.append(spike_list);
            }
            return data_list;
        }
        return pybind11::none();
    }

private:
    std::vector<std::vector<sanafe::NeuronAddress>> data;
};

class PyPotentialTrace : public PyTrace
{
public:
    PyPotentialTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyTrace(chip, trace_obj, overwrite_trace)
    {
    }
    ~PyPotentialTrace() override = default;
    void get_header_string(std::ostringstream &ss) override
    {
        parent_chip->sim_trace_write_potential_header(ss);
    }
    void get_trace_string(
            std::ostringstream &ss, const sanafe::Timestep &ts) override
    {
        parent_chip->sim_trace_record_potentials(ss, ts.timestep);
    }
    void append_trace_data(const sanafe::Timestep & /*ts*/) override
    {
        auto potentials = parent_chip->get_potentials();
        data.emplace_back(std::move(potentials));
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            return pybind11::cast(data);
        }
        return pybind11::none();
    }

private:
    std::vector<std::vector<double>> data;
};

class PyPerfTrace : public PyTrace
{
public:
    ~PyPerfTrace() override = default;
    PyPerfTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyTrace(chip, trace_obj, overwrite_trace)
    {
    }
    void get_header_string(std::ostringstream &ss) override
    {
        parent_chip->sim_trace_write_perf_header(ss);
    }
    void get_trace_string(
            std::ostringstream &ss, const sanafe::Timestep &ts) override
    {
        parent_chip->sim_trace_record_perf(ss, ts);
    }
    void append_trace_data(const sanafe::Timestep &ts) override
    {
        std::map<std::string, PerfStatistic> stats = timestep_data_to_map(ts);

        // Get any optional traces and cast double values to Python objects
        const std::map<std::string, double> optional_perf_traces =
                parent_chip->sim_trace_get_optional_traces();
        for (const auto &[name, value] : optional_perf_traces)
        {
            stats[name] = value;
        }

        for (const auto &[key, value] : stats)
        {
            data[key].emplace_back(value);
        }
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            return pybind11::cast(data);
        }
        return pybind11::none();
    }

private:
    std::map<std::string, std::vector<PerfStatistic>> data;
};

class PyMessageTrace : public PyTrace
{
public:
    PyMessageTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyTrace(chip, trace_obj, overwrite_trace)
    {
    }
    ~PyMessageTrace() override = default;
    void get_header_string(std::ostringstream &ss) override
    {
        sanafe::SpikingChip::sim_trace_write_message_header(ss);
    }
    void get_trace_string(
            std::ostringstream &ss, const sanafe::Timestep &ts) override
    {
        std::vector<std::reference_wrapper<const sanafe::Message>> all_messages;
        for (const sanafe::MessageFifo &q : ts.messages)
        {
            for (const sanafe::Message &m : q)
            {
                all_messages.emplace_back(m);
            }
        }
        // Sort messages in message ID order (with placeholders last)
        std::sort(all_messages.begin(), all_messages.end(),
                sanafe::CompareMessagesByID{});
        // Save the messages in sorted order
        for (const sanafe::Message &m : all_messages)
        {
            sim_trace_record_message(ss, m);
        }
    }
    void append_trace_data(const sanafe::Timestep &ts) override
    {
        std::vector<sanafe::Message> timestep_messages;

        // Not crucial, but its nice to print messages in ID order.
        //  Copy all the messages into a single vector and sort
        for (const sanafe::MessageFifo &q : ts.messages)
        {
            for (const sanafe::Message &m : q)
            {
                timestep_messages.emplace_back(m);
            }
        }
        // Sort messages in message order, starting at lowest mid first
        std::sort(timestep_messages.begin(), timestep_messages.end(),
                [](const sanafe::Message &left, const sanafe::Message &right) {
                    return left.mid < right.mid;
                });
        data.emplace_back(std::move(timestep_messages));
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            // Convert each message to a dict, meaning we need to manually
            //  iterate over every message in data
            pybind11::list data_list;
            for (const auto &ts_messages : data)
            {
                pybind11::list ts_list;
                for (const auto &m : ts_messages)
                {
                    ts_list.append(message_to_dict(m));
                }
                data_list.append(ts_list);
            }
            return data_list;
        }
        return pybind11::none();
    }

private:
    std::vector<std::vector<sanafe::Message>> data;
};

// A similar implementation to SpikingChip::flush_timestep_data, but supporting
//  Python interfaces and the Python file system
void pyflush_timestep_data(sanafe::SpikingChip *self, sanafe::RunData &rd,
        PyPerfTrace &perf, sanafe::Scheduler &scheduler,
        PyMessageTrace &message_trace)
{
    while (!scheduler.timesteps_to_write.empty())
    {
        sanafe::TimestepHandle timestep_handle;
        scheduler.timesteps_to_write.pop(timestep_handle);
        sanafe::Timestep &ts = timestep_handle.get();
        TRACE1(CHIP, "retiring ts:%ld\n", ts.timestep);

        perf.update(ts);
        message_trace.update(ts);
        sanafe::SpikingChip::update_run_data(rd, ts);
        self->retire_timestep(ts);
    }
}

pybind11::dict pysim(sanafe::SpikingChip *self, const long int timesteps,
        std::string timing_model_str, const int processing_threads,
        const int scheduler_threads, pybind11::object spike_trace,
        pybind11::object potential_trace, pybind11::object perf_trace,
        pybind11::object message_trace, const bool write_trace_headers)
{
    sanafe::TimingModel timing_model{sanafe::timing_model_detailed};
    timing_model = sanafe::parse_timing_model(timing_model_str);

#ifndef HAVE_OPENMP
    pybind11::print("Warning: multiple threads not supported; flag ignored");
#else
    const int total_threads_available = omp_get_num_procs();
    omp_set_num_threads(std::min(total_threads_available, processing_threads));
#endif

    const bool overwrite_if_trace_exists = write_trace_headers;
    PySpikeTrace spikes(self, spike_trace, overwrite_if_trace_exists);
    PyPotentialTrace potentials(
            self, potential_trace, overwrite_if_trace_exists);
    PyMessageTrace messages(self, message_trace, overwrite_if_trace_exists);
    PyPerfTrace perf(self, perf_trace, overwrite_if_trace_exists);

    if (write_trace_headers)
    {
        spikes.write_header();
        potentials.write_header();
        messages.write_header();
        perf.write_header();
    }

    sanafe::RunData rd(self->get_total_timesteps() + 1);
    rd.timesteps_executed += timesteps;

    sanafe::Scheduler scheduler;
    scheduler.noc_width_in_tiles = self->noc_width_in_tiles;
    scheduler.noc_height_in_tiles = self->noc_height_in_tiles;
    scheduler.buffer_size = self->noc_buffer_size;
    scheduler.core_count = self->core_count;
    scheduler.max_cores_per_tile = self->max_cores_per_tile;
    scheduler.timing_model = timing_model;
    schedule_create_threads(scheduler, scheduler_threads);

    auto last_check = std::chrono::steady_clock::now();
    auto last_print = std::chrono::steady_clock::now();

    const std::string first_message =
            "Executed steps: [0/" + std::to_string(timesteps) + "]";
    pybind11::print(first_message, pybind11::arg("end") = "");

    const pybind11::gil_scoped_release release;
    for (long int timestep = 1; timestep <= timesteps; timestep++)
    {
        self->step(scheduler);
        spikes.update(timestep);
        potentials.update(timestep);

        const auto now = std::chrono::steady_clock::now();
        constexpr std::chrono::milliseconds check_interval{100};
        constexpr std::chrono::seconds print_interval{1};
        if ((now - last_check) >= check_interval)
        {
            const pybind11::gil_scoped_acquire acquire;
            if (PyErr_CheckSignals() != 0)
            {
                throw pybind11::error_already_set();
            }
            last_check = now;
        }
        if ((now - last_print) >= print_interval)
        {
            const std::string message = "\rExecuted steps: [" +
                    std::to_string(timestep) + "/" + std::to_string(timesteps) +
                    "]";
            const pybind11::gil_scoped_acquire acquire;
            pybind11::print(message, pybind11::arg("end") = "");
            last_print = now;
        }
        // Retire and trace messages as we go along. Not strictly required, but
        //  avoids the message write queue growing too large
        pyflush_timestep_data(self, rd, perf, scheduler, messages);
    }

    const std::string last_message = "\rExecuted steps: [" +
            std::to_string(timesteps) + "/" + std::to_string(timesteps) + "]";
    const pybind11::gil_scoped_acquire acquire;
    pybind11::print(last_message);

    schedule_stop_all_threads(scheduler);
    pyflush_timestep_data(self, rd, perf, scheduler, messages);

    const auto spike_data = spikes.get_python_object();
    const auto potential_data = potentials.get_python_object();
    const auto perf_data = perf.get_python_object();
    const auto message_data = messages.get_python_object();

    pybind11::dict result = run_data_to_dict(rd);
    result["spike_trace"] = spike_data;
    result["potential_trace"] = potential_data;
    result["perf_trace"] = perf_data;
    result["message_trace"] = message_data;

    return result;
}

} // end of anonymous namespace

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

    [[nodiscard]] sanafe::Neuron *get() const
    {
        return neuron_;
    }

    [[nodiscard]] const std::vector<sanafe::Connection> &edges_out() const
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
    size_t current_{0};
    size_t size_;

public:
    NeuronGroupIterator(pybind11::object group_ref, sanafe::NeuronGroup *group)
            : group_ref_(std::move(group_ref))
            , group_(group)
            , size_(group->neurons.size())
    {
    }

    PyNeuronRef next()
    {
        if (current_ >= size_)
        {
            throw pybind11::stop_iteration();
        }
        return {group_ref_, &group_->neurons[current_++]};
    }
};

// NOLINTBEGIN(readability-function-cognitive-complexity)
PYBIND11_MODULE(sanafecpp, m)
{
    m.doc() = docstrings::module_doc;

    m.def("load_arch", &sanafe::load_arch, docstrings::load_arch_doc);
    m.def("load_net", &sanafe::load_net, docstrings::load_net_doc,
            pybind11::arg("path"), pybind11::arg("arch"),
            pybind11::arg("use_netlist_format") = false);

    pybind11::enum_<sanafe::BufferPosition>(m, "BufferPosition")
            .value("buffer_before_dendrite_unit",
                    sanafe::buffer_before_dendrite_unit)
            .value("buffer_before_soma_unit", sanafe::buffer_before_soma_unit)
            .value("buffer_before_axon_out_unit",
                    sanafe::buffer_before_axon_out_unit)
            .value("buffer_positions", sanafe::buffer_positions);

    pybind11::class_<sanafe::SpikingNetwork>(
            m, "Network", docstrings::network_doc)
            .def(pybind11::init<>())
            .def("__repr__", &sanafe::SpikingNetwork::info)
            .def("create_neuron_group", &pycreate_neuron_group,
                    docstrings::network_create_neuron_group_doc,
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
            .def("save", &sanafe::SpikingNetwork::save,
                    docstrings::network_save_doc, pybind11::arg("path"),
                    pybind11::arg("use_netlist_format") = false)
            .def_property_readonly(
                    "groups",
                    [](sanafe::SpikingNetwork &self)
                            -> std::map<std::string, sanafe::NeuronGroup> & {
                        return self.groups;
                    },
                    pybind11::return_value_policy::reference_internal)
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

    pybind11::class_<sanafe::NeuronGroup>(
            m, "NeuronGroup", docstrings::neuron_group_doc)
            .def("__repr__", &sanafe::NeuronGroup::info)
            .def("get_name", &sanafe::NeuronGroup::get_name)
            .def("connect_neurons_dense", &pyconnect_neurons_dense,
                    docstrings::neuron_group_connect_dense_doc)
            .def("connect_neurons_sparse", &pyconnect_neurons_sparse,
                    docstrings::neuron_group_connect_sparse_doc)
            .def("connect_neurons_conv2d", &pyconnect_neurons_conv2d,
                    docstrings::neuron_group_connect_conv2d_doc,
                    pybind11::arg("dest_group"), pybind11::arg("attributes"),
                    pybind11::arg("input_width"), pybind11::arg("input_height"),
                    pybind11::arg("input_channels"),
                    pybind11::arg("kernel_width"),
                    pybind11::arg("kernel_height"),
                    pybind11::arg("kernel_count") = 1,
                    pybind11::arg("stride_width") = 1,
                    pybind11::arg("stride_height") = 1)
            .def_property_readonly(
                    "neurons",
                    [](pybind11::object &self_obj) -> pybind11::object {
                        auto &self =
                                pybind11::cast<sanafe::NeuronGroup &>(self_obj);
                        std::vector<pybind11::object> neuron_refs;
                        for (sanafe::Neuron &neuron : self.neurons)
                        {
                            neuron_refs.emplace_back(pybind11::cast(
                                    PyNeuronRef(self_obj, &neuron)));
                        }
                        return pybind11::cast(neuron_refs);
                    },
                    pybind11::return_value_policy::reference_internal)
            .def("__getitem__",
                    [](pybind11::object self_obj,
                            pybind11::object index) -> pybind11::object {
                        auto &self =
                                pybind11::cast<sanafe::NeuronGroup &>(self_obj);
                        if (pybind11::isinstance<pybind11::int_>(index))
                        {
                            // Integer indexing
                            const auto i = pybind11::cast<size_t>(index);
                            if (i >= self.neurons.size())
                            {
                                throw pybind11::index_error();
                            }
                            // Return a NeuronRef that keeps the group reference alive
                            return pybind11::cast(
                                    PyNeuronRef(self_obj, &self.neurons[i]));
                        }
                        if (pybind11::isinstance<pybind11::slice>(index))
                        {
                            // Handle slice indexing
                            const auto slice =
                                    pybind11::cast<pybind11::slice>(index);
                            size_t start = 0;
                            size_t stop = 0;
                            size_t step = 0;
                            size_t slice_length = 0;
                            if (!slice.compute(self.neurons.size(), &start,
                                        &stop, &step, &slice_length))
                            {
                                throw pybind11::error_already_set();
                            }

                            // Create a Python list of NeuronRef objects
                            pybind11::list result;
                            for (size_t i = 0; i < slice_length; ++i)
                            {
                                const size_t idx = start + (i * step);
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
                auto &self = pybind11::cast<sanafe::NeuronGroup &>(self_obj);
                return NeuronGroupIterator(self_obj, &self);
            });
    // Neuron bindings are a little more complicated, because we use a wrapper
    //  reference class to help with managing object lifetimes. This is partly
    //  due to the fact PyBind11 struggles to keep the parent alive when
    //  returning lists of them e.g., when we want to return a slice of Neuron objects
    pybind11::class_<PyNeuronRef>(m, "Neuron", docstrings::neuron_doc)
            // Forward all methods from Neuron to PyNeuronRef
            .def("__repr__",
                    [](const PyNeuronRef &ref) { return ref.get()->info(); })
            .def("get_id",
                    [](const PyNeuronRef &ref) { return ref.get()->get_id(); })
            .def(
                    "map_to_core",
                    [](const PyNeuronRef &ref,
                            const sanafe::CoreConfiguration &core_configuration) {
                        ref.get()->map_to_core(core_configuration);
                        return;
                    },
                    docstrings::neuron_map_to_core_doc)
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
                        pyset_attributes(ref.get(), std::move(soma_hw_name),
                                std::move(default_synapse_hw_name),
                                std::move(dendrite_hw_name), log_spikes,
                                log_potential, force_synapse_update,
                                force_dendrite_update, force_soma_update,
                                model_attributes, dendrite_specific_attributes,
                                soma_specific_attributes);
                    },
                    docstrings::neuron_set_attributes_doc,
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
            .def(
                    "connect_to_neuron",
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
                        for (const auto &[key, attribute] : attributes)
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
                    },
                    docstrings::neuron_connect_to_neuron_doc)
            // Expose edges_out as a property
            .def_property_readonly("edges_out",
                    [](const PyNeuronRef &ref)
                            -> const std::vector<sanafe::Connection> & {
                        return ref.edges_out();
                    });
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
            .def_property_readonly("synapse_attributes",
                    [](const sanafe::Connection &self) -> pybind11::object {
                        return pymodel_attributes_to_pydict(
                                self.synapse_attributes);
                    })
            .def("__repr__", &sanafe::Connection::info);

    pybind11::class_<sanafe::NeuronAddress>(m, "NeuronAddress")
            .def_readonly("group_name", &sanafe::NeuronAddress::group_name)
            .def_readonly(
                    "neuron_offset", &sanafe::NeuronAddress::neuron_offset)
            .def("__repr__", &sanafe::NeuronAddress::info)
            // Make this tiny struct picklable so that we can easily dump it to
            //  YAML in Python
            .def(pybind11::pickle(
                    [](const sanafe::NeuronAddress &n) { // __getstate__
                        // All members except 'in_noc', which is only used by the
                        //  detailed scheduler
                        return pybind11::make_tuple(
                                n.group_name, n.neuron_offset);
                    },
                    [](const pybind11::tuple &t)
                            -> sanafe::NeuronAddress { // __setstate__
                        sanafe::NeuronAddress n;
                        n.group_name = pybind11::cast<std::string>(t[0]);
                        n.neuron_offset = pybind11::cast<long int>(t[1]);
                        return n;
                    }));
    pybind11::class_<sanafe::Architecture>(
            m, "Architecture", docstrings::architecture_doc)
            .def(pybind11::init<std::string,
                    sanafe::NetworkOnChipConfiguration>())
            .def("__repr__", &sanafe::Architecture::info)
            .def("create_tile", &pycreate_tile,
                    docstrings::architecture_create_tile_doc,
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
            .def("create_core", &pycreate_core,
                    docstrings::architecture_create_core_doc,
                    pybind11::arg("name"), pybind11::arg("parent_tile_id"),
                    pybind11::arg("buffer_position") =
                            sanafe::buffer_before_soma_unit,
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
                            sanafe::buffer_before_soma_unit,
                    pybind11::arg("buffer_inside_unit") = false,
                    pybind11::arg("max_neurons_supported") =
                            sanafe::default_max_neurons,
                    pybind11::arg("log_energy") = false,
                    pybind11::arg("log_latency") = false);
    pybind11::class_<sanafe::MappedNeuron>(m, "MappedNeuron")
            .def("set_model_attributes", &pyset_model_attributes,
                    pybind11::arg("model_attributes") = pybind11::dict(),
                    pybind11::arg("soma_attributes") = pybind11::dict(),
                    pybind11::arg("dendrite_attributes") = pybind11::dict());
    pybind11::class_<sanafe::SpikingChip>(
            m, "SpikingChip", docstrings::spiking_chip_doc)
            .def_property(
                    "mapped_neuron_groups",
                    [](sanafe::SpikingChip &self)
                            -> std::map<std::string,
                                    std::vector<std::reference_wrapper<
                                            sanafe::MappedNeuron>>> & {
                        return self.mapped_neuron_groups;
                    },
                    nullptr, pybind11::return_value_policy::reference_internal)
            .def(pybind11::init<sanafe::Architecture &>(),
                    pybind11::arg("arch"))
            .def("load", &sanafe::SpikingChip::load,
                    docstrings::spiking_chip_load_doc)
            .def("sim", &pysim, docstrings::spiking_chip_sim_doc,
                    pybind11::arg("timesteps") = 1,
                    pybind11::arg("timing_model") = "detailed",
                    pybind11::arg("processing_threads") = 0,
                    pybind11::arg("scheduler_threads") = 0,
                    pybind11::arg("spike_trace") = pybind11::none(),
                    pybind11::arg("potential_trace") = pybind11::none(),
                    pybind11::arg("perf_trace") = pybind11::none(),
                    pybind11::arg("message_trace") = pybind11::none(),
                    pybind11::arg("write_trace_headers") = true)
            .def("get_power", &sanafe::SpikingChip::get_power,
                    docstrings::spiking_chip_get_power_doc)
            .def("reset", &sanafe::SpikingChip::reset,
                    docstrings::spiking_chip_reset_doc);
}
// NOLINTEND(readability-function-cognitive-complexity)
// NOLINTEND(bugprone-easily-swappable-parameters,readability-function-size)
