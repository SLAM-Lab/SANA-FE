// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include "arch.hpp"
#include "chip.hpp"
#include "mapped.hpp"
#include "pipeline.hpp"
#include "print.hpp"

void sanafe::PipelineUnit::check_implemented(
        const bool check_implements_synapse,
        const bool check_implements_dendrite,
        const bool check_implements_soma) const
{
    if (check_implements_synapse != implements_synapse)
    {
        const std::string error = check_implements_synapse ?
                "Unit does not support synapse functionality, remove from arch "
                "description's 'synapse' section" :
                "Unit supports synapse functionality but was not added to "
                "'synapse' hardware section. Add it or pick a different model.";
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }

    if (check_implements_dendrite != implements_dendrite)
    {
        const std::string error = check_implements_dendrite ?
                "Unit does not support dendrite functionality, remove from arch "
                "description's 'dendrite' section" :
                "Unit supports dendrite functionality but was not added to "
                "'dendrite' hardware section. Add it or pick a different model.";
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }

    if (check_implements_soma != implements_soma)
    {
        const std::string error = check_implements_soma ?
                "Unit does not support soma functionality, remove from arch "
                "description's 'soma' section" :
                "Unit supports soma functionality but was not added to "
                "'soma' hardware section. Add it or pick a different model.";
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }

}

void sanafe::PipelineUnit::check_outputs(
        const MappedNeuron & /*n*/, const PipelineResult &result) const
{
    // Check the hw returns a valid value for the next unit to process
    if (implements_soma && result.status == invalid_neuron_state)
    {
        throw std::runtime_error("Soma output; should return valid "
                                 "neuron state.");
    }
    if (!implements_soma && (implements_synapse || implements_dendrite) &&
            !result.current.has_value())
    {
        throw std::runtime_error("Synaptic or dendritic output; should return "
                                 "synaptic/dendritic current");
    }
}

void sanafe::PipelineUnit::register_attributes(
        const std::set<std::string> &attribute_names)
{
    for (const auto &attr : attribute_names)
    {
        supported_attribute_names.insert(attr);
    }
}

void sanafe::PipelineUnit::check_attribute(const std::string attribute_name)
{
    if (supported_attribute_names.find(attribute_name) ==
                    supported_attribute_names.end() &&
            (attribute_warnings <= max_attribute_warnings))
    {
        INFO("Warning: Attribute (%s) not supported by model: %s, will be "
             "ignored.\nEither remove this attribute from the SNN/Architecture "
             "description file, be more specific which hardware requires it "
             "(using synapse/dendrite/soma sections) or register the attribute "
             "in the constructor to suppress this warning.\n",
                attribute_name.c_str(), name.c_str());
        if (attribute_warnings == max_attribute_warnings)
        {
            INFO("Warning: Reached max attribute warnings (%ld) for this h/w "
                 "unit and future warnings will be suppressed. To print all, "
                 "rebuild and increase PipelineUnit::max_attribute_warnings\n",
                    max_attribute_warnings);
        }
        attribute_warnings++;
    }
}

void sanafe::PipelineUnit::set_attributes(
        std::string unit_name, const ModelInfo &model)
{
    model_attributes = model.model_attributes;
    plugin_lib = model.plugin_library_path;
    name = std::move(unit_name);
    log_energy = model.log_energy;
    log_latency = model.log_latency;

    synapse_set_default_attributes();
    dendrite_set_default_attributes();
    soma_set_default_attributes();

    // Finally, forward all attributes from the architecture description to the
    //  model. This might be useful if you want to define any additional
    //  model-specific attributes here, e.g., fault-rate or maximum memory size.
    for (auto &[key, attribute] : model_attributes)
    {
        check_attribute(key);
        set_attribute_hw(key, attribute);
    }
}

void sanafe::PipelineUnit::synapse_set_default_attributes()
{
    if (model_attributes.find("energy_process_spike") != model_attributes.end())
    {
        default_energy_process_spike =
                static_cast<double>(model_attributes["energy_process_spike"]);
    }
    if (model_attributes.find("latency_process_spike") !=
            model_attributes.end())
    {
        default_latency_process_spike =
                static_cast<double>(model_attributes["latency_process_spike"]);
    }
}

void sanafe::PipelineUnit::dendrite_set_default_attributes()
{
    if (model_attributes.find("energy_update") != model_attributes.end())
    {
        default_energy_update =
                static_cast<double>(model_attributes["energy_update"]);
    }
    if (model_attributes.find("latency_update") != model_attributes.end())
    {
        default_latency_update =
                static_cast<double>(model_attributes["latency_update"]);
    }
}

void sanafe::PipelineUnit::soma_set_default_attributes()
{
    // Lambda for finding keys later
    auto key_exists = [this](const std::string &key) {
        return model_attributes.find(key) != model_attributes.end();
    };

    // Soma energy attributes
    const std::set<std::string> energy_metric_names{
            "energy_access_neuron", "energy_update_neuron", "energy_spike_out"};
    const bool parse_energy_metrics = std::any_of(
            energy_metric_names.begin(), energy_metric_names.end(), key_exists);
    if (parse_energy_metrics)
    {
        for (const auto &metric : energy_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Metric not defined: " + metric;
                INFO("Error: %s\n", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaEnergyMetrics energy_metrics;
        energy_metrics.energy_access_neuron =
                static_cast<double>(model_attributes["energy_access_neuron"]);
        energy_metrics.energy_update_neuron =
                static_cast<double>(model_attributes["energy_update_neuron"]);
        energy_metrics.energy_spike_out =
                static_cast<double>(model_attributes["energy_spike_out"]);
        default_soma_energy_metrics = energy_metrics;
    }

    // Soma latency attributes
    const std::set<std::string> latency_metric_names{"latency_access_neuron",
            "latency_update_neuron", "latency_spike_out"};
    const bool parse_latency_metrics = std::any_of(latency_metric_names.begin(),
            latency_metric_names.end(), key_exists);
    if (parse_latency_metrics)
    {
        for (const auto &metric : latency_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Missing metric: " + metric;
                INFO("Error: %s\n", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaLatencyMetrics latency_metrics;
        latency_metrics.latency_access_neuron =
                static_cast<double>(model_attributes["latency_access_neuron"]);
        latency_metrics.latency_update_neuron =
                static_cast<double>(model_attributes["latency_update_neuron"]);
        latency_metrics.latency_spike_out =
                static_cast<double>(model_attributes["latency_spike_out"]);
        default_soma_latency_metrics = latency_metrics;
    }
}

sanafe::BufferPosition sanafe::pipeline_parse_buffer_pos_str(
        const std::string &buffer_pos_str, const bool buffer_inside_unit)
{
    BufferPosition buffer_pos{buffer_before_dendrite_unit};
    // H/w either has SANA-FE insert and manage a time-step buffer on the h/w
    //  unit's inputs (before the unit), or can assume that the unit manages its
    //  own state that is internally buffered over consecutive time-steps e.g.,
    //  using a double buffer

    if (buffer_pos_str == "dendrite")
    {
        if (buffer_inside_unit)
        {
            buffer_pos = buffer_inside_dendrite_unit;
        }
        else
        {
            buffer_pos = buffer_before_dendrite_unit;
        }
    }
    else if (buffer_pos_str == "soma")
    {
        if (buffer_inside_unit)
        {
            buffer_pos = buffer_inside_soma_unit;
        }
        else
        {
            buffer_pos = buffer_before_soma_unit;
        }
    }
    else if (buffer_pos_str == "axon_out")
    {
        buffer_pos = buffer_before_axon_out_unit;
    }
    else
    {
        INFO("Error: Buffer position %s not supported", buffer_pos_str.c_str());
        throw std::invalid_argument("Error: Buffer position not supported");
    }

    return buffer_pos;
}
