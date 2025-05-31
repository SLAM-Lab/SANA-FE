// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <cstddef>
#include <filesystem>
#include <sstream>
#include <stdexcept>

#include "arch.hpp"
#include "core.hpp"
#include "models.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "plugins.hpp"
#include "print.hpp"

sanafe::AxonInUnit::AxonInUnit(const AxonInConfiguration &config)
        : name(config.name)
        , energy_spike_message(config.metrics.energy_message_in)
        , latency_spike_message(config.metrics.latency_message_in)
{
}

sanafe::AxonOutUnit::AxonOutUnit(const AxonOutConfiguration &config)
        : name(config.name)
        , energy_access(config.metrics.energy_message_out)
        , latency_access(config.metrics.latency_message_out)
{
}

sanafe::Core::Core(const CoreConfiguration &config)
        : pipeline_config(config.pipeline)
        , name(config.name)
        , id(config.address.id)
        , offset(config.address.offset_within_tile)
        , parent_tile_id(config.address.parent_tile_id)
        , log_energy(config.pipeline.log_energy)
        , log_latency(config.pipeline.log_latency)

{
    timestep_buffer.resize(pipeline_config.max_neurons_supported);
}

sanafe::PipelineUnit *sanafe::Core::map_neuron_to_dendrite(
        const Neuron &neuron_to_map)
{
    PipelineUnit *mapped_dendrite{nullptr};

    const bool choose_first_dendrite_by_default =
            (neuron_to_map.dendrite_hw_name.empty());
    bool dendrite_found = false;
    for (auto &hw : pipeline_hw)
    {
        if (hw->implements_dendrite &&
                (choose_first_dendrite_by_default ||
                        neuron_to_map.dendrite_hw_name == hw->name))
        {
            mapped_dendrite = hw.get();
            dendrite_found = true;
            break;
        }
    }
    if (!dendrite_found)
    {
        INFO("Error: Could not map neuron nid:%zu (hw:%s) "
             "to any dendrite h/w.\n",
                neuron_to_map.offset, neuron_to_map.dendrite_hw_name.c_str());
        throw std::runtime_error("Error: Could not map neuron to dendrite h/w");
    }

    return mapped_dendrite;
}

sanafe::PipelineUnit *sanafe::Core::map_neuron_to_soma(
        const Neuron &neuron_to_map)
{
    PipelineUnit *mapped_soma{nullptr};
    const bool choose_first_soma_by_default =
            (neuron_to_map.soma_hw_name.empty());
    bool soma_found = false;
    for (auto &hw : pipeline_hw)
    {
        if (hw->implements_soma &&
                (choose_first_soma_by_default ||
                        neuron_to_map.soma_hw_name == hw->name))
        {
            mapped_soma = hw.get();
            soma_found = true;
            break;
        }
    }
    if (!soma_found)
    {
        INFO("Error: Could not map neuron nid:%zu (hw:%s) "
             "to any soma h/w.\n",
                neuron_to_map.offset, neuron_to_map.soma_hw_name.c_str());
        throw std::runtime_error("Error: Could not map neuron to soma h/w");
    }
    mapped_soma->neuron_count++;

    return mapped_soma;
}

void sanafe::Core::map_neuron(
        const Neuron &neuron_to_map, const size_t neuron_id)
{
    TRACE1(CHIP, "Mapping nid:%s.%zu to core: %zu\n",
            neuron_to_map.parent_group_name.c_str(), neuron_to_map.offset, id);

    if (neurons.size() >= pipeline_config.max_neurons_supported)
    {
        INFO("Error: Exceeded maximum neurons per core (%zu)",
                pipeline_config.max_neurons_supported);
        throw std::runtime_error("Error: Exceeded maximum neurons per core.");
    }

    // Map neuron model to dendrite and soma hardware units in this core.
    //  Search through all models implemented by this core and return the
    //  one that matches. If no dendrite / soma hardware is specified,
    //  default to the first one defined
    if (pipeline_hw.empty())
    {
        INFO("Error: No pipeline units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No units defined");
    }
    PipelineUnit *mapped_dendrite = map_neuron_to_dendrite(neuron_to_map);
    PipelineUnit *mapped_soma = map_neuron_to_soma(neuron_to_map);

    if (axon_out_hw.empty())
    {
        INFO("Error: No axon out units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No axon out units defined");
    }
    AxonOutUnit *mapped_axon_out = axon_out_hw.data();

    // Map the neuron to the core and its hardware units
    const size_t address = neurons.size();
    neurons.emplace_back(neuron_to_map, neuron_id, this,
            mapped_soma, address, mapped_axon_out, mapped_dendrite);
}

sanafe::AxonInUnit &sanafe::Core::create_axon_in(
        const AxonInConfiguration &config)
{
    axon_in_hw.emplace_back(config);
    TRACE1(CHIP, "New axon in h/w unit created (%zu.%zu)\n", parent_tile_id,
            id);

    return axon_in_hw.back();
}

sanafe::PipelineUnit &sanafe::Core::create_pipeline_unit(
        const PipelineUnitConfiguration &config)
{
    // Create the synapse model
    if (config.model_info.plugin_library_path.has_value())
    {
        const std::filesystem::path plugin_lib_path =
                config.model_info.plugin_library_path.value();
        TRACE1(CHIP, "Creating unit from plugin: %s.\n",
                plugin_lib_path.c_str());
        pipeline_hw.emplace_back(
                plugin_get_hw(config.model_info.name, plugin_lib_path));
    }
    else
    {
        // Use built in models
        TRACE1(CHIP, "Creating built-in model %s.\n",
                config.model_info.name.c_str());
        pipeline_hw.emplace_back(
                model_get_pipeline_unit(config.model_info.name));
    }

    auto &new_unit = pipeline_hw.back();
    // Forward all attributes onto the new h/w unit
    new_unit->set_attributes(config.name, config.model_info);
    // Set the input/output interface of the pipeline unit and in doing so we
    //  configure which functionality the h/w unit supports
    new_unit->check_implemented(config.implements_synapse,
            config.implements_dendrite, config.implements_soma);
    TRACE1(CHIP, "New h/w unit created (%s) in core:%zu\n", config.name.c_str(),
            id);

    return *new_unit;
}

sanafe::AxonOutUnit &sanafe::Core::create_axon_out(
        const AxonOutConfiguration &config)
{
    axon_out_hw.emplace_back(config);
    TRACE1(CHIP, "New axon out h/w unit created: (%zu.%zu)\n", parent_tile_id,
            id);

    return axon_out_hw.back();
}

std::string sanafe::Core::info() const noexcept
{
    std::ostringstream ss;
    ss << "sanafe::Core(name= " << name << " tile=" << parent_tile_id << ")";
    return ss.str();
}
