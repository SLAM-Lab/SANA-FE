// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <cstddef>
#include <filesystem>
#include <functional> // For std::reference_wrapper
#include <sstream>
#include <stdexcept>
#include <vector>

#include "arch.hpp"
#include "core.hpp"
#include "mapped.hpp"
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
{
    timestep_buffer.resize(pipeline_config.max_neurons_supported);
}

void sanafe::Core::update_hw_in_use()
{
    std::vector<std::reference_wrapper<PipelineUnit>> hw_list;
    for (const auto &hw : pipeline_hw)
    {
        if (hw->is_used)
        {
            hw_list.emplace_back(*hw);
        }
    }

    pipeline_hw_in_use = hw_list;
}

sanafe::PipelineUnit *sanafe::Core::get_hw(const std::string &hw_name,
        const bool is_synapse, const bool is_dendrite, const bool is_soma)
{
    PipelineUnit *hw_unit{nullptr};
    const bool choose_first_available_by_default = (hw_name.empty());

    bool hw_found{false};
    for (auto &hw : pipeline_hw)
    {
        if ((is_synapse && !hw->implements_synapse) ||
                (is_dendrite && !hw->implements_dendrite) ||
                (is_soma && !hw->implements_soma))
        {
            // H/w does not support the required operations, skip
            continue;
        }
        if (choose_first_available_by_default || (hw_name == hw->name))
        {
            hw_unit = hw.get();
            hw_found = true;
            break;
        }
    }

    if (!hw_found)
    {
        INFO("Error: Could not find h/w with name:%s (that implements "
             "synapse:%d, dendrite:%d, soma:%d)\n ",
                hw_name.c_str(), is_synapse, is_dendrite, is_soma);
        throw std::runtime_error("Error: Could not find dendrite h/w");
    }

    return hw_unit;
}

sanafe::PipelineUnit *sanafe::Core::get_synapse_hw(
        const std::string &synapse_hw_name)
{
    return get_hw(synapse_hw_name, true, false, false);
}

sanafe::PipelineUnit *sanafe::Core::get_dendrite_hw(
        const std::string &dendrite_hw_name)
{
    return get_hw(dendrite_hw_name, false, true, false);
}

sanafe::PipelineUnit *sanafe::Core::get_soma_hw(const std::string &soma_hw_name)
{
    return get_hw(soma_hw_name, false, false, true);
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

    // Map neuron model to dendrite and soma h/w units in this core.
    //  Search through all models implemented by this core and return the
    //  one that matches. If no dendrite / soma h/w is specified,
    //  default to the first one defined
    if (pipeline_hw.empty())
    {
        INFO("Error: No pipeline units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No units defined");
    }
    PipelineUnit *mapped_dendrite_hw =
            get_dendrite_hw(neuron_to_map.dendrite_hw_name);
    PipelineUnit *mapped_soma_hw = get_soma_hw(neuron_to_map.soma_hw_name);

    if (axon_out_hw.empty())
    {
        INFO("Error: No axon out units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No axon out units defined");
    }
    AxonOutUnit *mapped_axon_out = axon_out_hw.data();

    // Map the neuron to the core and its h/w units
    const size_t neuron_offset_within_core = neurons.size();
    MappedNeuron &mapped_neuron = neurons.emplace_back(neuron_id, neuron_to_map,
            neuron_offset_within_core, this, mapped_soma_hw, mapped_axon_out,
            mapped_dendrite_hw);

    mapped_neuron.mapped_dendrite_hw_address = mapped_dendrite_hw->add_neuron();
    if (mapped_soma_hw != mapped_dendrite_hw)
    {
        mapped_neuron.mapped_soma_hw_address = mapped_soma_hw->add_neuron();
    }
    else // is a combined dendrite/soma h/w unit
    {
        // For combined units, do not map a neuron to the same h/w unit twice
        mapped_neuron.mapped_soma_hw_address =
                mapped_neuron.mapped_dendrite_hw_address;
    }

    mapped_neuron.set_model_attributes(neuron_to_map.model_attributes);
}

void sanafe::Core::map_connection(const Connection &con,
        MappedNeuron &pre_neuron, MappedNeuron &post_neuron,
        std::string synapse_hw_name)
{
    pre_neuron.connections_out.emplace_back(pre_neuron, post_neuron);
    MappedConnection &mapped_con = pre_neuron.connections_out.back();

    mapped_con.synapse_hw = get_synapse_hw(synapse_hw_name);
    mapped_con.synapse_hw->add_connection(mapped_con);

    mapped_con.build_message_processing_pipeline();
    mapped_con.set_model_attributes(con.synapse_attributes);
    TRACE2(CHIP, "Mapped connection to hw: %s\n",
            mapped_con.synapse_hw->name.c_str());
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
        // Use external plug-in model
        const std::filesystem::path plugin_lib_path =
                config.model_info.plugin_library_path.value();
        TRACE1(CHIP, "Creating unit from plugin: %s.\n",
                plugin_lib_path.c_str());
        pipeline_hw.emplace_back(
                plugin_get_hw(config.model_info.name, plugin_lib_path));
    }
    else
    {
        // Use built in model
        TRACE1(CHIP, "Creating built-in model %s.\n",
                config.model_info.name.c_str());
        pipeline_hw.emplace_back(
                model_get_pipeline_unit(config.model_info.name));
    }

    auto &new_unit = pipeline_hw.back();
    // Forward relevant h/w attributes onto the new pipeline unit, e.g., whether
    //  to log energy/latency or force updates every time-step
    new_unit->set_attributes_hw(config.name, config.model_info);
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
