// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include "core.hpp"
#include "models.hpp"
#include "network.hpp"
#include "plugins.hpp"

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
    PipelineUnit *mapped_dendrite;

    bool choose_first_dendrite_by_default =
            (neuron_to_map.dendrite_hw_name.length() == 0);
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

    PipelineUnit *mapped_soma;
    bool choose_first_soma_by_default =
            (neuron_to_map.soma_hw_name.length() == 0);
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

    if (axon_out_hw.empty())
    {
        INFO("Error: No axon out units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No axon out units defined");
    }
    AxonOutUnit *mapped_axon_out = &(axon_out_hw[0]);

    // Map the neuron to the core and its hardware units
    const size_t address = neurons.size();
    neurons.emplace_back(neuron_to_map, this, neuron_id, address,
            mapped_dendrite, mapped_soma, mapped_axon_out);

    return;
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
    TRACE1(CHIP, "implements synapse:%d dendrite:%d soma:%d\n",
            new_unit->implements_synapse, new_unit->implements_dendrite,
            new_unit->implements_soma);
    new_unit->set_attributes(config.name, config.model_info);
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
