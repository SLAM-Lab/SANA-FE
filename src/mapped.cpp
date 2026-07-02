
#include <cstddef>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>

#include "arch.hpp"
#include "attribute.hpp"
#include "core.hpp"
#include "mapped.hpp"
#include "network.hpp"
#include "pipeline.hpp"
#include "print.hpp"

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
sanafe::MappedConnection::MappedConnection(
        MappedNeuron &pre_neuron, MappedNeuron &post_neuron)
        : pre_neuron(&pre_neuron)
        , post_neuron(&post_neuron)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
}

void sanafe::MappedConnection::build_message_processing_pipeline()
{
    MappedNeuron &n = *post_neuron;
    const Core &mapped_core = *(n.core);

    // Cache whether *any* of the post-synaptic neuron's connections requires
    //  forcing updates every time-step. The connections might be mapped to
    //  different synapse units, so if this is true we still need to check each
    //  connection before updating. However, if false, it avoids us iterating
    //  over every connection every time-step
    n.check_for_synapse_updates_every_timestep |=
            synapse_hw->update_every_timestep;

    // We don't support putting the buffer inside or before the synapse unit, so
    //  unconditionally push the synapse h/w. This is because putting the buffer
    //  here could result in a spike being sent before all computation for the
    //  time-step has being completed (i.e. the order of spike arrival can
    //  affect the results)
    message_processing_pipeline.push_back(synapse_hw);
    if ((mapped_core.pipeline_config.buffer_position >
                buffer_before_dendrite_unit) &&
            (n.dendrite_hw != synapse_hw))
    {
        message_processing_pipeline.push_back(n.dendrite_hw);
    }
    if ((mapped_core.pipeline_config.buffer_position >
                buffer_before_soma_unit) &&
            (n.soma_hw != n.dendrite_hw))
    {
        message_processing_pipeline.push_back(n.soma_hw);
    }
}

void sanafe::MappedConnection::set_attributes(
        const std::map<std::string, sanafe::ModelAttribute> &model_attributes)
        const
{
    for (const auto &[key, value] : model_attributes)
    {
        bool attribute_supported = false;
        if (value.forward_to_synapse)
        {
            attribute_supported |= synapse_hw->check_attribute(key);
            synapse_hw->set_attribute_edge(
                    mapped_synapse_hw_address, key, value);
        }
        if (value.forward_to_dendrite)
        {
            MappedNeuron &n = *post_neuron;
            attribute_supported |= n.dendrite_hw->check_attribute(key);
            n.dendrite_hw->set_attribute_edge(
                    mapped_synapse_hw_address, key, value);
        }
        if (!attribute_supported)
        {
            const std::string error("Attribute '" + key +
                    "' not supported by any message processing h/w unit. "
                    "Mapping to h/w failed.");
            INFO("Error: %s\n", error.c_str());
            throw HardwareMappingError(error);
        }
    }
}

sanafe::MappedNeuron::MappedNeuron(const size_t nid,
        const Neuron &neuron_to_map, const size_t mapped_offset_within_core,
        Core *mapped_core, PipelineUnit *mapped_soma,
        AxonOutUnit *mapped_axon_out, PipelineUnit *mapped_dendrite)
        : parent_group_name(neuron_to_map.parent_group_name)
        , offset(neuron_to_map.offset)
        , id(nid)
        , core(mapped_core)
        , dendrite_hw(mapped_dendrite)
        , soma_hw(mapped_soma)
        , axon_out_hw(mapped_axon_out)
        , mapped_offset_within_core(mapped_offset_within_core)
        , mapping_order(neuron_to_map.mapping_order)
        , log_spikes(neuron_to_map.log_spikes)
        , log_potential(neuron_to_map.log_potential)
{
    // The neuron processing pipeline is a sequence of hardware to update the
    //  neuron, the message processing pipeline is built later when mapping
    //  connections
    build_neuron_processing_pipeline();
}

void sanafe::MappedNeuron::set_attributes(
        const std::map<std::string, sanafe::ModelAttribute> &model_attributes,
        std::optional<bool> set_log_spikes)
{
    // Note that neuron attributes related to h/w mapping (i.e.,
    //  default_synapse_hw_name, dendrite_hw_name) are not arguments here.
    //  At this point the neuron has been mapped, and so h/w cannot be changed

    if (set_log_spikes.has_value())
    {
        log_spikes = set_log_spikes.value();
    }
    // Note that setting potential traces (log_potential) is also not supported
    //  in mapped neurons; if and once a simulation has started, changing which
    //  neurons record their potentials would mess up the trace file, which is
    //  represented as a (CSV) table of recorded values

    for (const auto &[key, attribute] : model_attributes)
    {
        if (is_reserved_neuron_attribute(key))
        {
            const std::string error = ("Reserved neuron attribute '" + key +
                    "' cannot be used as a model attribute. " +
                    "Pass it as a direct argument instead (if supported).");
            throw std::invalid_argument(error);
        }

        TRACE2(CHIP, "Forwarding attribute: %s (dendrite:%d soma:%d)\n",
                key.c_str(), attribute.forward_to_dendrite,
                attribute.forward_to_soma);

        bool attribute_supported = false;
        if (attribute.forward_to_dendrite && (dendrite_hw != nullptr))
        {
            attribute_supported |= dendrite_hw->check_attribute(key);
            dendrite_hw->set_attribute_neuron(
                    mapped_dendrite_hw_address, key, attribute);
        }
        if (attribute.forward_to_soma && (soma_hw != nullptr))
        {
            attribute_supported |= soma_hw->check_attribute(key);
            soma_hw->set_attribute_neuron(
                    mapped_soma_hw_address, key, attribute);
        }
        if (!attribute_supported)
        {
            const std::string error("Attribute '" + key +
                    "' not supported by any neuron processing h/w unit. "
                    "Mapping to h/w failed.");
            INFO("Error: %s\n", error.c_str());
            throw HardwareMappingError(error);
        }
    }
}

void sanafe::MappedNeuron::build_neuron_processing_pipeline()
{
    bool dendrite_hw_added = false;
    if (core->pipeline_config.buffer_position < buffer_before_dendrite_unit)
    {
        throw std::runtime_error("Error: Buffer must be after synaptic h/w");
    }
    if (core->pipeline_config.buffer_position <= buffer_inside_dendrite_unit)
    {
        neuron_processing_pipeline.push_back(dendrite_hw);
        dendrite_hw_added = true;
    }
    if (core->pipeline_config.buffer_position <= buffer_inside_soma_unit)
    {
        // Avoid pushing the same h/w unit twice
        if ((soma_hw != dendrite_hw) || !dendrite_hw_added)
        {
            neuron_processing_pipeline.push_back(soma_hw);
        }
    }
}
