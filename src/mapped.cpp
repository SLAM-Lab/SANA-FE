
#include <cstddef>
#include <functional>
#include <map>
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
        std::reference_wrapper<MappedNeuron> pre_neuron,
        std::reference_wrapper<MappedNeuron> post_neuron)
        : pre_neuron_ref(pre_neuron)
        , post_neuron_ref(post_neuron)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
}

void sanafe::MappedConnection::build_message_processing_pipeline()
{
    MappedNeuron &n = post_neuron_ref;
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

void sanafe::MappedConnection::set_model_attributes(
        const std::map<std::string, sanafe::ModelAttribute> &model_attributes)
        const
{
    for (const auto &[key, value] : model_attributes)
    {
        if (value.forward_to_synapse)
        {
            synapse_hw->check_attribute(key);
            synapse_hw->set_attribute_edge(
                    mapped_synapse_hw_address, key, value);
        }
        if (value.forward_to_dendrite)
        {
            MappedNeuron &n = post_neuron_ref;
            n.dendrite_hw->check_attribute(key);
            n.dendrite_hw->set_attribute_edge(
                    mapped_synapse_hw_address, key, value);
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

void sanafe::MappedNeuron::set_model_attributes(
        const std::map<std::string, sanafe::ModelAttribute> &model_attributes)
        const
{
    for (const auto &[key, attribute] : model_attributes)
    {
        TRACE2(CHIP, "Forwarding attribute: %s (dendrite:%d soma:%d)\n",
                key.c_str(), attribute.forward_to_dendrite,
                attribute.forward_to_soma);
        if (attribute.forward_to_dendrite && (dendrite_hw != nullptr))
        {
            dendrite_hw->check_attribute(key);
            dendrite_hw->set_attribute_neuron(
                    mapped_dendrite_hw_address, key, attribute);
        }
        if (attribute.forward_to_soma && (soma_hw != nullptr))
        {
            soma_hw->check_attribute(key);
            soma_hw->set_attribute_neuron(
                    mapped_soma_hw_address, key, attribute);
        }
    }
}

void sanafe::MappedNeuron::build_neuron_processing_pipeline()
{
    if (core->pipeline_config.buffer_position < buffer_before_dendrite_unit)
    {
        throw std::runtime_error("Error: Buffer must be after synaptic h/w");
    }
    if (core->pipeline_config.buffer_position == buffer_inside_dendrite_unit)
    {
        neuron_processing_pipeline.push_back(dendrite_hw);
    }
    if ((core->pipeline_config.buffer_position <= buffer_inside_soma_unit) &&
            (soma_hw != dendrite_hw))
    {
        neuron_processing_pipeline.push_back(soma_hw);
    }
}
