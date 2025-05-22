#include <map>

#include "network.hpp"
#include "mapped.hpp"

sanafe::MappedConnection::MappedConnection(
        std::reference_wrapper<MappedNeuron> pre_neuron,
        std::reference_wrapper<MappedNeuron> post_neuron)
        : pre_neuron(pre_neuron)
        , post_neuron(post_neuron)
        , synapse_hw(nullptr)
{
}

void sanafe::MappedConnection::build_message_processing_pipeline()
{
    MappedNeuron &n = post_neuron;
    Core &mapped_core = *(n.core);

    // We don't support putting the buffer inside or before the synapse unit, so
    //  unconditionally push the synapse h/w. This is because putting the buffer
    //  here could cause a spike sent that shouldn't be
    message_processing_pipeline.push_back(synapse_hw);
    if ((mapped_core.pipeline_config.buffer_position >
                BUFFER_BEFORE_DENDRITE_UNIT) &&
            (n.dendrite_hw != synapse_hw))
    {
        message_processing_pipeline.push_back(n.dendrite_hw);
    }
    if ((mapped_core.pipeline_config.buffer_position >
                BUFFER_BEFORE_SOMA_UNIT) &&
            (n.soma_hw != n.dendrite_hw))
    {
        message_processing_pipeline.push_back(n.soma_hw);
    }

    return;
}

sanafe::MappedNeuron::MappedNeuron(const Neuron &neuron_to_map,
        Core *mapped_core, const size_t nid, const size_t mapped_address,
        PipelineUnit *mapped_dendrite, PipelineUnit *mapped_soma,
        AxonOutUnit *mapped_axon_out)
        : parent_group_name(neuron_to_map.parent_group_name)
        , offset(neuron_to_map.offset)
        , id(nid)
        , core(mapped_core)
        , dendrite_hw(mapped_dendrite)
        , soma_hw(mapped_soma)
        , axon_out_hw(mapped_axon_out)
        , mapped_address(mapped_address)
        , mapping_order(neuron_to_map.mapping_order)
        , force_synapse_update(neuron_to_map.force_synapse_update)
        , force_dendrite_update(neuron_to_map.force_dendrite_update)
        , force_soma_update(neuron_to_map.force_soma_update)
        , log_spikes(neuron_to_map.log_spikes)
        , log_potential(neuron_to_map.log_potential)

{
    set_model_attributes(neuron_to_map.model_attributes);
    build_neuron_processing_pipeline();
}

void sanafe::MappedNeuron::set_model_attributes(
        const std::map<std::string, sanafe::ModelAttribute> &model_attributes)
{
    for (auto &[key, attribute] : model_attributes)
    {
        TRACE2(CHIP, "Forwarding attribute: %s (dendrite:%d soma:%d)\n",
                key.c_str(), attribute.forward_to_dendrite,
                attribute.forward_to_soma);
        if (attribute.forward_to_dendrite && (dendrite_hw != nullptr))
        {
            dendrite_hw->check_attribute(key);
            dendrite_hw->set_attribute_neuron(mapped_address, key, attribute);
        }
        if (attribute.forward_to_soma && (soma_hw != nullptr))
        {
            soma_hw->check_attribute(key);
            soma_hw->set_attribute_neuron(mapped_address, key, attribute);
        }
    }
}

void sanafe::MappedNeuron::build_neuron_processing_pipeline()
{
    if (core->pipeline_config.buffer_position < BUFFER_BEFORE_DENDRITE_UNIT)
    {
        throw std::runtime_error("Error: Buffer must be after synaptic h/w");
    }
    if (core->pipeline_config.buffer_position == BUFFER_INSIDE_DENDRITE_UNIT)
    {
        neuron_processing_pipeline.push_back(dendrite_hw);
    }
    if ((core->pipeline_config.buffer_position <= BUFFER_INSIDE_SOMA_UNIT) &&
            (soma_hw != dendrite_hw))
    {
        neuron_processing_pipeline.push_back(soma_hw);
    }
}
