// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef MAPPED_HEADER_INCLUDED_
#define MAPPED_HEADER_INCLUDED_

#include <string>
#include <vector>

#include "attribute.hpp"
#include "fwd.hpp"

namespace sanafe
{

enum NeuronStatus : int
{
    INVALID_NEURON_STATE,
    IDLE,
    UPDATED,
    FIRED
};


class MappedConnection
{
public:
    MappedNeuron &pre_neuron;
    MappedNeuron &post_neuron;
    PipelineUnit *synapse_hw{nullptr};
    std::vector<PipelineUnit *> message_processing_pipeline{};
    size_t synapse_address{0UL};

    explicit MappedConnection(MappedNeuron &pre_neuron, MappedNeuron &post_neuron);
    void build_message_processing_pipeline();
};

class MappedNeuron
{
public:
    std::vector<MappedConnection> connections_out;
    std::vector<int> axon_out_addresses;
    std::string parent_group_name;
    size_t offset;
    size_t id;

    // Internal pointers to mapped hardware
    Core *core{nullptr};
    Core *post_synaptic_cores{nullptr};
    PipelineUnit *dendrite_hw{nullptr};
    PipelineUnit *soma_hw{nullptr};
    AxonOutUnit *axon_out_hw{nullptr};
    std::vector<PipelineUnit *> neuron_processing_pipeline{};

    size_t mapped_address{-1ULL};
    size_t mapping_order;
    int spike_count{0};
    int maps_in_count{0};
    int maps_out_count{0};
    NeuronStatus status{INVALID_NEURON_STATE};

    // Flags and traces
    bool force_synapse_update{false};
    bool force_dendrite_update{false};
    bool force_soma_update{false};
    bool log_spikes{false};
    bool log_potential{false};

    // Track spikes
    bool axon_out_input_spike{false};

    MappedNeuron(const Neuron &neuron_to_map, Core *mapped_core, const size_t nid, const size_t mapped_address, PipelineUnit *mapped_dendrite, PipelineUnit *mapped_soma, AxonOutUnit *mapped_axon_out);
    MappedNeuron(const MappedNeuron &copy) = default;
    MappedNeuron& operator=(const MappedNeuron& other) = default;
    MappedNeuron(MappedNeuron&& other) = default;
    MappedNeuron& operator=(MappedNeuron&& other) = default;
    void set_model_attributes(const std::map<std::string, ModelAttribute> &model_attributes);

private:
    void build_neuron_processing_pipeline();
};

}

#endif