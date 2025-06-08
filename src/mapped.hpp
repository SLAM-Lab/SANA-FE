// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef MAPPED_HEADER_INCLUDED_
#define MAPPED_HEADER_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "attribute.hpp"
#include "fwd.hpp"

namespace sanafe
{

enum NeuronStatus : uint8_t
{
    invalid_neuron_state = 0U,
    idle = 1U,
    updated = 2U,
    fired = 3U
};

class MappedConnection
{
public:
    std::reference_wrapper<MappedNeuron> pre_neuron_ref;
    std::reference_wrapper<MappedNeuron> post_neuron_ref;
    PipelineUnit *synapse_hw{nullptr};
    std::vector<PipelineUnit *> message_processing_pipeline;
    size_t synapse_address{0UL};

    explicit MappedConnection(std::reference_wrapper<MappedNeuron> pre_neuron, std::reference_wrapper<MappedNeuron> post_neuron);
    void build_message_processing_pipeline();
};

class MappedNeuron
{
public:
    std::vector<MappedConnection> connections_out;
    std::vector<size_t> axon_out_addresses;
    std::string parent_group_name;
    size_t offset;
    size_t id;

    // Internal pointers to mapped hardware
    Core *core{nullptr};
    Core *post_synaptic_cores{nullptr};
    PipelineUnit *dendrite_hw{nullptr};
    PipelineUnit *soma_hw{nullptr};
    AxonOutUnit *axon_out_hw{nullptr};
    std::vector<PipelineUnit *> neuron_processing_pipeline;

    size_t mapped_address{-1ULL};
    size_t mapping_order;
    int spike_count{0};
    int maps_in_count{0};
    int maps_out_count{0};
    NeuronStatus status{invalid_neuron_state};

    // Flags and traces
    bool force_synapse_update{false};
    bool force_dendrite_update{false};
    bool force_soma_update{false};
    bool log_spikes{false};
    bool log_potential{false};

    MappedNeuron(const Neuron &neuron_to_map, size_t nid, Core *mapped_core, PipelineUnit *mapped_soma, size_t mapped_address, AxonOutUnit *mapped_axon_out, PipelineUnit *mapped_dendrite);
    MappedNeuron(const MappedNeuron &copy) = default;
    ~MappedNeuron() = default;
    MappedNeuron& operator=(const MappedNeuron& other) = default;
    MappedNeuron(MappedNeuron&& other) = default;
    MappedNeuron& operator=(MappedNeuron&& other) = default;
    void set_model_attributes(const std::map<std::string, ModelAttribute> &model_attributes) const;

private:
    void build_neuron_processing_pipeline();
};

}

#endif