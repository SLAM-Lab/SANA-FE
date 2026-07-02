// Copyright (c) 2026 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef MAPPED_HEADER_INCLUDED_
#define MAPPED_HEADER_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_set>

#include "attribute.hpp"
#include "fwd.hpp"

namespace sanafe
{
enum NeuronStatus : uint8_t
{
    neuron_state_unset = 0U,
    idle = 1U,
    updated = 2U,
    fired = 3U
};

class HardwareMappingError : public std::runtime_error
{
public:
    // Pass the specific message to the base class constructor
    explicit HardwareMappingError(const std::string &message)
            : std::runtime_error(message)
    {
    }
};

class MappedConnection
{
public:
    // Note raw pointers are required here, since there's a circular dependency
    //  between MappedConnection and MappedNeuron
    MappedNeuron *const pre_neuron;
    MappedNeuron *const post_neuron;
    PipelineUnit *synapse_hw{nullptr};
    std::vector<PipelineUnit *> message_processing_pipeline;
    size_t connection_offset{0UL};
    size_t mapped_dendrite_hw_address{0UL};
    size_t mapped_synapse_hw_address{0UL};

    explicit MappedConnection(MappedNeuron &pre_neuron, MappedNeuron &post_neuron);
    void set_attributes(const std::map<std::string, sanafe::ModelAttribute> &attributes) const;
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

    size_t mapped_offset_within_core{std::numeric_limits<size_t>::max()};
    size_t mapped_dendrite_hw_address{std::numeric_limits<size_t>::max()};
    size_t mapped_soma_hw_address{std::numeric_limits<size_t>::max()};
    size_t mapping_order;
    int spike_count{0};
    int maps_in_count{0};
    int maps_out_count{0};
    NeuronStatus status{neuron_state_unset};

    // Flags and traces
    std::unordered_set<std::string> trace_names{};
    bool log_spikes{false};
    bool log_potential{false};
    bool check_for_synapse_updates_every_timestep{false};

    MappedNeuron(size_t nid, const Neuron &neuron_to_map, size_t mapped_offset_within_core, Core *mapped_core, PipelineUnit *mapped_soma, AxonOutUnit *mapped_axon_out, PipelineUnit *mapped_dendrite);
    MappedNeuron(const MappedNeuron &copy) = default;
    ~MappedNeuron() = default;
    MappedNeuron& operator=(const MappedNeuron& other) = default;
    MappedNeuron(MappedNeuron&& other) = default;

    MappedNeuron& operator=(MappedNeuron&& other) = default;
    void set_attributes(const std::map<std::string, ModelAttribute> &attributes, std::optional<bool> set_log_spikes=std::nullopt);

private:
    void build_neuron_processing_pipeline();
};

}

#endif