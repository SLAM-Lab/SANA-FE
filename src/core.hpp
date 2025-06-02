// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// tile.hpp
#ifndef CORE_HEADER_INCLUDED_
#define CORE_HEADER_INCLUDED_

#include <cstddef>
#include <functional>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "arch.hpp"
#include "fwd.hpp"
#include "mapped.hpp"
#include "message.hpp"
#include "pipeline.hpp"

namespace sanafe
{
struct AxonInModel
{
    // List of all neuron connections to send spikes to
    std::vector<size_t> synapse_addresses;
    Message *message{nullptr};
    int spikes_received{0};
    int active_synapses{0};
};

struct AxonOutModel
{
    size_t dest_axon_id{};
    size_t dest_tile_id{};
    size_t dest_core_offset{};
    size_t src_neuron_offset{};
};

class AxonInUnit
{
public:
    std::string name;
    long int spike_messages_in{0L};
    double energy{0.0};
    double time{0.0};
    double energy_spike_message;
    double latency_spike_message;

    explicit AxonInUnit(const AxonInConfiguration &config);
};

class AxonOutUnit
{
public:
    // The axon output points to a number of axons, stored at the
    //  post-synaptic core. A neuron can point to a number of these
    std::string name;
    long int packets_out{0L};
    double energy{0.0};
    double time{0.0};
    double energy_access;
    double latency_access;

    explicit AxonOutUnit(const AxonOutConfiguration &config);
};

class Core
{
public:
    std::vector<AxonInUnit> axon_in_hw;
    std::vector<std::shared_ptr<PipelineUnit>> pipeline_hw;
    std::vector<AxonOutUnit> axon_out_hw;

    std::vector<std::reference_wrapper<Message>> messages_in;
    std::vector<AxonInModel> axons_in;
    std::vector<MappedNeuron> neurons;
    std::vector<MappedConnection *> connections_in;
    std::vector<AxonOutModel> axons_out;
    std::vector<PipelineResult> timestep_buffer;

    std::list<BufferPosition> neuron_processing_units;
    std::list<BufferPosition> message_processing_units;
    CorePipelineConfiguration pipeline_config{};
    std::string name;
    double energy{0.0};
    double next_message_generation_delay{0.0};
    size_t id;
    size_t offset;
    size_t parent_tile_id;
    int message_count{0};
    bool log_energy{false};
    bool log_latency{false};

    explicit Core(const CoreConfiguration &config);
    void map_neuron(const Neuron &n, size_t neuron_id);
    AxonInUnit &create_axon_in(const AxonInConfiguration &config);
    PipelineUnit &create_pipeline_unit(const PipelineUnitConfiguration &config);
    AxonOutUnit &create_axon_out(const AxonOutConfiguration &config);
    [[nodiscard]] size_t get_id() const { return id; }
    [[nodiscard]] size_t get_offset() const { return offset; }
    [[nodiscard]] std::string info() const noexcept;

private:
    PipelineUnit *map_neuron_to_dendrite(const Neuron &neuron_to_map);
    PipelineUnit *map_neuron_to_soma(const Neuron &neuron_to_map);
};

}

#endif