// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// message.hpp
#ifndef MESSAGE_HEADER_INCLUDED_
#define MESSAGE_HEADER_INCLUDED_

#include <cstddef>
#include <limits>
#include <string>

#include "fwd.hpp"

namespace sanafe
{

constexpr long int placeholder_mid = -1L; // An invalid message id for placeholders
struct Message
{
    double generation_delay{0.0};
    double receive_delay{0.0};
    double network_delay{0.0};
    // Specific hardware delays
    double min_hop_delay{0.0};
    double blocked_delay{0.0};
    double sent_timestamp{-std::numeric_limits<double>::infinity()};
    double received_timestamp{-std::numeric_limits<double>::infinity()};
    double processed_timestamp{-std::numeric_limits<double>::infinity()};
    long int timestep;
    long int mid;
    size_t spikes{0UL};
    size_t hops{0UL};
    size_t src_neuron_offset;
    std::string src_neuron_group_id;
    size_t src_x{0UL};
    size_t dest_x{0UL};
    size_t src_y{0UL};
    size_t dest_y{0UL};
    size_t src_tile_id;
    size_t src_core_id;
    size_t src_core_offset;
    size_t dest_tile_id{0UL};
    size_t dest_core_id{0UL};
    size_t dest_core_offset{0UL};
    int dest_axon_hw{0};
    size_t dest_axon_id{0UL};
    bool placeholder{true};
    bool in_noc{false};

    explicit Message();
    explicit Message(long int id, const SpikingChip &hw, const MappedNeuron &n,
            long int timestep);
    explicit Message(long int id, const SpikingChip &hw, size_t axon_address,
            const MappedNeuron &n, long int timestep);
    Message(const Message &copy) = default;
    Message(Message &&move) = default;
    ~Message() = default;
    Message& operator=(const Message &copy) = default;
    Message& operator=(Message &&move) = default;
};

class CompareMessagesBySentTime // Used by scheduler
{
public:
    bool operator()(const Message &first, const Message &second) const noexcept;
};

class CompareMessagesByID // Used by chip to sort message traces
{
public:
    bool operator()(const Message &first, const Message &second) const noexcept;
};

}

#endif