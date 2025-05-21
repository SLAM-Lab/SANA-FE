// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef MESSAGE_HEADER_INCLUDED_
#define MESSAGE_HEADER_INCLUDED_

#include <limits>
#include <string>

namespace sanafe
{

struct MappedNeuron;
class SpikingChip;

constexpr long int placeholder_mid = -1L; // An invalid message id for placeholders
struct Message
{
    double generation_delay{0.0};
    double network_delay{0.0};
    double receive_delay{0.0};
    double blocked_delay{0.0};
    double sent_timestamp{-std::numeric_limits<double>::infinity()};
    double received_timestamp{-std::numeric_limits<double>::infinity()};
    double processed_timestamp{-std::numeric_limits<double>::infinity()};
    long int timestep;
    long int mid;
    int spikes{0};
    size_t hops{0UL};
    size_t src_neuron_id;
    std::string src_neuron_group_id;
    int src_x;
    int dest_x{0};
    int src_y;
    int dest_y{0};
    int src_tile_id;
    int src_core_id;
    int src_core_offset;
    int dest_tile_id{0};
    int dest_core_id{0};
    int dest_core_offset{0};
    int dest_axon_hw{0};
    int dest_axon_id{0};
    bool placeholder{true};
    bool in_noc{false};

    explicit Message(const long int id, const SpikingChip &hw, const MappedNeuron &n, long int timestep);
    explicit Message(const long int id, const SpikingChip &hw, const MappedNeuron &n, long int timestep, int axon_address);
    Message(const Message &copy) = default;
};

}

#endif