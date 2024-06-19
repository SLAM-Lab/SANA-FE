// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// schedule.hpp: Schedule global order of messages on a neuromorphic chip
#ifndef SCHEDULE_HEADER_INCLUDED_
#define SCHEDULE_HEADER_INCLUDED_

#include <functional> // For std::reference_wrapper
#include <queue>
#include <vector>

namespace sanafe
{

struct Message;
struct RouterMessageStats;

enum Direction
{
    // Router link direction
    north = 0,
    east,
    south,
    west,
    ndirections,
};

class CompareMessages
{
public:
    bool operator()(const Message &first, const Message &second) const { return first.sent_timestamp > second.sent_timestamp; }
};

using MessageRef = std::reference_wrapper<Message>;
using MessagePriorityQueue = std::priority_queue<MessageRef, std::vector<MessageRef>, CompareMessages>;
using MessageFifo = std::list<MessageRef>;

struct Scheduler
{
    int noc_width, noc_height, buffer_size, core_count, max_cores_per_tile;
};

// NocInfo is used by the scheduler to track the high-level state of the NoC
//  at a given time
struct NocInfo
{
    std::vector<MessageFifo> messages_received;
    const int noc_width, noc_height, core_count;
    const size_t max_cores_per_tile;
    // Message density is the distribution of messages buffered in different
    //  links across the NoC. This vector is flattened row-major order
    //  where indexes are (x, y, links) and the links idx change fastest
    std::vector<double> message_density;
    // The time that cores finish receiving their last message, i.e. the
    //  time at which they become free again to process more messages
    std::vector<double> core_finished_receiving;
    double mean_in_flight_receive_delay;
    long int messages_in_noc;

    NocInfo(const int width, const int height, const int core_count, const size_t max_cores_per_tile);
    int idx(const int x, const int y, const int link)
    {
        const int links_per_router = max_cores_per_tile + ndirections;
        return (x * noc_height * links_per_router) + (y * links_per_router) +
                link;
    }
};

MessagePriorityQueue schedule_init_timing_priority(std::vector<MessageFifo> &message_queues_per_core);
double schedule_messages(std::vector<std::list<Message>> &messages_sent_per_core, const Scheduler &scheduler);
void schedule_update_noc_message_counts(const Message &m, NocInfo &noc, const bool message_in);
double schedule_calculate_messages_along_route(const Message &m, NocInfo &noc);
void schedule_update_noc(const double t, NocInfo &noc);

}

#endif