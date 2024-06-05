
// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// schedule.hpp: Schedule global order of messages on a neuromorphic chip
//  The schedule then determines on-chip timing and predicts run-time.
//  The scheduler maintains a priority queue of messages and accounts for
//  message generation delays, receiving delays and network delays. Generation
//  and receive delays are calculated earlier in sim.cpp. To calculate
//  network delays accurately, keep track of predicted NoC state by counting
//  the numbers of messages that are travelling through routes / flows at
//  different times
#include <cassert>
#include <cmath>
#include <functional> // For std::reference_wrapper
#include <iterator>
#include <vector>

#include "arch.hpp"
#include "network.hpp"
#include "print.hpp"
#include "schedule.hpp"
#include "sim.hpp" // For sim_message_fifo

sanafe::NocInfo::NocInfo(const int width, const int height,
        const int total_cores, const size_t cores_per_tile)
        : noc_width(width)
        , noc_height(height)
        , core_count(total_cores)
        , max_cores_per_tile(cores_per_tile)
{
    messages_in_noc = 0UL;
    mean_in_flight_receive_delay = 0.0;
    return;
}

double sanafe::schedule_messages(
        std::vector<std::list<Message>> &messages, const Scheduler &scheduler)
{
    // Schedule the global order of messages. Takes a vector containing
    //  a list of messages per core, and scheduler parameters (mostly
    //  NoC configuration parameters). Returns the timestamp of the last
    //  scheduled event, i.e., the total time-step delay
    MessagePriorityQueue priority;
    NocInfo noc(scheduler.noc_width, scheduler.noc_height, scheduler.core_count,
            scheduler.max_cores_per_tile);
    double last_timestamp;
    const int total_links = noc.noc_height * noc.noc_width *
            (sanafe::ndirections + noc.max_cores_per_tile);
    noc.message_density = std::vector<double>(total_links);

    std::vector<MessageFifo> messages_sent_per_core(noc.core_count);
    for (size_t core = 0; core < messages.size(); core++)
    {
        auto &q = messages[core];
        for (Message &m : q)
        {
            messages_sent_per_core[core].push_back(m);
        }
    }

    noc.messages_received = std::vector<MessageFifo>(noc.core_count);
    noc.core_finished_receiving = std::vector<double>(noc.core_count);

    priority = schedule_init_timing_priority(messages_sent_per_core);
    last_timestamp = 0.0;
    SIM_TRACE1("Scheduling global order of messages.\n");

    // Each core has a queue of received messages. A structure tracks how
    //  many in-flight messages are in the NoC and occupy each tile. We
    //  track the number of messages passing through each tile at
    //  the point of sending, and the average processing delay of
    //  all of those messages. When a message is added or removed from the
    //  NoC we update the average counts.
    while (priority.size() > 0)
    {
        // Get the core's queue with the earliest simulation time
        Message &m = priority.top();
        priority.pop();
        last_timestamp = fmax(last_timestamp, m.sent_timestamp);

        // Update the Network-on-Chip state
        schedule_update_noc(m.sent_timestamp, noc);

        // Messages without a destination (neuron) are dummy messages.
        //  Dummy messages account for processing time that does not
        //  result in any spike messages. Otherwise, messages are sent
        //  from a src neuron to a dest neuron
        if (!m.placeholder)
        {
            TRACE1("Processing message for nid:%d.%d\n",
                    m.src_neuron->parent_group_id, m.src_neuron->id);
            TRACE1("Send delay:%e\n", m.generation_delay);
            TRACE1("Receive delay:%e\n", m.receive_delay);
            const int dest_core = m.dest_core_id;
            // Figure out if we are able to send a message into the
            //  network i.e., is the route to the dest core
            //  saturated and likely to block? Sum along the route
            //  and see the density of messages along all links.
            const double messages_along_route =
                    schedule_calculate_messages_along_route(m, noc);

            const int path_capacity = (m.hops + 1) * scheduler.buffer_size;
            if (messages_along_route > path_capacity)
            {
                m.sent_timestamp += (messages_along_route - path_capacity) *
                        noc.mean_in_flight_receive_delay;
            }

            // Now, push the message into the right receiving queue
            //  Calculate the network delay and when the message
            //  is received
            m.in_noc = true;
            noc.messages_received[dest_core].push_back(m);

            // Update the rolling average for message
            //  receiving times in-flight in the network
            schedule_update_noc_message_counts(m, noc, true);

            double network_delay = messages_along_route *
                    noc.mean_in_flight_receive_delay / (m.hops + 1.0);
            TRACE1("Path capacity:%d messages:%lf delay:%e\n", path_capacity,
                    messages_along_route, network_delay);

            double earliest_received_time =
                    m.sent_timestamp + fmax(m.network_delay, network_delay);
            m.received_timestamp = fmax(noc.core_finished_receiving[dest_core],
                    earliest_received_time);
            noc.core_finished_receiving[dest_core] = fmax(
                    (noc.core_finished_receiving[dest_core] + m.receive_delay),
                    (earliest_received_time + m.receive_delay));
            m.processed_timestamp = noc.core_finished_receiving[dest_core];
            last_timestamp = fmax(last_timestamp, m.processed_timestamp);
        }

        // Get the next message for this core
        const size_t src_core = m.src_core_id;
        if (messages_sent_per_core[src_core].size() > 0)
        {
            auto &q = messages_sent_per_core[src_core];
            Message &next_message = q.front();
            // If applicable, schedule this next message immediately
            //  after the current message finishes sending
            next_message.sent_timestamp =
                    m.sent_timestamp + next_message.generation_delay;
            last_timestamp = fmax(last_timestamp, next_message.sent_timestamp);
            priority.push(next_message);
            q.pop_front();
        }
        else
        {
            TRACE1("\tCore finished simulating\n");
        }

#ifdef DEBUG
        auto view = priority;
        INFO("***\n");
        while (view.size())
        {
            Message &m = view.top();
            INFO("m:%e\n", m.sent_timestamp);
            view.pop();
        }
        INFO("***\n");
#endif

        TRACE1("Priority size:%lu\n", priority.size());
    }
    TRACE1("Scheduler finished.\n");

    return last_timestamp;
}

void sanafe::schedule_update_noc_message_counts(
        const Message &m, NocInfo &noc, const bool message_in)
{
    // Update the tracked state of the NoC, accounting for a single message
    //  either entering or leaving the NoC (message_in)
    // Go along x path, then y path (dimension order routing), and increment
    //  or decrement counter depending on whether a message is coming in or
    //  out
    int x_increment, y_increment;
    // Adjust by dividing by the total number of links along the path, also
    //  including the output link at the sending core and input link at the
    //  receiving core, i.e. the hops plus 2. The total sum of the added
    //  densities along the path should equal one for one new message.
    double adjust = (1.0 / (2.0 + m.hops));

    if (!message_in)
    {
        adjust *= -1.0;
    }

    if (m.src_x < m.dest_x)
    {
        x_increment = 1;
    }
    else
    {
        x_increment = -1;
    }
    if (m.src_y < m.dest_y)
    {
        y_increment = 1;
    }
    else
    {
        y_increment = -1;
    }
    assert(m.src_core_offset == m.src_core_offset);
    int prev_direction = sanafe::ndirections + (m.src_core_offset);
    for (int x = m.src_x; x != m.dest_x; x += x_increment)
    {
        int direction;
        if (x_increment > 0)
        {
            direction = sanafe::east;
        }
        else
        {
            direction = sanafe::west;
        }
        if (x == m.src_x)
        {
            assert(m.src_core_offset == m.src_core_offset);
            const int link = sanafe::ndirections + (m.src_core_offset);
            noc.message_density[noc.idx(x, m.src_y, link)] += adjust;
        }
        else
        {
            noc.message_density[noc.idx(x, m.src_y, direction)] += adjust;
        }
        prev_direction = direction;
    }
    for (int y = m.src_y; y != m.dest_y; y += y_increment)
    {
        int direction;
        if (y_increment > 0)
        {
            direction = sanafe::north;
        }
        else
        {
            direction = sanafe::south;
        }
        if ((m.src_x == m.dest_x) && (y == m.src_y))
        {
            const int link = sanafe::ndirections + m.src_core_offset;
            noc.message_density[noc.idx(m.dest_x, y, link)] += adjust;
        }
        else
        {
            noc.message_density[noc.idx(m.dest_x, y, prev_direction)] += adjust;
        }

        prev_direction = direction;
    }

    if ((m.src_x == m.dest_x) && (m.src_y == m.dest_y))
    {
        assert(m.src_core_offset == m.src_core_offset);
        const int link = sanafe::ndirections + (m.src_core_offset);
        noc.message_density[noc.idx(m.dest_x, m.dest_y, link)] += adjust;
    }
    else
    {
        noc.message_density[noc.idx(m.dest_x, m.dest_y, prev_direction)] +=
                adjust;
    }

    // Update rolling averages and message counts
    if (message_in)
    {
        // Message entering NoC
        noc.mean_in_flight_receive_delay +=
                (m.receive_delay - noc.mean_in_flight_receive_delay) /
                (static_cast<double>(noc.messages_in_noc) + 1);
        noc.messages_in_noc++;
    }
    else
    {
        // Message leaving the NoC
        if ((noc.messages_in_noc) > 1)
        {
            noc.mean_in_flight_receive_delay +=
                    (noc.mean_in_flight_receive_delay - m.receive_delay) /
                    (static_cast<double>(noc.messages_in_noc) - 1);
        }
        else
        {
            noc.mean_in_flight_receive_delay = 0.0;
        }

        noc.messages_in_noc--;
    }

    return;
}

double sanafe::schedule_calculate_messages_along_route(Message &m, NocInfo &noc)
{
    // Calculate the total flow density along a spike message route.
    //  Sum the densities for all links the message will travel
    double flow_density;
    int x_increment, y_increment;
    int direction;

    flow_density = 0.0;
    if (m.src_x < m.dest_x)
    {
        x_increment = 1;
    }
    else
    {
        x_increment = -1;
    }
    if (m.src_y < m.dest_y)
    {
        y_increment = 1;
    }
    else
    {
        y_increment = -1;
    }
    int prev_direction = sanafe::ndirections + (m.src_core_offset);
    for (int x = m.src_x; x != m.dest_x; x += x_increment)
    {
        int direction = 0;
        if (x_increment > 0)
        {
            direction = sanafe::east;
        }
        else
        {
            direction = sanafe::west;
        }
        if (x == m.src_x)
        {
            const int link = sanafe::ndirections + m.src_core_offset;
            flow_density += noc.message_density[noc.idx(x, m.src_y, link)];
        }
        else
        {
            flow_density += noc.message_density[noc.idx(x, m.src_y, direction)];
        }
        prev_direction = direction;
    }

    for (int y = m.src_y; y != m.dest_y; y += y_increment)
    {
        if (y_increment > 0)
        {
            direction = sanafe::north;
        }
        else
        {
            direction = sanafe::south;
        }
        if (m.src_x == m.dest_x && y == m.src_y)
        {
            const int link = sanafe::ndirections + m.src_core_offset;
            flow_density += noc.message_density[noc.idx(m.dest_x, y, link)];
        }
        else
        {
            flow_density +=
                    noc.message_density[noc.idx(m.dest_x, y, prev_direction)];
        }
        prev_direction = direction;
    }
    // Handle the last (destination) tile
    if ((m.src_x == m.dest_x) && (m.src_y == m.dest_y))
    {
        const int link = sanafe::ndirections + m.src_core_offset;
        flow_density += noc.message_density[noc.idx(m.dest_x, m.dest_y, link)];
    }
    else
    {
        flow_density +=
                noc.message_density[noc.idx(m.dest_x, m.dest_y, prev_direction)];
    }

    assert(flow_density >= -0.1);
    return flow_density;
}

void sanafe::schedule_update_noc(const double t, NocInfo &noc)
{
    // Update the tracked state of the NoC at a given time, t
    for (auto &q : noc.messages_received)
    {
        // Go through all messages in the NoC and check to see if that
        //  message has been fully received by time t. If so, remove it
        //  from the NoC.
        auto it = q.begin();
        while (it != q.end())
        {
            Message &m = (*it);
            bool remove_message = false;
            auto curr = it;
            if (m.in_noc && (t >= m.received_timestamp))
            {
                m.in_noc = false;
                // Go along the message path and decrement tile
                //  message counters for all links traversed
                schedule_update_noc_message_counts(m, noc, false);
                remove_message = true;
            }
            std::advance(it, 1);
            if (remove_message)
            {
                q.erase(curr);
            }
        }
    }

    return;
}

sanafe::MessagePriorityQueue sanafe::schedule_init_timing_priority(
        std::vector<MessageFifo> &message_queues_per_core)
{
    // Create the priority queue of messages
    MessagePriorityQueue priority;

    TRACE1("Initializing priority queue.\n");
    for (auto &q : message_queues_per_core)
    {
        // For each per-core queue, get the first message for that core
        //  and add it to the corresponding priority queue
        if (q.size() > 0) // messages
        {
            Message &m = q.front();
            q.pop_front();
            m.sent_timestamp = m.generation_delay;
            priority.push(m);
        }
        else
        {
            TRACE1("No messages for core %d\n", i);
        }
    }

#ifdef DEBUG2
    int i = 0;
    for (Core *curr = priority_queue; curr != NULL; curr = curr->next)
    {
        // TODO
    }
#endif

    return priority;
}
