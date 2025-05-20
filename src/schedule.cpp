// Copyright (c) 2025 - The University of Texas at Austin
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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional> // For std::reference_wrapper
#include <iterator>
#include <vector>

#include <booksim_lib.hpp>

#include "chip.hpp"
#include "network.hpp"
#include "print.hpp"
#include "schedule.hpp"

sanafe::NocInfo::NocInfo(const int width, const int height,
        const int core_count, const size_t max_cores_per_tile)
        : noc_width_in_tiles(width)
        , noc_height_in_tiles(height)
        , core_count(core_count)
        , max_cores_per_tile(max_cores_per_tile)
{
}

double sanafe::schedule_messages_simple(Timestep &ts, Scheduler &scheduler)
{
    // Simple analytical model, that takes the maximum of either neuron or
    //  message processing for each core, and takes the maximum latency of
    //  any core in the design.
    //
    // This scheduler is extremely simple and therefore
    //  fast (also easily done in parallel), however, it won't capture any
    //  network interactions or resource contention.
    const size_t cores = ts.messages->size();
    std::vector<double> neuron_processing_latencies(cores, 0.0);
    std::vector<double> message_processing_latencies(cores, 0.0);
    for (size_t sending_core = 0; sending_core < ts.messages->size();
            ++sending_core)
    {
        std::list<Message> &q = ts.messages->at(sending_core);
        for (Message &m : q)
        {
            neuron_processing_latencies[sending_core] += m.generation_delay;
            message_processing_latencies[m.dest_core_id] += m.receive_delay;
        }
    }

    double max_message_processing =
            *std::max_element(message_processing_latencies.begin(),
                    message_processing_latencies.end());
    double max_neuron_processing =
            *std::max_element(neuron_processing_latencies.begin(),
                    neuron_processing_latencies.end());

    ts.sim_time = std::max(max_message_processing, max_neuron_processing);
    scheduler.timesteps_to_write.push(ts);
    return ts.sim_time;
}

double sanafe::schedule_messages_cycle_accurate(
        Timestep &ts, const BookSimConfig &config, Scheduler &scheduler)
{
    // Cycle-accurate (NoC)-based timing model for highly accurate network
    //  simulation, using external simulator Booksim 2. This is the most
    //  accurate model, but is probably too slow for meaningful design-space
    //  exploration (about 1000x/100x slower than the simple/detailed models).
    //  The timings within each core's pipeline is predictible, but network
    //  effects are more complex. Booksim 2 has been developed in Stanford and
    //  accurately models the transactions happening within an NoC [Jiang 2013]
    //
    // The version of Booksim 2 used here has had substantial modifications,
    //  including a new static library interface (instead of being standalone)
    // TODO: support running across multiple threads like the detailed model
    for (auto &core_messages : *(ts.messages))
    {
        for (auto &message : core_messages)
        {
            if (message.mid == placeholder_mid)
            {
                booksim_create_processing_event(message.timestep,
                        std::make_pair(message.src_neuron_group_id,
                                message.src_neuron_id),
                        std::make_pair(
                                message.src_tile_id, message.src_core_offset),
                        message.generation_delay);
            }
            else
            {
                booksim_create_spike_event(message.timestep,
                        std::make_pair(message.src_neuron_group_id,
                                message.src_neuron_id),
                        std::make_pair(
                                message.src_tile_id, message.src_core_offset),
                        std::make_pair(
                                message.dest_tile_id, message.dest_core_offset),
                        message.generation_delay, message.receive_delay);
            }
        }
    }

    // Messages have been sent to the library, so now just execute the
    //  simulation and return simulated time
    TRACE1(SCHEDULER, "Running Booksim2 simulation\n");
    double booksim_time = booksim_run(config);

    ts.sim_time = booksim_time;
    scheduler.timesteps_to_write.push(ts);
    return booksim_time;
}

void sanafe::schedule_messages_detailed(Timestep &ts, Scheduler &scheduler)
{
    if (scheduler.scheduler_threads.empty())
    {
        schedule_messages_timestep(ts, scheduler);
    }
    else
    {
        TRACE1(SCHEDULER, "Pushing timestep:%ld to be scheduled\n",
                ts.timestep);
        scheduler.timesteps_to_schedule.push(ts);
        return;
    }
}

void sanafe::schedule_messages_thread(
        sanafe::Scheduler &scheduler, const int tid)
{
    while (!scheduler.should_stop)
    {
        if (scheduler.should_stop)
        {
            break;
        }

        Timestep ts;
        bool got_ts = scheduler.timesteps_to_schedule.pop(ts);
        if (got_ts)
        {
            TRACE1(SCHEDULER, "tid:%d Scheduling ts:%ld\n", tid, ts.timestep);
            schedule_messages_timestep(ts, scheduler);
        }
    }

    TRACE1(SCHEDULER, "Scheduler thread tid:%d terminating gracefully\n", tid);
    return;
}

void sanafe::schedule_stop_all_threads(sanafe::Scheduler &scheduler,
        std::ofstream &message_trace, sanafe::RunData &rd)
{
    TRACE1(SCHEDULER, "Stopping all scheduling threads.\n");
    scheduler.timesteps_to_schedule.wait_until_empty();
    TRACE1(SCHEDULER, "All messages scheduled so terminate threads.\n");
    scheduler.should_stop = true;
    scheduler.timesteps_to_schedule.set_terminate();
    scheduler.timesteps_to_write.set_terminate();

    // This function is blocking i.e., waits for all workers to finish
    for (auto &thread : scheduler.scheduler_threads)
    {
        if (thread.joinable())
        {
            thread.join();
        }
    }

    TRACE1(SCHEDULER, "All threads stopped successfully.\n");
    return;
}

double sanafe::schedule_messages_timestep(Timestep &ts, Scheduler &scheduler)
{
    // Schedule the global order of messages using a semi-analytical timing
    //  model. This sits in between the complexity of the simple and
    //  cycle-accurate model, capturing some but not all network effects. This
    //  takes a vector containing a list of messages per core, and scheduler
    //  parameters (mostly NoC configuration parameters). Returns the timestamp
    //  of the last scheduled event, i.e., the total time-step delay
    MessagePriorityQueue priority;
    NocInfo noc(scheduler.noc_width, scheduler.noc_height, scheduler.core_count,
            scheduler.max_cores_per_tile);
    double last_timestamp;
    const size_t total_links = noc.noc_height_in_tiles *
            noc.noc_width_in_tiles *
            (sanafe::ndirections + noc.max_cores_per_tile);
    noc.message_density = std::vector<double>(total_links, 0.0);

    std::vector<MessageFifo> messages_sent_per_core(noc.core_count);
    for (size_t core = 0; core < ts.messages->size(); core++)
    {
        auto &q = ts.messages->at(core);
        for (Message &m : q)
        {
            messages_sent_per_core[core].push_back(m);
        }
    }

    noc.messages_received = std::vector<MessageFifo>(noc.core_count);
    noc.core_finished_receiving = std::vector<double>(noc.core_count);

    priority = schedule_init_timing_priority(messages_sent_per_core);
    last_timestamp = 0.0;
    TRACE1(SCHEDULER, "Scheduling global order of messages.\n");

    // Each core has a queue of received messages. A structure tracks how
    //  many in-flight messages are in the NoC and occupy each tile. We
    //  track the number of messages passing through each tile at
    //  the point of sending, and the average processing delay of
    //  all of those messages. When a message is added or removed from the
    //  NoC we update the average counts.
    while (!priority.empty())
    {
        // Get the core's queue with the earliest simulation time
        Message m = priority.top();
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
            TRACE1(SCHEDULER, "Processing message for nid:%s.%zu\n",
                    m.src_neuron_group_id.c_str(), m.src_neuron_id);
            TRACE1(SCHEDULER, "Send delay:%e\n", m.generation_delay);
            TRACE1(SCHEDULER, "Receive delay:%e\n", m.receive_delay);
            const int dest_core = m.dest_core_id;
            // Figure out if we are able to send a message into the
            //  network i.e., is the route to the dest core
            //  saturated and likely to block? Sum along the route
            //  and see the density of messages along all links.
            const double messages_along_route =
                    schedule_calculate_messages_along_route(m, noc);

            const size_t path_capacity = (m.hops + 1) * scheduler.buffer_size;
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
            TRACE1(SCHEDULER, "Path capacity:%zu messages:%lf delay:%e\n",
                    path_capacity, messages_along_route, network_delay);

            const double earliest_received_time =
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
        if (!messages_sent_per_core[src_core].empty())
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
            TRACE1(SCHEDULER, "\tCore finished simulating\n");
        }

#ifdef DEBUG
        // Print contents of priority queue. Because of how the queue works,
        //  to view all elements we need to copy and pop off elements
        auto view = priority;
        INFO("***\n");
        while (const Message &curr = view.pop())
        {
            INFO("m:%e\n", curr.sent_timestamp);
        }
        INFO("***\n");
#endif

        TRACE1(SCHEDULER, "Priority size:%zu\n", priority.size());
    }
    TRACE1(SCHEDULER, "Scheduler finished.\n");

    ts.sim_time = last_timestamp;
    scheduler.timesteps_to_write.push(ts);

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
    int x_increment;
    int y_increment;
    // Adjust by dividing by the total number of links along the path, also
    //  including the output link at the sending core and input link at the
    //  receiving core, i.e. the hops plus 2. The total sum of the added
    //  densities along the path should equal one for one new message.
    constexpr double input_plus_output_link = 2.0;
    double adjust = (1.0 / (input_plus_output_link + m.hops));

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
}

double sanafe::schedule_calculate_messages_along_route(
        const Message &m, NocInfo &noc)
{
    // Calculate the total flow density along a spike message route to act as a
    //  link congestion metric. This is given by the sum of the densities, for
    //  all links the message will travel i.e., the message path
    int x_increment;
    int y_increment;

    double flow_density = 0.0;
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
        const int direction = (x_increment > 0) ? sanafe::east : sanafe::west;
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
        const int direction = (y_increment > 0) ? sanafe::north : sanafe::south;
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

#ifndef NDEBUG
    constexpr double epsilon = 0.1; // In case density is very slightly below 0
#endif
    assert(flow_density >= (-epsilon));
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
        // Use remove_if to identify messages that should be removed
        q.remove_if([&](Message &m) -> bool {
            if (m.in_noc && (t >= m.received_timestamp))
            {
                m.in_noc = false;
                // Update message counts before removing
                schedule_update_noc_message_counts(m, noc, false);
                return true; // Remove this message
            }
            return false; // Keep this message
        });
    }
}

sanafe::MessagePriorityQueue sanafe::schedule_init_timing_priority(
        std::vector<MessageFifo> &message_queues_per_core)
{
    MessagePriorityQueue priority;

    TRACE1(SCHEDULER, "Initializing priority queue.\n");
    for (auto &message_queue : message_queues_per_core)
    {
        // Initialize each per-core queue with the first message for that core
        //  We use priority queues so that we can sort by simulation time and
        //  always pop the next timing event
        if (!message_queue.empty()) // messages
        {
            Message message = message_queue.front();
            message_queue.pop_front();
            message.sent_timestamp = message.generation_delay;
            priority.push(message);
        }
        else
        {
            TRACE1(SCHEDULER, "No messages for core\n");
        }
    }

    return priority;
}
