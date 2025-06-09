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
#include <cstddef>
#include <functional> // For std::reference_wrapper
#include <list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <booksim_lib.hpp>

#include "chip.hpp"
#include "message.hpp"
#include "print.hpp"
#include "schedule.hpp"

sanafe::NocInfo::NocInfo(const Scheduler &scheduler)
        : noc_width_in_tiles(scheduler.noc_width_in_tiles)
        , noc_height_in_tiles(scheduler.noc_height_in_tiles)
        , core_count(scheduler.core_count)
        , max_cores_per_tile(scheduler.max_cores_per_tile)
{
    messages_received = std::vector<sanafe::MessageFifo>(core_count);
    core_finished_receiving = std::vector<double>(core_count);
}

void sanafe::schedule_messages(TimestepHandle &ts, Scheduler &scheduler,
        const BookSimConfig &booksim_config)
{
    if (scheduler.timing_model == timing_model_simple)
    {
        TRACE1(CHIP, "Running simple timing model\n");
        schedule_messages_simple(ts, scheduler);
    }
    else if (scheduler.timing_model == timing_model_detailed)
    {
        TRACE1(CHIP, "Running detailed timing model\n");
        schedule_messages_detailed(ts, scheduler);
    }
    else if (scheduler.timing_model == timing_model_cycle_accurate)
    {
        TRACE1(CHIP, "Running cycle-accurate timing model\n");
        schedule_messages_cycle_accurate(ts, booksim_config, scheduler);
    }
    else
    {
        INFO("Error: Timing model:%d not recognized\n", scheduler.timing_model);
        throw std::invalid_argument("Timing model not recognized");
    }
}

void sanafe::schedule_messages_simple(
        TimestepHandle &ts, Scheduler &scheduler)
{
    // Simple analytical model, that takes the maximum of either neuron or
    //  message processing for each core, and takes the maximum latency of
    //  any core in the design.
    //
    // This scheduler is extremely simple and therefore
    //  fast (also easily done in parallel), however, it won't capture any
    //  network interactions or resource contention.
    Timestep &ts_data = ts.get();
    const size_t cores = ts_data.messages.size();
    std::vector<double> neuron_processing_latencies(cores, 0.0);
    std::vector<double> message_processing_latencies(cores, 0.0);
    for (size_t sending_core = 0; sending_core < ts_data.messages.size();
            ++sending_core)
    {
        std::list<Message> &q = ts_data.messages.at(sending_core);
        for (Message &m : q)
        {
            neuron_processing_latencies[sending_core] += m.generation_delay;
            message_processing_latencies[m.dest_core_id] += m.receive_delay;
            // Update message delays using very simple timing model
            m.blocked_delay = 0.0; // No blocking modeled
            m.network_delay = m.min_hop_delay;
        }
    }

    const double max_message_processing =
            *std::max_element(message_processing_latencies.begin(),
                    message_processing_latencies.end());
    const double max_neuron_processing =
            *std::max_element(neuron_processing_latencies.begin(),
                    neuron_processing_latencies.end());

    ts_data.sim_time = std::max(max_message_processing, max_neuron_processing);
    // Account for fixed costs per timestep e.g., house-keeping or global sync
    ts_data.sim_time += scheduler.timestep_sync_delay;
    scheduler.timesteps_to_write.push(ts);
}

void sanafe::schedule_messages_cycle_accurate(
        TimestepHandle &ts, const BookSimConfig &config, Scheduler &scheduler)
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
    Timestep &ts_data = ts.get();
    for (auto &core_messages : ts_data.messages)
    {
        for (auto &message : core_messages)
        {
            if (message.mid == placeholder_mid)
            {
                std::pair<std::string, int> src_neuron =
                        std::make_pair(message.src_neuron_group_id,
                                static_cast<int>(message.src_neuron_offset));
                std::pair<int, int> src_hw =
                        std::make_pair(static_cast<int>(message.src_tile_id),
                                static_cast<int>(message.src_core_offset));

                booksim_create_processing_event(
                        static_cast<int>(message.timestep),
                        std::move(src_neuron), std::move(src_hw),
                        message.generation_delay);
            }
            else
            {
                std::pair<std::string, int> src_neuron =
                        std::make_pair(message.src_neuron_group_id,
                                static_cast<int>(message.src_neuron_offset));
                std::pair<int, int> src_hw =
                        std::make_pair(static_cast<int>(message.src_tile_id),
                                static_cast<int>(message.src_core_offset));
                std::pair<int, int> dest_hw =
                        std::make_pair(static_cast<int>(message.dest_tile_id),
                                static_cast<int>(message.dest_core_offset));

                booksim_create_spike_event(static_cast<int>(message.timestep),
                        std::move(src_neuron), std::move(src_hw),
                        std::move(dest_hw), message.generation_delay,
                        message.receive_delay);
            }
        }
    }

    // Messages have been sent to the library, so now just execute the
    //  simulation and return simulated time
    TRACE1(SCHEDULER, "Running Booksim2 simulation\n");
    const double booksim_time = booksim_run(config);

    ts_data.sim_time = booksim_time;
    // Account for fixed costs per timestep e.g., house-keeping or global sync
    ts_data.sim_time += scheduler.timestep_sync_delay;
    scheduler.timesteps_to_write.push(ts);
}

void sanafe::schedule_create_threads(
        Scheduler &scheduler, const int scheduler_thread_count)
{
    TRACE1(CHIP, "Creating %d scheduler threads\n", scheduler_thread_count);
    for (int thread_id = 0; thread_id < scheduler_thread_count; thread_id++)
    {
        TRACE1(CHIP, "Created scheduler thread:%d\n", thread_id);
        scheduler.scheduler_threads.emplace_back(
                &schedule_messages_thread, std::ref(scheduler), thread_id);
    }
}

// **** Detailed scheduler implementation ****

void sanafe::schedule_messages_detailed(TimestepHandle &ts, Scheduler &scheduler)
{
    if (scheduler.scheduler_threads.empty())
    {
        schedule_messages_timestep(ts, scheduler);
    }
    else
    {
        TRACE1(SCHEDULER, "Pushing timestep:%ld to be scheduled\n",
                ts->timestep);
        scheduler.timesteps_to_schedule.push(ts);
        return;
    }
}

double sanafe::schedule_messages_timestep(
        TimestepHandle &ts, Scheduler &scheduler)
{
    // Schedule the global order of messages using a semi-analytical timing
    //  model. This sits in between the complexity of the simple and
    //  cycle-accurate model, capturing some but not all network effects. This
    //  takes a vector containing a list of messages per core, and scheduler
    //  parameters (mostly NoC configuration parameters). Returns the timestamp
    //  of the last scheduled event, i.e., the total time-step delay
    MessagePriorityQueue priority;
    NocInfo noc(scheduler);
    double last_timestamp = 0.0;
    Timestep &ts_data = ts.get();

    auto messages_sent_per_core = schedule_init_message_queues(ts_data, noc);
    std::vector<MessageFifo> scheduled_messages_per_core(noc.core_count);

    priority = schedule_init_timing_priority(messages_sent_per_core);

    // Each core has a queue of received messages. A structure tracks how
    //  many in-flight messages are in the NoC and occupy each tile. We
    //  track the number of messages passing through each tile at
    //  the point of sending, and the average processing delay of
    //  all of those messages. When a message is added or removed from the
    //  NoC we update the average counts.
    TRACE1(SCHEDULER, "Scheduling global order of messages.\n");
    while (!priority.empty())
    {
        // Get the core's queue with the earliest simulation time
        Message m = priority.top();
        priority.pop();
        last_timestamp = std::max(last_timestamp, m.sent_timestamp);

        // Update the Network-on-Chip state
        noc_update_all_tracked_messages(m.sent_timestamp, noc);

        if (!m.placeholder)
        {
            schedule_handle_message(m, scheduler, noc);
            last_timestamp = std::max(last_timestamp, m.processed_timestamp);
        }
        // Else messages without a destination (neuron) are placeholders.
        //  Dummy messages account for processing time that does not
        //  result in any spike messages.

        // Get the next message for this core
        const size_t src_core = m.src_core_id;
        if (!messages_sent_per_core[src_core].empty())
        {
            const double next_sent_timestamp = schedule_push_next_message(
                    messages_sent_per_core, priority, m);
            last_timestamp = std::max(last_timestamp, next_sent_timestamp);
        }
        else
        {
            TRACE1(SCHEDULER, "\tCore finished simulating\n");
        }

        // Push the updated message with its delays and timestamps set
        scheduled_messages_per_core[src_core].emplace_back(m);

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

    ts_data.sim_time = last_timestamp;
    // Account for fixed costs per timestep e.g., house-keeping or global sync
    ts_data.sim_time += scheduler.timestep_sync_delay;

    // Update all messages with delays and timestamps set
    ts_data.messages = std::move(scheduled_messages_per_core);
    scheduler.timesteps_to_write.push(ts);

    return ts_data.sim_time;
}

std::vector<sanafe::MessageFifo>
sanafe::schedule_init_message_queues(const Timestep &ts, NocInfo &noc)
{
    const size_t total_links = noc.noc_height_in_tiles *
            noc.noc_width_in_tiles *
            (sanafe::ndirections + noc.max_cores_per_tile);
    noc.message_density = std::vector<double>(total_links, 0.0);

    // Return a local copy of the per-core messages
    return ts.messages;
}

void sanafe::schedule_handle_message(
        Message &m, Scheduler &scheduler, NocInfo &noc)
{
    TRACE1(SCHEDULER, "Processing message for nid:%s.%zu\n",
            m.src_neuron_group_id.c_str(), m.src_neuron_offset);
    TRACE1(SCHEDULER, "Send delay:%e\n", m.generation_delay);
    TRACE1(SCHEDULER, "Receive delay:%e\n", m.receive_delay);
    const size_t dest_core = m.dest_core_id;
    // Figure out if we are able to send a message into the
    //  network i.e., is the route to the dest core
    //  saturated and likely to block? Sum along the route
    //  and see the density of messages along all links.
    const double messages_along_route = noc.calculate_route_congestion(m);
    const auto path_capacity =
            static_cast<double>((m.hops + 1UL) * scheduler.buffer_size);

    if (messages_along_route > path_capacity)
    {
        // Use heuristic for estimating delay based on route congestion,
        //  path capacity in messages and the mean delay per message
        m.blocked_delay = (messages_along_route - path_capacity) *
                noc.mean_in_flight_receive_delay;
        m.sent_timestamp += m.blocked_delay;
    }
    else
    {
        m.blocked_delay = 0.0; // Path isn't at capacity; no blocking
    }

    const double congestion_delay = messages_along_route *
            noc.mean_in_flight_receive_delay /
            (static_cast<double>(m.hops) + 1.0);
    TRACE1(SCHEDULER, "Path capacity:%lf messages:%lf congestion delay:%e\n",
            path_capacity, messages_along_route, congestion_delay);

    // Update the messages timestamps, both when the message is received and
    //  when the receiving core has finished processing it
    m.network_delay = std::max(m.min_hop_delay, congestion_delay);
    const double earliest_received_time = m.sent_timestamp + m.network_delay;
    m.received_timestamp = std::max(
            noc.core_finished_receiving[dest_core], earliest_received_time);

    noc.core_finished_receiving[dest_core] =
            std::max((noc.core_finished_receiving[dest_core] + m.receive_delay),
                    (earliest_received_time + m.receive_delay));
    m.processed_timestamp = noc.core_finished_receiving[dest_core];

    // Now, push the message into the right receiving queue. Calculate the
    //  network delay and when the message is received
    m.in_noc = true;
    noc.messages_received[dest_core].push_back(m);
    noc_update_message_tracking(m, noc, true);
}

double sanafe::schedule_push_next_message(
        std::vector<MessageFifo> &messages_sent_per_core,
        MessagePriorityQueue &priority, const Message &current_message)
{
    auto &q = messages_sent_per_core[current_message.src_core_id];
    Message &next_message = q.front();

    // If applicable, schedule this next message immediately
    //  after the current message finishes sending
    next_message.sent_timestamp =
            current_message.sent_timestamp + next_message.generation_delay;
    priority.push(next_message);

    // Record the latest timestamp before deleting the copy of the message
    const double last_timestamp = next_message.sent_timestamp;
    q.pop_front();

    return last_timestamp;
}

void sanafe::noc_update_all_tracked_messages(const double t, NocInfo &noc)
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
                noc_update_message_tracking(m, noc, false);
                TRACE1(SCHEDULER, "Removing message mid:%zu\n", m.mid);
                return true; // Remove this message
            }
            return false; // Keep this message
        });
    }
}

void sanafe::noc_update_message_tracking(
        const Message &m, NocInfo &noc, const bool entering_noc)
{
    // Update the tracked state of the NoC, accounting for a single message
    //  either entering or leaving the NoC (message_in)
    noc.update_message_density(m, entering_noc);
    noc.update_rolling_averages(m, entering_noc);
#if (DEBUG_LEVEL_SCHEDULER > 0)
    INFO("Message density:");
    for (auto &val : noc.message_density)
    {
        printf("%.2lf,", val);
    }
    printf("\n");
#endif
    TRACE1(SCHEDULER, "Mean receive delay:%e\n",
            noc.mean_in_flight_receive_delay);
    TRACE1(SCHEDULER, "Messages:%ld\n", noc.messages_in_noc);
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

void sanafe::NocInfo::update_rolling_averages(
        const Message &message, const bool entering_noc)
{
    if (entering_noc)
    {
        // Message entering NoC
        mean_in_flight_receive_delay +=
                (message.receive_delay - mean_in_flight_receive_delay) /
                (static_cast<double>(messages_in_noc) + 1.0);
        messages_in_noc++;
    }
    else
    {
        // Message leaving the NoC
        if (messages_in_noc > 1)
        {
            mean_in_flight_receive_delay +=
                    (mean_in_flight_receive_delay - message.receive_delay) /
                    (static_cast<double>(messages_in_noc) - 1.0);
        }
        else
        {
            mean_in_flight_receive_delay = 0.0;
        }

        messages_in_noc--;
    }
}

void sanafe::NocInfo::update_message_density(
        const Message &message, const bool entering_noc)
{
    // Adjust by dividing by the total number of links along the path, also
    //  including the output link at the sending core and input link at the
    //  receiving core, i.e. the hops plus 2. The total sum of the added
    //  densities along the path should equal one for one new message.
    constexpr double input_plus_output_link = 2.0;
    double adjust = (1.0 /
            (input_plus_output_link + static_cast<double>(message.hops)));

    if (!entering_noc)
    {
        // Remove message from NoC
        adjust *= -1.0;
    }

    auto [x_increment, y_increment] = get_route_xy_increments(message);
    size_t prev_direction = sanafe::ndirections + (message.src_core_offset);
    for (size_t x = message.src_x; x != message.dest_x; x += x_increment)
    {
        const int direction = (x_increment > 0) ? sanafe::east : sanafe::west;
        if (x == message.src_x)
        {
            const size_t link = sanafe::ndirections + (message.src_core_offset);
            message_density[idx(x, message.src_y, link)] += adjust;
        }
        else
        {
            message_density[idx(x, message.src_y, direction)] += adjust;
        }
        prev_direction = direction;
    }
    for (size_t y = message.src_y; y != message.dest_y; y += y_increment)
    {
        const int direction = (y_increment > 0) ? sanafe::north : sanafe::south;
        if ((message.src_x == message.dest_x) && (y == message.src_y))
        {
            const size_t link = sanafe::ndirections + message.src_core_offset;
            message_density[idx(message.dest_x, y, link)] += adjust;
        }
        else
        {
            message_density[idx(message.dest_x, y, prev_direction)] += adjust;
        }

        prev_direction = direction;
    }

    if ((message.src_x == message.dest_x) && (message.src_y == message.dest_y))
    {
        const size_t link = sanafe::ndirections + (message.src_core_offset);
        message_density[idx(message.dest_x, message.dest_y, link)] += adjust;
    }
    else
    {
        message_density[idx(message.dest_x, message.dest_y, prev_direction)] +=
                adjust;
    }
}

double sanafe::NocInfo::calculate_route_congestion(const Message &m) const
{
    // Calculate the total flow density as a metric for route congestion along a
    //  spike message's route. This is given by the sum of the densities, for
    //  all links the message will travel i.e., the message path. Note that we
    //  calculat ethe path assuming a dimension-order routing scheme.
    // TODO: extend this to generalize to different routing schemes.
    auto [x_increment, y_increment] = get_route_xy_increments(m);
    double flow_density = 0.0;

    size_t prev_direction = sanafe::ndirections + (m.src_core_offset);
    for (size_t x = m.src_x; x != m.dest_x; x += x_increment)
    {
        const int direction = (x_increment > 0) ? sanafe::east : sanafe::west;
        if (x == m.src_x)
        {
            const size_t link = sanafe::ndirections + m.src_core_offset;
            flow_density += message_density[idx(x, m.src_y, link)];
        }
        else
        {
            flow_density += message_density[idx(x, m.src_y, direction)];
        }
        prev_direction = direction;
    }

    for (size_t y = m.src_y; y != m.dest_y; y += y_increment)
    {
        const int direction = (y_increment > 0) ? sanafe::north : sanafe::south;
        if (m.src_x == m.dest_x && y == m.src_y)
        {
            const size_t link = sanafe::ndirections + m.src_core_offset;
            flow_density += message_density[idx(m.dest_x, y, link)];
        }
        else
        {
            flow_density += message_density[idx(m.dest_x, y, prev_direction)];
        }
        prev_direction = direction;
    }
    // Handle the last (destination) tile
    if ((m.src_x == m.dest_x) && (m.src_y == m.dest_y))
    {
        const size_t link = sanafe::ndirections + m.src_core_offset;
        flow_density += message_density[idx(m.dest_x, m.dest_y, link)];
    }
    else
    {
        flow_density +=
                message_density[idx(m.dest_x, m.dest_y, prev_direction)];
    }

#ifndef NDEBUG
    constexpr double epsilon = 0.1; // In case density is very slightly below 0
#endif
    assert(flow_density >= (-epsilon));
    return flow_density;
}


std::pair<int, int> sanafe::NocInfo::get_route_xy_increments(
        const Message &m) noexcept
{
    const int x_increment = (m.src_x < m.dest_x) ? 1 : -1;
    const int y_increment = (m.src_y < m.dest_y) ? 1 : -1;

    return std::make_pair(x_increment, y_increment);
}


// **** Thread management ****
// TODO: make this agnostic to scheduling algorithm, so it can be applied to
//  the simple and cycle accurate models in future

void sanafe::schedule_messages_thread(Scheduler &scheduler, const int thread_id)
{
    while (!scheduler.should_stop)
    {
        if (scheduler.should_stop)
        {
            break;
        }

        TimestepHandle ts;
        const bool got_ts = scheduler.timesteps_to_schedule.pop(ts);
        if (got_ts)
        {
            TRACE1(SCHEDULER, "tid:%d Scheduling ts:%ld\n", thread_id,
                    ts->timestep);
            schedule_messages_timestep(ts, scheduler);
        }
    }

    TRACE1(SCHEDULER, "Scheduler thread tid:%d terminating gracefully\n",
            thread_id);
}

void sanafe::schedule_stop_all_threads(Scheduler &scheduler)
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
}