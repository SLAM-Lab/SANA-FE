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
#include <unordered_set>
#include <utility>
#include <vector>

#include <booksim_lib.hpp>

#include "chip.hpp"
#include "message.hpp"
#include "print.hpp"
#include "schedule.hpp"

std::vector<size_t> sanafe::Flow::generate_path(const sanafe::NewNocInfo &noc) const
{
    std::vector<size_t> path;

    // Assuming dimension-order (XY) routing
    auto [x_increment, y_increment] = get_path_xy_increments(*this, noc);
    const size_t src_core_offset = src % noc.max_cores_per_tile;
    const size_t src_tile = src / noc.max_cores_per_tile;
    const size_t src_x = src_tile / noc.noc_height_in_tiles;
    const size_t src_y = src_tile % noc.noc_height_in_tiles;

    const size_t dest_core_offset = dest % noc.max_cores_per_tile;
    const size_t dest_tile = dest / noc.max_cores_per_tile;
    const size_t dest_x = dest_tile / noc.noc_height_in_tiles;
    const size_t dest_y = dest_tile % noc.noc_height_in_tiles;

    size_t prev_direction = sanafe::ndirections + src_core_offset;
    for (size_t x = src_x; x != dest_x; x += x_increment)
    {
        const int direction = (x_increment > 0) ? sanafe::east : sanafe::west;
        if (x == src_x)
        {
            const size_t link = sanafe::ndirections + (src_core_offset);
            path.push_back(noc.idx(x, src_y, link));
        }
        else
        {
            path.push_back(noc.idx(x, src_y, direction));
        }
        prev_direction = direction;
    }
    for (size_t y = src_y; y != dest_y; y += y_increment)
    {
        const int direction = (y_increment > 0) ? sanafe::north : sanafe::south;
        if ((src_x == dest_x) && (y == src_y))
        {
            const size_t link = sanafe::ndirections + src_core_offset;
            path.push_back(noc.idx(dest_x, y, link));
        }
        else
        {
            path.push_back(noc.idx(dest_x, y, prev_direction));
        }

        prev_direction = direction;
    }

    if ((src_x == dest_x) && (src_y == dest_y))
    {
        const size_t link = sanafe::ndirections + (src_core_offset);
        path.push_back(noc.idx(dest_x, dest_y, link));
    }
    else
    {
        path.push_back(noc.idx(dest_x, dest_y, prev_direction));
    }

    // Finally push back the link to the desination core
    path.push_back(noc.idx(dest_x, dest_y,
            sanafe::ndirections + noc.max_cores_per_tile + dest_core_offset));

    return path;
}

// New file for noc.cpp or something like this.. 1 class per file generally
sanafe::NewNocInfo::NewNocInfo(const Scheduler &scheduler)
        : noc_width_in_tiles(scheduler.noc_width_in_tiles)
        , noc_height_in_tiles(scheduler.noc_height_in_tiles)
        , core_count(scheduler.core_count)
        , max_cores_per_tile(scheduler.max_cores_per_tile)
{
    messages_received = std::vector<sanafe::MessageFifo>(core_count);
    core_finished_receiving = std::vector<double>(core_count);
}

std::pair<std::optional<size_t>, size_t> sanafe::NewNocInfo::find_overlapping_links(const sanafe::Flow &flow1, const sanafe::Flow &flow2) const
{
    auto path1 = flow1.generate_path(*this);
    auto path2 = flow2.generate_path(*this);

    // Store elements of path2 in an unordered_set for efficient lookups
    const std::unordered_set<size_t> path2_set(path2.begin(), path2.end());

    std::optional<size_t> first_overlapping_link;
    size_t overlapping_links = 0UL;
    // Iterate through path1 to find common links with path2
    for (size_t i = 0; i < path1.size(); ++i)
    {
        if (path2_set.count(path1[i]))
        {
            if (!first_overlapping_link.has_value())
            {
                first_overlapping_link = i;
            }
            ++overlapping_links;
        }
    }

    return std::make_pair(first_overlapping_link, overlapping_links);
}

void sanafe::NewNocInfo::update_message_density(const Message &message, bool entering_noc)
{
    // Using the following psuedocode:
    //  if entering_noc:
    //    if flow does not exist
    //      create and initialize new flow: fnew
    //      for all existing flows fexist:
    //        overlap_pair = find_overlapping_links(fnew, fexist):
    //        add overlap fnew to list of overlaps
    //    push message to flow for given src/dst pair
    //  else (if leaving noc):
    //    subtract one from flow (assert >= 0)
    //    if flow is empty:
    //      for now do nothing, leave it as an empty flow in case used later

    // TODO: in the future, change this code to build a reverse index within
    //  each link that tracks all flows over that link. This way we can
    //  efficiently search for overlapping flows when adding a new flow without
    //  having to scan through all possible flows. As we track the overlapping
    //  flow, we would need to consider the relative link position of the
    //  overlap, and the number of overlapping links, both of which could be
    //  easily obtained as we follow the path and its reverse indexes!
    std::pair<size_t, size_t> src_dest_pair = {
            message.src_core_id, message.dest_core_id};
    if (entering_noc)
    {
        auto flow_it = flows.find(src_dest_pair);

        bool flow_exists = (flow_it != flows.end());
        if (!flow_exists)
        {
            Flow new_flow(message.src_core_id, message.dest_core_id);
            // Find all overlapping links. TODO: we'll do this more efficiently in future
            //  Associate links with their flows and just go over all links to
            //  avoid searching through every possible flow...
            for (auto &[s, f_exist] : flows)
            {
                auto overlap_with_new = find_overlapping_links(new_flow, f_exist);
                // Simply pushing back is useless, we need to associate the flow f_exist
                //  with the information...
                std::pair<size_t, size_t> existing_src_dest_pair = {f_exist.src, f_exist.dest};
                new_flow.flows_sharing_links[existing_src_dest_pair] = overlap_with_new;

                auto overlap_with_existing = find_overlapping_links(f_exist, new_flow);
                f_exist.flows_sharing_links[src_dest_pair] = overlap_with_existing;
            }
            // Add new flow
            new_flow.messages_in_flight = 1;
            flows.try_emplace(src_dest_pair, new_flow);

            // INFO("Added new flows src:%zu dest:%zu (total=%zu)\n",
            //         src_dest_pair.first, src_dest_pair.second, flows.size());
        }
        else
        {
            ++flow_it->second.messages_in_flight;
            // INFO("Flow exists src:%zu dest:%zu (in flight=%zu)\n",
            //     flow_it->second.src,
            //     flow_it->second.dest,
            //     flow_it->second.messages_in_flight);
        }
    }
    else // removing message from noc
    {
        if (flows.at(src_dest_pair).messages_in_flight <= 0)
        {
            INFO("Error: no messages in flow, can't pop\n");
            throw std::logic_error("No messages in flow");
        }
        --flows.at(src_dest_pair).messages_in_flight;
        // INFO("Removing from flow src:%zu dest:%zu (in flight:%zu)\n",
        //         flows.at(src_dest_pair).src, flows.at(src_dest_pair).dest,
        //         flows.at(src_dest_pair).messages_in_flight);
    }
}

void sanafe::NewNocInfo::update_rolling_averages(const Message& message, bool entering_noc)
{
    // TODO: same as before
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

void sanafe::schedule_messages_simple(TimestepHandle &ts, Scheduler &scheduler)
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

// TODO: we should decouple the parallel job dispatcher from the scheduling
//  algorithm better. In future, we may want to make the simple or cycle-accurate
//  models dispatched in parallel in the same way as the semi-analytical model
// TODO: I think the current way of implementing the scheduler is confusing
//  and inflexible, when considering different scheduling algorithms.
//  It's a mix of c++ classes and c-like scheduling functions. A better approach
//  might be to have some base Scheduler class, with some virtual functions
//  that are required. New schedulers will derive from this class, and can have
//  their own internal functions and state if needed.
//  In future, we may want to break models and schedule into their own directories
//  in addition to the files in src/
void sanafe::schedule_messages_detailed(
        TimestepHandle &ts, Scheduler &scheduler)
{
    if (scheduler.scheduler_threads.empty())
    {
        // For 1000 tests
        // TODO: old function real runtime 0m5.9s
        //schedule_messages_timestep(ts, scheduler);
        // TODO: new function real runtime 12s about 2x slower with optimization
        //  if the algorithm is accurate, we can optimize it to hopefully get it around
        //  the same ballpark
        schedule_messages_timestep_new(ts, scheduler);
    }
    else
    {
        TRACE1(SCHEDULER, "Pushing timestep:%ld to be scheduled\n",
                ts->timestep);
        scheduler.timesteps_to_schedule.push(ts);
        return;
    }
}

// ************************************************************

double sanafe::schedule_messages_timestep_new(
        TimestepHandle &ts, Scheduler &scheduler)
{
    // TODO: HACK: Trying out a new model that better accounts for shared
    //  paths for different flows
    //
    // Schedule the global order of messages using a semi-analytical timing
    //  model. This sits in between the complexity of the simple and
    //  cycle-accurate model, capturing some but not all network effects. This
    //  takes a vector containing a list of messages per core, and scheduler
    //  parameters (mostly NoC configuration parameters). Returns the timestamp
    //  of the last scheduled event, i.e., the total time-step delay
    MessagePriorityQueue priority;
    NewNocInfo noc(scheduler);
    double last_timestamp = 0.0;
    Timestep &ts_data = ts.get();

    auto messages_sent_per_core = schedule_init_message_queues_new(ts_data, noc);
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
        noc_update_all_tracked_messages_new(m.sent_timestamp, noc);

        if (!m.placeholder)
        {
            schedule_handle_message_new(m, scheduler, noc);
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

void sanafe::schedule_handle_message_new(
        Message &m, Scheduler &scheduler, NewNocInfo &noc)
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
    m.messages_along_route = noc.calculate_route_congestion(m, scheduler);
    const auto path_capacity =
            static_cast<double>((m.hops + 1UL) * scheduler.buffer_size);

    if (m.messages_along_route > path_capacity)
    {
        // Use heuristic for estimating delay based on route congestion,
        //  path capacity in messages and the mean delay per message
        m.blocked_delay = (m.messages_along_route - path_capacity) *
                noc.mean_in_flight_receive_delay;
        m.sent_timestamp += m.blocked_delay;
    }
    else
    {
        m.blocked_delay = 0.0; // Path isn't at capacity; no blocking
    }
    INFO("mid:%zu path_capacity:%lf messages_along_route:%lf blocked:%e\n", m.mid,
            path_capacity, m.messages_along_route, m.blocked_delay);

    // HACK: TODO: this is my major change to algorithm. Also need to develop this
    //  further to consider when there are multiple destinations handling
    //  messages along the route i.e. its no longer a simple product but could
    //  take into account the parallel servicing of messages by different
    //  receiving cores.. maybe this would move logic into the route congestion
    //  calculation
    // TODO: we may need this additional check
    const double max_messages_within_path =
            std::min(static_cast<double>((m.hops + 1) * scheduler.buffer_size),
                    m.messages_along_route);
     const double congestion_delay = max_messages_within_path *
             noc.mean_in_flight_receive_delay;
    //const double congestion_delay = 0.0;
    TRACE1(SCHEDULER, "Path capacity:%lf messages:%lf congestion delay:%e\n",
            path_capacity, m.messages_along_route, congestion_delay);

    // Update the messages timestamps, both when the message is received and
    //  when the receiving core has finished processing it
    m.network_delay = std::max(m.min_hop_delay, congestion_delay);
    const double earliest_received_time = m.sent_timestamp + m.network_delay;
    m.received_timestamp = std::max(
            noc.core_finished_receiving.at(dest_core), earliest_received_time);

    noc.core_finished_receiving[dest_core] =
            std::max((noc.core_finished_receiving[dest_core] + m.receive_delay),
                    (earliest_received_time + m.receive_delay));
    m.processed_timestamp = noc.core_finished_receiving[dest_core];

    // Now, push the message into the right receiving queue. Calculate the
    //  network delay and when the message is received
    m.in_noc = true;
    noc.messages_received[dest_core].push_back(m);
    noc_update_message_tracking_new(m, noc, true);
}

std::vector<sanafe::MessageFifo> sanafe::schedule_init_message_queues_new(
        const Timestep &ts, NewNocInfo &noc)
{
    // const size_t total_links = noc.noc_height_in_tiles *
    //         noc.noc_width_in_tiles *
    //         (sanafe::ndirections + noc.max_cores_per_tile);
    //noc.message_density = std::vector<double>(total_links, 0.0);

    // Return a local copy of the per-core messages
    return ts.messages;
}

void sanafe::noc_update_all_tracked_messages_new(const double t, NewNocInfo &noc)
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
                noc_update_message_tracking_new(m, noc, false);
                TRACE1(SCHEDULER, "Removing message mid:%zu\n", m.mid);
                return true; // Remove this message
            }
            return false; // Keep this message
        });
    }
}

void sanafe::noc_update_message_tracking_new(
        const Message &m, NewNocInfo &noc, const bool entering_noc)
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

double sanafe::NewNocInfo::calculate_route_congestion(
        const Message &m, const Scheduler &scheduler) const
{
    // Currently tracks messages in the direct flow, and then also messages
    //  in overlapping flows (that may impact congestion)
    // TODO: we don't consider when multiple destination cores are processing
    //   messages, since this will affect the service rate, possibly significantly
    // TODO: we don't consider multiple subnetworks. Although this seems like
    //  a small detail, it has a large impact on the network latency, especially
    //  for small or large delays.
    //
    // Currently, simulations are much more accurate when modeling buffer size of
    //  8 rather than 16 (so basically ignoring 1 of the subnets, 1/2 of the total
    //  buffer capacity). However, the destination core Rx buffers should be twice
    //  as large, which I think would hide some of the latency for smaller delays
    //  (in some cases there would be less backpressure congestion in the NoC)
    //  However, I think effects have compounded since there should be many more
    //  messages buffered within the NoC, and this should ultimately halve the
    //  service rate (two subnets competing over destination cores). Maybe some
    //  inaccuracies are cancelling one another out though.
    const auto src_dest = std::make_pair(m.src_core_id, m.dest_core_id);

    double messages = 0.0;
    const bool flow_exists = flows.count(src_dest) > 0;
    if (flow_exists)
    {
        const auto &f_existing = flows.at(src_dest);
        messages = static_cast<double>(f_existing.messages_in_flight);

        // INFO("flow src:%zu dst:%zu messages:%lf\n", src_dest.first, src_dest.second, messages);
        // Work out any messages on links shared with this flow
        const auto total_links = m.hops + 2;
        for (auto &[shared_src_dest, shared_links_start_end_pair] :
                f_existing.flows_sharing_links)
        {
            const auto [optional_start_link, shared_links] =
                    shared_links_start_end_pair;
            const bool flow_sharing_links_exists =
                    flows.count(shared_src_dest) > 0;
            if (flow_sharing_links_exists && optional_start_link.has_value())
            {
                // If there is overlapping paths between the flows, figure out
                //  if there might be messages buffered along those shared links
                auto &flow_sharing_links = flows.at(shared_src_dest);
                const size_t maximum_messages_sharing_buffers =
                        shared_links * scheduler.buffer_size;

                // We need to track the last link that overlaps,
                //  because buffers will fill up from the back, last link first!
                const size_t last_shared_link = optional_start_link.value() + shared_links;

                const size_t minimum_buffered_for_contention =
                        (total_links - last_shared_link) * scheduler.buffer_size;
                if (flow_sharing_links.messages_in_flight >
                        minimum_buffered_for_contention)
                {
                    const size_t messages_sharing_buffers =
                            std::min(maximum_messages_sharing_buffers,
                                    flow_sharing_links.messages_in_flight -
                                            minimum_buffered_for_contention);
                    // INFO("shared flow src:%zu dst:%zu msg:%zu start_link:%zu shared_links:%zu end_link:%zu total_links:%zu shared_min:%zu shared_max:%zu sharing:%zu hops:%zu\n",
                    //         shared_src_dest.first, shared_src_dest.second, flow_sharing_links.messages_in_flight,
                    //         optional_start_link.value_or(-1), shared_links,
                    //         last_shared_link, total_links,
                    //         minimum_buffered_for_contention,
                    //         maximum_messages_sharing_buffers,
                    //         messages_sharing_buffers,
                    //         m.hops);
                    messages += messages_sharing_buffers;
                }
                else
                {
                    // INFO("shared flow src:%zu dst:%zu msg:%zu start_link:%zu shared_links:%zu end_link:%zu total_links:%zu shared_min:%zu shared_max:%zu sharing:%zu hops:%zu\n",
                    //         shared_src_dest.first, shared_src_dest.second,
                    //         flow_sharing_links.messages_in_flight,
                    //         optional_start_link.value_or(-1), shared_links,
                    //         last_shared_link, total_links,
                    //         minimum_buffered_for_contention,
                    //         maximum_messages_sharing_buffers, 0UL, m.hops);
                }
            }
            // TODO: removed the check limiting the messages in the path
        }
    }

    //INFO("src:%zu dest:%zu messages:%lf\n", src_dest.first, src_dest.second, messages);
    return messages;
}

std::pair<int, int> sanafe::get_path_xy_increments(
        const Flow &flow, const NewNocInfo &noc) noexcept
{
    const size_t src_tile = flow.src / noc.max_cores_per_tile;
    const size_t src_x = src_tile / noc.noc_height_in_tiles;
    const size_t src_y = src_tile % noc.noc_height_in_tiles;

    const size_t dest_tile = flow.dest / noc.max_cores_per_tile;
    const size_t dest_x = dest_tile / noc.noc_height_in_tiles;
    const size_t dest_y = dest_tile % noc.noc_height_in_tiles;

    int x_increment;
    if (src_x == dest_x)
    {
        x_increment = 0;
    }
    else
    {
        x_increment = (src_x < dest_x) ? 1 : -1;
    }

    int y_increment;
    if (src_y == dest_y)
    {
        y_increment = 0;
    }
    else
    {
        y_increment = (src_y < dest_y) ? 1 : -1;
    }

    return std::make_pair(x_increment, y_increment);
}

// ***************************************

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

std::vector<sanafe::MessageFifo> sanafe::schedule_init_message_queues(
        const Timestep &ts, NocInfo &noc)
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
    m.messages_along_route = noc.calculate_route_congestion(m);
    const auto path_capacity =
            static_cast<double>((m.hops + 1UL) * scheduler.buffer_size);

    if (m.messages_along_route > path_capacity)
    {
        // Use heuristic for estimating delay based on route congestion,
        //  path capacity in messages and the mean delay per message
        m.blocked_delay = (m.messages_along_route - path_capacity) *
                noc.mean_in_flight_receive_delay;
        m.sent_timestamp += m.blocked_delay;
    }
    else
    {
        m.blocked_delay = 0.0; // Path isn't at capacity; no blocking
    }

    const double congestion_delay = m.messages_along_route *
            noc.mean_in_flight_receive_delay /
            (static_cast<double>(m.hops) + 1.0);
    TRACE1(SCHEDULER, "Path capacity:%lf messages:%lf congestion delay:%e\n",
            path_capacity, m.messages_along_route, congestion_delay);

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

    auto [x_increment, y_increment] = get_path_xy_increments(message);
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
    auto [x_increment, y_increment] = get_path_xy_increments(m);
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

std::pair<int, int> sanafe::get_path_xy_increments(
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
            schedule_messages_timestep_new(ts, scheduler);
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