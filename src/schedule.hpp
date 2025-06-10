// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// schedule.hpp: Schedule global order of messages on a neuromorphic chip
#ifndef SCHEDULE_HEADER_INCLUDED_
#define SCHEDULE_HEADER_INCLUDED_

#include <array>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <functional> // For std::reference_wrapper
#include <list>
#include <mutex>
#include <queue>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <booksim_config.hpp>

#include "chip.hpp"
#include "fwd.hpp"
#include "timestep.hpp"

namespace sanafe
{

enum Direction : uint8_t
{
    // Router link direction
    north = 0U,
    east = 1U,
    south = 2U,
    west = 3U,
    ndirections = 4U,
};

using MessageFifo = std::list<Message>;
using MessagePriorityQueue = std::priority_queue<Message, std::vector<Message>, CompareMessagesBySentTime>;

template <typename T, typename Container = std::vector<T>,
        typename Compare = std::less<T>>
class ThreadSafePriorityQueue
{
public:
    ThreadSafePriorityQueue()
            : pq_()
    {
    }

    // Wait until the queue is empty or terminated
    void wait_until_empty()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        empty_cv_.wait(lock, [this] { return pq_.empty() || terminated_; });
    }

    template <typename... Args> void push(Args &&...args)
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            pq_.push(std::forward<Args>(args)...);
        }
        cv_.notify_one();
    }

    bool pop(T &value)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this] { return !pq_.empty() || terminated_; });
        if (pq_.empty())
        {
            return false;
        }

        value = pq_.top();
        pq_.pop();
        // If queue is now empty, notify any threads waiting for empty state
        if (pq_.empty())
        {
            empty_cv_.notify_all();
        }

        return true;
    }

    bool try_pop(T &value)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (pq_.empty())
        {
            return false;
        }
        value = pq_.top();
        pq_.pop();
        if (pq_.empty())
        {
            empty_cv_.notify_all();
        }

        return true;
    }

    bool empty() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return pq_.empty();
    }

    size_t size() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return pq_.size();
    }

    void set_terminate()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            terminated_ = true;
        }
        cv_.notify_all();
    }

private:
    std::priority_queue<T, Container, Compare> pq_;
    mutable std::mutex mutex_;
    std::condition_variable cv_;
    std::condition_variable empty_cv_;
    bool terminated_ = false;
};

// For comparing timesteps when maintaining priority queue of timesteps
struct CompareTimesteps
{
    bool operator()(const TimestepHandle &lhs, const TimestepHandle &rhs) const
    {
        return lhs->timestep > rhs->timestep;
    }
};

struct Scheduler
{
    std::vector<std::thread> scheduler_threads;
    std::atomic<bool> should_stop{false};
    ThreadSafePriorityQueue<TimestepHandle, std::vector<TimestepHandle>,
            CompareTimesteps>
            timesteps_to_schedule;
    ThreadSafePriorityQueue<TimestepHandle, std::vector<TimestepHandle>,
            CompareTimesteps>
            timesteps_to_write;

    TimingModel timing_model{timing_model_detailed};
    size_t noc_width_in_tiles;
    size_t noc_height_in_tiles;
    size_t buffer_size;
    size_t core_count;
    size_t max_cores_per_tile;
    double timestep_sync_delay;
};

// NocInfo is used by the scheduler to track the high-level state of the NoC
//  at a given time
class NocInfo
{
public:
    std::vector<MessageFifo> messages_received;
    size_t noc_width_in_tiles;
    size_t noc_height_in_tiles;
    size_t core_count;
    size_t max_cores_per_tile;
    // Message density is the distribution of messages buffered in different
    //  links across the NoC. This vector is flattened row-major order
    //  where indexes are (x, y, links) and the links idx change fastest
    std::vector<double> message_density;
    // The time that cores finish receiving their last message, i.e. the
    //  time at which they become free again to process more messages
    std::vector<double> core_finished_receiving;
    double mean_in_flight_receive_delay{0.0};
    long int messages_in_noc{0L};

    NocInfo(const Scheduler &scheduler);
    [[nodiscard]] size_t idx(const size_t x, const size_t y, const size_t link) const
    {
        const size_t links_per_router = max_cores_per_tile + ndirections;
        return (x * noc_height_in_tiles * links_per_router) + (y * links_per_router) +
                link;
    }

    void update_message_density(const Message &message, bool entering_noc);
    void update_rolling_averages(const Message& message, bool entering_noc);
    [[nodiscard]] double calculate_route_congestion(
            const Message &message) const;

private:
    static std::pair<int, int> get_route_xy_increments(const Message &m) noexcept;
};

MessagePriorityQueue schedule_init_timing_priority(std::vector<MessageFifo> &message_queues_per_core);
void schedule_messages(TimestepHandle &timestep_handle, Scheduler &scheduler, const BookSimConfig &booksim_config);
void schedule_messages_simple(TimestepHandle &timestep_handle, Scheduler &scheduler);
void schedule_messages_detailed(TimestepHandle &timestep_handle, Scheduler &scheduler);
void schedule_messages_cycle_accurate(TimestepHandle &timestep_handle, const BookSimConfig &config, Scheduler &scheduler);

void schedule_create_threads(Scheduler &scheduler, int scheduler_thread_count);
void schedule_messages_thread(Scheduler &scheduler, int thread_id);
void schedule_stop_all_threads(Scheduler &scheduler);

std::vector<MessageFifo> schedule_init_message_queues(const Timestep &ts, NocInfo &noc);
double schedule_messages_timestep(TimestepHandle &timestep_handle, Scheduler &scheduler);
void schedule_handle_message(Message &m, Scheduler &scheduler, NocInfo &noc);
double schedule_push_next_message(std::vector<MessageFifo> &messages_sent_per_core, MessagePriorityQueue &priority, const Message &current_message);
void noc_update_message_tracking(const Message &m, NocInfo &noc, bool entering_noc);
double schedule_calculate_messages_along_route(const Message &m, NocInfo &noc);
void noc_update_all_tracked_messages(double t, NocInfo &noc);

// TODO: in the future, allow these strings to be setup dynamically. That way
//  the user can specify the NoC parameters within the architecture
//  description file.
constexpr std::array<std::string_view, 27> booksim_config_str = {{
        "topology = cmesh",
        "subnets = 2",
        "k = 8",
        "n = 2",
        "x = 8",
        "y = 4",
        "c = 4",
        "xr = 2",
        "yr = 2",
        "routing_function = dor_no_express",
        "use_noc_latency = 0",
        "num_vcs = 1",
        "vc_buf_size = 8",
        "wait_for_tail_credit = 1",
        "vc_allocator = islip",
        "sw_allocator = islip",
        "alloc_iters = 1",
        "credit_delay = 0",
        "routing_delay = 0",
        "vc_alloc_delay = 1",
        "sw_alloc_delay = 1",
        "input_speedup = 1",
        "output_speedup = 1",
        "internal_speedup = 1.0",
        "packet_size = 1",
        "injection_rate = 0.0",
        "clock_period = 1e-9"
}};

}

#endif