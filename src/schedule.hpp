// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// schedule.hpp: Schedule global order of messages on a neuromorphic chip
#ifndef SCHEDULE_HEADER_INCLUDED_
#define SCHEDULE_HEADER_INCLUDED_

#include <condition_variable>
#include <functional> // For std::reference_wrapper
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include <booksim_config.hpp>

#include "chip.hpp"

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

using MessagePtr = std::shared_ptr<Message>;
using MessageFifo = std::list<MessagePtr>;

class CompareMessages
{
public:
    bool operator()(const MessagePtr &first, const MessagePtr &second) const { return first->sent_timestamp > second->sent_timestamp; }
};

using MessagePriorityQueue = std::priority_queue<MessagePtr, std::vector<MessagePtr>, CompareMessages>;


template <typename T, typename Container = std::vector<T>,
        typename Compare = std::less<T>>
class ThreadSafePriorityQueue
{
public:
    ThreadSafePriorityQueue()
            : pq_()
            , mutex_()
            , cv_()
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
        else
        {
            value = pq_.top();
            pq_.pop();
            // If queue is now empty, notify any threads waiting for empty state
            if (pq_.empty())
            {
                empty_cv_.notify_all();
            }

            return true;
        }
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
struct TimestepComparator
{
    bool operator()(const Timestep &lhs, const Timestep &rhs) const
    {
        return lhs.timestep > rhs.timestep;
    }
};

struct Scheduler
{
    std::vector<std::thread> scheduler_threads{};
    // TODO: allow for threads to be individually stopped, not just globally
    std::atomic<bool> should_stop{false};
    ThreadSafePriorityQueue<Timestep, std::vector<Timestep>, TimestepComparator> timesteps_to_schedule{};
    ThreadSafePriorityQueue<Timestep, std::vector<Timestep>, TimestepComparator> timesteps_to_write{};

    TimingModel timing_model{TIMING_MODEL_DETAILED};
    int noc_width;
    int noc_height;
    int buffer_size;
    int core_count;
    int max_cores_per_tile;
};

// NocInfo is used by the scheduler to track the high-level state of the NoC
//  at a given time
struct NocInfo
{
    std::vector<MessageFifo> messages_received;
    const size_t noc_width;
    const size_t noc_height;
    const size_t core_count;
    const size_t max_cores_per_tile;
    // Message density is the distribution of messages buffered in different
    //  links across the NoC. This vector is flattened row-major order
    //  where indexes are (x, y, links) and the links idx change fastest
    std::vector<double> message_density{};
    // The time that cores finish receiving their last message, i.e. the
    //  time at which they become free again to process more messages
    std::vector<double> core_finished_receiving{};
    double mean_in_flight_receive_delay{0.0};
    long int messages_in_noc{0L};

    NocInfo(int width, int height, int core_count, size_t max_cores_per_tile);
    [[nodiscard]] size_t idx(const size_t x, const size_t y, const size_t link) const
    {
        const size_t links_per_router = max_cores_per_tile + ndirections;
        return (x * noc_height * links_per_router) + (y * links_per_router) +
                link;
    }
};

MessagePriorityQueue schedule_init_timing_priority(std::vector<MessageFifo> &message_queues_per_core);
double schedule_messages_simple(Timestep &ts, const Scheduler &scheduler);
void schedule_messages_detailed(Timestep &ts, Scheduler &scheduler);
double schedule_messages_cycle_accurate(Timestep &ts, const BookSimConfig &config);


void schedule_messages_task(Scheduler &scheduler, const int tid);
void schedule_stop_all_threads(Scheduler &scheduler);

double schedule_messages_timestep(Timestep &ts, const Scheduler &scheduler);
void schedule_update_noc_message_counts(const Message &m, NocInfo &noc, bool message_in);
double schedule_calculate_messages_along_route(const Message &m, NocInfo &noc);
void schedule_update_noc(double t, NocInfo &noc);

// TODO: in the future, allow these strings to be setup dynamically. That way
//  the user can specify the NoC parameters within the architecture
//  description file.
constexpr const char booksim_config_str[][128] = {"topology = cmesh",
        "subnets = 2", "k = 8", "n = 2", "x = 8", "y = 4", "c = 4", "xr = 2",
        "yr = 2", "routing_function = dor_no_express", "use_noc_latency = 0",
        "num_vcs = 1", "vc_buf_size = 8", "wait_for_tail_credit = 1",
        "vc_allocator = islip", "sw_allocator = islip", "alloc_iters = 1",
        "credit_delay = 0", "routing_delay = 0", "vc_alloc_delay = 1",
        "sw_alloc_delay = 1", "input_speedup = 1", "output_speedup = 1",
        "internal_speedup = 1.0", "packet_size = 1", "injection_rate = 0.0",
        "clock_period = 1e-9"};
}

#endif