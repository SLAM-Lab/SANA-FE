// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// schedule.hpp: Schedule global order of messages on a neuromorphic chip
#ifndef TIMESTEP_HEADER_INCLUDED_
#define TIMESTEP_HEADER_INCLUDED_

#include <list>
#include <memory>
#include <vector>

#include "fwd.hpp"

namespace sanafe
{
constexpr long int invalid_timestep = -1L;

// Timestep data, including exchanged messages
struct Timestep
{
    std::vector<std::list<Message>> messages; // Per-sending core
    long int timestep{invalid_timestep};
    long int spike_count{0L};
    size_t total_hops{0UL};
    long int packets_sent{0L};
    long int neurons_updated{0L};
    long int neurons_fired{0L};

    double total_energy{0.0};
    double synapse_energy{0.0};
    double dendrite_energy{0.0};
    double soma_energy{0.0};
    double network_energy{0.0};
    double sim_time{0.0};
    double wall_time{0.0};

    Timestep() = default;
    Timestep(long int ts);
    void set_cores(size_t core_count);
};

// A wrapper class to efficiently share Timestep objects. The Timestep struct
//  contains info on runtime, energy and latency, but importantly also holds all
//  messages exchanged during the timestep. This may be in the order of
//  millions of messages or more. To avoid copying a huge amount of data, make
//  sure to store and pass handles, with the message data living on the heap.
class TimestepHandle
{
public:
    explicit TimestepHandle();
    explicit TimestepHandle(long int timestep_num);

    Timestep &operator*() { return *handle; }
    const Timestep &operator*() const { return *handle; }

    Timestep *operator->() { return handle.get(); }
    const Timestep *operator->() const { return handle.get(); }
    Timestep &get() { return *handle; }
    const Timestep &get() const { return *handle; }
    bool valid() const { return handle != nullptr; }

private:
    std::shared_ptr<Timestep> handle{nullptr};
};
}

#endif