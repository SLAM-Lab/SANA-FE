// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// tile.hpp
#ifndef TILE_HEADER_INCLUDED_
#define TILE_HEADER_INCLUDED_

#include <cstddef>
#include <string>
#include <vector>

#include "core.hpp"
#include "fwd.hpp"

namespace sanafe
{
class Tile
{
public:
    std::vector<Core> cores;
    std::string name;
    double energy{0.0};
    double energy_north_hop;
    double latency_north_hop;
    double energy_east_hop;
    double latency_east_hop;
    double energy_south_hop;
    double latency_south_hop;
    double energy_west_hop;
    double latency_west_hop;
    size_t hops{0UL};
    long int messages_received{0L};
    long int total_neurons_fired{0L};
    size_t north_hops{0UL};
    size_t east_hops{0UL};
    size_t south_hops{0UL};
    size_t west_hops{0UL};
    size_t id;
    size_t x{0};
    size_t y{0};
    bool log_energy{false};
    bool log_latency{false};

    explicit Tile(const TileConfiguration &config);
    [[nodiscard]] size_t get_id() const { return id; }
    [[nodiscard]] std::string info() const;
};
}

#endif
