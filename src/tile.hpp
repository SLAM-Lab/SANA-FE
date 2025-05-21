// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// tile.hpp
#ifndef TILE_HEADER_INCLUDED_
#define TILE_HEADER_INCLUDED_

#include <string>
#include <vector>

#include "core.hpp"
#include "fwd.hpp"

namespace sanafe
{
class Tile
{
public:
    std::vector<Core> cores{};
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
    long int hops{0L};
    long int messages_received{0L};
    long int total_neurons_fired{0L};
    long int north_hops{0L};
    long int east_hops{0L};
    long int south_hops{0L};
    long int west_hops{0L};
    size_t id;
    size_t x{0};
    size_t y{0};
    bool log_energy{false};
    bool log_latency{false};

    explicit Tile(const TileConfiguration &config);
    [[nodiscard]] int get_id() const { return id; }
    [[nodiscard]] std::string info() const;
};
}

#endif
