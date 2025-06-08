// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <sstream>

#include "arch.hpp"
#include "tile.hpp"

sanafe::Tile::Tile(const TileConfiguration &config)
        : name(config.name)
        , energy_north_hop(config.power_metrics.energy_north_hop)
        , latency_north_hop(config.power_metrics.latency_north_hop)
        , energy_east_hop(config.power_metrics.energy_east_hop)
        , latency_east_hop(config.power_metrics.latency_east_hop)
        , energy_south_hop(config.power_metrics.energy_south_hop)
        , latency_south_hop(config.power_metrics.latency_south_hop)
        , energy_west_hop(config.power_metrics.energy_west_hop)
        , latency_west_hop(config.power_metrics.latency_west_hop)
        , id(config.id)
        , x(config.x)
        , y(config.y)
        , log_energy(config.power_metrics.log_energy)
        , log_latency(config.power_metrics.log_latency)
{
}

std::string sanafe::Tile::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Tile(tile=" << id << " cores=";
    ss << cores.size() << ")";

    return ss.str();
}
