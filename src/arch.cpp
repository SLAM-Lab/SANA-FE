// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.cpp
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem> // For std::filesystem::path
#include <fstream>
#include <functional> // For std::reference_wrapper
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

#include "arch.hpp"
#include "print.hpp"
#include "yaml_arch.hpp"

sanafe::Architecture::Architecture(
        std::string name, const NetworkOnChipConfiguration &noc)
        : ts_sync_delay_table(noc.ts_sync_delay_table)
        , name(std::move(name))
        , noc_width_in_tiles(noc.width_in_tiles)
        , noc_height_in_tiles(noc.height_in_tiles)
        , noc_buffer_size(noc.link_buffer_size)
        , timestep_delay(noc.timestep_delay)
{
}

std::vector<std::reference_wrapper<sanafe::CoreConfiguration>>
sanafe::Architecture::cores()
{
    std::vector<std::reference_wrapper<CoreConfiguration>> all_cores_in_arch;

    for (TileConfiguration &tile : tiles)
    {
        std::copy(tile.cores.begin(), tile.cores.end(),
                std::back_inserter(all_cores_in_arch));
    }

    return all_cores_in_arch;
}

std::string sanafe::Architecture::info() const noexcept
{
    std::ostringstream ss;
    ss << "sanafe::Architecture(tiles=" << tiles.size();
    ss << ", cores=" << core_count << ")";

    return ss.str();
}

sanafe::TileConfiguration::TileConfiguration(
        std::string name, const size_t id, const TilePowerMetrics &metrics)
        : power_metrics(metrics)
        , name(std::move(name))
        , id(id)
{
}

sanafe::CoreConfiguration::CoreConfiguration(std::string name,
        const CoreAddress &address, const CorePipelineConfiguration &pipeline)
        : pipeline(pipeline)
        , name(std::move(name))
        , address(address)
{
}

std::pair<int, int> sanafe::Architecture::calculate_tile_coordinates(
        const size_t tile_id) const
{
    // Map linear tile IDs to 2D coordinates for physical layout representation.
    //  This conversion assumes a row-major NoC grid arrangement where
    //  consecutive IDs are placed vertically before moving to the next column.
    const size_t x = tile_id / noc_height_in_tiles; // floor(id/height)
    const size_t y = tile_id % noc_height_in_tiles;
    assert(x < noc_width_in_tiles);
    return std::make_pair(x, y);
}

sanafe::TileConfiguration &sanafe::Architecture::create_tile(
        std::string name, const TilePowerMetrics &power_metrics)
{
    // Tiles are assigned sequential IDs based on creation order to ensure
    //  deterministic addressing in the network-on-chip topology
    const size_t new_tile_id = tiles.size();

    // Tile IDs serve as both array indices and as network addresses in the NoC
    //  topology, enabling O(1) tile lookups and efficient message routing
    tiles.emplace_back(std::move(name), new_tile_id, power_metrics);
    TileConfiguration &new_tile = tiles[new_tile_id];
    std::tie(new_tile.x, new_tile.y) = calculate_tile_coordinates(new_tile.id);

    return new_tile;
}

sanafe::Architecture sanafe::load_arch(const std::filesystem::path &path)
{
    std::ifstream arch_fp_stream(path);

    if (arch_fp_stream.fail())
    {
        throw std::system_error(std::make_error_code(std::errc::io_error),
                "Failed to open architecture file: " + path.string());
    }
    INFO("Loading architecture from file: %s\n", path.c_str());
    return description_parse_arch_file_yaml(arch_fp_stream);
}

sanafe::CoreConfiguration &sanafe::Architecture::create_core(std::string name,
        const size_t parent_tile_id,
        const CorePipelineConfiguration &pipeline_config)
{
    if (parent_tile_id >= tiles.size())
    {
        INFO("Error: Tile ID (%zu) out of range (>=%zu)\n", parent_tile_id,
                tiles.size());
        throw std::invalid_argument("Tile ID out of range");
    }
    // Cores must be attached to existing tiles representing on-chip hierarchy
    TileConfiguration &parent_tile = tiles.at(parent_tile_id);
    const size_t offset_within_tile = parent_tile.cores.size();
    const size_t new_core_id = core_count++;
    const CoreAddress new_core_address = {
            parent_tile_id, offset_within_tile, new_core_id};
    // Cores have dual referencing: within their parent tile's local space and
    //  globally within the architecture to support both local and cross-tile
    //  routing
    parent_tile.cores.emplace_back(
            std::move(name), new_core_address, pipeline_config);

    // The architecture tracks the maximum cores in *any* of its tiles.
    //  This information is needed later by the scheduler, when creating
    //  structures to track spike messages and congestion in the NoC
    max_cores_per_tile =
            std::max<size_t>(max_cores_per_tile, offset_within_tile + 1);
    TRACE1(ARCH, "Core created id:%zu.%zu.\n", parent_tile_id, new_core_id);
    CoreConfiguration &new_core = parent_tile.cores[offset_within_tile];

    return new_core;
}

sanafe::AxonInConfiguration &sanafe::CoreConfiguration::create_axon_in(
        std::string name, const AxonInPowerMetrics &power_metrics)
{
    axon_in.emplace_back(power_metrics, name);
    AxonInConfiguration &new_axon = axon_in.back();

    // Return a reference to the newly created component to allow caller to
    //  perform further configuration using method chaining
    return new_axon;
}

sanafe::PipelineUnitConfiguration &
sanafe::CoreConfiguration::create_hardware_unit(
        std::string name, const ModelInfo &model_details)
{
    pipeline_hw.emplace_back(model_details, name);
    PipelineUnitConfiguration &new_hw = pipeline_hw.back();

    return new_hw;
}

sanafe::AxonOutConfiguration &sanafe::CoreConfiguration::create_axon_out(
        std::string name, const AxonOutPowerMetrics &power_metrics)
{
    axon_out.emplace_back(power_metrics, name);
    AxonOutConfiguration &new_axon_out = axon_out.back();

    return new_axon_out;
}
