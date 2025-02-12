// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.cpp
#include <algorithm>
#include <any>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem> // For std::filesystem::path
#include <fstream>
#include <functional> // For std::reference_wrapper
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <sstream>

#include "arch.hpp"
#include "description.hpp"
#include "models.hpp"
#include "network.hpp"
#include "plugins.hpp"
#include "print.hpp"
#include "chip.hpp"

sanafe::Architecture::Architecture(
        std::string name, const NetworkOnChipConfiguration &noc)
        : name(std::move(name))
        , noc_width(noc.width_in_tiles)
        , noc_height(noc.height_in_tiles)
        , noc_buffer_size(noc.link_buffer_size)
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

std::string sanafe::Architecture::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Architecture(tiles=" << tiles.size();
    ss << ", cores=" << core_count << ")";

    return ss.str();
}

sanafe::TileConfiguration::TileConfiguration(
        std::string name, const size_t id, const TilePowerMetrics &metrics)
        : power_metrics(metrics)
        , name(name)
        , id(id)
{
}

sanafe::CoreConfiguration::CoreConfiguration(std::string name,
        const CoreAddress &address,
        const CorePipelineConfiguration &pipeline)
        : pipeline(pipeline)
        , name(name)
        , address(address)
{
}

sanafe::TileConfiguration &sanafe::Architecture::create_tile(
        std::string name, const TilePowerMetrics &power_metrics)
{
    // Initialize a new tile given metrics by the user, and push it into the
    //  Architecture's list of tiles
    const size_t new_tile_id = tiles.size();
    // The tile id is a unique global value that be used to index into the
    //  Architecture's tile array
    tiles.emplace_back(std::move(name), new_tile_id, power_metrics);
    TileConfiguration &new_tile = tiles[new_tile_id];
    new_tile.x = new_tile.id / noc_height;
    new_tile.y = new_tile.id % noc_height;

    return new_tile;
}

sanafe::Architecture sanafe::load_arch(const std::filesystem::path &path)
{
    std::ifstream arch_fp(path);

    if (arch_fp.fail())
    {
        const std::string error = "Error: Architecture file: " +
                std::string(path) + "failed to open.";
        throw std::invalid_argument(error);
    }
    INFO("Loading architecture from file: %s\n", path.c_str());
    Architecture arch = description_parse_arch_file_yaml(arch_fp);
    arch_fp.close();

    return arch;
}

sanafe::CoreConfiguration &sanafe::Architecture::create_core(std::string name,
        const size_t parent_tile_id,
        const CorePipelineConfiguration &pipeline_config)
{
    if (parent_tile_id > tiles.size())
    {
        throw std::invalid_argument("Error: Tile ID > total tiles");
    }
    // Lookup the parent tile to assign a new core to
    TileConfiguration &parent_tile = tiles[parent_tile_id];
    const size_t offset_within_tile = parent_tile.cores.size();
    const size_t new_core_id = core_count++;
    CoreAddress new_core_address = {
            parent_tile_id, offset_within_tile, new_core_id};

    // Initialize the new core and refer to it at both tile and arch levels
    parent_tile.cores.emplace_back(CoreConfiguration(
            std::move(name), new_core_address, pipeline_config));

    // The architecture tracks the maximum cores in *any* of its tiles.
    //  This information is needed later by the scheduler, when creating
    //  structures to track spike messages and congestion in the NoC
    max_cores_per_tile =
            std::max<size_t>(max_cores_per_tile, offset_within_tile + 1);
    TRACE1(ARCH, "Core created id:%zu.%zu.\n", parent_tile_id, new_core_id);
    CoreConfiguration &new_core = parent_tile.cores[offset_within_tile];

    return new_core;
}

std::string sanafe::Core::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Core(name= " << name << " tile=" << parent_tile_id << ")";
    return ss.str();
}

sanafe::AxonInConfiguration &sanafe::CoreConfiguration::create_axon_in(
        const std::string &name, const AxonInPowerMetrics &power_metrics)
{
    axon_in.emplace_back(power_metrics, name);
    AxonInConfiguration &new_axon = axon_in.back();

    return new_axon;
}

sanafe::PipelineUnitConfiguration &sanafe::CoreConfiguration::create_hw(
        const std::string &name, const ModelInfo &model_details)
{
    pipeline_hw.emplace_back(model_details, name);
    PipelineUnitConfiguration &new_hw = pipeline_hw.back();

    return new_hw;
}

sanafe::AxonOutConfiguration &sanafe::CoreConfiguration::create_axon_out(
        const std::string &name, const AxonOutPowerMetrics &power_metrics)
{
    axon_out.emplace_back(power_metrics, name);
    AxonOutConfiguration &new_axon_out = axon_out.back();

    return new_axon_out;
}
