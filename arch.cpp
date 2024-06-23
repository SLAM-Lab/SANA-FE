// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// arch.cpp
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
//#include "pipeline.hpp"
#include "plugins.hpp"
#include "print.hpp"

sanafe::Architecture::Architecture(
        std::string name, const NetworkOnChipConfiguration &noc)
        : name(std::move(name))
        , noc_width(noc.width_in_tiles)
        , noc_height(noc.height_in_tiles)
        , noc_buffer_size(noc.link_buffer_size)
{
}

std::vector<std::reference_wrapper<sanafe::Core>> sanafe::Architecture::cores()
{
    std::vector<std::reference_wrapper<Core>> all_cores_in_arch;

    for (auto &tile : tiles)
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

/*
std::string sanafe::Architecture::description() const
{
    std::map<std::string, std::string> attributes;
    attributes["link_buffer_size"] = print_int(noc_buffer_size);
    attributes["topology"] = "mesh";
    attributes["width"] = print_int(noc_width);
    attributes["height"] = print_int(noc_height);

    std::ostringstream ss;
    ss << '@' << print_format_attributes(attributes) << std::endl;
    return ss.str();
}
*/

sanafe::Message::Message(
        const Architecture &arch, const Neuron &n, const long int timestep)
        : timestep(timestep)
        , src_neuron_id(n.id)
        , src_neuron_group_id(n.parent_group_id)
{
    // If no axon was given create a message with no destination. By
    //  default, messages without destinations act as a placeholder for neuron
    //  processing
    const Core &src_core = *(n.core);
    const Tile &src_tile = arch.tiles[src_core.parent_tile_id];
    src_x = src_tile.x;
    src_y = src_tile.y;
    src_tile_id = src_tile.id;
    src_core_id = src_core.id;
    src_core_offset = src_core.offset;
}

sanafe::Message::Message(const Architecture &arch, const Neuron &n,
        const long int timestep, const int axon_address)
        : Message(arch, n, timestep)
{
    const Core &src_core = *(n.core);
    const AxonOutModel &src_axon = src_core.axons_out[axon_address];
    const Tile &dest_tile = arch.tiles[src_axon.dest_tile_id];
    const Core &dest_core = dest_tile.cores[src_axon.dest_core_offset];
    const AxonInModel &dest_axon = dest_core.axons_in[src_axon.dest_axon_id];

    placeholder = false;
    spikes = dest_axon.synapse_addresses.size();
    dest_x = dest_tile.x;
    dest_y = dest_tile.y;
    dest_tile_id = dest_tile.id;
    dest_core_id = dest_core.id;
    dest_core_offset = dest_core.offset;
    dest_axon_id = src_axon.dest_axon_id;
    // TODO: support multiple axon output units, included in the synapse
    dest_axon_hw = 0;
}

sanafe::Tile::Tile(std::string name, const size_t tile_id,
        const TilePowerMetrics &power_metrics)
        : name(std::move(name))
        , energy_north_hop(power_metrics.energy_north_hop)
        , latency_north_hop(power_metrics.latency_north_hop)
        , energy_east_hop(power_metrics.energy_east_hop)
        , latency_east_hop(power_metrics.latency_east_hop)
        , energy_south_hop(power_metrics.energy_south_hop)
        , latency_south_hop(power_metrics.latency_south_hop)
        , energy_west_hop(power_metrics.energy_south_hop)
        , latency_west_hop(power_metrics.latency_east_hop)
        , id(tile_id)
{
}

std::string sanafe::Tile::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Tile(tile=" << id << " cores=";
    ss << cores.size() << ")";

    return ss.str();
}

/*
std::string sanafe::Tile::description() const
{
    std::map<std::string, std::string> attributes;
    attributes["energy_east"] = print_float(energy_east_hop);
    attributes["energy_west"] = print_float(energy_west_hop);
    attributes["energy_north"] = print_float(energy_north_hop);
    attributes["energy_south"] = print_float(energy_south_hop);
    attributes["latency_east"] = print_float(energy_east_hop);
    attributes["latency_west"] = print_float(energy_west_hop);
    attributes["latency_north"] = print_float(energy_north_hop);
    attributes["latency_south"] = print_float(energy_south_hop);

    std::ostringstream ss;
    ss << "t " << name << print_format_attributes(attributes) << std::endl;
    return ss.str();
}
*/

sanafe::Core::Core(std::string name, const CoreAddress &address,
        const CorePipelineConfiguration &pipeline)
        : pipeline_config(pipeline)
        , name(std::move(name))
        , id(address.id)
        , offset(address.offset_within_tile)
        , parent_tile_id(address.parent_tile_id)
{
}

sanafe::Tile &sanafe::Architecture::create_tile(
        std::string name, const TilePowerMetrics &power_metrics)
{
    // Initialize a new tile given metrics by the user, and push it into the
    //  Architecture's list of tiles
    const size_t new_tile_id = tiles.size();
    // The tile id is a unique global value that be used to index into the
    //  Architecture's tile array
    tiles.emplace_back(std::move(name), new_tile_id, power_metrics);
    Tile &new_tile = tiles[new_tile_id];
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

sanafe::AxonInUnit::AxonInUnit(std::string axon_in_name,
        const CoreAddress &parent_core, const AxonInPowerMetrics &power_metrics)
        : name(std::move(axon_in_name))
        , parent_core_address(parent_core)
        , energy_spike_message(power_metrics.energy_message_in)
        , latency_spike_message(power_metrics.latency_message_in)
{
}

sanafe::SynapseUnit::SynapseUnit(std::string synapse_name,
        const CoreAddress &parent_core,
        const SynapsePowerMetrics &power_metrics, const ModelInfo &model)
        : plugin_lib(model.plugin_library_path)
        , name(std::move(synapse_name))
        , model(model.name)
        , parent_core_address(parent_core)
        , energy_spike_op(power_metrics.energy_process_spike)
        , latency_spike_op(power_metrics.latency_process_spike)
{
}

sanafe::DendriteUnit::DendriteUnit(std::string dendrite_name,
        const CoreAddress &parent_core,
        const DendritePowerMetrics &power_metrics, const ModelInfo &model)
        : plugin_lib(model.plugin_library_path)
        , name(std::move(dendrite_name))
        , model(model.name)
        , parent_core_address(parent_core)
        , energy_access(power_metrics.energy_access)
        , latency_access(power_metrics.latency_access)
{
}

sanafe::SomaUnit::SomaUnit(std::string soma_name,
        const CoreAddress &parent_core, const SomaPowerMetrics &power_metrics,
        const ModelInfo &model)
        : plugin_lib(model.plugin_library_path)
        , name(std::move(soma_name))
        , model(model.name)
        , parent_core_address(parent_core)
        , energy_update_neuron(power_metrics.energy_update_neuron)
        , latency_update_neuron(power_metrics.latency_update_neuron)
        , energy_access_neuron(power_metrics.energy_access_neuron)
        , latency_access_neuron(power_metrics.latency_access_neuron)
        , energy_spiking(power_metrics.energy_spike_out)
        , latency_spiking(power_metrics.latency_spike_out)
{
}

sanafe::AxonOutUnit::AxonOutUnit(std::string axon_out_name,
        const CoreAddress &parent_core,
        const AxonOutPowerMetrics &power_metrics)
        : name(std::move(axon_out_name))
        , parent_core_address(parent_core)
        , energy_access(power_metrics.energy_message_out)
        , latency_access(power_metrics.latency_message_out)
{
}

sanafe::Core &sanafe::Architecture::create_core(std::string name,
        const size_t parent_tile_id,
        const CorePipelineConfiguration &pipeline_config)
{
    if (parent_tile_id > tiles.size())
    {
        throw std::invalid_argument("Error: Tile ID > total tiles");
    }
    // Lookup the parent tile to assign a new core to
    Tile &parent_tile = tiles[parent_tile_id];
    const size_t offset_within_tile = parent_tile.cores.size();
    const size_t new_core_id = core_count++;
    CoreAddress new_core_address = {
            parent_tile_id, offset_within_tile, new_core_id};

    // Initialize the new core and refer to it at both tile and arch levels
    parent_tile.cores.emplace_back(
            Core(std::move(name), new_core_address, pipeline_config));

    // The architecture tracks the maximum cores in *any* of its tiles.
    //  This information is needed later by the scheduler, when creating
    //  structures to track spike messages and congestion in the NoC
    max_cores_per_tile =
            std::max<size_t>(max_cores_per_tile, offset_within_tile + 1);
    TRACE1("Core created id:%zu.%zu.\n", parent_tile_id, new_core_id);
    Core &new_core = parent_tile.cores[offset_within_tile];

    return new_core;
}

std::string sanafe::Core::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Core(name= " << name << " tile=" << parent_tile_id << ")";
    return ss.str();
}

/*
std::string sanafe::Core::description() const
{
    std::ostringstream ss;
    ss << "c " << name << ' ' << parent_tile_id << std::endl;
    return ss.str();
}

std::string sanafe::AxonInUnit::description() const
{
    std::map<std::string, std::string> attributes;
    attributes["energy_message"] = print_float(energy_spike_message);
    attributes["latency_message"] = print_float(latency_spike_message);
    std::ostringstream ss;
    ss << "i " << name << ' ' << parent_core_address.parent_tile_id;
    ss << ' ' << parent_core_address.offset_within_tile;
    ss << print_format_attributes(attributes) << std::endl;
    return ss.str();
}

std::string sanafe::SynapseUnit::description() const
{
    std::map<std::string, std::string> attributes;
    attributes["energy_spike"] = print_float(energy_spike_op);
    attributes["latency_spike"] = print_float(latency_spike_op);
    attributes["model"] = model;
    std::ostringstream ss;
    ss << "s " << name << ' ' << parent_core_address.parent_tile_id << ' ';
    ss << parent_core_address.offset_within_tile;
    ss << print_format_attributes(attributes) << std::endl;
    return ss.str();
}

std::string sanafe::DendriteUnit::description() const
{
    std::map<std::string, std::string> attributes;
    attributes["model"] = model;
    attributes["energy"] = print_float(energy_access);
    attributes["latency"] = print_float(latency_access);
    std::ostringstream ss;
    ss << "d " << name << ' ' << parent_core_address.parent_tile_id;
    ss << ' ' << parent_core_address.offset_within_tile;
    ss << print_format_attributes(attributes) << std::endl;
    return ss.str();
}

std::string sanafe::SomaUnit::description() const
{
    std::map<std::string, std::string> attributes;
    attributes["model"] = model;
    attributes["energy_access_neuron"] = print_float(energy_access_neuron);
    attributes["latency_access_neuron"] = print_float(latency_access_neuron);
    attributes["energy_update_neuron"] = print_float(energy_update_neuron);
    attributes["latency_update_neuron"] = print_float(latency_update_neuron);
    attributes["energy_spike_out"] = print_float(energy_access_neuron);
    attributes["latency_spike_out"] = print_float(latency_access_neuron);
    std::ostringstream ss;
    ss << "+ " << name << ' ' << parent_core_address.parent_tile_id;
    ss << ' ' << parent_core_address.offset_within_tile;
    ss << print_format_attributes(attributes) << std::endl;
    return ss.str();
}

std::string sanafe::AxonOutUnit::description() const
{
    std::map<std::string, std::string> attributes;
    attributes["energy"] = print_float(energy_access);
    attributes["latency"] = print_float(latency_access);
    std::ostringstream ss;
    ss << "o " << name << ' ' << parent_core_address.parent_tile_id;
    ss << ' ' << parent_core_address.offset_within_tile;
    ss << print_format_attributes(attributes) << std::endl;
    return ss.str();
}
*/

sanafe::AxonInUnit &sanafe::Core::create_axon_in(
        const std::string &name, const AxonInPowerMetrics &power_metrics)
{
    const CoreAddress parent_core_address = {parent_tile_id, offset, id};
    axon_in_hw.emplace_back(
            AxonInUnit(name, parent_core_address, power_metrics));
    AxonInUnit &new_axon_in_hw_unit = axon_in_hw.back();

    return new_axon_in_hw_unit;
}

sanafe::SynapseUnit &sanafe::Core::create_synapse(const std::string &name,
        const SynapsePowerMetrics &power_metrics, const ModelInfo &model)
{
    const CoreAddress parent_core_address = {parent_tile_id, offset, id};
    synapse.emplace_back(
            SynapseUnit(name, parent_core_address, power_metrics, model));
    SynapseUnit &new_synapse_hw_unit = synapse.back();
    TRACE1("New synapse h/w unit created (cid:%d.%d)\n",
            parent_core_address.parent_tile_id,
            parent_core_address.offset_within_tile);

    return new_synapse_hw_unit;
}

sanafe::DendriteUnit &sanafe::Core::create_dendrite(const std::string &name,
        const DendritePowerMetrics &power_metrics, const ModelInfo &model)
{
    const CoreAddress parent_core_address = {parent_tile_id, offset, id};
    dendrite.emplace_back(
            DendriteUnit(name, parent_core_address, power_metrics, model));
    TRACE1("New dendrite h/w unit created (c:%d.%d)\n",
            parent_core_address.parent_tile_id,
            parent_core_address.offset_within_tile);

    return dendrite.back();
}

sanafe::SomaUnit &sanafe::Core::create_soma(std::string name,
        const SomaPowerMetrics &power_metrics, const ModelInfo &model)
{
    const CoreAddress parent_core_address = {parent_tile_id, offset, id};
    soma.emplace_back(SomaUnit(
            std::move(name), parent_core_address, power_metrics, model));
    TRACE1("New soma h/w unit created (c:%d.%d)\n",
            parent_core_address.parent_tile_id,
            parent_core_address.offset_within_tile);

    return soma.back();
}

sanafe::AxonOutUnit &sanafe::Core::create_axon_out(const std::string &name,
        const AxonOutPowerMetrics &power_metrics)
{
    const CoreAddress parent_core_address = {parent_tile_id, offset, id};
    axon_out_hw.emplace_back(AxonOutUnit(
            name, parent_core_address, power_metrics));
    TRACE1("New axon out h/w unit created (c:%d.%d)\n",
            parent_core_address.parent_tile_id,
            parent_core_address.offset_within_tile);

    return axon_out_hw.back();
}

void sanafe::arch_create_axons(Architecture &arch)
{
    TRACE1("Creating all connection maps.\n");
    for (Tile &tile : arch.tiles)
    {
        for (Core &c : tile.cores)
        {
            for (auto *n_ptr : c.neurons)
            {
                arch_map_neuron_connections(*n_ptr);
            }
        }
    }

    TRACE1("Finished creating connection maps.\n");
    arch_print_axon_summary(arch);
}

void sanafe::arch_print_axon_summary(Architecture &arch)
{
    int in_count = 0;
    int out_count = 0;

    INFO("** Mapping summary **\n");
    for (Tile &tile : arch.tiles)
    {
        // For debug only, print the axon maps
        for (Core &c : tile.cores)
        {
            bool core_used = false;
            for (size_t k = 0; k < c.neurons.size(); k++)
            {
#ifdef DEBUG
                Neuron *n = c->neurons[k];
                TRACE2("\tnid:%d.%d ", n->group->id, n->id);
                TRACE2("i:%d o:%d\n", n->maps_in_count, n->maps_out_count);
#endif
                core_used = true;
            }

            if (core_used)
            {
                // TODO: update these counts of axons
                //in_count += c.axon_in_hw.axons;
                //out_count += c.axon_out.map_count;
            }
        }
    }
    INFO("Total cores: %zu\n", arch.core_count);
    INFO("Average in map count: %lf\n",
            static_cast<double>(in_count) / arch.core_count);
    INFO("Average out map count: %lf\n",
            static_cast<double>(out_count) / arch.core_count);
}

void sanafe::arch_map_neuron_connections(Neuron &pre_neuron)
{
    // Setup the connections between neurons and map them to hardware
    assert(pre_neuron.core != nullptr);

    // Figure out the unique set of cores that this neuron broadcasts to
    TRACE1("Counting connections for neuron nid:%d\n", pre_neuron.id);
    std::set<Core *> cores_out;
    for (Connection &curr_connection : pre_neuron.connections_out)
    {
        TRACE1("Looking at connection id: %d\n", curr_connection.id);
        Core *dest_core = curr_connection.post_neuron->core;
        cores_out.insert(dest_core);
        TRACE1("Connected to dest core: %zu\n", dest_core->id);
    }

    TRACE1("Creating connections for neuron nid:%d to %zu core(s)\n",
            pre_neuron.id, cores_out.size());
    for (Core *dest_core : cores_out)
    {
        // Create the axon, and add it to both the destination and
        //  source cores
        arch_allocate_axon(pre_neuron, *dest_core);
    }
    TRACE3("Counted all maps for nid:%d count: %d\n", pre_neuron.id);

    for (Connection &curr_connection : pre_neuron.connections_out)
    {
        // Add every connection to the axon. Also link to the map in the
        //  post synaptic core / neuron
        Core &post_core = *(curr_connection.post_neuron->core);
        TRACE1("Adding connection:%d\n", curr_connection.id);
        arch_add_connection_to_axon(curr_connection, post_core);

        // Map to synapse hardware unit
        curr_connection.synapse_hw = &(post_core.synapse[0]);
        if (curr_connection.synapse_hw_name.length() > 0)
        {
            bool synapse_found = false;
            for (auto &synapse_hw : post_core.synapse)
            {
                if (curr_connection.synapse_hw_name == synapse_hw.name)
                {
                    curr_connection.synapse_hw = &synapse_hw;
                    synapse_found = true;
                }
            }
            if (!synapse_found)
            {
                INFO("Error: Could not map connection (hw:%s) "
                     "to any dendrite h/w.\n",
                        curr_connection.synapse_hw_name.c_str());
                throw std::runtime_error(
                        "Error: Could not map connection to synapse h/w");
            }
        }
        // Create the synapse model
        if (curr_connection.synapse_model == nullptr)
        {
            if (curr_connection.synapse_hw->plugin_lib.has_value())
            {
                const std::filesystem::path plugin_lib_path =
                        curr_connection.synapse_hw->plugin_lib.value();
                INFO("Creating synapse from plugin: %s.\n",
                        plugin_lib_path.c_str());
                curr_connection.synapse_model = plugin_get_synapse(
                        curr_connection.synapse_hw->model, plugin_lib_path);
            }
            else
            {
                // Use built in models
                TRACE1("Creating synapse built-in model %s.\n",
                        curr_connection.synapse_hw->model.c_str());
                curr_connection.synapse_model = sanafe::model_get_synapse(
                        curr_connection.synapse_hw->model);
            }
        }
        curr_connection.synapse_model->set_attributes(
                curr_connection.synapse_params);
    }
    TRACE1("Finished mapping connections to hardware for nid:%d.%d.\n",
            pre_neuron.parent_group_id, pre_neuron.id);
}

void sanafe::Core::map_neuron(Neuron &n)
{
    TRACE1("Mapping nid:%d to core: %zu\n", n.id, id);
    // Map the neuron to hardware units
    if (n.core != nullptr)
    {
        throw std::runtime_error(
                "Error: Neuron already mapped, dynamic remapping not "
                "supported.");
    }

    n.core = this;
    TRACE1("Mapping neuron %d to core %zu\n", n.id, id);
    neurons.push_back(&n);

    if (neurons.size() > pipeline_config.max_neurons_supported)
    {
        INFO("Error: Exceeded maximum neurons per core (%zu)",
                pipeline_config.max_neurons_supported);
        throw std::runtime_error("Error: Exceeded maximum neurons per core.");
    }

    // Map neuron model to dendrite and soma hardware units in this core.
    //  Search through all models implemented by this core and return the
    //  one that matches. If no dendrite / soma hardware is specified,
    //  default to the first one defined
    if (dendrite.empty())
    {
        INFO("Error: No dendrite units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No dendrite units defined");
    }
    n.dendrite_hw = &(dendrite[0]);
    if (n.dendrite_hw_name.length() > 0)
    {
        bool dendrite_found = false;
        for (auto &dendrite_hw : dendrite)
        {
            if (n.dendrite_hw_name == dendrite_hw.name)
            {
                n.dendrite_hw = &dendrite_hw;
                dendrite_found = true;
            }
        }
        if (!dendrite_found)
        {
            INFO("Error: Could not map neuron nid:%d (hw:%s) "
                 "to any dendrite h/w.\n",
                    n.id, n.dendrite_hw_name.c_str());
            throw std::runtime_error(
                    "Error: Could not map neuron to dendrite h/w");
        }
    }

    if (soma.empty())
    {
        INFO("Error: No soma units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No soma units defined");
    }
    n.soma_hw = &(soma[0]);
    if (n.soma_hw_name.length() > 0)
    {
        bool soma_found = false;
        for (auto &soma_hw : soma)
        {
            if (n.soma_hw_name == soma_hw.name)
            {
                n.soma_hw = &soma_hw;
                soma_found = true;
            }
        }
        if (!soma_found)
        {
            INFO("Error: Could not map neuron nid:%d (hw:%s) "
                 "to any soma h/w.\n",
                    n.id, n.soma_hw_name.c_str());
            throw std::runtime_error("Error: Could not map neuron to soma h/w");
        }
    }
    n.soma_hw->neuron_count++;

    // TODO: support multiple axon outputs
    if (axon_out_hw.empty())
    {
        INFO("Error: No axon out units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No axon out units defined");
    }
    n.axon_out_hw = &(axon_out_hw[0]);

    // Pass all the model specific arguments
    if (n.soma_model == nullptr)
    {
        // Setup the soma model
        TRACE1("Soma hw name: %s", soma_hw_name.c_str());
        assert(n.soma_hw != nullptr);
        if (n.soma_hw->plugin_lib.has_value())
        {
            // Use external plug-in
            const std::filesystem::path &plugin_lib_path =
                    n.soma_hw->plugin_lib.value();
            INFO("Creating soma from plugin %s.\n", plugin_lib_path.c_str());
            n.soma_model = plugin_get_soma(
                    n.soma_hw->model, n.parent_group_id, n.id, plugin_lib_path);
        }
        else
        {
            // Use built in models
            INFO("Creating soma built-in model %s.\n",
                    n.soma_hw->model.c_str());
            n.soma_model =
                    model_get_soma(n.soma_hw->model, n.parent_group_id, n.id);
        }
        const NeuronGroup &group = *(n.parent_net->groups[n.parent_group_id]);
        assert(n.soma_model != nullptr);
        // First set the group's default attribute values, and then
        //  any defined by the neuron
        n.soma_model->set_attributes(
                group.default_neuron_config.soma_model_params);
        n.soma_model->set_attributes(n.soma_model_params);
    }

    if (n.dendrite_model == nullptr)
    {
        // Setup the dendrite model
        TRACE1("Dendrite hw name: %s", dendrite_hw_name.c_str());
        if (n.dendrite_hw->plugin_lib.has_value())
        {
            const std::filesystem::path &plugin_lib_path =
                    n.dendrite_hw->plugin_lib.value();
            INFO("Creating dendrite from plugin %s.\n",
                    plugin_lib_path.c_str());
            n.dendrite_model =
                    plugin_get_dendrite(n.dendrite_hw->model, plugin_lib_path);
        }
        else
        {
            // Use built in models
            TRACE1("Creating dendrite built-in model %s.\n",
                    n.dendrite_hw->model.c_str());
            n.dendrite_model = model_get_dendrite(n.dendrite_hw->model);
        }
        const NeuronGroup &group = *(n.parent_net->groups[n.parent_group_id]);
        assert(n.dendrite_model != nullptr);
        // Global attributes for all compartments in the dendrite
        n.dendrite_model->set_attributes(
                group.default_neuron_config.dendrite_model_params);
        n.dendrite_model->set_attributes(n.dendrite_model_params);
    }
}

void sanafe::arch_allocate_axon(Neuron &pre_neuron, Core &post_core)
{
    // Create a new input axon at a receiving (destination) core
    //  Then create the output axon at the sending core. Finally
    //  update the presynaptic neuron and postsynaptic neuron

    Core &pre_core = *(pre_neuron.core);

    TRACE3("Adding connection to core.\n");
    // Allocate the axon and its connections at the post-synaptic core
    post_core.axons_in.emplace_back(AxonInModel());
    const size_t new_axon_in_address = post_core.axons_in.size() - 1;

    // Add the axon at the sending, pre-synaptic core
    TRACE1("Axon in address:%zu for core:%zu.%zu\n", new_axon_in_address,
            post_core.parent_tile_id, post_core.id);
    AxonOutModel out;
    out.dest_axon_id = new_axon_in_address;
    out.dest_core_offset = post_core.offset;
    out.dest_tile_id = post_core.parent_tile_id;
    out.src_neuron_id = pre_neuron.id;
    pre_core.axons_out.push_back(out);
    const size_t new_axon_out_address = pre_core.axons_out.size() - 1;

    // Then add the output axon to the sending pre-synaptic neuron
    pre_neuron.axon_out_addresses.push_back(new_axon_out_address);
    TRACE1("nid:%d.%d cid:%zu.%zu added one output axon address %zu.\n",
            pre_neuron.parent_group_id, pre_neuron.id, pre_core.parent_tile_id,
            pre_core.offset, new_axon_out_address);
}

void sanafe::arch_add_connection_to_axon(Connection &con, Core &post_core)
{
    // Add a given connection to the axon in the post-synaptic
    //  (destination) core
    TRACE3("Adding to connection to axon:%zu\n",
            post_core.axons_out.size() - 1);

    post_core.connections_in.push_back(&con);
    const size_t synapse_address = post_core.connections_in.size() - 1;

    // Access the most recently created axon in for the post-synaptic core
    AxonInModel &last_added_target_axon = post_core.axons_in.back();
    last_added_target_axon.synapse_addresses.push_back(synapse_address);

    // Map the connection to the synapse hardware. Default to the first defined
    //  hardware unit (there must be at least one hardware unit defined)
    const size_t default_hw_id = 0;
    con.synapse_hw = &(post_core.synapse[default_hw_id]);
    if (!con.synapse_hw_name.empty())
    {
        const auto mapped_hw = std::find_if(post_core.synapse.begin(),
                post_core.synapse.end(), [&](const SynapseUnit &syn) {
                    return syn.name == con.synapse_hw_name;
                });
        const bool mapped_successfully = (mapped_hw != post_core.synapse.end());
        if (!mapped_successfully)
        {
            throw std::runtime_error(
                    "Error: Could not map connection to synapse h/w.\n");
        }
    }
}

/*
void sanafe::Architecture::save_arch_description(
        const std::filesystem::path &path)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        INFO("Error: Couldn't open arch file to save to: %s\n", path.c_str());
        throw std::invalid_argument(
                "Error: Couldn't open arch file to save to.");
    }

    INFO("tile count:%zu\n", tiles.size());

    for (const Tile &tile : tiles)
    {
        out << tile.description();
        for (const Core &core : tile.cores)
        {
            out << core.description();
            for (const AxonInUnit &in : core.axon_in_hw)
            {
                out << in.description();
            }
            for (const SynapseUnit &s : core.synapse)
            {
                out << s.description();
            }
            for (const DendriteUnit &d : core.dendrite)
            {
                out << d.description();
            }
            for (const SomaUnit &s : core.soma)
            {
                out << s.description();
            }
            for (const AxonOutUnit &o : core.axon_out_hw)
            {
                out << o.description();
            }
        }
    }
    out << description();
}
*/
