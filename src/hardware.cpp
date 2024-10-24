
#include "hardware.hpp"
#include "models.hpp"
#include "plugins.hpp"

sanafe::AxonInUnit::AxonInUnit(const AxonInConfiguration &config)
        : name(config.name)
        , energy_spike_message(config.metrics.energy_message_in)
        , latency_spike_message(config.metrics.latency_message_in)
{
}

void sanafe::SynapseUnit::configure(
        std::string synapse_name, const ModelInfo &model)
{
    model_parameters = model.model_parameters;
    plugin_lib = model.plugin_library_path;
    name = synapse_name;

    if (model_parameters.find("energy_process_spike") != model_parameters.end())
    {
        default_energy_process_spike =
                static_cast<double>(model_parameters["energy_process_spike"]);
    }
    if (model_parameters.find("latency_process_spike") !=
            model_parameters.end())
    {
        default_latency_process_spike =
                static_cast<double>(model_parameters["latency_process_spike"]);
    }
}

void sanafe::DendriteUnit::configure(
        std::string dendrite_name, const ModelInfo &model_details)
{
    model_parameters = model_details.model_parameters;
    plugin_lib = model_details.plugin_library_path;
    name = dendrite_name;
    model = model_details.name;
    if (model_parameters.find("energy_update") != model_parameters.end())
    {
        default_energy_update =
                static_cast<double>(model_parameters["energy_update"]);
    }
    if (model_parameters.find("latency_update") != model_parameters.end())
    {
        default_latency_update =
                static_cast<double>(model_parameters["latency_update"]);
    }
}

void sanafe::SomaUnit::configure(
        const std::string &soma_name, const ModelInfo &model_details)
{
    model_parameters = model_details.model_parameters;
    plugin_lib = model_details.plugin_library_path;
    name = soma_name;
    model = model_details.name;
    auto key_exists = [this](const std::string &key) {
        return model_parameters.find(key) != model_parameters.end();
    };

    const std::set<std::string> energy_metric_names{
            "energy_access_neuron", "energy_update_neuron", "energy_spike_out"};
    bool parse_energy_metrics = std::any_of(
            energy_metric_names.begin(), energy_metric_names.end(), key_exists);
    if (parse_energy_metrics)
    {
        for (const auto &metric : energy_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Error: Metric not defined: " + metric + "\n";
                INFO("%s", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaEnergyMetrics energy_metrics;
        energy_metrics.energy_access_neuron =
                static_cast<double>(model_parameters["energy_access_neuron"]);
        energy_metrics.energy_update_neuron =
                static_cast<double>(model_parameters["energy_update_neuron"]);
        energy_metrics.energy_spike_out =
                static_cast<double>(model_parameters["energy_spike_out"]);
        default_energy_metrics = energy_metrics;
    }

    const std::set<std::string> latency_metric_names{"latency_access_neuron",
            "latency_update_neuron", "latency_spike_out"};
    bool parse_latency_metrics = std::any_of(latency_metric_names.begin(),
            latency_metric_names.end(), key_exists);
    if (parse_latency_metrics)
    {
        for (const auto &metric : latency_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Error: Missing metric: " + metric + "\n";
                INFO("%s", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaLatencyMetrics latency_metrics;
        latency_metrics.latency_access_neuron =
                static_cast<double>(model_parameters["latency_access_neuron"]);
        latency_metrics.latency_update_neuron =
                static_cast<double>(model_parameters["latency_update_neuron"]);
        latency_metrics.latency_spike_out =
                static_cast<double>(model_parameters["latency_spike_out"]);
        default_latency_metrics = latency_metrics;
    }
}

sanafe::AxonOutUnit::AxonOutUnit(const AxonOutConfiguration &config)
        : name(std::move(config.name))
        , energy_access(config.metrics.energy_message_out)
        , latency_access(config.metrics.latency_message_out)
{
}

sanafe::MappedNeuron &sanafe::Core::map_neuron(const Neuron &neuron)
{
    // TODO: we need to insert the neuron according to its mapped order...
    // TODO: we need to do this stuff last... we need to insert all neurons first
    //  the sort on mapped order
    TRACE1(HW, "Mapping nid:%s.%zu to core: %zu\n",
            neuron.parent_group_id.c_str(), neuron.id, id);

    if (neurons.size() >= pipeline_config.max_neurons_supported)
    {
        INFO("Error: Exceeded maximum neurons per core (%zu)",
                pipeline_config.max_neurons_supported);
        throw std::runtime_error("Error: Exceeded maximum neurons per core.");
    }

    // Map the neuron to hardware units
    neurons.emplace_back(neuron.id, neuron.mapping_order);
    MappedNeuron &mapped = neurons.back();

    mapped.parent_group_name = neuron.parent_group_id;
    mapped.core = this;

    // Map neuron model to dendrite and soma hardware units in this core.
    //  Search through all models implemented by this core and return the
    //  one that matches. If no dendrite / soma hardware is specified,
    //  default to the first one defined
    if (dendrite.empty())
    {
        INFO("Error: No dendrite units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No dendrite units defined");
    }
    mapped.dendrite_hw = dendrite[0].get();
    if (neuron.dendrite_hw_name.length() > 0)
    {
        bool dendrite_found = false;
        for (auto &dendrite_hw : dendrite)
        {
            if (neuron.dendrite_hw_name == dendrite_hw->name)
            {
                mapped.dendrite_hw = dendrite_hw.get();
                dendrite_found = true;
            }
        }
        if (!dendrite_found)
        {
            INFO("Error: Could not map neuron nid:%zu (hw:%s) "
                 "to any dendrite h/w.\n",
                    neuron.id, neuron.dendrite_hw_name.c_str());
            throw std::runtime_error(
                    "Error: Could not map neuron to dendrite h/w");
        }
    }

    if (soma.empty())
    {
        INFO("Error: No soma units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No soma units defined");
    }
    mapped.soma_hw = soma[0].get();
    if (neuron.soma_hw_name.length() > 0)
    {
        bool soma_found = false;
        for (auto &soma_hw : soma)
        {
            if (neuron.soma_hw_name == soma_hw->name)
            {
                mapped.soma_hw = soma_hw.get();
                soma_found = true;
            }
        }
        if (!soma_found)
        {
            INFO("Error: Could not map neuron nid:%zu (hw:%s) "
                 "to any soma h/w.\n",
                    neuron.id, neuron.soma_hw_name.c_str());
            throw std::runtime_error("Error: Could not map neuron to soma h/w");
        }
    }

    if (axon_out_hw.empty())
    {
        INFO("Error: No axon out units defined for cid:%zu\n", id);
        throw std::runtime_error("Error: No axon out units defined");
    }
    mapped.axon_out_hw = &(axon_out_hw[0]);
    mapped.soma_hw->neuron_count++;

    return mapped;
}

sanafe::AxonInUnit &sanafe::Core::create_axon_in(
        const AxonInConfiguration &config)
{
    axon_in_hw.emplace_back(config);
    TRACE1(HW, "New axon in h/w unit created (%zu.%zu)\n", parent_tile_id, id);

    return axon_in_hw.back();
}

sanafe::SynapseUnit &sanafe::Core::create_synapse(
        const SynapseConfiguration &config)
{
    // Create the synapse model
    if (config.model_info.plugin_library_path.has_value())
    {
        const std::filesystem::path plugin_lib_path =
                config.model_info.plugin_library_path.value();
        TRACE1(HW, "Creating synapse from plugin: %s.\n",
                plugin_lib_path.c_str());
        synapse.emplace_back(
                plugin_get_synapse(config.model_info.name, plugin_lib_path));
    }
    else
    {
        // Use built in models
        TRACE1(HW, "Creating synapse built-in model %s.\n",
                config.model_info.name.c_str());
        synapse.emplace_back(model_get_synapse(config.model_info.name));
    }

    auto &new_unit = synapse.back();
    new_unit->configure(config.name, config.model_info);
    TRACE1(HW, "New synapse h/w unit created\n");

    return *new_unit;
}

sanafe::DendriteUnit &sanafe::Core::create_dendrite(
        const DendriteConfiguration &config)
{
    TRACE1(HW, "New dendrite h/w unit created\n");

    if (config.model_info.plugin_library_path.has_value())
    {
        const std::filesystem::path &plugin_library_path =
                config.model_info.plugin_library_path.value();
        TRACE1(HW, "Creating dendrite from plugin %s.\n",
                plugin_library_path.c_str());
        dendrite.emplace_back(plugin_get_dendrite(
                config.model_info.name, plugin_library_path));
    }
    else
    {
        // Use built in models
        TRACE1(HW, "Creating dendrite built-in model %s.\n",
                config.model_info.name.c_str());
        dendrite.emplace_back(model_get_dendrite(config.model_info.name));
    }

    auto &unit = dendrite.back();
    unit->configure(config.name, config.model_info);
    return *unit;
}

sanafe::SomaUnit &sanafe::Core::create_soma(const SomaConfiguration &config)
{
    TRACE1(HW, "New soma h/w unit created (%s)\n", config.name.c_str());

    if (config.model_info.plugin_library_path.has_value())
    {
        // Use external plug-in
        auto &plugin_library_path =
                config.model_info.plugin_library_path.value();
        TRACE1(HW, "Creating soma from plugin %s.\n",
                plugin_library_path.c_str());
        soma.emplace_back(
                plugin_get_soma(config.model_info.name, plugin_library_path));
    }
    else
    {
        // Use built-in unit models
        TRACE1(HW, "Creating soma built-in model %s.\n",
                config.model_info.name.c_str());
        soma.emplace_back(model_get_soma(config.model_info.name));
    }
    auto &unit = soma.back();
    unit->configure(config.name, config.model_info);

    return *unit;
}

sanafe::AxonOutUnit &sanafe::Core::create_axon_out(
        const AxonOutConfiguration &config)
{
    axon_out_hw.emplace_back(config);
    TRACE1(HW, "New axon out h/w unit created: (%zu.%zu)\n", parent_tile_id,
            id);

    return axon_out_hw.back();
}

sanafe::Message::Message(const SpikingHardware &hw, const MappedNeuron &n,
        const long int timestep)
        : timestep(timestep)
        , src_neuron_id(n.id)
        , src_neuron_group_id(n.parent_group_name)
{
    // If no axon was given create a message with no destination. By
    //  default, messages without destinations act as a placeholder for neuron
    //  processing
    const Core &src_core = *(n.core);
    const Tile &src_tile = hw.tiles[src_core.parent_tile_id];
    src_x = src_tile.x;
    src_y = src_tile.y;
    src_tile_id = src_tile.id;
    src_core_id = src_core.id;
    src_core_offset = src_core.offset;
}

sanafe::Message::Message(const SpikingHardware &hw, const MappedNeuron &n,
        const long int timestep, const int axon_address)
        : Message(hw, n, timestep)
{
    const Core &src_core = *(n.core);
    const AxonOutModel &src_axon = src_core.axons_out[axon_address];
    const Tile &dest_tile = hw.tiles[src_axon.dest_tile_id];
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
{
}

std::string sanafe::Tile::info() const
{
    std::ostringstream ss;
    ss << "sanafe::Tile(tile=" << id << " cores=";
    ss << cores.size() << ")";

    return ss.str();
}

sanafe::Core::Core(const CoreConfiguration &config)
        : pipeline_config(config.pipeline)
        , name(std::move(config.name))
        , id(config.address.id)
        , offset(config.address.offset_within_tile)
        , parent_tile_id(config.address.parent_tile_id)
{
}

sanafe::MappedConnection::MappedConnection(const int connection_id)
        : post_neuron(nullptr)
        , pre_neuron(nullptr)
        , synapse_hw(nullptr)
        , id(connection_id)
{
}

void sanafe::MappedNeuron::set_attributes(const NeuronTemplate &attributes)
{
    if (attributes.log_spikes.has_value())
    {
        log_spikes = attributes.log_spikes.value();
    }
    if (attributes.log_potential.has_value())
    {
        log_potential = attributes.log_potential.value();
    }
    if (attributes.force_dendrite_update)
    {
        force_dendrite_update = attributes.force_dendrite_update.value();
    }
    if (attributes.force_soma_update.has_value())
    {
        force_soma_update = attributes.force_soma_update.value();
    }
    if (attributes.force_synapse_update)
    {
        force_synapse_update = attributes.force_synapse_update.value();
    }

    for (auto &[key, param] : attributes.model_parameters)
    {
        TRACE2(HW, "Forwarding param: %s (dendrite:%d soma:%d)\n", key.c_str(),
                param.forward_to_dendrite, param.forward_to_soma);
        if (param.forward_to_dendrite && (dendrite_hw != nullptr))
        {
            dendrite_hw->set_attribute(mapped_address, key, param);
        }
        if (param.forward_to_soma && (soma_hw != nullptr))
        {
            soma_hw->set_attribute(mapped_address, key, param);
        }
    }
}
