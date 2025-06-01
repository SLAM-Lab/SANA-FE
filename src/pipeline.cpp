// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

#include "arch.hpp"
#include "chip.hpp"
#include "mapped.hpp"
#include "pipeline.hpp"
#include "print.hpp"

sanafe::PipelineUnit::PipelineUnit(const bool implements_synapse,
        const bool implements_dendrite, const bool implements_soma)
        : implements_synapse(implements_synapse)
        , implements_dendrite(implements_dendrite)
        , implements_soma(implements_soma)
{
    // When constructing the pipeline h/w unit, setup input/output interfaces.
    //  A hardware unit may implement synaptic, dendritic
    //  or somatic functionality, or a combination of the three. The
    //  functionality that a derived class implements will determine the correct
    //  input and output interfaces. The input interface is then selected based
    //  on the first supported functional unit (synapse/dendrite/soma), whereas
    //  the output interface is selected for the last supported functionality.
    //  The input and output interfaces are dynamically assigned using
    //  std::function references (similar conceptually to function pointers)
    //
    // This approach using std::function was chosen over other approaches e.g.,
    //  using a templated class or a class containing templated interfaces. The
    //  main reason for this is we require pipelines of multiple h/w units
    //  stored in a single container. This implies each pipeline unit should be
    //  the same type, which is not possible / very complex with a templated
    //  approach. Additionally, using templates would require us creating a
    //  different type for every possible combination. The std::function
    //  approach is simpler and probably more flexible for the future e.g., if
    //  we expand the pipeline.
    //
    // Note that supporting synaptic and somatic, but not dendritic operations
    //  is an invalid combination and will raise an exception
    TRACE1(CHIP, "Setting interfaces %d %d %d\n", implements_synapse,
            implements_dendrite, implements_soma);
    if (implements_synapse && implements_soma && !implements_dendrite)
    {
        throw std::logic_error(
                "Invalid pipeline configuration: h/w supports synapse and soma "
                "but not dendrite functionality. To fix this, either add this "
                "to the core's dendrite section, or remove from either the "
                "synapse or soma sections.");
    }

    // Set input interface
    if (implements_synapse)
    {
        process_input_fn = [this](auto &&ts, auto &&n, auto &&con,
                                   auto &&input) {
            return this->process_synapse_input(ts, n, con, input);
        };
    }
    else if (implements_dendrite)
    {
        process_input_fn = [this](auto &&ts, auto &&n, auto &&con,
                                   auto &&input) {
            return this->process_dendrite_input(ts, n, con, input);
        };
    }
    else if (implements_soma)
    {
        process_input_fn = [this](auto &&ts, auto &&n, auto &&con,
                                   auto &&input) {
            return this->process_soma_input(ts, n, con, input);
        };
    }
    else
    {
        throw std::logic_error(
                "H/w must implement at least one functional unit out of "
                "synapse/dendrite/soma");
    }

    // Set output interface
    if (implements_soma)
    {
        process_output_fn = [](auto &&n, auto &&con, auto &&output) {
            process_soma_output(n, con, output);
        };
    }
    else if (implements_dendrite)
    {
        process_output_fn = [](auto &&n, auto &&con, auto &&output) {
            process_dendrite_output(n, con, output);
        };
    }
    else if (implements_synapse)
    {
        process_output_fn = [](auto &&n, auto &&con, auto &&output) {
            process_synapse_output(n, con, output);
        };
    }
    // else the fallthrough case where nothing is implemented should have
    //  already been handled above
}

void sanafe::PipelineUnit::check_implemented(
        const bool check_implements_synapse,
        const bool check_implements_dendrite,
        const bool check_implements_soma) const
{
    if (check_implements_synapse != implements_synapse)
    {
        const std::string error = check_implements_synapse ?
                "Unit does not support synapse functionality, remove from arch "
                "description's 'synapse' section" :
                "Unit supports synapse functionality but was not added to "
                "'synapse' hardware section. Add it or pick a different model.";
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }

    if (check_implements_dendrite != implements_dendrite)
    {
        const std::string error = check_implements_dendrite ?
                "Unit does not support dendrite functionality, remove from arch "
                "description's 'dendrite' section" :
                "Unit supports dendrite functionality but was not added to "
                "'dendrite' hardware section. Add it or pick a different model.";
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }

    if (check_implements_soma != implements_soma)
    {
        const std::string error = check_implements_soma ?
                "Unit does not support soma functionality, remove from arch "
                "description's 'soma' section" :
                "Unit supports soma functionality but was not added to "
                "'soma' hardware section. Add it or pick a different model.";
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }

}

void sanafe::PipelineUnit::check_outputs(
        const MappedNeuron & /*n*/, const PipelineResult &result) const
{
    // Check the hw returns a valid value for the next unit to process
    if (implements_soma && result.status == INVALID_NEURON_STATE)
    {
        throw std::runtime_error("Soma output; should return valid "
                                 "neuron state.");
    }
    if (!implements_soma && (implements_synapse || implements_dendrite) &&
            !result.current.has_value())
    {
        throw std::runtime_error("Synaptic or dendritic output; should return "
                                 "synaptic/dendritic current");
    }
}

void sanafe::PipelineUnit::register_attributes(
        const std::set<std::string> &attribute_names)
{
    for (const auto &attr : attribute_names)
    {
        supported_attribute_names.insert(attr);
    }
}

void sanafe::PipelineUnit::check_attribute(const std::string attribute_name)
{
    if (supported_attribute_names.find(attribute_name) ==
            supported_attribute_names.end())
    {
        INFO("Warning: Attribute (%s) not supported by model: %s, will be "
             "ignored.\nEither remove this attribute from the SNN/Architecture "
             "description file, be more specific which hardware requires it "
             "(using synapse/dendrite/soma sections) or register the attribute "
             "in the constructor to suppress this warning.\n",
                attribute_name.c_str(), name.c_str());
    }
}

sanafe::PipelineResult sanafe::PipelineUnit::process_dendrite_input(
        Timestep & /*ts*/, MappedNeuron &n,
        std::optional<MappedConnection *> con, const PipelineResult &input)
{
    PipelineResult output{};

    std::optional<size_t> synapse_address{std::nullopt};
    if (con.has_value() && (con.value() != nullptr))
    {
        synapse_address = con.value()->synapse_address;
    }
    output = update(n.mapped_address, input.current, synapse_address);

    return output;
}

sanafe::PipelineResult sanafe::PipelineUnit::process_soma_input(
        Timestep & /*ts*/, MappedNeuron &n,
        std::optional<MappedConnection *> /*con*/, const PipelineResult &input)
{
    PipelineResult output = update(n.mapped_address, input.current);
    return output;
}

void sanafe::PipelineUnit::process_synapse_output(MappedNeuron & /*n*/,
        std::optional<MappedConnection *> con, PipelineResult &output)
{
    calculate_synapse_default_energy_latency(*(con.value_or(nullptr)), output);
}

void sanafe::PipelineUnit::process_dendrite_output(MappedNeuron &n,
        std::optional<MappedConnection *> /*con*/, PipelineResult &output)
{
    calculate_dendrite_default_energy_latency(n, output);
}

void sanafe::PipelineUnit::process_soma_output(MappedNeuron &n,
        std::optional<MappedConnection *> /*con*/, PipelineResult &output)
{
    calculate_soma_default_energy_latency(n, output);
    update_soma_activity(n, output);
}

sanafe::PipelineResult sanafe::PipelineUnit::process_synapse_input(
        Timestep & /*ts*/, MappedNeuron & /*n*/,
        std::optional<MappedConnection *> con, const PipelineResult & /*input*/)
{
    if (!con.has_value())
    {
        throw std::runtime_error(
                "Pipeline error, didn't receive "
                "synaptic connection info. Check that no h/w unit is being "
                "invoked before this one in the pipeline.");
    }
    PipelineResult output = update(con.value()->synapse_address, true);
    ++spikes_processed;

    return output;
}

void sanafe::PipelineUnit::calculate_synapse_default_energy_latency(
        MappedConnection &con, PipelineResult &simulation_result)
{
    const bool energy_simulated = simulation_result.energy.has_value();
    const bool latency_simulated = simulation_result.latency.has_value();

    const bool default_synapse_energy_metrics_set =
            con.synapse_hw->default_energy_process_spike.has_value();
    if (energy_simulated && default_synapse_energy_metrics_set)
    {
        const std::string error(
                "Synapse unit simulates energy and also has "
                "default energy metrics set.");
        throw std::runtime_error(error);
    }
    if (default_synapse_energy_metrics_set)
    {
        simulation_result.energy = con.synapse_hw->default_energy_process_spike;
    }

    const bool default_synapse_latency_metrics_set =
            con.synapse_hw->default_latency_process_spike.has_value();
    if (latency_simulated && default_synapse_latency_metrics_set)
    {
        const std::string error(
                "Synapse unit simulates latency and also has "
                "default latency metrics set. Remove the default "
                "metric from the architecture description.");
        throw std::runtime_error(error);
    }

    if (default_synapse_latency_metrics_set)
    {
        if (simulation_result.latency.has_value())
        {
            const std::string error(
                    "Synapse unit simulates latency and also has "
                    "default latency metrics set. Remove the default "
                    "metric from the architecture description.");
            throw std::runtime_error(error);
        }
        simulation_result.latency =
                con.synapse_hw->default_latency_process_spike;
    }

    if (!simulation_result.energy.has_value())
    {
        const std::string error(
                "Synapse unit does not simulate energy or provide "
                "a default energy cost in the architecture "
                "description.");
        throw std::runtime_error(error);
    }
    if (!simulation_result.latency.has_value())
    {
        const std::string error(
                "Synapse unit does not simulate latency or "
                "provide a default latency cost in the architecture "
                "description.");
        throw std::runtime_error(error);
    }
}

void sanafe::PipelineUnit::calculate_dendrite_default_energy_latency(
        MappedNeuron &n, PipelineResult &simulation_result)
{
    const bool energy_simulated = simulation_result.energy.has_value();
    const bool latency_simulated = simulation_result.latency.has_value();

    const bool default_dendrite_energy_metrics_set =
            n.dendrite_hw->default_energy_update.has_value();
    if (energy_simulated && default_dendrite_energy_metrics_set)
    {
        const std::string error(
                "Dendrite unit simulates energy and also has "
                "default energy metrics set.");
        throw std::runtime_error(error);
    }
    if (default_dendrite_energy_metrics_set)
    {
        simulation_result.energy = n.dendrite_hw->default_energy_update;
    }

    const bool default_dendrite_latency_metrics_set =
            n.dendrite_hw->default_latency_update.has_value();
    if (latency_simulated && default_dendrite_latency_metrics_set)
    {
        const std::string error(
                "Dendrite unit simulates latency and also has "
                "default latency metrics set.");
        throw std::runtime_error(error);
    }

    if (default_dendrite_latency_metrics_set)
    {
        if (simulation_result.latency.has_value())
        {
            const std::string error(
                    "Error: Dendrite unit simulates latency and also has "
                    "default latency metrics set.");
            throw std::runtime_error(error);
        }
        simulation_result.latency = n.dendrite_hw->default_latency_update;
    }

    if (!simulation_result.energy.has_value())
    {
        const std::string error(
                "Error: Dendrite unit does not simulate energy or provide "
                "a default energy cost in the architecture "
                "description.");
        throw std::runtime_error(error);
    }
    if (!simulation_result.latency.has_value())
    {
        const std::string error(
                "Error: Dendrite unit does not simulate latency or "
                "provide a default latency cost in the architecture "
                "description.");
        throw std::runtime_error(error);
    }
}

void sanafe::PipelineUnit::calculate_soma_default_energy_latency(
        MappedNeuron &n, PipelineResult &simulation_result)
{
    const bool energy_simulated = simulation_result.energy.has_value();
    const bool latency_simulated = simulation_result.latency.has_value();

    const bool soma_energy_metrics_set =
            n.soma_hw->default_soma_energy_metrics.has_value();
    if (energy_simulated && soma_energy_metrics_set)
    {
        const std::string error(
                "Error: Soma unit simulates energy and also has "
                "default energy metrics set. Remove the default energy metrics "
                "from the architecture description.");
        throw std::runtime_error(error);
    }
    if (soma_energy_metrics_set)
    {
        simulation_result.energy =
                n.soma_hw->default_soma_energy_metrics->energy_access_neuron;
    }

    const bool soma_latency_metrics_set =
            n.soma_hw->default_soma_latency_metrics.has_value();
    if (latency_simulated && soma_latency_metrics_set)
    {
        const std::string error(
                "Error: Soma unit simulates latency and also has "
                "default latency costs set. Remove the default latency metrics "
                "from the architecture description");
        throw std::runtime_error(error);
    }
    if (soma_latency_metrics_set)
    {
        simulation_result.latency =
                n.soma_hw->default_soma_latency_metrics->latency_access_neuron;
    }

    if ((simulation_result.status == sanafe::UPDATED) ||
            (simulation_result.status == sanafe::FIRED))
    {
        if (soma_energy_metrics_set)
        {
            simulation_result.energy.value() +=
                    n.soma_hw->default_soma_energy_metrics->energy_update_neuron;
        }
        if (soma_latency_metrics_set)
        {
            simulation_result.latency.value() +=
                    n.soma_hw->default_soma_latency_metrics
                            ->latency_update_neuron;
        }
    }
    if (simulation_result.status == sanafe::FIRED)
    {
        if (soma_energy_metrics_set)
        {
            simulation_result.energy.value() +=
                    n.soma_hw->default_soma_energy_metrics->energy_spike_out;
        }
        if (soma_latency_metrics_set)
        {
            simulation_result.latency.value() +=
                    n.soma_hw->default_soma_latency_metrics->latency_spike_out;
        }
    }

    if (!simulation_result.energy.has_value())
    {
        const std::string error(
                "Soma unit does not simulate energy or "
                "provide default energy costs in the architecture "
                "description.");
        throw std::runtime_error(error);
    }
    if (!simulation_result.latency.has_value())
    {
        const std::string error(
                "Soma unit does not simulate latency or "
                "provide default latency costs in the architecture "
                "description.");
        throw std::runtime_error(error);
    }
}

void sanafe::PipelineUnit::update_soma_activity(
        MappedNeuron &n, const PipelineResult &simulation_result)
{
    if ((simulation_result.status == sanafe::UPDATED) ||
            (simulation_result.status == sanafe::FIRED))
    {
        n.soma_hw->neurons_updated++;

        if (simulation_result.status == sanafe::FIRED)
        {
            n.soma_hw->neurons_fired++;
            TRACE1(CHIP, "Neuron %s.%zu fired\n", n.parent_group_name.c_str(),
                    n.id);
        }
    }
}

void sanafe::PipelineUnit::set_attributes(
        std::string unit_name, const ModelInfo &model)
{
    model_attributes = model.model_attributes;
    plugin_lib = model.plugin_library_path;
    name = std::move(unit_name);
    log_energy = model.log_energy;
    log_latency = model.log_latency;

    synapse_set_default_attributes();
    dendrite_set_default_attributes();
    soma_set_default_attributes();

    // Finally, forward all attributes from the architecture description to the
    //  model. This might be useful if you want to define any additional
    //  model-specific attributes here, e.g., fault-rate or maximum memory size.
    for (auto &[key, attribute] : model_attributes)
    {
        check_attribute(key);
        set_attribute_hw(key, attribute);
    }
}

void sanafe::PipelineUnit::synapse_set_default_attributes()
{
    if (model_attributes.find("energy_process_spike") != model_attributes.end())
    {
        default_energy_process_spike =
                static_cast<double>(model_attributes["energy_process_spike"]);
    }
    if (model_attributes.find("latency_process_spike") !=
            model_attributes.end())
    {
        default_latency_process_spike =
                static_cast<double>(model_attributes["latency_process_spike"]);
    }
}

void sanafe::PipelineUnit::dendrite_set_default_attributes()
{
    if (model_attributes.find("energy_update") != model_attributes.end())
    {
        default_energy_update =
                static_cast<double>(model_attributes["energy_update"]);
    }
    if (model_attributes.find("latency_update") != model_attributes.end())
    {
        default_latency_update =
                static_cast<double>(model_attributes["latency_update"]);
    }
}

void sanafe::PipelineUnit::soma_set_default_attributes()
{
    // Lambda for finding keys later
    auto key_exists = [this](const std::string &key) {
        return model_attributes.find(key) != model_attributes.end();
    };

    // Soma energy attributes
    const std::set<std::string> energy_metric_names{
            "energy_access_neuron", "energy_update_neuron", "energy_spike_out"};
    const bool parse_energy_metrics = std::any_of(
            energy_metric_names.begin(), energy_metric_names.end(), key_exists);
    if (parse_energy_metrics)
    {
        for (const auto &metric : energy_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Metric not defined: " + metric;
                INFO("Error: %s\n", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaEnergyMetrics energy_metrics;
        energy_metrics.energy_access_neuron =
                static_cast<double>(model_attributes["energy_access_neuron"]);
        energy_metrics.energy_update_neuron =
                static_cast<double>(model_attributes["energy_update_neuron"]);
        energy_metrics.energy_spike_out =
                static_cast<double>(model_attributes["energy_spike_out"]);
        default_soma_energy_metrics = energy_metrics;
    }

    // Soma latency attributes
    const std::set<std::string> latency_metric_names{"latency_access_neuron",
            "latency_update_neuron", "latency_spike_out"};
    const bool parse_latency_metrics = std::any_of(latency_metric_names.begin(),
            latency_metric_names.end(), key_exists);
    if (parse_latency_metrics)
    {
        for (const auto &metric : latency_metric_names)
        {
            if (!key_exists(metric))
            {
                const std::string error =
                        "Missing metric: " + metric;
                INFO("Error: %s\n", error.c_str());
                throw std::invalid_argument(error);
            }
        }
        SomaLatencyMetrics latency_metrics;
        latency_metrics.latency_access_neuron =
                static_cast<double>(model_attributes["latency_access_neuron"]);
        latency_metrics.latency_update_neuron =
                static_cast<double>(model_attributes["latency_update_neuron"]);
        latency_metrics.latency_spike_out =
                static_cast<double>(model_attributes["latency_spike_out"]);
        default_soma_latency_metrics = latency_metrics;
    }
}

sanafe::BufferPosition sanafe::pipeline_parse_buffer_pos_str(
        const std::string &buffer_pos_str, const bool buffer_inside_unit)
{
    BufferPosition buffer_pos{BUFFER_BEFORE_DENDRITE_UNIT};
    // H/w either has SANA-FE insert and manage a time-step buffer on the h/w
    //  unit's inputs (before the unit), or can assume that the unit manages its
    //  own state that is internally buffered over consecutive time-steps e.g.,
    //  using a double buffer

    if (buffer_pos_str == "dendrite")
    {
        if (buffer_inside_unit)
        {
            buffer_pos = BUFFER_INSIDE_DENDRITE_UNIT;
        }
        else
        {
            buffer_pos = BUFFER_BEFORE_DENDRITE_UNIT;
        }
    }
    else if (buffer_pos_str == "soma")
    {
        if (buffer_inside_unit)
        {
            buffer_pos = BUFFER_INSIDE_SOMA_UNIT;
        }
        else
        {
            buffer_pos = BUFFER_BEFORE_SOMA_UNIT;
        }
    }
    else if (buffer_pos_str == "axon_out")
    {
        buffer_pos = BUFFER_BEFORE_AXON_OUT_UNIT;
    }
    else
    {
        INFO("Error: Buffer position %s not supported", buffer_pos_str.c_str());
        throw std::invalid_argument("Error: Buffer position not supported");
    }

    return buffer_pos;
}
