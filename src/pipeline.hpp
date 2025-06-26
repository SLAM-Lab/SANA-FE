// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef PIPELINE_HEADER_INCLUDED_
#define PIPELINE_HEADER_INCLUDED_

#include <cstddef>
#include <filesystem>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>

#include "attribute.hpp"
#include "fwd.hpp"
#include "mapped.hpp"

namespace sanafe
{
struct SomaEnergyMetrics
{
    double energy_update_neuron{0.0};
    double energy_access_neuron{0.0};
    double energy_spike_out{0.0};
};

struct SomaLatencyMetrics
{
    double latency_update_neuron{0.0};
    double latency_access_neuron{0.0};
    double latency_spike_out{0.0};
};

struct PipelineResult
{
    // Hardware outputs
    std::optional<double> current{std::nullopt};
    NeuronStatus status{invalid_neuron_state};
    // Optionally simulate energy and/or latency
    std::optional<double> energy{std::nullopt};
    std::optional<double> latency{std::nullopt};
};

class PipelineUnit
{
public:
    PipelineUnit(bool implements_synapse, bool implements_dendrite, bool implements_soma);
    PipelineUnit(const PipelineUnit &copy) = default;
    PipelineUnit(PipelineUnit &&other) = default;
    virtual ~PipelineUnit() = default;
    PipelineUnit &operator=(const PipelineUnit &other) = delete;
    PipelineUnit &operator=(PipelineUnit &&other) = delete;

    // Virtual member functions
    virtual void set_attribute_hw(const std::string &attribute_name, const ModelAttribute &param) = 0;
    virtual void set_attribute_neuron(size_t neuron_address, const std::string &attribute_name, const ModelAttribute &param) = 0;
    virtual void set_attribute_edge(size_t synapse_address, const std::string &attribute_name, const ModelAttribute &param) = 0;
    virtual void reset() = 0;

    // Optional virtual functions
    // The user of this class must implement the interfaces they wish to support
    //  Depending on whether you want to support Synapse, Dendrite, Soma
    //  operations, or a combination of the three.
    // If using synaptic inputs (address and read vs. synapse update without read)
    virtual PipelineResult update(
            size_t /*synapse_address*/, bool /*read*/ = false)
    {
        throw std::logic_error("Error: Synapse input not implemented");
    }
    // If using dendritic inputs (neuron address, synaptic current and synaptic address for additional info)
    virtual PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/,
            std::optional<size_t> /*synaptic_address*/)
    {
        throw std::logic_error("Error: Dendrite input not implemented");
    }
    // If using somatic inputs (address and current in)
    virtual PipelineResult update(
            size_t /*neuron_address*/, std::optional<double> /*current_in*/)
    {
        throw std::logic_error("Error: Soma input not implemented");
    }
    virtual void map_connection(MappedConnection &con) {}
    virtual double get_potential(size_t /*neuron_address*/)
    {
        return 0.0;
    }

    // Normal member functions and function pointers
    void set_time(long int timestep) { simulation_time = timestep; }
    void set_attributes(std::string unit_name, const ModelInfo &model);
    PipelineResult process(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    void register_attributes(const std::set<std::string> &attribute_names);
    void check_attribute(std::string attribute_name);
    void check_implemented(bool check_implements_synapse, bool check_implements_dendrite, bool check_implements_soma) const;

    // using InputInterfaceFunc = PipelineResult (PipelineUnit:: *)(Timestep &, MappedNeuron &, std::optional<MappedConnection*>, const PipelineResult &);
    // using OutputInterfaceFunc = void (*)(MappedNeuron &, std::optional<MappedConnection *>, PipelineResult &);

    using InputInterfaceFunc = std::function<PipelineResult(Timestep&, MappedNeuron&, std::optional<MappedConnection*>, const PipelineResult&)>;
    using OutputInterfaceFunc = std::function<void(MappedNeuron&, std::optional<MappedConnection*>, PipelineResult&)>;

    // Model information
    std::map<std::string, ModelAttribute> model_attributes;
    std::optional<std::filesystem::path> plugin_lib{std::nullopt};
    std::string name;
    std::string model;

    // Performance metrics
    std::optional<double> default_energy_process_spike{std::nullopt};
    std::optional<double> default_latency_process_spike{std::nullopt};
    std::optional<double> default_energy_update{std::nullopt};
    std::optional<double> default_latency_update{std::nullopt};
    std::optional<SomaEnergyMetrics> default_soma_energy_metrics;
    std::optional<SomaLatencyMetrics> default_soma_latency_metrics;
    double energy{0.0};
    double time{0.0};

    // Performance counters
    long int spikes_processed{0L};
    long int neurons_updated{0L};
    long int neurons_fired{0L};
    long int neuron_count{0L};

    // Warning counter (in case we get attributes we don't recognize)
    long int attribute_warnings{0L};
    static constexpr long int max_attribute_warnings{100L}; // Stops spam

    // Implementation flags, set whichever operations your derived unit supports
    //  to 'true'. Note that a hardware unit must support one or more of these
    bool implements_synapse;
    bool implements_dendrite;
    bool implements_soma;

    // Performance monitoring flags
    bool log_energy{false};
    bool log_latency{false};

protected:
    std::set<std::string> supported_attribute_names{
            "log_spikes",
            "log_potential",
            "force_update",
            "force_synapse_update",
            "force_dendrite_update",
            "force_soma_update",
            "synapse_hw_name",
            "dendrite_hw_name",
            "soma_hw_name",
            "model",
            "plugin",
            "energy_message_in",
            "latency_message_in",
            "energy_access_neuron",
            "latency_access_neuron",
            "energy_update_neuron",
            "latency_update_neuron",
            "energy_spike_out",
            "latency_spike_out",
            "energy_process_spike",
            "latency_process_spike",
            "energy_update",
            "latency_update",
            "energy_message_out",
            "latency_message_out",
            // Legacy attributes (v1)
            "log_v",
            "connections_out",
    };
    InputInterfaceFunc process_input_fn;
    OutputInterfaceFunc process_output_fn;
    long int simulation_time{0L};

    PipelineResult process_synapse_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_dendrite_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_soma_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);

    static void process_synapse_output(MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);
    static void process_dendrite_output(MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);
    static void process_soma_output(MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);

private:
    static void calculate_synapse_default_energy_latency(MappedConnection &con, PipelineResult &simulation_result);
    static void calculate_dendrite_default_energy_latency(MappedNeuron &n, PipelineResult &simulation_result);
    static void calculate_soma_default_energy_latency(MappedNeuron &n, PipelineResult &simulation_result);
    static void update_soma_activity(MappedNeuron &n, const PipelineResult &simulation_result);
    void check_outputs(const MappedNeuron &n, const PipelineResult &result) const;

    void synapse_set_default_attributes();
    void dendrite_set_default_attributes();
    void soma_set_default_attributes();
};

// Specific unit base classes, for the normal use case where the model
//  implements only synaptic, dendritic or somatic functionality (and not a)
//  combination
class SynapseUnit : public PipelineUnit
{
public:
    SynapseUnit() : PipelineUnit(true, false, false) {}
    PipelineResult update(size_t synapse_address, bool read = false) override = 0;
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/,
            std::optional<size_t> /*synaptic_address*/) final
    {
        throw std::logic_error(
                "Error: Synapse H/W called with dendrite inputs");
    }
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/) final
    {
        throw std::logic_error("Error: Synapse H/W called with soma inputs");
    }
    void set_attribute_neuron(size_t neuron_address,  const std::string &attribute_name, const ModelAttribute &param) final {};
};

class DendriteUnit : public PipelineUnit
{
public:
    DendriteUnit() : PipelineUnit(false, true, false) {}
    PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) override = 0;
    PipelineResult update(
            size_t /*synapse_address*/, bool /*read*/ = false) final
    {
        throw std::logic_error(
                "Error: Dendrite H/W called with synapse inputs");
    }
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/) final
    {
        throw std::logic_error("Error: Dendrite H/W called with soma inputs");
    }
};

class SomaUnit : public PipelineUnit
{
public:
    SomaUnit() : PipelineUnit(false, false, true) {}
    PipelineResult update(size_t neuron_address, std::optional<double> current_in) override = 0;
    PipelineResult update(
            size_t /*synapse_address*/, bool /*read*/ = false) final
    {
        throw std::logic_error("Error: Soma H/W called with synapse inputs");
    }
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/,
            std::optional<size_t> /*synaptic_address*/) final
    {
        throw std::logic_error("Error: Soma H/W called with dendrite inputs");
    }
    void set_attribute_edge(size_t synapse_address, const std::string &attribute_name, const ModelAttribute &param) final {};
    void map_connection(MappedConnection &con) final {};
};

BufferPosition pipeline_parse_buffer_pos_str(
        const std::string &buffer_pos_str, bool buffer_inside_unit);
}

// Include these member function definition in-line rather than in pipeline.cpp
//  to avoid linkage issues with plug-ins. Generally, header-only approaches
//  are useful in this situation. In this case, I've inlined only public and
//  protected members. Private members can be defined externally.
inline sanafe::PipelineUnit::PipelineUnit(const bool implements_synapse,
            const bool implements_dendrite, const bool implements_soma)
            : implements_synapse(implements_synapse)
            , implements_dendrite(implements_dendrite)
            , implements_soma(implements_soma)
{
    // When constructing the pipeline h/w unit, setup input/output interfaces.
    //  A hardware unit may implement synaptic, dendritic
    //  or somatic functionality, or a combination of the three. The
    //  functionality that a derived class implements will determine the
    //  correct input and output interfaces. The input interface is then
    //  selected based on the first supported functional unit
    //  (synapse/dendrite/soma), whereas the output interface is selected
    //  for the last supported functionality. The input and output
    //  interfaces are dynamically assigned using std::function references
    //  (similar conceptually to function pointers)
    //
    // This approach using std::function was chosen over other approaches
    //  e.g., using a templated class or a class containing templated
    //  interfaces. The main reason for this is we require pipelines of
    //  multiple h/w units stored in a single container. This implies each
    //  pipeline unit should be the same type, which is not possible / very
    //  complex with a templated approach. Additionally, using templates
    //  would require us creating a different type for every possible
    //  combination. The std::function approach is simpler and probably more
    //  flexible for the future e.g., if we expand the pipeline.
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

inline void sanafe::PipelineUnit::process_synapse_output(MappedNeuron & /*n*/,
        std::optional<MappedConnection *> con, PipelineResult &output)
{
    calculate_synapse_default_energy_latency(*(con.value_or(nullptr)), output);
}

inline void sanafe::PipelineUnit::process_dendrite_output(MappedNeuron &n,
        std::optional<MappedConnection *> /*con*/, PipelineResult &output)
{
    calculate_dendrite_default_energy_latency(n, output);
}

inline void sanafe::PipelineUnit::process_soma_output(MappedNeuron &n,
        std::optional<MappedConnection *> /*con*/, PipelineResult &output)
{
    calculate_soma_default_energy_latency(n, output);
    update_soma_activity(n, output);
}

inline sanafe::PipelineResult sanafe::PipelineUnit::process_synapse_input(
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

inline sanafe::PipelineResult sanafe::PipelineUnit::process_dendrite_input(
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

inline sanafe::PipelineResult sanafe::PipelineUnit::process_soma_input(
        Timestep & /*ts*/, MappedNeuron &n,
        std::optional<MappedConnection *> /*con*/, const PipelineResult &input)
{
    PipelineResult output = update(n.mapped_address, input.current);
    return output;
}


inline void sanafe::PipelineUnit::calculate_synapse_default_energy_latency(
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

inline void sanafe::PipelineUnit::calculate_dendrite_default_energy_latency(
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

inline void sanafe::PipelineUnit::calculate_soma_default_energy_latency(
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

    if ((simulation_result.status == sanafe::updated) ||
            (simulation_result.status == sanafe::fired))
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
    if (simulation_result.status == sanafe::fired)
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

inline void sanafe::PipelineUnit::update_soma_activity(
        MappedNeuron &n, const PipelineResult &simulation_result)
{
    if ((simulation_result.status == sanafe::updated) ||
            (simulation_result.status == sanafe::fired))
    {
        n.soma_hw->neurons_updated++;

        if (simulation_result.status == sanafe::fired)
        {
            n.soma_hw->neurons_fired++;
            TRACE1(CHIP, "Neuron %s.%zu fired\n", n.parent_group_name.c_str(),
                    n.id);
        }
    }
}

#endif
