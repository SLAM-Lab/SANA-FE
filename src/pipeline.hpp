// Copyright (c) 2026 - The University of Texas at Austin
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
#include <unordered_map>

#include "attribute.hpp"
#include "fwd.hpp"
#include "mapped.hpp"
#include "timestep.hpp"  // For Timestep reference

namespace sanafe
{
enum HardwareBitfield : uint8_t
{
    implements_invalid = 0U,
    implements_synapse = 1U << 0,
    implements_dendrite = 1U << 1,
    implements_soma = 1U << 2
};

inline HardwareBitfield operator|(HardwareBitfield a, HardwareBitfield b)
{
    return static_cast<HardwareBitfield>(
        static_cast<uint8_t>(a) | static_cast<uint8_t>(b));
}

inline HardwareBitfield operator&(HardwareBitfield a, HardwareBitfield b)
{
    return static_cast<HardwareBitfield>(
        static_cast<uint8_t>(a) & static_cast<uint8_t>(b));
}

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
    NeuronStatus status{neuron_state_unset};
    // Optionally simulate energy and/or latency
    std::optional<double> energy{std::nullopt};
    std::optional<double> latency{std::nullopt};
};

class PipelineUnit
{
public:
    PipelineUnit(const HardwareBitfield implemented);
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
    virtual PipelineResult update(size_t /*synapse_address*/,
            bool /*read*/, long int /*timestep*/)
    {
        throw std::logic_error("Error: Synapse input not implemented");
    }
    // If using dendritic inputs (neuron address, synaptic current and synaptic address for additional info)
    virtual PipelineResult update(
            size_t /*neuron_address*/, std::optional<double> /*current_in*/,
            std::optional<size_t> /*synaptic_address*/, long int /*timestep*/)
    {
        throw std::logic_error("Error: Dendrite input not implemented");
    }
    // If using somatic inputs (address and current in)
    virtual PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/, long int /*timestep*/)
    {
        throw std::logic_error("Error: Soma input not implemented");
    }
    virtual void track_connection(size_t /*synapse_address*/, size_t /*src_neuron_id*/, size_t /*dest_neuron_id*/) {}
    virtual double get_potential(size_t /*neuron_address*/)
    {
        return 0.0;
    }
    virtual std::map<std::string, double> get_neuron_traces(size_t /*neuron_address*/)
    {
        // Default is to have no neuron traces in addition to potentials
        return {};
    }

    // Normal member functions and function pointers
    void set_attributes_hw(std::string unit_name, const ModelInfo &model_info);
    PipelineResult process(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    void register_attributes(const std::set<std::string> &attribute_names);
    void register_attributes(const std::unordered_map<std::string, std::string> &attributes_with_descriptions);
    bool check_attribute(const std::string &attribute_name);
    std::vector<std::string> get_attributes() const;
    std::string get_attribute_help(const std::string &attribute_name) const;

    void check_implemented(bool check_implements_synapse, bool check_implements_dendrite, bool check_implements_soma) const;
    size_t add_neuron();
    size_t add_connection(MappedConnection &con);

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
    double latency{0.0};

    // Performance counters
    long int spikes_processed{0L};
    long int neurons_updated{0L};
    long int neurons_fired{0L};

    // Mapped neuron and connection counts
    long int neuron_count{0L};
    long int connection_count{0L};

    // Warning counter (in case we get attributes we don't recognize)
    long int attribute_warnings{0L};
    static constexpr long int max_attribute_warnings{0L}; // Disabled by default

    // Implementation flags, set whichever operations your derived unit supports
    //  to 'true'. Note that a hardware unit must support one or more of these
    bool implements_synapse;
    bool implements_dendrite;
    bool implements_soma;

    // Performance monitoring flags
    bool log_energy{false};
    bool log_latency{false};

    // Track whether the hardware is used by any mapped neurons or connections
    bool is_used{false};

    // Track whether to update the hardware every time-step, regardless of
    //  whether it received input from the network/preceding hardware unit.
    //  This may be needed in some cases e.g., for an analog hardware circuit
    //  that has a small leakage.
    bool update_every_timestep{false};

    static inline const std::unordered_map<std::string, std::string> framework_attributes{
            {"force_update", "(bool) Force updates every time-step."},
            {"synapse_hw_name", "(str) Unique name of the synapse H/W unit."},
            {"dendrite_hw_name", "(str) Unique name of the dendrite H/W unit."},
            {"soma_hw_name", "(str) Unique name of the soma H/W unit."},
            {"model", "(str) Unique model name, either built-in or plugin."},
            {"plugin", "(str) Plug-in library path." },
            {"energy_message_in", "(float) Energy cost of receiving a spike message (J)."},
            {"latency_message_in", "(float) Latency cost of receiving a spike message (s)."},
            {"energy_access_neuron", "(float) Energy cost for a soma to access a neuron (J)."},
            {"latency_access_neuron", "(float) Latency cost for a soma to access a neuron (s)."},
            {"energy_update_neuron", "(float) Energy cost for a soma to update (J)."},
            {"latency_update_neuron", "(float) Energy cost for a soma to update (s)."},
            {"energy_spike_out", "(float) Energy cost for a soma to spike (J)."},
            {"latency_spike_out", "(float) Latency cost for a soma to spike (s)."},
            {"energy_process_spike", "(float) Energy cost for one synapse look-up/access (J)."},
            {"latency_process_spike", "(float) Latency cost for one synapse look-up/access (s)."},
            {"energy_update", "(float) Energy cost of updating a dendrite (s)"},
            {"latency_update", "(float) Latency cost of updating a dendrite (s)"},
            {"energy_message_out", "(float) Energy cost of sending a spike message (J)"},
            {"latency_message_out", "(float) Latency cost of sending a spike message (s)"},
            // Legacy attributes (v1)
            {"connections_out", "(int) Connections outgoing from a neuron (deprecated)"},
    };

protected:
    std::unordered_map<std::string, std::string> supported_attributes{framework_attributes};

    InputInterfaceFunc process_input_fn;
    OutputInterfaceFunc process_output_fn;

    PipelineResult process_synapse_input(const Timestep &ts, const MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_dendrite_input(const Timestep &ts, const MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_soma_input(const Timestep &ts, const MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);

    static void process_synapse_output(const MappedNeuron &/*n*/, std::optional<MappedConnection *> con, PipelineResult &output);
    static void process_dendrite_output(const MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);
    static void process_soma_output(MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);

private:
    static void calculate_synapse_default_energy_latency(const MappedConnection &con, PipelineResult &simulation_result);
    static void calculate_dendrite_default_energy_latency(const MappedNeuron &n, PipelineResult &simulation_result);
    static void calculate_soma_default_energy_latency(const MappedNeuron &n, PipelineResult &simulation_result);
    static void update_soma_activity(MappedNeuron &n, const PipelineResult &simulation_result);
    void check_outputs(const MappedNeuron &n, const PipelineResult &result) const;

    void synapse_set_default_attributes();
    void dendrite_set_default_attributes();
    void soma_set_default_attributes();

    //long int simulation_time{0L};
};

// Specific unit base classes, for the normal use case where the model
//  implements only synaptic, dendritic or somatic functionality (and not a)
//  combination
class SynapseUnit : public PipelineUnit
{
public:
    SynapseUnit() : PipelineUnit(HardwareBitfield::implements_synapse) {}
    PipelineResult update(size_t synapse_address,
            bool read, long int /*timestep*/) override = 0;
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/,
            std::optional<size_t> /*synaptic_address*/,
            long int /*timestep*/) final
    {
        throw std::logic_error(
                "Error: Synapse H/W called with dendrite inputs");
    }
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/, long int /*timestep*/) final
    {
        throw std::logic_error("Error: Synapse H/W called with soma inputs");
    }
    void set_attribute_neuron(size_t /*neuron_address*/,  const std::string &/*attribute_name*/, const ModelAttribute &/*param*/) final {};
};

class DendriteUnit : public PipelineUnit
{
public:
    DendriteUnit()
            : PipelineUnit(HardwareBitfield::implements_dendrite)
    {
    }
    PipelineResult update(size_t neuron_address,
            std::optional<double> current_in,
            std::optional<size_t> synaptic_address,
            long int timestep) override = 0;
    PipelineResult update(size_t /*synapse_address*/, bool /*read*/,
            long int /*timestep*/) final
    {
        throw std::logic_error(
                "Error: Dendrite H/W called with synapse inputs");
    }
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/, long int /*timestep*/) final
    {
        throw std::logic_error("Error: Dendrite H/W called with soma inputs");
    }
};

class SomaUnit : public PipelineUnit
{
public:
    SomaUnit() : PipelineUnit(HardwareBitfield::implements_soma) {}
    PipelineResult update(size_t neuron_address, std::optional<double> current_in, long int timestep) override = 0;
    PipelineResult update(size_t /*synapse_address*/, bool /*read*/, long int /*timestep*/) final
    {
        throw std::logic_error("Error: Soma H/W called with synapse inputs");
    }
    PipelineResult update(size_t /*neuron_address*/,
            std::optional<double> /*current_in*/,
            std::optional<size_t> /*synaptic_address*/, long int /*timestep*/) final
    {
        throw std::logic_error("Error: Soma H/W called with dendrite inputs");
    }
    void set_attribute_edge(size_t /*synapse_address*/, const std::string &/*attribute_name*/, const ModelAttribute &/*param*/) final {};
    void track_connection(size_t /*synapse_address*/, size_t /*src_neuron_id*/, size_t /*dest_neuron_id*/) final {};
};

BufferPosition pipeline_parse_buffer_pos_str(
        const std::string &buffer_pos_str, bool buffer_inside_unit);
}

// Include these member function definition inline rather than in pipeline.cpp
//  to avoid linkage issues with plug-ins. Generally, header-only approaches
//  are useful in this situation. In this case, I've inlined only public and
//  protected members. Private members can be defined externally.
inline sanafe::PipelineUnit::PipelineUnit(const HardwareBitfield hw_implemented)
        : implements_synapse(
                  hw_implemented & HardwareBitfield::implements_synapse)
        , implements_dendrite(
                  hw_implemented & HardwareBitfield::implements_dendrite)
        , implements_soma(hw_implemented & HardwareBitfield::implements_soma)
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
        process_input_fn = [this](const Timestep &ts, const MappedNeuron &n,
                                   std::optional<MappedConnection *> con,
                                   const PipelineResult &input) {
            return this->process_synapse_input(ts, n, con, input);
        };
    }
    else if (implements_dendrite)
    {
        process_input_fn = [this](const Timestep &ts, const MappedNeuron &n,
                                   std::optional<MappedConnection *> con,
                                   const PipelineResult &input) {
            return this->process_dendrite_input(ts, n, con, input);
        };
    }
    else if (implements_soma)
    {
        process_input_fn = [this](const Timestep &ts, const MappedNeuron &n,
                                   std::optional<MappedConnection *> con,
                                   const PipelineResult &input) {
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

inline std::vector<std::string> sanafe::PipelineUnit::get_attributes() const
{
    std::vector<std::string> keys;
    keys.reserve(supported_attributes.size());

    for (const auto &[attribute_name, description] : supported_attributes)
    {
        keys.push_back(attribute_name);
    }
    return keys;
}

inline void sanafe::PipelineUnit::register_attributes(
        const std::unordered_map<std::string, std::string> &attributes_with_descriptions)
{
    for (const auto &[attr, description] : attributes_with_descriptions)
    {
        supported_attributes[attr] = description;
    }
}

inline void sanafe::PipelineUnit::register_attributes(
        const std::set<std::string> &attribute_names)
{
    for (const auto &attr : attribute_names)
    {
        // Add supported attribute without any helpful description. If possible,
        //  users should always add a helpful description that lets API users
        //  introspect the model features
        supported_attributes[attr] = "";
    }
}

inline void sanafe::PipelineUnit::process_synapse_output(
        const MappedNeuron & /*n*/, std::optional<MappedConnection *> con,
        PipelineResult &output)
{
    calculate_synapse_default_energy_latency(*(con.value_or(nullptr)), output);
}

inline void sanafe::PipelineUnit::process_dendrite_output(const MappedNeuron &n,
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
        const Timestep &ts, const MappedNeuron & /*n*/,
        std::optional<MappedConnection *> con, const PipelineResult & /*input*/)
{
    bool read = false;
    size_t synapse_address = 0UL;

    if (!con.has_value())
    {
        INFO("Pipeline warning, didn't receive "
             "synaptic connection info. Check that no h/w unit is being "
             "invoked before this one in the pipeline.");
    }
    else
    {
        synapse_address = con.value()->mapped_synapse_hw_address;
        read = true;
    }
    PipelineResult output = update(synapse_address, read, ts.timestep);
    ++spikes_processed;

    return output;
}

inline sanafe::PipelineResult sanafe::PipelineUnit::process_dendrite_input(
        const Timestep &ts, const MappedNeuron &n,
        std::optional<MappedConnection *> con, const PipelineResult &input)
{
    PipelineResult output{};

    std::optional<size_t> synapse_address{std::nullopt};
    if (con.has_value() && (con.value() != nullptr))
    {
        synapse_address = con.value()->mapped_synapse_hw_address;
    }
    output = update(n.mapped_dendrite_hw_address,
            input.current, synapse_address, ts.timestep);

    return output;
}

inline sanafe::PipelineResult sanafe::PipelineUnit::process_soma_input(
        const Timestep &ts, const MappedNeuron &n,
        std::optional<MappedConnection *> /*con*/, const PipelineResult &input)
{
    PipelineResult output =
            update(n.mapped_soma_hw_address, input.current, ts.timestep);
    return output;
}


inline void sanafe::PipelineUnit::calculate_synapse_default_energy_latency(
        const MappedConnection &con, PipelineResult &simulation_result)
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
        const MappedNeuron &n, PipelineResult &simulation_result)
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
        simulation_result.energy = n.dendrite_hw->default_energy_update.value();
    }

    const bool default_dendrite_latency_metrics_set =
            n.dendrite_hw->default_latency_update.has_value();
    //INFO("default dendrite latency metric set:%d\n", default_dendrite_latency_metrics_set);
    if (latency_simulated && default_dendrite_latency_metrics_set)
    {
        const std::string error(
                "Dendrite unit simulates latency and also has "
                "default latency metrics set.");
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }
    if (default_dendrite_latency_metrics_set)
    {
        simulation_result.latency =
                n.dendrite_hw->default_latency_update.value();
    }

    if (!simulation_result.energy.has_value())
    {
        const std::string error(
                "Dendrite unit does not simulate energy or provide "
                "a default energy cost in the architecture "
                "description.");
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }
    if (!simulation_result.latency.has_value())
    {
        const std::string error(
                "Dendrite unit does not simulate latency or "
                "provide a default latency cost in the architecture "
                "description.");
        INFO("Error: %s\n", error.c_str());
        throw std::runtime_error(error);
    }
}

inline void sanafe::PipelineUnit::calculate_soma_default_energy_latency(
        const MappedNeuron &n, PipelineResult &simulation_result)
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
