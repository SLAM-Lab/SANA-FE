// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
#ifndef PIPELINE_HEADER_INCLUDED_
#define PIPELINE_HEADER_INCLUDED_

#include <filesystem>
#include <optional>
#include <set>
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
    NeuronStatus status{INVALID_NEURON_STATE};
    // Optionally simulate energy and/or latency
    std::optional<double> energy{std::nullopt};
    std::optional<double> latency{std::nullopt};
};

class PipelineUnit
{
public:
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
    virtual PipelineResult update(size_t synapse_address, bool read = false) { throw std::logic_error("Error: Synapse input not implemented"); }
    // If using dendritic inputs (neuron address, synaptic current and synaptic address for additional info)
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) { throw std::logic_error("Error: Dendrite input not implemented"); }
    // If using somatic inputs (address and current in)
    virtual PipelineResult update(size_t neuron_address, std::optional<double> current_in) { throw std::logic_error("Error: Soma input not implemented"); }
    virtual void map_connection(MappedConnection &con) {}
    virtual double get_potential(size_t neuron_address) { return 0.0; }

    // Normal member functions and function pointers
    void set_time(long int timestep) { simulation_time = timestep; }
    void set_attributes(std::string unit_name, const ModelInfo &model);
    PipelineResult process(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    void register_attributes(const std::set<std::string> &attribute_names);
    void check_attribute(std::string attribute_name);

    using InputInterfaceFunc = PipelineResult (PipelineUnit:: *)(Timestep &, MappedNeuron &, std::optional<MappedConnection*>, const PipelineResult &);
    using OutputInterfaceFunc = void (PipelineUnit:: *)(MappedNeuron &, std::optional<MappedConnection *>, PipelineResult &);
    InputInterfaceFunc process_input_fn{nullptr};
    OutputInterfaceFunc process_output_fn{nullptr};

    // Model information
    std::map<std::string, ModelAttribute> model_attributes{};
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
    long int neuron_updates{0L};
    long int neurons_fired{0L};
    long int neuron_count{0L};

    // Implementation flags, set whichever operations your derived unit supports
    //  to 'true'. Note that a hardware unit must support one or more of these
    const bool implements_synapse;
    const bool implements_dendrite;
    const bool implements_soma;

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
    long int simulation_time{0L};
    PipelineUnit(bool implements_synapse, bool implements_dendrite, bool implements_soma);

    PipelineResult process_synapse_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_dendrite_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);
    PipelineResult process_soma_input(Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);

    void process_synapse_output(MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);
    void process_dendrite_output(MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);
    void process_soma_output(MappedNeuron &n, std::optional<MappedConnection *> con, PipelineResult &output);

private:
    static void calculate_synapse_default_energy_latency(
            MappedConnection &con, PipelineResult &simulation_result);
    static void calculate_dendrite_default_energy_latency(
            MappedNeuron &n, PipelineResult &simulation_result);
    static void calculate_soma_default_energy_latency(
            MappedNeuron &n, PipelineResult &simulation_result);
    static void update_soma_activity(
            MappedNeuron &n, const PipelineResult &simulation_result);
    void check_outputs(
            const MappedNeuron &n, const PipelineResult &result) const;
};

// Specific unit base classes, for the normal use case where the model
//  implements only synaptic, dendritic or somatic functionality (and not a)
//  combination
class SynapseUnit : public PipelineUnit
{
public:
    SynapseUnit() : PipelineUnit(true, false, false) {};
    PipelineResult update(size_t synapse_address, bool read = false) override = 0;
    PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) final { throw std::logic_error("Error: Synapse H/W called with dendrite inputs"); }
    PipelineResult update(size_t neuron_address, std::optional<double> current_in) final { throw std::logic_error("Error: Synapse H/W called with soma inputs"); }
    void set_attribute_neuron(size_t neuron_address,  const std::string &attribute_name, const ModelAttribute &param) final {};
};

class DendriteUnit : public PipelineUnit
{
public:
    DendriteUnit() : PipelineUnit(false, true, false) {};
    PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) override = 0;
    PipelineResult update(size_t synapse_address, bool read = false) final { throw std::logic_error("Error: Dendrite H/W called with synapse inputs"); }
    PipelineResult update(size_t neuron_address, std::optional<double> current_in) final { throw std::logic_error("Error: Dendrite H/W called with soma inputs"); }
};

class SomaUnit : public PipelineUnit
{
public:
    SomaUnit() : PipelineUnit(false, false, true) {};
    PipelineResult update(size_t neuron_address, std::optional<double> current_in) override = 0;
    PipelineResult update(size_t synapse_address, bool read = false) final { throw std::logic_error("Error: Soma H/W called with synapse inputs"); }
    PipelineResult update(size_t neuron_address, std::optional<double> current_in, std::optional<size_t> synaptic_address) final { throw std::logic_error("Error: Soma H/W called with dendrite inputs"); }
    void set_attribute_edge(size_t synapse_address, const std::string &attribute_name, const ModelAttribute &param) final {};
    void map_connection(MappedConnection &con) final {};
};

BufferPosition pipeline_parse_buffer_pos_str(const std::string &buffer_pos_str, const bool buffer_inside_unit);

}
#endif
