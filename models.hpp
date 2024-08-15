// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// models.hpp
#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include <map>
#include <memory>  // For shared_ptr<T>
#include <optional>
#include <string>
#include <vector>

#include "network.hpp"

namespace sanafe
{
enum NeuronResetModes
{
    NEURON_NO_RESET,
    NEURON_RESET_SOFT,
    NEURON_RESET_HARD,
    NEURON_RESET_SATURATE,
    NEURON_RESET_MODE_COUNT,
};

class SynapseModel
{
public:
    SynapseModel() = default;
    SynapseModel(const SynapseModel &copy) = default;
    SynapseModel(SynapseModel &&other) = default;
    virtual ~SynapseModel() = default;
    SynapseModel &operator=(const SynapseModel &other) = default;
    SynapseModel &operator=(SynapseModel &&other) = default;

    virtual double update(std::optional<int> synapse_address = std::nullopt, bool step = true) = 0;
    // Set synapse attributes
    virtual void set_attributes(const std::map<std::string, ModelParam> &attr) = 0;
};

class DendriteModel
{
public:
    DendriteModel() = default;
    DendriteModel(const DendriteModel &copy) = default;
    DendriteModel(DendriteModel &&other) = default;
    virtual ~DendriteModel() = default;
    DendriteModel &operator=(const DendriteModel &other) = default;
    DendriteModel &operator=(DendriteModel &&other) = default;

    virtual double update(std::optional<Synapse> synapse_in = std::nullopt, bool step = true) = 0;
    virtual void set_attributes(const std::map<std::string, ModelParam> &attr) = 0;
};

constexpr int default_weight_bits = 8;  // Based on real-world H/W e.g., Loihi
class CurrentBasedSynapseModel : public SynapseModel
{
public:
    CurrentBasedSynapseModel() = default;
    CurrentBasedSynapseModel(const CurrentBasedSynapseModel &copy) = default;
    CurrentBasedSynapseModel(CurrentBasedSynapseModel &&other) = default;
    ~CurrentBasedSynapseModel() override = default;
    CurrentBasedSynapseModel &operator=(const CurrentBasedSynapseModel &other) = default;
    CurrentBasedSynapseModel &operator=(CurrentBasedSynapseModel &&other) = default;

    double update(std::optional<int> synapse_address = std::nullopt,  bool step = true) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;

private:
    double weight{0.0};
    double min_synaptic_resolution{0.0};
    double current{0.0};
    double synaptic_current_decay{0.0};
    int weight_bits{default_weight_bits};
};

class SingleCompartmentModel : public DendriteModel
{
public:
    SingleCompartmentModel();
    SingleCompartmentModel(const SingleCompartmentModel &copy) = default;
    SingleCompartmentModel(SingleCompartmentModel &&other) = default;
    ~SingleCompartmentModel() override = default;
    SingleCompartmentModel &operator=(const SingleCompartmentModel &other) = default;
    SingleCompartmentModel &operator=(SingleCompartmentModel &&other) = default;

    double update(std::optional<Synapse> synapse_in = std::nullopt, bool step = true) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
private:
    double accumulated_charge{0.0};
    double leak_decay{0.0};
};

class MultiTapModel1D : public DendriteModel
{
public:
    MultiTapModel1D();
    MultiTapModel1D(const MultiTapModel1D &copy) = default;
    MultiTapModel1D(MultiTapModel1D &&other) = default;
    ~MultiTapModel1D() override = default;
    MultiTapModel1D &operator=(const MultiTapModel1D &other) = default;
    MultiTapModel1D &operator=(MultiTapModel1D &&other) = default;

    double update(std::optional<Synapse> synapse_in = std::nullopt, bool step = true) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
private:
    // Assuming a 1D tap dendrite
    std::vector<double> tap_voltages, next_voltages;
    std::vector<double> space_constants, time_constants;
};

class SomaModel
{
public:
    SomaModel(std::string gid, size_t nid) : group_id(std::move(gid)), neuron_id(nid) {}
    SomaModel(const SomaModel &copy) = default;
    SomaModel(SomaModel &&other) = default;
    virtual ~SomaModel() = default;
    SomaModel &operator=(const SomaModel &other) = delete;
    SomaModel &operator=(SomaModel &&other) = delete;

    virtual NeuronStatus update(std::optional<double> current_in = std::nullopt, bool step = true) = 0;
    virtual void set_attributes(const std::map<std::string, ModelParam> &attr) = 0;
    virtual double get_potential() { return 0.0; }
protected:
    const std::string group_id;
    const size_t neuron_id;
};

class LoihiLifModel : public SomaModel
{
public:
    LoihiLifModel(std::string gid, size_t nid);
    LoihiLifModel(const LoihiLifModel &copy) = default;
    LoihiLifModel(LoihiLifModel &&other) = default;
    ~LoihiLifModel() override = default;
    LoihiLifModel &operator=(const LoihiLifModel &other) = delete;
    LoihiLifModel &operator=(LoihiLifModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    NeuronStatus update(std::optional<double> current_in,  bool step = true) override;
    double get_potential() override { return potential; }
private:
    bool force_update{false};
    int reset_mode{NEURON_RESET_HARD};
    int reverse_reset_mode{NEURON_NO_RESET};
    //int noise_type;
    double potential{0.0};
    double leak_decay{1.0};
    double bias{0.0};
    double threshold{0.0};
    double reverse_threshold{0.0};
    double reset{0.0};
    double reverse_reset{0.0};
};

class TrueNorthModel : public SomaModel
{
public:
    TrueNorthModel(const std::string &gid, size_t nid);
    TrueNorthModel(const TrueNorthModel &copy) = default;
    TrueNorthModel(TrueNorthModel &&other) = default;
    ~TrueNorthModel() override = default;
    TrueNorthModel &operator=(const TrueNorthModel &other) = delete;
    TrueNorthModel &operator=(TrueNorthModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    NeuronStatus update(std::optional<double> current_in = std::nullopt, bool step = true) override;
    double get_potential() override { return potential; }
private:
    bool force_update{false};
    unsigned int random_range_mask{0U};
    int reset_mode{NEURON_RESET_HARD};
    int reverse_reset_mode{NEURON_NO_RESET};
    bool leak_towards_zero{true};
    // int noise_type;
    double potential{0.0};
    double leak{0.0};  // Default is no leak (potential decay)
    double bias{0.0};
    double threshold{0.0};
    double reverse_threshold{0.0};
    double reset{0.0};
    double reverse_reset{0.0};
};

class InputModel: public SomaModel
{
public:
    InputModel(const std::string &gid, size_t nid) : SomaModel(gid, nid) {}
    InputModel(const InputModel &copy) = default;
    InputModel(InputModel &&other) = default;
    ~InputModel() override = default;
    InputModel &operator=(const InputModel &other) = delete;
    InputModel &operator=(InputModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    NeuronStatus update(std::optional<double> current_in = std::nullopt, bool step = true) override;

private:
    std::vector<bool> spikes{};
    std::vector<bool>::const_iterator curr_spike{spikes.begin()};
    bool send_spike{false};
};

NeuronResetModes model_parse_reset_mode(const std::string &str);
std::shared_ptr<SynapseModel> model_get_synapse(const std::string &model_name);
std::shared_ptr<DendriteModel> model_get_dendrite(const std::string &model_name);
std::shared_ptr<SomaModel> model_get_soma(const std::string &model_name, const std::string &group_id, size_t id);

}

#endif
