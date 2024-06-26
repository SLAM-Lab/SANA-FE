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

class CurrentBasedSynapseModel : public SynapseModel
{
public:
    CurrentBasedSynapseModel();
    CurrentBasedSynapseModel(const CurrentBasedSynapseModel &copy) = default;
    CurrentBasedSynapseModel(CurrentBasedSynapseModel &&other) = default;
    virtual ~CurrentBasedSynapseModel() = default;
    CurrentBasedSynapseModel &operator=(const CurrentBasedSynapseModel &other) = default;
    CurrentBasedSynapseModel &operator=(CurrentBasedSynapseModel &&other) = default;

    double update(std::optional<int> synapse_address = std::nullopt,  bool step = true) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;

private:
    double weight{0.0};
    double min_synaptic_resolution{0.0};
    double current{0.0};
    double synaptic_current_decay{1.0};
    int weight_bits{8};
};

class SingleCompartmentModel : public DendriteModel
{
public:
    SingleCompartmentModel();
    SingleCompartmentModel(const SingleCompartmentModel &copy) = default;
    SingleCompartmentModel(SingleCompartmentModel &&other) = default;
    virtual ~SingleCompartmentModel() = default;
    SingleCompartmentModel &operator=(const SingleCompartmentModel &other) = default;
    SingleCompartmentModel &operator=(SingleCompartmentModel &&other) = default;

    double update(std::optional<Synapse> current_in = std::nullopt, const bool step = true) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
private:
    double accumulated_charge, leak_decay;
};

class MultiTapModel1D : public DendriteModel
{
public:
    MultiTapModel1D();
    MultiTapModel1D(const MultiTapModel1D &copy) = default;
    MultiTapModel1D(MultiTapModel1D &&other) = default;
    virtual ~MultiTapModel1D() = default;
    MultiTapModel1D &operator=(const MultiTapModel1D &other) = default;
    MultiTapModel1D &operator=(MultiTapModel1D &&other) = default;

    double update(std::optional<Synapse> current_in = std::nullopt, const bool step = true) override;
    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
private:
    // Assuming a 1D tap dendrite
    std::vector<double> tap_voltages, next_voltages;
    std::vector<double> space_constants, time_constants;
};

class SomaModel
{
public:
    SomaModel(const std::string &gid, const std::string &nid) : group_id(gid), neuron_id(nid) {}
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
    const std::string neuron_id;
};

class LoihiLifModel : public SomaModel
{
public:
    LoihiLifModel(const std::string &gid, const std::string &nid);
    LoihiLifModel(const LoihiLifModel &copy) = default;
    LoihiLifModel(LoihiLifModel &&other) = default;
    virtual ~LoihiLifModel() = default;
    LoihiLifModel &operator=(const LoihiLifModel &other) = delete;
    LoihiLifModel &operator=(LoihiLifModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    NeuronStatus update(const std::optional<double> current_in,  const bool step = true) override;
    double get_potential() override { return potential; }
private:
    bool force_update;
    int reset_mode, reverse_reset_mode;
    //int noise_type;
    double potential, leak_decay, bias, threshold, reverse_threshold;
    double reset, reverse_reset;
};

class TrueNorthModel : public SomaModel
{
public:
    TrueNorthModel(const std::string &gid, const std::string &nid);
    TrueNorthModel(const TrueNorthModel &copy) = default;
    TrueNorthModel(TrueNorthModel &&other) = default;
    virtual ~TrueNorthModel() = default;
    TrueNorthModel &operator=(const TrueNorthModel &other) = delete;
    TrueNorthModel &operator=(TrueNorthModel &&other) = delete;

    void set_attributes(const std::map<std::string, ModelParam> &attr) override;
    NeuronStatus update(const std::optional<double> current_in = std::nullopt,  const bool step = true) override;
    double get_potential() override { return potential; }
private:
    bool force_update;
    unsigned int random_range_mask;
    int reset_mode, reverse_reset_mode;
    bool leak_towards_zero;
    // int noise_type;
    double potential, leak, bias, threshold, reverse_threshold;
    double reset, reverse_reset;
};

NeuronResetModes model_parse_reset_mode(const std::string &str);
std::shared_ptr<SynapseModel> model_get_synapse(const std::string &model_name);
std::shared_ptr<DendriteModel> model_get_dendrite(const std::string &model_name);
std::shared_ptr<SomaModel> model_get_soma(const std::string &model_name, const std::string &group_id, const std::string &id);

}

#endif
