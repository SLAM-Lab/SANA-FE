// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// models.hpp
#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include "network.hpp"

#include <map>
#include <memory> // For shared_ptr<T>
#include <optional>
#include <vector>

namespace sanafe
{
enum NeuronStatus: int { IDLE, UPDATED, FIRED};
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
    SynapseModel() {}
    virtual ~SynapseModel() {}
    virtual double update(std::optional<int> synapse_address=std::nullopt, const bool step=true) = 0;
    // Set synapse attributes
    virtual void set_attributes(const std::map<std::string, NeuronAttribute> &attr) = 0;
};

class DendriteModel
{
public:
    DendriteModel() {}
    virtual ~DendriteModel() {}
    // TODO: Can't forward declare synapses... find a better way?
    //  Might be easier to just make it two optional arguments... unpack it from the Synapse struct
    //  If it's optional, it gets harder to pass a reference, maybe should pass a pointer instead...
    virtual double update(std::optional<Synapse> synapse_in=std::nullopt, const bool step=true) = 0;
    virtual void set_attributes(const std::map<std::string, NeuronAttribute> &attr) = 0;
};

class CurrentBasedSynapseModel: public SynapseModel
{
public:
    CurrentBasedSynapseModel();
    ~CurrentBasedSynapseModel() {}
    virtual double update(std::optional<int> synapse_address=std::nullopt,  const bool step=true);
    virtual void set_attributes(const std::map<std::string, NeuronAttribute> &attr);

private:
    double weight, min_synaptic_resolution, current, synaptic_current_decay;
    int weight_bits;
};

class SingleCompartmentModel: public DendriteModel
{
public:
    SingleCompartmentModel();
    ~SingleCompartmentModel() {}
    virtual double update(std::optional<Synapse> current_in=std::nullopt, const bool step=true);
    virtual void set_attributes(const std::map<std::string, NeuronAttribute> &attr);
private:
    double accumulated_charge, leak_decay;
};

class MultiTapModel1D: public DendriteModel
{
public:
    MultiTapModel1D();
    virtual double update(std::optional<Synapse> current_in=std::nullopt, const bool step=true);
    virtual void set_attributes(const std::map<std::string, NeuronAttribute> &attr);
private:
    // Assuming a 1D tap dendrite
    std::vector<double> tap_voltages, next_voltages;
    std::vector<double> space_constants, time_constants;
};

class SomaModel
{
public:
    SomaModel(const int gid, const int nid) : group_id(gid), neuron_id(nid) {}
    virtual ~SomaModel() {}
    virtual NeuronStatus update(const std::optional<double> current_in=std::nullopt, const bool step=true) = 0;
    virtual void set_attributes(const std::map<std::string, NeuronAttribute2> &attr) = 0;
    virtual double get_potential() { return 0.0; }
protected:
    const int group_id, neuron_id;
};

class LoihiLifModel: public SomaModel
{
public:
    LoihiLifModel(const int gid, const int nid);
    ~LoihiLifModel() {}
    void set_attributes(const std::map<std::string, NeuronAttribute2> &attr);
    NeuronStatus update(const std::optional<double> current_in,  const bool step=true);
    double get_potential() { return potential; }
private:
    bool force_update;
    int reset_mode, reverse_reset_mode;
    //int noise_type;
    double potential, leak_decay, bias, threshold, reverse_threshold;
        double reset, reverse_reset;
};

class TrueNorthModel: public SomaModel
{
public:
    TrueNorthModel(const int gid, const int nid);
    ~TrueNorthModel() {}
    void set_attributes(const std::map<std::string, NeuronAttribute2> &attr);
    NeuronStatus update(const std::optional<double> current_in=std::nullopt,  const bool step=true);
    double get_potential() { return potential; }
private:
    bool force_update;
    unsigned int random_range_mask;
    int reset_mode, reverse_reset_mode;
    bool leak_towards_zero;
    //int noise_type;
    double potential, leak, bias, threshold, reverse_threshold;
    double reset, reverse_reset;
};

NeuronResetModes model_parse_reset_mode(const std::string &str);
std::shared_ptr<SynapseModel> model_get_synapse(const std::string &model_name);
std::shared_ptr<DendriteModel> model_get_dendrite(const std::string &model_name);
std::shared_ptr<SomaModel> model_get_soma(const std::string &model_name, const int group_id, const int id);

}

#endif
