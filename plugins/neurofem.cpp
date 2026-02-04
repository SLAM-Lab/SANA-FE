// Copyright (c) 2026 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.

// TODO: this is a placeholder model for the NeuroFEM neuron in Loihi 2.
//  In reality, Loihi 2 uses a similar compartment to Loihi 1 with micro-ops.
//  Since I don't know the inner workings of Loihi 2, mainly how micro-ops work
//  and which ones are supported, stick with a dummy hardcoded model for now.

#include <cassert>
#include <stack>
#include <random>
#include <variant>

#include "arch.hpp"
#include "models.hpp"
#include "plugins.hpp"
#include "print.hpp"

//#include <iostream>

class NeuroFEMModel : public sanafe::PipelineUnit
{
public:
    std::random_device rd;
    //std::mt19937 gen{42};
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d{};

    static inline const std::set<std::string> supported_attributes{"weight",
            "w", "lambda_v", "lambda_d", "sigma_v", "ki", "kp", "bias",
            "threshold", "reset", "dt", "compartment"};

    NeuroFEMModel()
            : PipelineUnit(false, true, true) // implements dendrite+soma
    {
        register_attributes(supported_attributes);
    }
    NeuroFEMModel(const NeuroFEMModel &copy) = default;
    NeuroFEMModel(NeuroFEMModel &&other) = default;
    ~NeuroFEMModel() override = default;
    NeuroFEMModel &operator=(const NeuroFEMModel &other) = delete;
    NeuroFEMModel &operator=(NeuroFEMModel &&other) = delete;

    void set_attribute_hw(const std::string &param_name,
            const sanafe::ModelAttribute &param) override {};
    void set_attribute_neuron(size_t neuron_address,
            const std::string &param_name,
            const sanafe::ModelAttribute &param) override;
    void set_attribute_edge(size_t synapse_address,
            const std::string &param_name,
            const sanafe::ModelAttribute &param) override;
    sanafe::PipelineResult update(size_t neuron_address,
            std::optional<double> current_in,
            std::optional<size_t> synaptic_address,
            long int timestep) override;
    void reset() override;
    double get_potential(const size_t neuron_address) override
    {
        return neurons[neuron_address].potential;
    }

    struct NeuroFEMNeuron
    {
        double potential{0.0};
        double u1{0.0};
        double u2{0.0};
        double u_error{0.0};
        double u_integrated{0.0};

        double lambda_v{0.0};
        double lambda_d{0.0};
        double sigma_v{0.0};
        double ki{0.0};
        double kp{0.0};
        double bias{0.0};
        double threshold{0.0};
        double reset{0.0};
        double dt{1.0e-3};

        // TODO: generalize
        std::optional<double> u1_dendritic_accumulator{std::nullopt};
        std::optional<double> u2_dendritic_accumulator{std::nullopt};
        std::optional<double> next_u1_dendritic_accumulator{std::nullopt};
        std::optional<double> next_u2_dendritic_accumulator{std::nullopt};

        long int timesteps_simulated{0L};
        bool force_update{false};
    };

private:
    std::vector<NeuroFEMNeuron> neurons{};
    std::map<int, int> synapse_to_compartment{};

    sanafe::NeuronStatus process_fem(NeuroFEMNeuron &n);
    void check_compartments();
};

constexpr size_t neurofem_max_compartments{1024UL};
void NeuroFEMModel::check_compartments()
{
    // Count the total mapped compartments on this H/W core
    if (neurons.size() > neurofem_max_compartments)
    {
        std::string error = "Error: Mapped too many neurons for NeuroFEM (" +
                std::to_string(neurons.size()) + "> " +
                std::to_string(neurofem_max_compartments) + ")";
        INFO("%s\n", error.c_str());
        throw std::runtime_error(error);
    }
}

void NeuroFEMModel::set_attribute_edge(const size_t synapse_address,
        const std::string &param_name, const sanafe::ModelAttribute &param)
{
    if (param_name == "compartment")
    {
        // Attributes for mapped connections/synapses
        int compartment = static_cast<int>(param);
        synapse_to_compartment[synapse_address] = compartment;
        TRACE1(PLUGINS, "mapping synapse:%zu to compartment:%d\n",
                synapse_address, compartment);
        if ((compartment < 0) || (compartment > 1))
        {
            INFO("Error: compartment must be 0 or 1\n");
            throw std::runtime_error("Error: compartment must be 0 or 1");
        }
    }
}

void NeuroFEMModel::set_attribute_neuron(const size_t neuron_address,
        const std::string &param_name, const sanafe::ModelAttribute &param)
{
    // Attributes for mapped neurons
    if (neuron_address >= neurons.size())
    {
        neurons.resize(neuron_address + 1);
        INFO("Increasing mapped neurons to:%zu\n", neuron_address + 1);
        check_compartments();
    }
    NeuroFEMNeuron &n = neurons.at(neuron_address);

    if (param_name == "threshold")
    {
        n.threshold = static_cast<double>(param);
    }
    else if (param_name == "reset")
    {
        n.reset = static_cast<double>(param);
    }
    else if (param_name == "lambda_d")
    {
        double leak = static_cast<double>(param);
        n.lambda_d = leak;
    }
    else if (param_name == "lambda_v")
    {
        double leak = static_cast<double>(param);
        n.lambda_v = leak;
    }
    else if (param_name == "bias")
    {
        n.bias = static_cast<double>(param);
        TRACE2(PLUGINS, "Setting bias of %zu=%lf\n", neuron_address,
                n.bias);
    }
    else if (param_name == "dt")
    {
        n.dt = static_cast<double>(param);
    }
    else if (param_name == "kp")
    {
        n.kp = static_cast<double>(param);
    }
    else if (param_name == "ki")
    {
        n.ki = static_cast<double>(param);
    }
    else if (param_name == "sigma_v")
    {
        n.sigma_v = static_cast<double>(param);
    }
    else if ((param_name == "force_update") ||
            (param_name == "force_soma_update"))
    {
        n.force_update = static_cast<bool>(param);
    }

    TRACE1(MODELS, "Set parameter: %s\n", param_name.c_str());
}

sanafe::PipelineResult NeuroFEMModel::update(size_t neuron_address,
        std::optional<double> current_in,
        std::optional<size_t> synaptic_address,
        const long int simulation_time)
{
    NeuroFEMNeuron &n = neurons[neuron_address];
    sanafe::NeuronStatus state{sanafe::invalid_neuron_state};

    // Update time-step double buffered dendritic accumulators
    if (n.timesteps_simulated < (simulation_time - 1))
    {
        std::string error("Error: Must update neurons every time-step");
        INFO("%s\n", error.c_str());
        throw std::runtime_error(error);
    }
    if (n.timesteps_simulated == (simulation_time - 1))
    {
        n.u1_dendritic_accumulator = n.next_u1_dendritic_accumulator;
        n.u2_dendritic_accumulator = n.next_u2_dendritic_accumulator;
        n.next_u1_dendritic_accumulator = std::nullopt;
        n.next_u2_dendritic_accumulator = std::nullopt;

        TRACE1(PLUGINS,
                "Updating potential (nid:%zu cx:%zu ts:%ld), before:%lf\n",
                neuron_address, cx_address, n.timesteps_simulated,
                cx.potential);

        state = process_fem(n);
        TRACE1(PLUGINS,
                "Updating potential (nid:%zu cx:%zu ts:%ld), after:%lf\n",
                neuron_address, cx_address, n.timesteps_simulated,
                cx.potential);

        ++n.timesteps_simulated;
    }

    if (current_in.has_value())
    {
        TRACE1(PLUGINS, "Received synaptic current:%e\n", current_in.value());
        int cx = 0;
        if (synaptic_address.has_value())
        {
            // Note if not stored explicitly, the default compartment is 0
            cx = synapse_to_compartment[synaptic_address.value()];
        }
        if (cx == 0)
        {
            n.next_u1_dendritic_accumulator =
                    n.next_u1_dendritic_accumulator.value_or(0.0) +
                    current_in.value();
        }
        else // cx == 1
        {
            n.next_u2_dendritic_accumulator =
                    n.next_u2_dendritic_accumulator.value_or(0.0) +
                    current_in.value();
        }
    }

    return {std::nullopt, state, std::nullopt, std::nullopt};
}


//  ######################
//         # compute update to u1
//         du1 = dt*(-lambda_d*u1) + gtAg @ spikes[:, 0]
//         u1 += du1
//         ######################

//         ######################
//         # comptue update to u2
//         du2 = dt*(-lambda_d*u2) + lambda_d * omega_f @ spikes[:, 0]
//         u2 += du2
//         ######################

//         #############
//         # compute err
//         u_err = u1 + c_in_fixed_gamma
//         #############

//         ##############
//         # update u_int
//         u_int += dt*u_err
//         ##############

//         ##########
//         # update v
//         dv = dt*(-lambda_v*V + kp*u_err + ki*u_int + u2) - omega_f.dot(spikes[:, 0]) + sigma_v*np.random.randn(n_neurons)
//         V += dv
//         ##########

#include <iostream>

sanafe::NeuronStatus NeuroFEMModel::process_fem(NeuroFEMNeuron &n)
{
    //n.potential = static_cast<int>(n.potential * 64.0) / 64.0; // TODO: quantization
    TRACE1(PLUGINS, "Adding bias:%lf\n", n.bias);
    //cx.potential += cx.bias;

    n.u1 -= n.lambda_d * n.dt * n.u1;
    n.u2 -= n.lambda_d * n.dt * n.u2;

    n.u1 += n.u1_dendritic_accumulator.value_or(0.0);
    n.u2 += n.lambda_d * n.u2_dendritic_accumulator.value_or(0.0);

    n.u_error = n.u1 + n.bias;
    n.u_integrated += n.dt * n.u_error;

    double noise = d(gen);

    n.potential = n.potential - (n.lambda_v * n.dt * n.potential);
    n.potential = n.potential + (n.dt * n.kp * n.u_error) +
            (n.dt * n.ki * n.u_integrated) + (n.dt * n.u2) +
            (n.sigma_v * noise) - n.u2_dendritic_accumulator.value_or(0.0);

    sanafe::NeuronStatus state{sanafe::updated};

    // Check against threshold potential (for spiking)
    if (n.potential > n.threshold)
    {
        //n.potential -= n.threshold;
        n.potential = n.reset;
        state = sanafe::fired;
        TRACE1(PLUGINS, "Compartment fired.\n");
    }

    return state;
}

void NeuroFEMModel::reset()
{
    for (NeuroFEMNeuron &n : neurons)
    {
        n.potential = 0.0;
        n.u1 = 0.0;
        n.u2 = 0.0;
        n.u_integrated = 0.0;
        n.u_error = 0.0;
        n.u1_dendritic_accumulator = std::nullopt;
        n.u2_dendritic_accumulator = std::nullopt;
        n.next_u1_dendritic_accumulator = std::nullopt;
        n.next_u2_dendritic_accumulator = std::nullopt;
    }

    return;
}

// the Class factories
extern "C" sanafe::PipelineUnit *create_neurofem()
{
    TRACE1(PLUGINS, "Creating NeuroFEM dendrite/soma instance\n");
    return (sanafe::PipelineUnit *) new NeuroFEMModel();
}
