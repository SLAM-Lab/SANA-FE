#include <list>
#include <cmath>
#include <cassert>
#include <sstream>

#include "plugins.hpp"
#include "models.hpp"
#include "description.hpp"
#include "print.hpp"

// **** Soma models ****
LoihiLifModel::LoihiLifModel()
{
	potential = 0.0;
	soma_last_updated = 0;
	// Default is no leak (potential decay), i.e., the potential for the
	//  next timestep is 100% of the previous timestep's
	leak_decay = 1.0;
	bias = 0.0;
	threshold = 0.0;
	reverse_threshold = 0.0;
	reset = 0.0;
	reverse_reset = 0.0;
	reset_mode = sanafe::NEURON_RESET_HARD;
	reverse_reset_mode = sanafe::NEURON_NO_RESET;
	force_update = 0;

	return;
}

LoihiLifModel::~LoihiLifModel()
{
	return;
}

void LoihiLifModel::set_attributes(const std::list<Attribute> &attr)
{
	for (auto a: attr)
	{
		std::istringstream ss(a.value_str);
		if (a.key == "threshold")
		{
			ss >> threshold;
		}
		else if (a.key == "reverse_threshold")
		{
			ss >> reverse_threshold;
		}
		else if (a.key == "reset")
		{
			ss >> reset;
		}
		else if (a.key == "reverse_reset")
		{
			ss >> reverse_reset;
		}
		else if (a.key == "reset_mode")
		{
			reset_mode = model_parse_reset_mode(a.value_str);
		}
		else if (a.key =="reverse_reset_mode")
		{
			reverse_reset_mode =
				model_parse_reset_mode(a.value_str);
		}
		else if (a.key =="leak_decay")
		{
			ss >> leak_decay;
		}
		else if (a.key =="bias")
		{
			ss >> bias;
		}
		else if (a.key == "force_update")
		{
			ss >> force_update;
		}
        }
}

sanafe::NeuronStatus LoihiLifModel::update(const double current_in)
{
	// Calculate the change in potential since the last update e.g.
	//  integate inputs and apply any potential leak
	INFO("Updating potential, before:%f\n", potential);
	sanafe::NeuronStatus state = sanafe::IDLE;
	potential *= leak_decay;

	// Add randomized noise to potential if enabled
	/*
	if (noise_type == NOISE_FILE_STREAM)
	{
		// TODO: fix noise generation. This depends on which core
		//  is simulating the neuron. So somehow the neuron can still
		//  need information about which core it is executing on
		double random_potential = sim_generate_noise(n);
		n->potential += random_potential;
	}
	*/
	// Add the synaptic / dendrite current to the potential
	//printf("n->bias:%lf n->potential before:%lf current_in:%lf\n", n->bias, n->potential, current_in);
	potential += current_in + bias;
	INFO("leak decay:%lf bias:%lf threshold:%lf potential after:%lf\n",
		leak_decay, bias, threshold, potential);
	INFO("Updating potential, after:%f\n", potential);

	// Update soma, if there are any received spikes, there is a non-zero
	//  bias or we force the neuron to update every time-step
	if ((fabs(potential) > 0.0) || (current_in != 0) ||
		(fabs(bias) > 0.0) || force_update)
	{
		// Neuron is turned on and potential write
		state = sanafe::UPDATED;
	}

	// Check against threshold potential (for spiking)
	if (((bias != 0.0) && (potential > threshold)) ||
		((bias == 0.0) && (potential >= threshold)))
	{
		if (reset_mode == sanafe::NEURON_RESET_HARD)
		{
			potential = reset;
		}
		else if (reset_mode == sanafe::NEURON_RESET_SOFT)
		{
			potential -= threshold;
		}
		state = sanafe::FIRED;
		INFO("Neuron fired!\n");
	}
	// Check against reverse threshold
	if (potential < reverse_threshold)
	{
		if (reverse_reset_mode == sanafe::NEURON_RESET_SOFT)
		{
			potential -= reverse_threshold;
		}
		else if (reverse_reset_mode == sanafe::NEURON_RESET_HARD)
		{
			potential = reverse_reset;
		}
		else if (reverse_reset_mode == sanafe::NEURON_RESET_SATURATE)
		{
			potential = reverse_threshold;
		}
	}
	return state;
}


sanafe::NeuronResetModes model_parse_reset_mode(const std::string &str)
{
	sanafe::NeuronResetModes reset_mode;

	if (str == "none")
	{
		reset_mode = sanafe::NEURON_NO_RESET;
	}
	else if (str == "soft")
	{
		reset_mode = sanafe::NEURON_RESET_SOFT;
	}
	else if (str == "hard")
	{
		reset_mode = sanafe::NEURON_RESET_HARD;
	}
	else if (str == "saturate")
	{
		reset_mode = sanafe::NEURON_RESET_SATURATE;
	}
	else
	{
		throw std::invalid_argument("Error: reset mode not recognized");
	}

	return reset_mode;
}

// double sim_update_soma_truenorth(
// 	Timestep &ts, Neuron *n, const double current_in)
// {
// 	struct SomaProcessor *soma = n->soma_hw;
// 	double v, latency = 0.0;
// 	int randomize_threshold;

// 	// Apply leak
// 	while (n->soma_last_updated <= ts->timestep)
// 	{
// 		// Linear leak
// 		if (soma->leak_towards_zero)
// 		{
// 			// TODO: what happens if we're above zero but by less
// 			//  than the leak amount (for convergent), will we
// 			//  oscillate between the two? Does it matter
// 			if (n->potential > 0.0)
// 			{
// 				n->potential -= n->leak_bias;
// 			}
// 			else if (n->potential < 0.0)
// 			{
// 				n->potential += n->leak_bias;
// 			}
// 			// else equals zero, so no leak is applied
// 		}
// 		else
// 		{
// 			n->potential += n->leak_decay;
// 		}
// 		n->soma_last_updated++;
// 	}

// 	// Add the synaptic currents, processed by the dendrite
// 	n->potential += current_in + n->bias;
// 	n->current = 0.0;
// 	n->charge = 0.0;

// 	// Apply thresholding and reset
// 	v = n->potential;
// 	randomize_threshold = (n->random_range_mask != 0);
// 	if (randomize_threshold)
// 	{
// 		unsigned int r = rand() & n->random_range_mask;
// 		v += (double) r;
// 	}

// 	TRACE2("v:%lf +vth:%lf mode:%d -vth:%lf mode:%d\n", v, n->threshold,
// 		n->group->reset_mode, n->reverse_threshold,
// 		n->group->reverse_reset_mode);
// 	if (v >= n->threshold)
// 	{
// 		struct SomaProcessor *soma = n->soma_hw;
// 		int reset_mode = n->group->reset_mode;

// 		if (reset_mode == NEURON_RESET_HARD)
// 		{
// 			n->potential = n->reset;
// 		}
// 		else if (reset_mode == NEURON_RESET_SOFT)
// 		{
// 			n->potential -= n->threshold;
// 		}
// 		else if (reset_mode == NEURON_RESET_SATURATE)
// 		{
// 			n->potential = n->threshold;
// 		}
// 		n->fired = 1;
// 		soma->neurons_fired++;
// 		latency += soma->latency_spiking;
// 		if (n->core->buffer_pos != BUFFER_AXON_OUT)
// 		{
// 			sim_neuron_send_spike_message(ts, n);
// 		}
// 	}
// 	else if (v <= n->reverse_threshold)
// 	{
// 		int reset_mode = n->group->reverse_reset_mode;
// 		if (reset_mode == NEURON_RESET_HARD)
// 		{
// 			n->potential = n->reverse_reset;
// 		}
// 		else if (reset_mode == NEURON_RESET_SOFT)
// 		{
// 			n->potential += n->reverse_threshold;
// 		}
// 		else if (reset_mode == NEURON_RESET_SATURATE)
// 		{
// 			n->potential = n->reverse_threshold;
// 		}
// 		// No spike is generated
// 	}
// 	TRACE2("potential:%lf threshold %lf\n", n->potential, n->threshold);

// 	return latency;
// }