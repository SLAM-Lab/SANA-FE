// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// hodgkin_huxley.cpp
// Plugin implementation of the Hodgkin-Huxley neuron model. Implemented by
//  Robin Sam.
// Model inspired by this paper: https://ieeexplore.ieee.org/document/9235538
//  and this textbook: https://mrgreene09.github.io/computational-neuroscience-textbook
#include <cstring>
#include <cmath>
#include <iostream>
#include <map>

#include "../print.hpp"
#include "../models.hpp"
#include "../plugins.hpp"
#include "../arch.hpp"

class HH: public sanafe::SomaModel
{
// HH specific
public:
	// system variables
	double C_m;
	double g_Na;
	double g_K;
	double g_L;
	double V_Na;
	double V_K;
	double V_L;
	double dt;

	// main parameters
	double V, prev_V;	// Membrane potential
	double I;		// Stimulation current per area
	double m;		// m, n, h are coeff
	double n;
	double h;

	// internal results of various differential equations
	double alpha_m;
	double alpha_n;
	double alpha_h;
	double beta_m;
	double beta_n;
	double beta_h;

	double tau_m, tau_n, tau_h;
	double pm, pn, ph;
	double denominator, tau_V, Vinf;

	HH(const int gid, const int nid) : SomaModel(gid, nid)
	{
		V = 0.0;
		C_m = 10.0;	// Effective capacitance per area of membrane; default is 1
		g_Na= 1200.0;	// Conductance of sodium
		g_K = 360.0;	// Conductance of potassium
		g_L = 3.0;	// Conductance of leak channel
		V_Na = 50.0;	// Reverse potential of sodium
		V_K = -77.0;	// Reverse potential of potassium
		V_L = -54.387;	// Reverse potential of leak channel
		dt = 0.1;
	}

	virtual double get_potential() { return V; }
	virtual void set_attributes(
		const std::map<std::string, std::string> &attr)
	{
		/*** Set attributes ***/
		for (auto a: attr)
		{
			int ret = 1;
			const std::string &key = a.first;
			const std::string &value_str = a.second;
			if (strncmp("m", key.c_str(), key.length()) == 0)
			{
				ret = sscanf(value_str.c_str(), "%lf", &m);
			}
			else if (strncmp("n", key.c_str(), key.length()) == 0)
			{
				ret = sscanf(value_str.c_str(), "%lf", &n);
			}
			else if (strncmp("h", key.c_str(), key.length()) == 0)
			{
				ret = sscanf(value_str.c_str(), "%lf", &h);
			}
			else if (strncmp(
				"current", key.c_str(), key.length()) == 0)
			{
				ret = sscanf(value_str.c_str(), "%lf", &I);
			}
			if (ret < 1)
			{
				INFO("Invalid attribute (%s:%s)\n",
					key.c_str(), value_str.c_str());
				exit(1);
			}
		}
	}

	sanafe::NeuronStatus update(const std::optional<double> current_in)
	{
		sanafe::NeuronStatus status = sanafe::IDLE;

		// Calculate the change in potential since the last update e.g.
		//  integate inputs and apply any potential leak
		TRACE1("Updating potential, before:%f\n", V);

		alpha_n = (0.01*(V+55)) / (1-exp(-0.1*(V+55)));
		alpha_m = (0.1*(V+40)) / (1-exp(-0.1*(V+40)));
		alpha_h = 0.07*exp(-0.05*(V+65));

		beta_n = 0.125*exp(-0.01125*(V+55));
		beta_m = 4*exp(-0.05556*(V+65));
		beta_h = 1/(1 + exp(-0.1*(V+35)));

		tau_n = 1 / (alpha_n + beta_n);
		tau_m = 1 / (alpha_m + beta_m);
		tau_h = 1 / (alpha_h + beta_h);

		pm = alpha_m / (alpha_m + beta_m);
		pn = alpha_n / (alpha_n + beta_n);
		ph = alpha_h / (alpha_h + beta_h);

		denominator = g_L + g_K*(pow(n,4)) + g_Na*(pow(m,3)*h);
		tau_V = C_m/denominator;
		Vinf = ((g_L)*V_L + g_K*(pow(n,4))*V_K + g_Na*(pow(m,3))*h*V_Na + I)/denominator;

		// update main parameters
		prev_V = V;
		V = Vinf + (V - Vinf)*exp(-1*dt/tau_V);
		m = pm + (m - pm)*exp(-1*dt/tau_m);
		n = pn + (n - pn)*exp(-1*dt/tau_n);
		h = ph + (h - ph)*exp(-1*dt/tau_h);

		// Check against threshold potential (for spiking)
		if ((prev_V < 25) && (V > 25))
		{
			// If voltage just crossed the 25 mV boundary, then
			//  spike
			status = sanafe::FIRED;
		}
		else
		{
			status = sanafe::UPDATED;
		}

		INFO("Updating potential, after:%f\n", V);

		return status;
	}
};

// the Class factories
extern "C" sanafe::SomaModel *create_HH(const int gid, const int nid)
{
	TRACE1("Creating HH soma instance\n");
	return (sanafe::SomaModel *) new HH(gid, nid);
}

// Memory Leak?
//extern "C" void destroy_HH(sanafe::SomaModel *HH)
//{
//	delete HH;
//}
