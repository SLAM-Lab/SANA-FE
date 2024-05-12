// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// models.hpp
#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include <map>
#include <vector>
#include "description.hpp"

#define MAX_NOISE_FILE_ENTRY 128

namespace sanafe
{
	enum NeuronStatus { IDLE, UPDATED, FIRED};
	enum NeuronResetModes
	{
		NEURON_NO_RESET,
		NEURON_RESET_SOFT,
		NEURON_RESET_HARD,
		NEURON_RESET_SATURATE,
		NEURON_RESET_MODE_COUNT,
	};

class SomaModel
{
public:
	SomaModel(const int gid, const int nid) : group_id(gid), neuron_id(nid) {}
	virtual ~SomaModel() {}
	virtual NeuronStatus update(const double current_in) = 0;
	virtual void set_attributes(const std::map<std::string, std::string> &attr) = 0;
	virtual double get_potential() { return 0.0; }
protected:
	const int group_id, neuron_id;
	// TODO: this might be useful context
	//const int mapped_tile_id, mapped_core_id, mapped_core_offset;
};

class LoihiLifModel: public SomaModel
{
public:
        LoihiLifModel(const int gid, const int nid);
	~LoihiLifModel();
	void set_attributes(const std::map<std::string, std::string> &attr);
	NeuronStatus update(const double current_in);
	double get_potential() { return potential; }
private:
	bool force_update;
	int soma_last_updated, reset_mode, reverse_reset_mode;
        int noise_type;
	double potential, leak_decay, bias, threshold, reverse_threshold;
        double reset, reverse_reset, leak_bias;
};

NeuronResetModes model_parse_reset_mode(const std::string &str);
}

#endif
