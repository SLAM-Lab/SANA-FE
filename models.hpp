#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include <list>
#include "plugins.hpp"
#include "arch.hpp"

#define MAX_NOISE_FILE_ENTRY 128

class LoihiLifModel: public SomaModel
{
public:
        LoihiLifModel();
	~LoihiLifModel();
	void set_attributes(const std::list<Attribute> &attr);
	NeuronStatus update(const double current_in);
private:
	int soma_last_updated, reset_mode, reverse_reset_mode;
        int noise_type;
        bool force_update;
	double potential, leak_decay, bias, threshold, reverse_threshold;
        double reset, reverse_reset, leak_bias;
};

int model_parse_reset_mode(const std::string &str);

#endif
