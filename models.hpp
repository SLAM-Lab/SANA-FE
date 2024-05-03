#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include <list>
struct Attribute;

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
}

class SomaModel
{
public:
	SomaModel(){}
	virtual ~SomaModel(){}
	virtual sanafe::NeuronStatus update(const double current_in) = 0;
	virtual void set_attributes(const std::list<Attribute> &attr) = 0;
	virtual double get_potential(){ return 0.0; }
};

class LoihiLifModel: public SomaModel
{
public:
        LoihiLifModel();
	~LoihiLifModel();
	void set_attributes(const std::list<Attribute> &attr);
	sanafe::NeuronStatus update(const double current_in);
	double get_potential() { return potential; }
private:
	int soma_last_updated, reset_mode, reverse_reset_mode;
        int noise_type;
        bool force_update;
	double potential, leak_decay, bias, threshold, reverse_threshold;
        double reset, reverse_reset, leak_bias;
};

sanafe::NeuronResetModes model_parse_reset_mode(const std::string &str);

#endif
