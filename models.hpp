#ifndef MODELS_HEADER_INCLUDED_
#define MODELS_HEADER_INCLUDED_

#include <list>
#include "plugins.hpp"
#include "arch.hpp"

#define MAX_NOISE_FILE_ENTRY 128

class Loihi_Lif_Model: public Soma_Model
{
public:
        Loihi_Lif_Model();
	~Loihi_Lif_Model();
	void set_attributes(const std::list<attribute> &attr);
	Neuron_Status update(const double current_in);
private:
	int soma_last_updated, reset_mode, reverse_reset_mode;
        int noise_type;
        bool force_update;
	double potential, leak_decay, bias, threshold, reverse_threshold;
        double reset, reverse_reset, leak_bias;
};

#endif
