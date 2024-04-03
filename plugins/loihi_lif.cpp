#include <string.h>
#include "../print.hpp"
#include "../plugins.hpp"
#include "../arch.hpp"
#include <iostream>
#include <math.h>

using namespace std;

int network_parse_reset_mode(const char *str)
{
	int reset_mode = -1;

	if (strcmp(str, "none") == 0)
	{
		reset_mode = NEURON_NO_RESET;
	}
	else if (strcmp(str, "soft") == 0)
	{
		reset_mode = NEURON_RESET_SOFT;
	}
	else if (strcmp(str, "hard") == 0)
	{
		reset_mode = NEURON_RESET_HARD;
	}
	else if (strcmp(str, "saturate") == 0)
	{
		reset_mode = NEURON_RESET_SATURATE;
	}
	else
	{
		INFO("Error: reset mode not recognized.");
		exit(1);
	}

	return reset_mode;
}


class loihi_lif: public Base_Soma {

    // LIF specific
    public:
        unsigned int random_range_mask;
        double potential, current, charge, bias;
        double reset, reverse_reset, threshold, reverse_threshold;
        double leak_decay, leak_bias, potential_time_const;
        // double dendritic_current_decay, processing_latency;
        int reset_mode, reverse_reset_mode;

        loihi_lif(){
            reset = 0.0;
            reverse_reset = 0.0;
            threshold = 1.0;
            reverse_threshold = -1.0;

            // Default is no leak (potential decay), i.e., the potential for the
            //  next timestep is 100% of the previous timestep's
            leak_decay = 1.0;
            leak_bias = 0.0;
            potential_time_const = 0.0;

            random_range_mask = 0;
            potential = 0.0;
            current = 0.0;
            // charge = 0.0;
            bias = 0.0;
            
            // dendritic_current_decay = 0.0;
            reset_mode = 0;
            reverse_reset_mode = 0;
        }

        virtual void parameters(struct attributes* attr, const int attribute_count) {
            /*** Set attributes ***/
            for (int i = 0; i < attribute_count; i++)
            {
                struct attributes *a = &(attr[i]);
                int ret = 1;

                if (strncmp("bias", a->key, MAX_FIELD_LEN) == 0)
                {
                    ret = sscanf(a->value_str, "%lf", &bias);
                }
                else if (strncmp("reset", a->key, MAX_FIELD_LEN) == 0)
                {
                    ret = sscanf(a->value_str, "%lf", &reset);
                }
                else if (strncmp("reverse_reset", a->key, MAX_FIELD_LEN) == 0)
                {
                    ret = sscanf(a->value_str, "%lf", &reverse_reset);
                }
                else if (strncmp("threshold", a->key, MAX_FIELD_LEN) == 0)
                {
                    ret = sscanf(a->value_str, "%lf", &threshold);
                }
                else if (strncmp("reverse_threshold", a->key, MAX_FIELD_LEN) ==
                    0)
                {
                    ret = sscanf(
                        a->value_str, "%lf", &reverse_threshold);
                }
                else if (strncmp("leak_decay", a->key, MAX_FIELD_LEN) == 0)
                {
                    ret = sscanf(a->value_str, "%lf",
                        &leak_decay);
                }
                else if (strncmp("leak_bias", a->key, MAX_FIELD_LEN) == 0)
                {
                    ret = sscanf(
                        a->value_str, "%lf", &leak_bias);
                }
                else if (strncmp("reset_mode", a->key, MAX_FIELD_LEN) == 0)
                {
                    reset_mode =
                        network_parse_reset_mode(a->value_str);
                    // Was parsed successfully if we got here
                    ret = 1;
                }
                else if (strncmp("reverse_reset_mode", a->key, MAX_FIELD_LEN) ==
                    0)
                {
                    reverse_reset_mode =
                        network_parse_reset_mode(a->value_str);
                    // Was parsed successfully if we got here
                    ret = 1;
                }
                else if (strncmp("input_spike", a->key, MAX_FIELD_LEN) == 0){
                    double res;
                    ret = sscanf(a->value_str, "%lf", &res);
                    potential += res;
                }
                if (ret < 1)
                {
                    INFO("Invalid attribute (%s:%s)\n", a->key,
                        a->value_str);
                    exit(1);
                }
            }
        }

        Neuron_Status update_soma(double input){

            Neuron_Status neuron_status = IDLE;

        	// Calculate the change in potential since the last update e.g.
        	//  integate inputs and apply any potential leak
        	TRACE1("Updating potential, before:%f\n", potential);

            potential *= leak_decay;

        	// Add the synaptic / dendrite current to the potential
        	//printf("n->bias:%lf n->potential before:%lf current_in:%lf\n", n->bias, n->potential, current_in);
        	potential += input + bias;
        	charge = 0.0;
        	//printf("n->bias:%lf n->potential after:%lf\n", n->bias, n->potential);

        	TRACE1("Updating potential, after:%f\n", potential);

            if ((fabs(potential) > 0.0) || (fabs(bias) > 0.0))
            {
                neuron_status = UPDATED;
            }

        	// Check against threshold potential (for spiking)
        	if (((bias != 0.0) && (potential > threshold)) ||
        		((bias == 0.0) && (potential >= threshold)))
        	{
        		if (reset_mode == NEURON_RESET_HARD)
        		{
        			potential = reset;
        		}
        		else if (reset_mode == NEURON_RESET_SOFT)
        		{
        			potential -= threshold;
        		}
                neuron_status = FIRED;
        	}

        	// Check against reverse threshold
        	if (potential < reverse_threshold)
        	{
        		if (reverse_reset_mode == NEURON_RESET_SOFT)
        		{
        			potential -= reverse_threshold;
        		}
        		else if (reverse_reset_mode == NEURON_RESET_HARD)
        		{
        			potential = reverse_reset;
        		}
        		else if (reverse_reset_mode == NEURON_RESET_SATURATE)
        		{
        			potential = reverse_threshold;
        		}
        	}

            return neuron_status;
        }
};

// the class factories

extern "C" Base_Soma* create_loihi_lif() {
    return new loihi_lif();
}

extern "C" void destroy_loihi_lif(Base_Soma* lif) {
    delete lif;
}