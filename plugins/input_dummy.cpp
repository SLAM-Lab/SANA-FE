#include <string.h>
#include "../print.hpp"
#include "../plugins.hpp"
#include "../arch.hpp"
#include <iostream>
#include <math.h>
#include <queue>

using namespace std;

void clear( queue<Neuron_Status> &q )
{
   queue<Neuron_Status> empty;
   swap( q, empty );
}

class input_dummy: public Base_Soma {
    public:
        int inputs_received, outputs_sent;
        queue<Neuron_Status> output_spikes;
        Neuron_Status default_out, out_spike;
        int hold_mode; // Whether to hold output or go back to default after
        // 0 for specified num of spikes mode, 1 for hold mode
        input_dummy(){
            inputs_received = 0;
            outputs_sent = 0;
            default_out = IDLE;
            out_spike = IDLE;
            hold_mode = 0;
        }
        virtual void parameters(struct attributes* attr, const int attribute_count) {
            /*** Set attributes ***/
            for (int i = 0; i < attribute_count; i++)
            {
                struct attributes *a = &(attr[i]);
                int ret = 1;

                if (strncmp("hold_mode", a->key, MAX_FIELD_LEN) == 0)
                {
                    ret = sscanf(a->value_str, "%d", &hold_mode);
                    clear(output_spikes);
                }
                else if (strncmp("default_out", a->key, MAX_FIELD_LEN) == 0){
                    int res;
                    ret = sscanf(a->value_str, "%d", &res);
                    default_out = Neuron_Status(res);
                }
                else if (strncmp("input_spike", a->key, MAX_FIELD_LEN) == 0){
                    ++inputs_received;
                    int res;
                    ret = sscanf(a->value_str, "%d", &res);
                    if (!hold_mode){
                        output_spikes.push(Neuron_Status(res));
                    }
                    out_spike = Neuron_Status(res);
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
            ++outputs_sent;
            Neuron_Status out;
            if (hold_mode){
                out = out_spike;
            }
            else if (output_spikes.size() > 0){
                out = output_spikes.front();
                output_spikes.pop();
            }
            else{
                out = default_out;
            }
            return out;
        }
};

// the class factories

extern "C" Base_Soma* create_input_dummy() {
    return new input_dummy();
}

extern "C" void destroy_input_dummy(Base_Soma* neuron) {
    delete neuron;
}