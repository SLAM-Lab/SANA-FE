#ifndef PLUGINS_HPP
#define PLUGINS_HPP
#include "description.hpp"

#define MAX_SOMA_CLASSES MAX_FIELD_LEN
#define MAX_SOMA_LEN MAX_FIELD_LEN

enum Neuron_Status { IDLE, UPDATED, FIRED};

class Base_Soma {
public:
	Base_Soma(){}
    virtual ~Base_Soma(){}
	virtual Neuron_Status update_soma(double input) = 0;
	virtual void parameters(struct attributes* attr, const int attribute_count) = 0;
};

typedef Base_Soma* _create_soma();
typedef void _destroy_soma(Base_Soma*);

void init_soma(char* name);
Base_Soma* get_soma(char* name);

#endif