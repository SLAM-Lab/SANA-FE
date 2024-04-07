#ifndef PLUGINS_HPP
#define PLUGINS_HPP
#include <string>
#include <list>

#include "description.hpp"

enum Neuron_Status { IDLE, UPDATED, FIRED};

class Soma_Model {
public:
	Soma_Model(){}
	virtual ~Soma_Model(){}
	virtual Neuron_Status update(double input) = 0;
	virtual void set_attributes(const std::list<attribute> &attr) = 0;
};

typedef Soma_Model *_create_soma(void);
typedef void _destroy_soma(Soma_Model *model);

void plugin_init_soma(char* name);
Soma_Model *plugin_get_soma(const std::string &name);

#endif