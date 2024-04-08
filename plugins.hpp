#ifndef PLUGINS_HPP
#define PLUGINS_HPP
#include <string>
#include <list>

#include "description.hpp"

enum NeuronStatus { IDLE, UPDATED, FIRED};

class SomaModel {
public:
	SomaModel(){}
	virtual ~SomaModel(){}
	virtual NeuronStatus update(double input) = 0;
	virtual void set_attributes(const std::list<Attribute> &attr) = 0;
};

typedef SomaModel *_create_soma(void);
typedef void _destroy_soma(SomaModel *model);

void plugin_init_soma(char* name);
SomaModel *plugin_get_soma(const std::string &name);

#endif