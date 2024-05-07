#ifndef PLUGINS_HEADER_INCLUDED_
#define PLUGINS_HEADER_INCLUDED_
#include <string>
#include <list>

namespace sanafe
{
struct Attribute;
class SomaModel;

void plugin_init_soma(char* name);
SomaModel *plugin_get_soma(const std::string &name);
}

#endif