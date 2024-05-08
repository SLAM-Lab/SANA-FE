#ifndef PLUGINS_HEADER_INCLUDED_
#define PLUGINS_HEADER_INCLUDED_
#include <string>

namespace sanafe
{
class SomaModel;

void plugin_init_soma(char* name);
SomaModel *plugin_get_soma(const std::string &name);
}

#endif