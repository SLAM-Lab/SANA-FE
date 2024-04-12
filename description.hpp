#ifndef PARSE_HEADER_INCLUDED_
#define PARSE_HEADER_INCLUDED_

#define DEFAULT_LINE_LEN  (4096)

#include <iostream>
#include <fstream>
#include <vector>

enum DescriptionRet
{
	RET_FAIL = -1,
	RET_OK = 0, // Anything >= 0 means successfully parsed
};

struct Attribute
{
	std::string key, value_str;
};

// Forward struct declarations
struct Architecture;
struct Network;
struct Simulation;
#include <stdio.h>

int description_parse_arch_file(std::fstream &fp, Architecture &arch);
int description_parse_net_file(std::fstream &fp, Network &net, Architecture &arch);
int description_read_line(const std::string line, std::vector<std::string> fields, Network &net, Architecture &arch);
std::vector<std::string> description_get_fields(const std::string &line);
int description_read_arch_entry(const std::vector<std::string> &fields, Architecture &arch, const int line_number);
int description_read_network_entry(const std::vector<std::string> &fields, Architecture &arch, Network &net, const int line_number);
//int description_parse_command(char fields[][MAX_FIELD_LEN], const int field_count, Network *net, Architecture *arch, struct simulation *sim);

#endif
