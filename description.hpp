#ifndef PARSE_HEADER_INCLUDED_
#define PARSE_HEADER_INCLUDED_

#define DEFAULT_LINE_LEN  (4096)

#include <iostream>
#include <fstream>
#include <vector>

enum description_ret
{
	RET_FAIL = -1,
	RET_OK = 0, // Anything >= 0 means successfully parsed
};

struct attribute
{
	std::string key, value_str;
};

// Forward struct declarations
struct architecture;
struct network;
struct simulation;
#include <stdio.h>

int description_parse_file(std::fstream &fp, struct architecture &arch);
int description_parse_file(std::fstream &fp, struct network &net, struct architecture &arch);
int description_read_line(const std::string line, std::vector<std::string> fields, struct network &net, struct architecture &arch);
std::vector<std::string> description_get_fields(const std::string &line);
int description_read_arch_entry(const std::vector<std::string> &fields, struct architecture &arch);
int description_read_network_entry(const std::vector<std::string> &fields, struct architecture &arch, struct network &net);
//int description_parse_command(char fields[][MAX_FIELD_LEN], const int field_count, struct network *net, struct architecture *arch, struct simulation *sim);

#endif
