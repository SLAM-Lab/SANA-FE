#ifndef DESCRIPTION_HEADER_INCLUDED_
#define DESCRIPTION_HEADER_INCLUDED_

#define DEFAULT_LINE_LEN  (4096)

#include <iostream>
#include <fstream>
#include <vector>
#include <string_view>
#include <cstdio>

namespace sanafe
{
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
class Simulation;
class Architecture;
struct Network;
struct NeuronGroup;
struct Neuron;
struct Tile;
struct Core;

int description_parse_arch_file(std::ifstream &fp, Architecture &arch);
int description_parse_net_file(std::ifstream &fp, Network &net, Architecture &arch);
int description_read_line(const std::string line, std::vector<std::string> fields, Network &net, Architecture &arch);
void description_get_fields(std::vector<std::string_view> &fields, const std::string &line);
int description_read_arch_entry(const std::vector<std::string_view> &fields, Architecture &arch, const int line_number);
int description_read_network_entry(const std::vector<std::string_view> &fields, Architecture &arch, Network &net, const int line_number);
//int description_parse_command(char fields[][MAX_FIELD_LEN], const int field_count, Network *net, Architecture *arch, struct simulation *sim);
void parse_neuron_field(const std::string_view &neuron_field, std::vector<NeuronGroup>::size_type &group_id, std::vector<Neuron>::size_type &neuron_id);
void parse_core_field(const std::string_view &core_field, std::vector<NeuronGroup>::size_type &tile_id, std::vector<Neuron>::size_type &core_offset);
void parse_edge_field(const std::string_view &edge_field, std::vector<NeuronGroup>::size_type &group_id, std::vector<Neuron>::size_type &neuron_id, std::vector<NeuronGroup>::size_type &dest_group_id, std::vector<Neuron>::size_type &dest_neuron_id);
void parse_mapping_field(const std::string_view &mapping_field, std::vector<NeuronGroup>::size_type &group_id, std::vector<Neuron>::size_type &neuron_id, std::vector<Tile>::size_type &tile_id, std::vector<Core>::size_type &core_offset);
size_t field_to_int(const std::string_view &field);

}

#endif
