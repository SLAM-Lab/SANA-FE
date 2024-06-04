#ifndef DESCRIPTION_HEADER_INCLUDED_
#define DESCRIPTION_HEADER_INCLUDED_

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string_view>
#include <vector>

namespace sanafe
{
enum DescriptionRet
{
    RET_FAIL = -1,
    RET_OK = 0, // Anything >= 0 means successfully parsed
};

// Forward struct declarations
class Simulation;
class Architecture;
class Network;
class NeuronGroup;
class Neuron;
class Tile;
class Core;

int description_parse_arch_file(std::ifstream &fp, Architecture &arch);
int description_parse_net_file(std::ifstream &fp, Network &net, Architecture &arch);
void description_get_fields(std::vector<std::string_view> &fields, const std::string &line);
void description_read_arch_entry(const std::vector<std::string_view> &fields, Architecture &arch, const int line_number);
void description_read_network_entry(const std::vector<std::string_view> &fields, Architecture &arch, Network &net, const int line_number);
void parse_neuron_field(const std::string_view &neuron_field, size_t &group_id, size_t &neuron_id);
void parse_core_field(const std::string_view &core_field, size_t &tile_id, size_t &core_offset);
void parse_edge_field(const std::string_view &edge_field, size_t &group_id, size_t &neuron_id, std::optional<size_t> &dendrite_id, size_t &dest_group_id, size_t &dest_neuron_id, std::optional<size_t> &dest_dendrite_id);
void parse_neuron_with_compartment_field(const std::string_view &compartment_field, size_t &group_id, size_t &neuron_id, std::optional<size_t> &compartment_id);
void parse_mapping_field(const std::string_view &mapping_field, size_t &group_id, size_t &neuron_id, size_t &tile_id, size_t &core_offset);
size_t field_to_int(const std::string_view &field);

}

#endif
