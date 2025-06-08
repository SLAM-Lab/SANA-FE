// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// netlist.hpp
#ifndef NETLIST_HEADER_INCLUDED_
#define NETLIST_HEADER_INCLUDED_

#include <cstddef>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "attribute.hpp"
#include "fwd.hpp"

namespace sanafe
{
SpikingNetwork netlist_parse_file(std::ifstream &fp, Architecture &arch);
void netlist_get_fields(std::vector<std::string_view> &fields, const std::string &line);
size_t field_to_int(const std::string_view &field);
std::pair<std::string, size_t> netlist_parse_neuron_field(const std::string_view &neuron_field);
std::pair<size_t, size_t> netlist_parse_core_field(const std::string_view &core_field);
std::tuple<std::string, size_t, std::string, size_t> netlist_parse_edge_field(const std::string_view &edge_field);
std::tuple<std::string, size_t, size_t, size_t> netlist_parse_mapping_field(const std::string_view &mapping_field);
void netlist_read_network_entry(const std::vector<std::string_view> &fields, Architecture &arch, SpikingNetwork &net, int line_number);

std::string netlist_group_to_netlist(const NeuronGroup &group);
std::string netlist_neuron_to_netlist(const Neuron &neuron, const SpikingNetwork &net, const std::map<std::string, size_t> &group_name_to_id);
std::string netlist_mapping_to_netlist(const Neuron &neuron, const std::map<std::string, size_t> &group_name_to_id);
std::string netlist_connection_to_netlist(const Connection &con, const std::map<std::string, size_t> &group_name_to_id);
std::string netlist_attributes_to_netlist(const std::map<std::string, ModelAttribute> &model_attributes, const std::map<std::string, ModelAttribute> &default_attributes);

void netlist_read_group(const std::vector<std::string_view> &fields, SpikingNetwork &net, int line_number);
void netlist_read_neuron(const std::vector<std::string_view> &fields, SpikingNetwork &net, int line_number);
void netlist_read_edge(const std::vector<std::string_view> &fields, SpikingNetwork &net, int line_number);
void netlist_read_mapping(const std::vector<std::string_view> &fields, Architecture &arch, SpikingNetwork &net, int line_number);

std::variant<bool, int, double, std::string, std::vector<ModelAttribute>> netlist_parse_attribute_value(std::string value_str);
void netlist_parse_attribute_field(const std::string_view &field, std::map<std::string, ModelAttribute> &attributes, int line_number);
std::map<std::string, ModelAttribute> netlist_parse_attributes(const std::vector<std::string_view> &attribute_fields, int line_number);
std::map<std::string, ModelAttribute> netlist_parse_embedded_json(const std::vector<std::string_view> &attribute_fields, int line_number);
char netlist_get_closing_char(char opening_char);
size_t netlist_embedded_json_end_pos(char opening_char, const std::string &all_fields, int line_number);

void add_string_attribute_if_unique(std::string &entry, const std::string &attr_name, const std::string &neuron_value, const std::optional<std::string> &group_default);
void add_bool_attribute_if_unique(std::string &entry, const std::string &attr_name, bool neuron_value, const std::optional<bool> &group_default);

template <typename T>
bool is_unique_attribute(const std::optional<T> &group_default, const T &neuron_value, bool is_empty = false)
{
    if (is_empty)
    {
        return false;
    }
    return !group_default.has_value() ||
            (group_default.value() != neuron_value);
}

}

#endif