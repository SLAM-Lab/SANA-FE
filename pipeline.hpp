// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  pipeline.hpp
namespace sanafe
{
class Architecture;
struct Timestep;

void pipeline_process_neurons(Timestep &ts, Architecture &arch);
void pipeline_process_messages(Timestep &ts, Architecture &arch);

void pipeline_process_neuron(Timestep &ts, const Architecture &arch, Neuron &n);
void pipeline_receive_message(Architecture &arch, Message &m);
double pipeline_process_message(const Timestep &ts, Core &c, Message &m);

double pipeline_process_axon_in(Core &core, const Message &m);
double pipeline_process_synapse(const Timestep &ts, Connection &con);
double pipeline_process_dendrite(const Timestep &ts, Neuron &n);
double pipeline_process_soma(const Timestep &ts, Neuron &n);
double pipeline_process_axon_out(Timestep &ts, const Architecture &arch, Neuron &n);
BufferPosition pipeline_parse_buffer_pos_str(const std::string &buffer_pos_str);

std::pair<double, double> pipeline_apply_default_dendrite_power_model(Neuron &n, std::optional<double> energy, std::optional<double> latency);
size_t abs_diff(size_t a, size_t b);
}
