// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  pipeline.hpp
namespace sanafe
{
class SpikingHardware;
class Core;

struct Timestep;
struct MappedNeuron;
struct Message;
struct MappedConnection;
enum BufferPosition;

void pipeline_process_neurons(Timestep &ts, SpikingHardware &hw);
void pipeline_process_messages(Timestep &ts, SpikingHardware &hw);

void pipeline_process_neuron(Timestep &ts, const SpikingHardware &arch, MappedNeuron &n);
void pipeline_receive_message(SpikingHardware &arch, Message &m);
double pipeline_process_message(const Timestep &ts, Core &c, Message &m);

double pipeline_process_axon_in(Core &core, const Message &m);
double pipeline_process_synapse(const Timestep &ts, MappedConnection &con);
double pipeline_process_dendrite(const Timestep &ts, MappedNeuron &n);
double pipeline_process_soma(const Timestep &ts, MappedNeuron &n);
double pipeline_process_axon_out(Timestep &ts, const SpikingHardware &arch, MappedNeuron &n);
BufferPosition pipeline_parse_buffer_pos_str(const std::string &buffer_pos_str);

std::pair<double, double> pipeline_apply_default_dendrite_power_model(MappedNeuron &n, std::optional<double> energy, std::optional<double> latency);
size_t abs_diff(size_t a, size_t b);
}
