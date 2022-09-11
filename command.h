// command.h - This is the low level interface between the user and the
//  simulator. It is essentially a text-based way to control the API function
//  calls that control simulation. Creating spiking neural networks, defining
//  the hardware architecture, mapping SNNs to arch, defining input spikes and
//  stepping through the simulation. These could be run from the command line
//  or from a file, like a primitive script. In this repo, I found it easy
//  to use a higher level script to generate commands and run the experiments.
// TODO: In future, the plan is to support some bindings for other languages
//  e.g. Python
#ifndef COMMAND_HEADER_INCLUDED_
#define COMMAND_HEADER_INCLUDED_

// TODO: better define the max fields
#define MAX_FIELDS (NEURON_FIELDS + (4096*CONNECTION_FIELDS))
#define MAX_NEURONS (128*1024)
#define MAX_FIELD_LEN 32
#define MAX_CSV_LINE (1024 + (4096*MAX_FIELD_LEN))
#define TOKEN_SEPERATORS " \t"

enum
{
	COMMAND_FAIL = -1,
	COMMAND_OK = 0, // Anything >= 0 means successfully parsed
};

int command_parse_file(FILE *fp, struct network *net, struct architecture *arch, FILE *probe_spike_fp, FILE *probe_potential_fp, FILE *perf_fp);
int command_parse_line(char *line, char fields[][MAX_FIELD_LEN], struct network *net, struct architecture *arch, FILE *probe_spike_fp, FILE *probe_potential_fp, FILE *perf_fp);
int command_parse_command(char fields[][MAX_FIELD_LEN], const int field_count, struct network *net, struct architecture *arch, FILE *probe_spike_fp, FILE *probe_potential_fp, FILE *perf_fp);

// (Spiking) neural network
int command_parse_neuron_group(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_neuron(struct network *const net, struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_extern_input_group(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_extern_input_node(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count);

// Neuromorphic hardware
int command_parse_noc(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_tile(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_core(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_axon_input(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_synapse(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_soma(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_axon_output(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);

// Mapping
int command_map_hardware(struct network *const net, struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);

// Simulation control
int command_parse_input_spikes(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_step_sim(struct network *const net, struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count, FILE *probe_spikes_fp, FILE *probe_potential_fp, FILE *perf_fp);
int command_parse_load_commands(struct network *net, struct architecture *arch, char fields[][MAX_FIELD_LEN], const int field_count, FILE *probe_spike_fp, FILE *probe_potential_fp, FILE *perf_fp);

#endif
