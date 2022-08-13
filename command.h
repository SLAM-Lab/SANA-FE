#ifndef COMMAND_HEADER_INCLUDED
#define COMMAND_HEADER_INCLUDED

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

int command_read_file(FILE *fp, struct network *net, struct architecture *arch);
int command_parse_line(char *line, char fields[][MAX_FIELD_LEN], struct network *net, struct architecture *arch);
int command_parse_command(char fields[][MAX_FIELD_LEN], const int field_count, struct network *net, struct architecture *arch);
int command_parse_neuron_group(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_neuron(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_noc(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_tile(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_core(struct architecture *const arch, char fields[][MAX_FIELD_LEN], const int field_count);

#endif