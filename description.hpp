#ifndef PARSE_HEADER_INCLUDED_
#define PARSE_HEADER_INCLUDED_

#define MAX_FIELD_LEN 64
#define MAX_FIELDS 128
#define MAX_LINE (MAX_FIELDS * MAX_FIELD_LEN)
#define TOKEN_SEPERATORS " \t"
#define DESCRIPTION_MAX_ATTRIBUTES 128

enum description_ret
{
	RET_FAIL = -1,
	RET_OK = 0, // Anything >= 0 means successfully parsed
};

struct attributes
{
	char key[MAX_FIELD_LEN], value_str[MAX_FIELD_LEN];
};

// Forward struct declarations
struct architecture;
struct network;
struct simulation;
#include <stdio.h>

int description_parse_file(FILE *fp, struct network *net, struct architecture *arch);
int description_read_line(char *line, char fields[][MAX_FIELD_LEN], struct network *net, struct architecture *arch);
int description_read_arch_entry(char fields[][MAX_FIELD_LEN], const int field_count, struct architecture *arch);
int description_read_network_entry(char fields[][MAX_FIELD_LEN], const int field_count, struct architecture *arch, struct network *net);
int description_parse_command(char fields[][MAX_FIELD_LEN], const int field_count, struct network *net, struct architecture *arch, struct simulation *sim);

#endif
