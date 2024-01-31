#include "arch.hpp"

struct description_block
{
	char name[MAX_FIELD_LEN];
	struct attributes attributes[ARCH_MAX_ATTRIBUTES];
	struct architecture *arch;
	struct tile *t;
	struct core *c;
	struct description_block *child, *parent, *next_sibling;
	int type, instances, indent, attribute_count;
};

int arch_is_field(const char *fieldname, char *str);
int arch_is_list(char *str);

int arch_get_description_line(FILE *fp, char *line);
int arch_parse_attributes(FILE *fp, struct description_block *block);
struct description_block *arch_parse_block(FILE *fp, const int block_type, struct description_block *const parent, struct description_block *const sibling);
struct description_block *arch_parse_list(FILE *fp, const int list_type, struct description_block *const parent);
int arch_parse_file(FILE *fp, struct architecture *arch, struct description_block *arch_description);
void arch_print_description(struct description_block *arch_description, const int level);
int arch_parse_field_int(char *str, int *val);
int arch_parse_field(char *str, char *field);
int arch_parse_name(char *field, char *name);
int arch_parse_neuron_model(char *model_str);

int arch_build_block(struct description_block *block);
int arch_build_arch(struct architecture *arch, struct description_block *arch_description);
