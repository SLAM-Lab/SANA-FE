#ifndef ARCH_HEADER_INCLUDED
#define ARCH_HEADER_INCLUDED

struct architecture
{
	int max_neurons, max_mem_blocks, max_routers, max_timers;
	int max_axon_inputs, max_axon_outputs, max_external_inputs;
        int initialized;
	struct neuron *neurons;
	struct synapse_mem *mem_blocks;
	struct router *routers;
	struct timer *timers;
	struct axon_input *axon_inputs;
	struct axon_output *axon_outputs;
	struct input *external_inputs;
};

struct range
{
	unsigned int min, max;
};

#define ARCH_LINE_LEN 512
#define ARCH_MAX_VALUES 128
#define ARCH_MAX_VALUE_DIGITS 8
#define ARCH_MAX_FIELDS 8
#define ARCH_MAX_FIELD_LEN ((ARCH_MAX_VALUE_DIGITS * 2) + 2) // min..max

struct architecture arch_read_file(FILE *fp);
void arch_free(struct architecture *arch);

//static void arch_read_router(struct architecture *arch, char *line);
//static void arch_read_timer(struct architecture *arch, char *line);

#endif
