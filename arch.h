#ifndef ARCH_HEADER_INCLUDED
#define ARCH_HEADER_INCLUDED

// TODO: for now hard define max numbers of hardware blocks. I did this so I
//  don't have quite as much dynamic allocation going on. If needed this can
//  fairly easily be removed and replaced with allocating arbitrary numbers of
//  elements off the heap e.g. in a linked list
#define ARCH_MAX_CORES 16
#define ARCH_MAX_LINKS 16
#define ARCH_MAX_PROCESSORS 4

#define ARCH_LINE_LEN 512
#define ARCH_MAX_VALUES 128
#define ARCH_MAX_VALUE_DIGITS 8
#define ARCH_MAX_FIELDS 8
#define ARCH_MAX_FIELD_LEN ((ARCH_MAX_VALUE_DIGITS * 2) + 2) // min..max

struct axon_input
{
	struct tile *t;
	int id, packet_size;
};

struct synapse_processor
{
	int id;
	double energy, time;
	double energy_spike_op, time_spike_op;
};

struct dendrite_processor
{
	int id;
	double energy, time;
};

struct soma_processor
{
	int id;
	double energy, time;
	double energy_active_neuron_update, time_active_neuron_update;
	double energy_inactive_neuron_update, time_inactive_neuron_update;
};

struct axon_output
{
	struct tile *t;
	int id;
	long int total_packets_sent;
	double energy, time;
	double energy_spike_within_tile, time_spike_within_tile;
};

// TODO: separate the spiking network from the architecture stuff i.e. SNN
//  from hardware
struct core
{
	struct axon_input axon_in[ARCH_MAX_PROCESSORS]; 
	struct synapse_processor synapse[ARCH_MAX_PROCESSORS];
	struct dendrite_processor dendrite[ARCH_MAX_PROCESSORS]; 
	struct soma_processor soma[ARCH_MAX_PROCESSORS];
	struct axon_output axon_out[ARCH_MAX_PROCESSORS];
	double energy, time;
	int core_count, id;
};

struct tile
{
	struct core cores[ARCH_MAX_CORES];
	// TODO: maybe can associate energy and latency with each link, that
	//  will be the most general way to implement this!
	struct tile *links[ARCH_MAX_LINKS];
	double energy, time;
	double energy_east_west_hop, time_east_west_hop;
	double energy_north_south_hop, time_north_south_hop;
	int id, x, y, core_count;
	int max_dimensions, width; // For now just support 2 dimensions
};

struct architecture
{
	struct tile *tiles;
	double time_barrier;
	int tile_count, initialized;
};

struct range
{
	unsigned int min, max;
};

//struct architecture arch_read_file(FILE *fp);
void arch_free(struct architecture *arch);
int arch_create_tile(void);
int arch_create_core(void);

//static void arch_read_router(struct architecture *arch, char *line);
//static void arch_read_timer(struct architecture *arch, char *line);

#endif
