#ifndef TECH_HEADER_INCLUDED
#define TECH_HEADER_INCLUDED

#define TECH_MAX_LINE 128

struct technology
{
	int max_cores, cores_per_tile, max_compartments, fan_out;
	double energy_active_neuron_update, energy_inactive_neuron_update;
	double energy_spike_op, energy_spike_within_tile;
	double energy_east_west_hop, energy_north_south_hop;
	double time_active_neuron_update, time_inactive_neuron_update;
	double time_spike_op, time_spike_within_tile;
	double time_east_west_hop, time_north_south_hop;
	double time_mesh_barrier;
};

void tech_init(struct technology *tech);
void tech_read_file(struct technology *tech, FILE *fp);
void tech_read_parameter(struct technology *tech, char *line);

#endif
