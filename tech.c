// tech.c
// Parse "technology" files, which describe key metrics and parameters of
//  a given chip design. Metrics describe the basic chip parameters e.g.
//  number of cores, plus energy and time metrics e.g. energy per spike op
#include <stdio.h>
#include <string.h>
#include "sim.h"  // For INFO()
#include "tech.h"

// Zero initialize all technology parameters
void tech_init(struct technology *tech)
{
	tech->energy_active_neuron_update = 0.0;
	tech->energy_inactive_neuron_update = 0.0;
	tech->energy_spike_op = 0.0;
	tech->energy_spike_within_tile = 0.0;
	tech->energy_east_west_hop = 0.0;
	tech->energy_north_south_hop = 0.0;
	tech->time_active_neuron_update = 0.0;
	tech->time_inactive_neuron_update = 0.0;
	tech->time_spike_op = 0.0;
	tech->time_spike_within_tile = 0.0;
	tech->time_east_west_hop = 0.0;
	tech->time_north_south_hop = 0.0;
	tech->time_mesh_barrier = 0.0;
}

// Read technology file, containing a simple list of parameters and values
void tech_read_file(struct technology *tech, FILE *fp)
{
	char line[TECH_MAX_LINE];

	if (fp == NULL)
	{
		INFO("Error: couldn't read tech file.\n");
		return;
	}

	while (fgets(line, TECH_MAX_LINE, fp))
	{
		tech_read_parameter(tech, line);
	}
}

// Read a line from the technology file, extract at most one parameter
//  Note that lines beginning with # are comments
void tech_read_parameter(struct technology *tech, char *line)
{
	char *metric_str, *value_str;

	metric_str = strtok(line, " \t");
	if ((metric_str == NULL) || (metric_str[0] == '#') ||
							(metric_str[0] == '\n'))
	{
		// Line is a comment or empty; ignore
		return;
	}
	value_str = strtok(NULL, " \t");
	if (value_str == NULL)
	{
		INFO("Error: No value given for %s\n", metric_str);
		return;
	}

	// Note: for now, the case does matter - the field must be all
	//  lowercase
	if (strcmp("energy_active_neuron_update", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->energy_active_neuron_update);
	}
	else if (strcmp("energy_inactive_neuron_update", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->energy_inactive_neuron_update);
	}
	else if (strcmp("energy_spike_op", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->energy_spike_op);
	}
	else if (strcmp("energy_spike_within_tile", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->energy_spike_within_tile);
	}
	else if (strcmp("energy_east_west_hop", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->energy_east_west_hop);
	}
	else if (strcmp("energy_north_south_hop", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->energy_north_south_hop);
	}
	else if (strcmp("time_active_neuron_update", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->time_active_neuron_update);
	}
	else if (strcmp("time_inactive_neuron_update", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->time_inactive_neuron_update);
	}
	else if (strcmp("time_spike_op", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->time_spike_op);
	}
	else if (strcmp("time_spike_within_tile", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->time_spike_within_tile);
	}
	else if (strcmp("time_east_west_hop", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->time_east_west_hop);
	}
	else if (strcmp("time_north_south_hop", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->time_north_south_hop);
	}
	else if (strcmp("time_mesh_barrier", metric_str) == 0)
	{
		sscanf(value_str, "%lf", &tech->time_mesh_barrier);
	}
	else
	{
		INFO("Error: field not recognized (%s).\n", metric_str);
	}

	return;
}

