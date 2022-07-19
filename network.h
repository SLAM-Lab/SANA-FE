#ifndef NETWORK_HEADER_INCLUDED
#define NETWORK_HEADER_INCLUDED

#define MAX_CSV_LINE (1024 + (4096*32))
#define MAX_FIELD_LEN 32

// This defines the structure of the csv format for specifying neurons in a
//  network for this simulator.   Each row in the file represents data for a
//  unique neuron in hardware.
//
// The format is a number of neuron fields, followed by a variable number of
//  synapses (each with SYNAPSE_FIELDS entries).  This means each line of the
//  csv file is variable length depending on the number of synapses for that
//  neuron.
enum neuron_config_format
{
	NEURON_ID = 0,
	COMPARTMENT_ID, // Which hardware compartment to map to
	THRESHOLD_VOLTAGE, // Spiking threshold value
	RESET_VOLTAGE,
	RECORD_SPIKES, // Record any spikes to a CSV file
	RECORD_VOLTAGE, // Write the voltage to a CSV file every timestep
	FORCE_UPDATE,
	NEURON_FIELDS,
};

enum synapse_config_format
{
	// This is the structure of the CSV format for specifying synaptic
	//  connections in a network for this simulator.  Each row represents
	//  a unique synapse.
	SYNAPSE_DEST_NID = 0,
	SYNAPSE_WEIGHT,
	SYNAPSE_FIELDS,
};

#include "sim.h"
#include "tech.h"
void network_create_empty(const struct technology *tech, struct architecture *arch);
void network_read_csv(FILE *fp, const struct technology *tech, struct architecture *arch, struct compartment **compartment_ptrs);
void network_init(const struct technology *tech, struct architecture *arch);

#endif
