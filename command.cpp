// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// command.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "print.hpp"
#include "arch.hpp"
#include "network.hpp"
#include "command.hpp"
#include "sim.hpp"

int command_parse_input_spikes(struct network *const net,
	char fields[][MAX_FIELD_LEN], const int field_count)
{
	double val;
	int ret, input_count;

	input_count = field_count - 1;
	if (input_count > net->external_input_count)
	{
		INFO("Error: Too many inputs given.\n");
		return RET_FAIL;
	}

	for (int i = 0; i < input_count; i++)
	{
		struct input *in = &(net->external_inputs[i]);
		ret = sscanf(fields[i + 1], "%lf", &val);
		if (ret < 1)
		{
			INFO("Error: Couldn't parse input value (%s)\n",
				fields[i + 1]);
			return RET_FAIL;
		}
		else if (val < 0)
		{
			INFO("Warning: id:%d input rate < 0 (%lf)\n", in->id,
				in->spike_val);
		}
		//else if (val > 1.0)
		//{
		//	INFO("Warning: id:%d input rate > 1 (%lf)\n",
		//					in->id, in->spike_val);
		//}
		network_set_input(net, i, val);
		INFO("Parsed input %d=%lf\n", in->id, in->spike_val);
	}

	return RET_OK;
}

int command_parse_step_sim(struct network *const net,
	struct architecture *const arch,
	//						char fields[][MAX_FIELD_LEN],
	//						const int field_count,
	struct simulation *sim)
{
	INFO("not implemented\n");
	exit(1);
	//sim_timestep(sim, net, arch, scheduler);
	return RET_OK;
}
