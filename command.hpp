// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
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

struct simulation;
struct architecture;
struct network;

int command_parse_input_spikes(struct network *const net, char fields[][MAX_FIELD_LEN], const int field_count);
int command_parse_step_sim(struct network *const net, struct architecture *const arch, struct simulation *sim);

#endif
