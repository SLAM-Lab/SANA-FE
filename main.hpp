#include "sim.hpp"
#include "network.hpp"
#include "arch.hpp"

enum program_args
{
	ARCH_FILENAME,
	NETWORK_FILENAME,
	TIMESTEPS,
	PROGRAM_NARGS,
};

void run(struct simulation *sim, struct network *net, struct architecture *arch);
struct timespec calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);