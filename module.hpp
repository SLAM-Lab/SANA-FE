#include "print.hpp"
#include "sim.hpp"
#include "network.hpp"
#include "arch.hpp"
#include "description.hpp"
#include "command.hpp"
#include "pybind11/pybind11.h"

enum program_args
{
	ARCH_FILENAME,
	NETWORK_FILENAME,
	TIMESTEPS,
	PROGRAM_NARGS,
};


class SANA_FE{
	public:
		struct simulation *sim;
		struct network net;
		struct architecture *arch;
		int timesteps;
		FILE *input_fp;

        SANA_FE();
		void init();
		int force_spike(int group_id, int n_id, int num_spikes);
		void run_timesteps(int timesteps = 1);
		void set_input(char *filename);
		void set_perf_flag(bool flag = true);
		void set_spike_flag(bool flag = true);
		void set_pot_flag(bool flag = true);
		void set_mess_flag(bool flag = true);
		void set_arch(char* filename);
		void set_net(char* filename);
        double get_power();
        void sim_summary();
		void clean_up(description_ret ret = RET_OK);
};

void run(struct simulation *sim, struct network *net, struct architecture *arch);
struct timespec calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);
