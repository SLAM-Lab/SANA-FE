#include "print.hpp"
#include "sim.hpp"
#include "network.hpp"
#include "arch.hpp"
#include "description.hpp"
#include "command.hpp"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

using namespace std;

enum program_args
{
	ARCH_FILENAME,
	NETWORK_FILENAME,
	TIMESTEPS,
	PROGRAM_NARGS,
};

struct run_ts_data{
	double energy;
	double time;
	double wall_time;
	long int spikes, packets, neurons;
	long int timestep_start, timesteps;
};


class SANA_FE{
	public:
		struct simulation *sim;
		struct network net;
		struct architecture *arch;
		int timesteps;
		FILE *input_fp;
		struct run_ts_data run_data;

        SANA_FE();
		void init();
		int update_neuron(int group_id, int n_id, vector<string> kwargs, int count);
		void run_timesteps(int timesteps = 1);
		void set_input(char *filename);
		void set_perf_flag(bool flag = true);
		void set_spike_flag(bool flag = true);
		void set_pot_flag(bool flag = true);
		void set_mess_flag(bool flag = true);
		void set_gui_flag(bool flag = true);
		void set_arch(char* filename);
		void set_net(char* filename);
        double get_power();
		vector<int> get_status(int gid);
        void sim_summary();
		vector<vector<int>> run_summary();
		void clean_up(int ret = RET_OK);
		~SANA_FE(){clean_up();};
};
class Vector_Cleanup_Class{
	public:
		struct architecture* arch;
		Vector_Cleanup_Class(){}
		Vector_Cleanup_Class(struct architecture* arch){
			this->arch = arch;
		}
		~Vector_Cleanup_Class(){
			arch->spike_vector.clear();
		}
};

void run(struct simulation *sim, struct network *net, struct architecture *arch);
struct timespec calculate_elapsed_time(struct timespec ts_start, struct timespec ts_end);
void store_data_init(run_ts_data* data, simulation* sim, int timesteps);
void store_data(run_ts_data* data, simulation* sim);
void print_run_data(FILE *fp, run_ts_data* data);
