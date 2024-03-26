// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// module.cpp - Pybind module interface
// Performance simulation of neuromorphic architectures
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "module.hpp"

namespace py = pybind11;

SANA_FE::SANA_FE(){
	init();
}
void SANA_FE::init(){
	arch = arch_init();
	network_init(&net);
	INFO("Initializing simulation.\n");
	sim = sim_init_sim();
}
int SANA_FE::update_neuron(int group_id, int n_id, vector<string> kwargs, int count){
	for (string item: kwargs){
		INFO("Kwarg: %s\n", item.c_str());
	}
	if (count < 1)
		return -1;
	if (group_id >= net.neuron_group_count)
		return -1;
	if (n_id >= net.groups[group_id].neuron_count)
		return -1;
	struct attributes attr[128];
	for (int i = 0; i < count; ++i){
		string s = kwargs[i];
		strcpy(&attr->key[i], s.substr(0, s.find('=')).c_str());
		strcpy(&attr->value_str[i], s.substr(s.find('=')+1).c_str());
		INFO("neuron: %d.%d updated with key: %s and val: %s\n", group_id, n_id, &attr->key[i], &attr->value_str[i]);
	}
	net.groups[group_id].neurons[n_id].soma_class->parameters(attr, count);
	return 1;
}
void SANA_FE::run_timesteps(int timesteps){
	store_data_init(&run_data, sim, timesteps);
	for (int i = 0; i < timesteps; ++i){
		++sim->ts.timestep;
		run(sim, &net, arch);
	}
	store_data(&run_data, sim);
}
void SANA_FE::set_input(char *filename){
	input_fp = fopen(filename, "r");
	if (input_fp == NULL)
	{
		INFO("Error: Couldn't open inputs %s.\n", filename);
		clean_up(RET_FAIL);
	}
}
void SANA_FE::set_perf_flag(bool flag){
	if (flag){
		sim->log_perf = 1;

		sim->perf_fp = fopen("perf.csv", "w");
		if (sim->perf_fp == NULL)
		{
			INFO("Error: Couldn't open perf file for writing.\n");
			clean_up(RET_FAIL);
		}
		sim_perf_write_header(sim->perf_fp);
	}
	else{
		sim->log_perf = 0;
		if (sim->perf_fp != NULL)
			fclose(sim->perf_fp);
	}
}
void SANA_FE::set_spike_flag(bool flag){
	if (flag){
		sim->log_spikes = 1;

		sim->spike_trace_fp = fopen("spikes.trace", "w");
		if (sim->spike_trace_fp == NULL)
		{
			INFO("Error: Couldn't open trace file for writing.\n");
			clean_up(RET_FAIL);
		}
		sim_spike_trace_write_header(sim);
	}
	else{
		sim->log_spikes = 0;
		if (sim->spike_trace_fp != NULL)
			fclose(sim->spike_trace_fp);
	}
}
void SANA_FE::set_pot_flag(bool flag){
	if (flag){
		sim->log_potential = 1;

		sim->potential_trace_fp = fopen("potential.trace", "w");
		if (sim->potential_trace_fp == NULL)
		{
			INFO("Error: Couldn't open trace file for writing.\n");
			clean_up(RET_FAIL);
		}
		sim_potential_trace_write_header(sim, &net);
	}
	else{
		sim->log_potential = 0;
		if (sim->potential_trace_fp != NULL)
			fclose(sim->potential_trace_fp);
	}
}
void SANA_FE::set_mess_flag(bool flag){
	if (flag){
		sim->log_messages = 1;

		sim->message_trace_fp = fopen("messages.trace", "w");
		if (sim->message_trace_fp == NULL)
		{
			INFO("Error: Couldn't open trace file for writing.\n");
			clean_up(RET_FAIL);
		}
		sim_message_trace_write_header(sim);
	}
	else{
		sim->log_messages = 0;
		if (sim->message_trace_fp != NULL)
			fclose(sim->message_trace_fp);
	}
}
void SANA_FE::set_arch(char* filename){
	FILE* arch_fp = fopen(filename, "r");
	if (arch_fp == NULL)
	{
		INFO("Error: Architecture file %s failed to open.\n", filename);
		clean_up(RET_FAIL);
	}
	int ret = description_parse_file(arch_fp, NULL, arch);
	//arch_print_description(&description, 0);
	fclose(arch_fp);
	if (ret == RET_FAIL)
	{
		clean_up(RET_FAIL);
	}
}
void SANA_FE::set_net(char* filename){
	FILE* network_fp = fopen(filename, "r");
	if (network_fp == NULL)
	{
		INFO("Network data (%s) failed to open.\n", filename);
		clean_up(RET_FAIL);
	}
	INFO("Reading network from file.\n");
	int ret = description_parse_file(network_fp, &net, arch);
	fclose(network_fp);
	if (ret == RET_FAIL)
	{
		clean_up(RET_FAIL);
	}
	network_check_mapped(&net);

	arch_create_connection_maps(arch);

	// Change Potential logging with new headers from net.
	if (sim->log_potential){
		set_pot_flag(true);
	}
}

double SANA_FE::get_power(){
	if (sim->total_sim_time > 0.0)
	{
		return sim->total_energy / sim->total_sim_time;
	}
	else
	{
		return 0.0;
	}
}

void SANA_FE::sim_summary(){
	sim_write_summary(stdout, sim);

	sim->stats_fp = fopen("run_summary.yaml", "w");
	if (sim->stats_fp != NULL)
	{
		sim_write_summary(sim->stats_fp, sim);
	}
}
void SANA_FE::run_summary(){
	print_run_data(stdout, &run_data);

	// Could write intermediate run data here
	// sim->stats_fp = fopen("run_summary.yaml", "w");
	// if (sim->stats_fp != NULL)
	// {
	// 	sim_write_summary(sim->stats_fp, sim);
	// }
}
void SANA_FE::clean_up(description_ret ret){
	// Free any larger structures here
	network_free(&net);
	arch_free(arch);

	// Close any open files here
	if (sim->potential_trace_fp != NULL)
	{
		fclose(sim->potential_trace_fp);
	}
	if (sim->spike_trace_fp != NULL)
	{
		fclose(sim->spike_trace_fp);
	}
	if (sim->message_trace_fp != NULL)
	{
		fclose(sim->message_trace_fp);
	}
	if (sim->perf_fp != NULL)
	{
		fclose(sim->perf_fp);
	}
	if (sim->stats_fp != NULL)
	{
		fclose(sim->stats_fp);
	}

	// Free the simulation structure only after we close all files
	free(sim);

	if (ret == RET_FAIL)
	{
		exit(1);
	}
	else
	{
		exit(0);
	}
}

void run(struct simulation *sim, struct network *net, struct architecture *arch)
{
	// TODO: remove the need to pass the network struct, only the arch
	//  should be needed (since it links back to the net anyway)
	// Run neuromorphic hardware simulation for one timestep
	//  Measure the CPU time it takes and accumulate the stats
	struct timespec ts_start, ts_end, ts_elapsed;
	struct timestep *ts = &(sim->ts);

	// Measure the wall-clock time taken to run the simulation
	//  on the host machine
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	sim_timestep(ts, net, arch);
	// Calculate elapsed time
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

	sim->total_energy += ts->energy;
	sim->total_sim_time += ts->sim_time;
	sim->total_spikes += ts->spike_count;
	sim->total_neurons_fired += ts->total_neurons_fired;
	sim->total_messages_sent += ts->packets_sent;
	if (sim->log_spikes)
	{
		sim_trace_record_spikes(sim, net);
	}
	if (sim->log_potential)
	{
		sim_trace_record_potentials(sim, net);
	}
	if (sim->log_perf)
	{
		sim_perf_log_timestep(ts, sim->perf_fp);
	}
	if (sim->log_messages)
	{
		for (int i = 0; i < ARCH_MAX_CORES; i++)
		{
			for (int j = 0; j < ts->message_queues[i].count; j++)
			{
				struct message *m = &(ts->messages[i][j]);
				if (m->dest_neuron != NULL)
				{
					// Ignore dummy messages (without a
					//  destination). These are inserted to
					//  account for processing that doesn't
					//  result in a spike being sent
					sim_trace_record_message(sim, m);
				}
			}
		}
	}
	sim->timesteps = ts->timestep;
	sim->wall_time += (double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9);
	TRACE1("Time-step took: %fs.\n",
		(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9));

}

struct timespec calculate_elapsed_time(struct timespec ts_start,
							struct timespec ts_end)
{
	// Calculate elapsed wall-clock time between ts_start and ts_end
	struct timespec ts_elapsed;

	ts_elapsed.tv_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	ts_elapsed.tv_sec = ts_end.tv_sec - ts_start.tv_sec;
	if (ts_end.tv_nsec < ts_start.tv_nsec)
	{
		ts_elapsed.tv_sec--;
		ts_elapsed.tv_nsec += 1000000000UL;
	}

	return ts_elapsed;
}

void store_data_init(run_ts_data* data, simulation* sim, int timesteps){
	data->energy = sim->total_energy;
	data->time = sim->total_sim_time;
	data->spikes = sim->total_spikes;
	data->packets = sim->total_messages_sent;
	data->neurons = sim->total_neurons_fired;
	data->wall_time = sim->wall_time;
	data->timestep_start = sim->timesteps;
	data->timesteps = timesteps;
}

void store_data(run_ts_data* data, simulation* sim){
	data->energy = sim->total_energy - data->energy;
	data->time = sim->total_sim_time - data->time;
	data->spikes = sim->total_spikes - data->spikes;
	data->packets = sim->total_messages_sent - data->packets;
	data->neurons = sim->total_neurons_fired - data->neurons;
	data->wall_time = sim->wall_time - data->wall_time;
}

void print_run_data(FILE *fp, run_ts_data* data){
	fprintf(fp, "energy: %e\n", data->energy);
	fprintf(fp, "time: %e\n", data->time);
	fprintf(fp, "total_spikes: %ld\n", data->spikes);
	fprintf(fp, "total_packets: %ld\n", data->packets);
	fprintf(fp, "total_neurons_fired: %ld\n", data->neurons);
	fprintf(fp, "wall_time: %lf\n", data->wall_time);
	fprintf(fp, "executed from: %ld to %ld timesteps\n", data->timestep_start,
	data->timestep_start+data->timesteps);
}

void test_pybind(void){
	INFO("Printing through Pybind!\n");
}

PYBIND11_MODULE(simcpp, m) {
    m.doc() = R"pbdoc(
        SANA-FE Cpp Module with Pybind11 
        --------------------------------

        .. currentmodule:: simcpp

        .. autosummary::
           :toctree: _generate

           test_pybind
		   SANA_FE
    )pbdoc";

    m.def("test_pybind", &test_pybind, R"pbdoc(
        test_pybind function from main.cpp

        Test pybind11 functionality.
    )pbdoc");

	py::class_<SANA_FE>(m, "SANA_FE")
		.def(py::init())
		.def("init", &SANA_FE::init)
        .def("update_neuron", &SANA_FE::update_neuron)
        .def("run_timesteps", &SANA_FE::run_timesteps, py::arg("timesteps") = 1)
		.def("set_input", &SANA_FE::set_input)
		.def("set_perf_flag", &SANA_FE::set_perf_flag, py::arg("flag") = true)
		.def("set_spike_flag", &SANA_FE::set_spike_flag, py::arg("flag") = true)
		.def("set_pot_flag", &SANA_FE::set_pot_flag, py::arg("flag") = true)
		.def("set_mess_flag", &SANA_FE::set_mess_flag, py::arg("flag") = true)
		.def("set_arch", &SANA_FE::set_arch)
		.def("set_net", &SANA_FE::set_net)
		.def("get_power", &SANA_FE::get_power)
		.def("sim_summary", &SANA_FE::sim_summary)
		.def("run_summary", &SANA_FE::run_summary)
		.def("clean_up", &SANA_FE::clean_up, py::arg("ret") = 1);
}
