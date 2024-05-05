// Copyright (c) 2023 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  sim.c
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
//#include <omp.h>
#include <vector>
#include <list>
#include <memory>

#include "print.hpp"
#include "sim.hpp"
#include "network.hpp"
#include "arch.hpp"

SanaFe::SanaFe()
{
	init();
}

void SanaFe::init()
{
	arch = new Architecture();
	net = Network();
	INFO("Initializing simulation.\n");
	sim = sim_init_sim();
}

int SanaFe::update_neuron(
	std::vector<NeuronGroup>::size_type group_id,
	std::vector<Neuron>::size_type n_id,
	std::vector<string> kwargs, int count)
{
	for (string item: kwargs)
	{
		INFO("Kwarg: %s\n", item.c_str());
	}
	if (count < 1)
	{
		return -1;
	}
	if (group_id >= net.groups.size())
	{
		return -1;
	}
	if (n_id >= net.groups[group_id].neurons.size())
	{
		return -1;
	}

	std::list<Attribute> attr;
	for (auto s: kwargs)
	{
		std::string key = s.substr(0, s.find('=')).c_str();
		std::string value_str = s.substr(s.find('=')+1).c_str();
		Attribute a = { key, value_str };
		attr.push_back(a);
		INFO("neuron: %lu.%lu updated with key: %s and val: %s\n",
			group_id, n_id, key.c_str(), value_str.c_str());
	}
	net.groups[group_id].neurons[n_id].model->set_attributes(attr);
	return 0;
}

void SanaFe::run_timesteps(int timesteps)
{
	store_data_init(&run_data, *sim, timesteps);
	for (int i = 0; i < timesteps; ++i)
	{
		sim->timesteps++;
		run(*sim, net, *arch);
	}
	store_data(&run_data, *sim);
}

void SanaFe::set_spike_trace(const bool enable)
{
	sim->log_spikes = enable;
}

void SanaFe::set_potential_trace(const bool enable)
{
	sim->log_potential = enable;
}

void SanaFe::set_perf_trace(const bool enable)
{
	sim->log_perf = enable;
}

void SanaFe::set_message_trace(const bool enable)
{
	sim->log_messages = enable;
}

void SanaFe::open_perf_trace(void)
{
	sim->log_perf = 1;

	// TODO: warning, this is specific to Linux
	// To be more portable, consider using the filesystem library in C++17
	std::string perf_path = out_dir + std::string("/perf.csv");
	sim->perf_fp = fopen(perf_path.c_str(), "w");
	if (sim->perf_fp == NULL)
	{
		INFO("Error: Couldn't open perf file for writing.\n");
		clean_up(RET_FAIL);
	}
	sim_perf_write_header(sim->perf_fp);

}

void SanaFe::open_spike_trace(void)
{
	sim->log_spikes = 1;

	std::string spike_path = out_dir + std::string("/spike.csv");
	sim->spike_trace_fp = fopen(spike_path.c_str(), "w");
	if (sim->spike_trace_fp == NULL)
	{
		INFO("Error: Couldn't open trace file for writing.\n");
		clean_up(RET_FAIL);
	}
	sim_spike_trace_write_header(*sim);
}

void SanaFe::open_potential_trace(void)
{
	sim->log_potential = 1;

	std::string potential_path = out_dir + "/potential.csv";
	sim->potential_trace_fp = fopen(potential_path.c_str(), "w");
	if (sim->potential_trace_fp == NULL)
	{
		INFO("Error: Couldn't open trace file for writing.\n");
		clean_up(RET_FAIL);
	}
	sim_potential_trace_write_header(*sim, net);
}

void SanaFe::open_message_trace(void)
{
	sim->log_messages = 1;

	std::string message_path = out_dir + "/messages.csv";
	sim->message_trace_fp = fopen(message_path.c_str(), "w");
	if (sim->message_trace_fp == NULL)
	{
		INFO("Error: Couldn't open trace file for writing.\n");
		clean_up(RET_FAIL);
	}
	sim_message_trace_write_header(*sim);
}

void SanaFe::set_gui_flag(bool flag)
{
	sim->gui_on = flag;
	arch->spike_vector_on = flag;
}

void SanaFe::set_out_dir(std::string dir)
{
	if (dir.length() == 0)
	{
		out_dir = ".";
	}
	else
	{
		out_dir = dir;
	}
	return;
}

void SanaFe::set_arch(const char *filename)
{
	//FILE* arch_fp = fopen(filename, "r");
	std::fstream arch_fp;
	arch_fp.open(filename);
	if (arch_fp.fail())
	{
		INFO("Error: Architecture file %s failed to open.\n", filename);
		clean_up(RET_FAIL);
	}
	int ret = description_parse_arch_file(arch_fp, *arch);
	arch_fp.close();
	if (ret == RET_FAIL)
	{
		clean_up(RET_FAIL);
	}
}

void SanaFe::set_net(const char *filename)
{
	std::fstream network_fp;
	network_fp.open(filename);
	if (network_fp.fail())
	{
		INFO("Network data (%s) failed to open.\n", filename);
		clean_up(RET_FAIL);
	}
	INFO("Reading network from file.\n");
	int ret = description_parse_net_file(network_fp, net, *arch);
	network_fp.close();
	if (ret == RET_FAIL)
	{
		clean_up(RET_FAIL);
	}
	network_check_mapped(net);
	arch_create_axons(*arch);

	// Change Potential logging with new headers from net.
	if (sim->log_potential)
	{
		open_potential_trace();
	}
}

double SanaFe::get_power()
{
	if (sim->total_sim_time > 0.0)
	{
		return sim->total_energy / sim->total_sim_time;
	}
	else
	{
		return 0.0;
	}
}

vector<int> SanaFe::get_status(const vector<NeuronGroup>::size_type gid){
	vector<int> statuses = vector<int>();

	if (gid >= net.groups.size())
	{
		INFO("Error: Got gid of %lu with only %lu groups in net.\n",
			gid, net.groups.size());
		return statuses;
	}

	for (std::vector<Neuron>::size_type i = 0;
		i < net.groups[gid].neurons.size(); i++)
	{
		statuses.push_back(net.groups[gid].neurons[i].neuron_status);
	}

	return statuses;
}

void SanaFe::sim_summary()
{
	sim_write_summary(stdout, *sim);

	sim->stats_fp = fopen("run_summary.yaml", "w");
	if (sim->stats_fp != NULL)
	{
		sim_write_summary(sim->stats_fp, *sim);
	}
}

vector<vector<int>> SanaFe::run_summary()
{
	print_run_data(stdout, &run_data);

	// Could write intermediate run data here
	// sim->stats_fp = fopen("run_summary.yaml", "w");
	// if (sim->stats_fp != NULL)
	// {
	// 	print_run_data(sim->stats_fp, &run_data);
	// }

	// Return 2D vector of spiking tiles, auto clears
	// vector after return.
	Vector_Cleanup_Class help_class(arch);
	return arch->spike_vector;
}

void SanaFe::clean_up(int ret)
{
	// Free any larger structures here
	//arch_free(arch);

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
	//free(sim);

	if (ret == RET_FAIL)
	{
		exit(1);
	}
	else
	{
		exit(0);
	}
}

Simulation::Simulation()
{
	total_energy = 0.0;   // Joules
	total_sim_time = 0.0; // Seconds
	wall_time = 0.0;	   // Seconds
	timesteps = 0;
	total_spikes = 0;
	total_messages_sent = 0;
	total_neurons_fired = 0;

	// All logging disabled by default
	log_perf = 0;
	log_potential = 0;
	log_spikes = 0;
	log_messages = 0;

	potential_trace_fp = NULL;
	spike_trace_fp = NULL;
	perf_fp = NULL;
	message_trace_fp = NULL;
	stats_fp = NULL;
}

void run(Simulation &sim, Network &net, Architecture &arch)
{
	// TODO: remove the need to pass the network struct, only the arch
	//  should be needed (since it links back to the net anyway)
	// Run neuromorphic hardware simulation for one timestep
	//  Measure the CPU time it takes and accumulate the stats
	struct Timestep ts = Timestep(sim.timesteps, arch.get_core_count());
	struct timespec ts_start, ts_end, ts_elapsed;

	// Measure the wall-clock time taken to run the simulation
	//  on the host machine
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	sim_timestep(ts, net, arch);
	// Calculate elapsed time
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

	sim.total_energy += ts.energy;
	sim.total_sim_time += ts.sim_time;
	sim.total_spikes += ts.spike_count;
	sim.total_neurons_fired += ts.total_neurons_fired;
	sim.total_messages_sent += ts.packets_sent;
	if (sim.log_spikes)
	{
		sim_trace_record_spikes(sim, net);
	}
	if (sim.log_potential)
	{
		sim_trace_record_potentials(sim, net);
	}
	if (sim.log_perf)
	{
		sim_perf_log_timestep(ts, sim.perf_fp);
	}
	if (sim.log_messages)
	{
		for (auto &q: ts.messages)
		{
			for (auto it = q.begin(); it != q.end(); ++it)
			{
				Message &m = *it;
				if (!m.dummy_message)
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
	sim.timesteps = ts.timestep;
	sim.wall_time += (double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9);
	TRACE1("Time-step took: %fs.\n",
		(double) ts_elapsed.tv_sec+(ts_elapsed.tv_nsec/1.0e9));

}

void store_data_init(run_ts_data *data, Simulation &sim, int timesteps)
{
	data->energy = sim.total_energy;
	data->time = sim.total_sim_time;
	data->spikes = sim.total_spikes;
	data->packets = sim.total_messages_sent;
	data->neurons = sim.total_neurons_fired;
	data->wall_time = sim.wall_time;
	data->timestep_start = sim.timesteps;
	data->timesteps = timesteps;
}

void store_data(run_ts_data *data, Simulation &sim)
{
	data->energy = sim.total_energy - data->energy;
	data->time = sim.total_sim_time - data->time;
	data->spikes = sim.total_spikes - data->spikes;
	data->packets = sim.total_messages_sent - data->packets;
	data->neurons = sim.total_neurons_fired - data->neurons;
	data->wall_time = sim.wall_time - data->wall_time;
}

void print_run_data(FILE *fp, run_ts_data* data)
{
	fprintf(fp, "energy: %e\n", data->energy);
	fprintf(fp, "time: %e\n", data->time);
	fprintf(fp, "total_spikes: %ld\n", data->spikes);
	fprintf(fp, "total_packets: %ld\n", data->packets);
	fprintf(fp, "total_neurons_fired: %ld\n", data->neurons);
	fprintf(fp, "wall_time: %lf\n", data->wall_time);
	fprintf(fp, "executed from: %ld to %ld timesteps\n", data->timestep_start,
	data->timestep_start+data->timesteps);
}

void sim_init_fifo(MessageFifo &f)
{
	f.count = 0;
	f.head = NULL;
	f.tail = NULL;
	f.next = NULL;
}

void sim_timestep(Timestep &ts, Network &net, Architecture &arch)
{
	Scheduler s;

	// Start the next time-step
	ts = Timestep(ts.timestep, arch.get_core_count());
	sim_reset_measurements(net, arch);

	// TODO: reimplement user input spikes again
	//sim_input_spikes(net);
	sim_process_neurons(ts, net, arch);
	sim_receive_messages(ts, arch);

	s.noc_width = arch.noc_width;
	s.noc_height = arch.noc_height;
	s.buffer_size = arch.noc_buffer_size;

	ts.sim_time = sim_schedule_messages(ts.message_queues, s);
	// Performance statistics for this time step
	ts.energy = sim_calculate_energy(arch);

	// Setup spike vector
	std::vector<int> spike_tile_vec(ARCH_MAX_TILES);

	for (auto &t: arch.tiles)
	{
		int tile_spike_count = 0;
		for (auto &c: t.cores)
		{
			for (std::vector<SynapseUnit>::size_type k = 0;
				k < c.synapse.size(); k++)
			{
				ts.spike_count +=
					c.synapse[k].spikes_processed;
				tile_spike_count +=
					c.synapse[k].spikes_processed;
			}
			for (std::vector<SomaUnit>::size_type k = 0;
				k < c.soma.size(); k++)
			{
				ts.total_neurons_fired +=
					c.soma[k].neurons_fired;
			}
			for (std::vector<AxonOutUnit>::size_type k = 0;
				k < c.axon_out_hw.size(); k++)
			{
				ts.packets_sent += c.axon_out_hw[k].packets_out;
			}
		}
		if (arch.spike_vector_on)
		{
			spike_tile_vec[t.id] = tile_spike_count;
		}
	}
	if (arch.spike_vector_on)
	{
		arch.spike_vector.push_back(spike_tile_vec);
	}

	TRACE1("Spikes sent: %ld\n", sim.total_spikes);
	return;
}

std::unique_ptr<Simulation> sim_init_sim(void)
{
	std::unique_ptr<Simulation> sim(new Simulation);

	return sim;
}

Timestep::Timestep(const long int ts, const int core_count)
{
	timestep = ts;
	spike_count = 0L;
	message_queues = std::vector<MessageFifo>(core_count);
	messages = std::vector<std::list<Message> >(core_count);
	for (auto &q: message_queues)
	{
		sim_init_fifo(q);
	}
	total_neurons_fired = 0L;
	spikes = 0L;
	total_hops = 0L;
	energy = 0.0;
	sim_time = 0.0;
	packets_sent = 0L;
}

void sim_process_neurons(Timestep &ts, Network &net, Architecture &arch)
{
// #pragma omp parallel for
	for (auto &t: arch.tiles)
	{
		for (auto &c: t.cores)
		{
			for (std::vector<Neuron>::size_type k = 0;
				k < c.neurons.size(); k++)
			{
				Neuron &n = *(c.neurons[k]);
				sim_process_neuron(ts, arch, n);
			}

			if (c.neurons.size() > 0)
			{
				// Add a dummy message to account for neuron
				//  processing that does not result in any sent
				//  messages. To do this, set the dest neuron
				//  set as invalid with a 0 receiving latency)
				Message dummy_message = c.next_message;
				dummy_message.dummy_message = true;
				dummy_message.src_neuron = c.neurons.back();
				dummy_message.receive_delay = 0.0;
				dummy_message.network_delay = 0.0;
				dummy_message.hops = 0;
				dummy_message.timestep = ts.timestep;

				assert(c.id >= 0);
				assert(static_cast<size_t>(c.id) <
					ts.messages.size());
				TRACE1("c.id:%d message q size:%lu\n", c.id,
					ts.messages.size());
				ts.messages[c.id].push_back(dummy_message);
				sim_message_fifo_push(
						ts.message_queues[c.id],
						ts.messages[c.id].back());
			}
		}
	}
}

void sim_receive_messages(Timestep &ts, Architecture &arch)
{
	// Assign outgoing spike messages to their respective destination
	//  cores, and calculate network costs
	for (auto &q: ts.messages)
	{
		for (Message &m: q)
		{
			if (!m.dummy_message)
			{
				const size_t src_tile_id =
					m.src_neuron->core->parent_tile_id;
				assert(src_tile_id < arch.tiles.size());
				Tile &src_tile = arch.tiles[src_tile_id];
				assert(m.dest_tile_id >= 0);
				assert(static_cast<size_t>(m.dest_tile_id) <
					arch.tiles.size());
				Tile &dest_tile = arch.tiles[m.dest_tile_id];
				m.network_delay = sim_estimate_network_costs(
					src_tile, dest_tile);
				m.hops = abs(src_tile.x - dest_tile.x) +
						abs(src_tile.y - dest_tile.y);

				Core &c = dest_tile.cores[m.dest_core_offset];
				c.messages_in.push_back(&m);
			}
		}
	}

	// Now process all messages at receiving cores
//#pragma omp parallel for
	for (Tile &tile: arch.tiles)
	{
		for (Core &c: tile.cores)
		{
			TRACE1("Processing %lu message(s) for cid:%d.%d\n",
				c.messages_in.size(), tile.id, c.offset);
			for (auto m: c.messages_in)
			{
				m->receive_delay = sim_pipeline_receive(
					ts, arch, c, *m);
			}
		}
	}
}

double sim_estimate_network_costs(Tile &src, Tile &dest)
{
	double network_latency;
	long int x_hops, y_hops;

	network_latency = 0.0;

	// Calculate the energy and time for sending spike packets
	x_hops = abs(src.x - dest.x);
	y_hops = abs(src.y - dest.y);
	// E-W hops

	if (src.x < dest.x)
	{
		dest.east_hops += x_hops;
		network_latency += (double) x_hops * src.latency_east_hop;
	}
	else
	{
		dest.west_hops += x_hops;
		network_latency += (double) x_hops * src.latency_west_hop;
	}

	// N-S hops
	if (src.y < dest.y)
	{
		dest.north_hops += y_hops;
		network_latency += (double) y_hops * src.latency_north_hop;
	}
	else
	{
		dest.south_hops += y_hops;
		network_latency += (double) y_hops * src.latency_south_hop;
	}

	dest.hops += (x_hops + y_hops);
	dest.messages_received++;
	TRACE1("xhops:%ld yhops%ld total hops:%ld latency:%e\n", x_hops, y_hops,
		x_hops + y_hops, network_latency);
	return network_latency;
}

Message *sim_message_fifo_pop(struct MessageFifo *queue)
{
	Message *m;

	assert(queue->count >= 0);
	if (queue->count == 0)
	{
		m = NULL;
	}
	else
	{
		assert(queue->tail != NULL);
		assert(queue->head != NULL);
		queue->count--;
		m = queue->tail;
		queue->tail = queue->tail->next;
		if (queue->count <= 0)
		{
			queue->head = NULL;
		}
	}

	return m;
}

void sim_message_fifo_push(struct MessageFifo &queue, Message &m)
{
	assert(queue.count >= 0);

	m.next = NULL;
	if (queue.count == 0)
	{
		queue.tail = &m;
	}
	else
	{
		queue.head->next = &m;
	}
	queue.head = &m;
	queue.count++;
}

void sim_update_noc_message_counts(
	const Message &m, NocInfo &noc, const int message_in)
{
	assert(m.src_neuron != nullptr);

	// Go along x path, then y path (dimension order routing), and increment
	//  or decrement counter depending on the operation
	int x_increment, y_increment;
	// Adjust by dividing by the total number of links along the path, also
	//  including the output link at the sending core and input link at the
	//  receiving core, i.e. the hops plus 2. The total sum of the densities
	//  along the path should equal one for the message
	double adjust = (1.0 / (2.0 + m.hops));

	if (!message_in)
	{
		adjust *= -1.0;
	}

	if (m.src_x < m.dest_x)
	{
		x_increment = 1;
	}
	else
	{
		x_increment = -1;
	}
	if (m.src_y < m.dest_y)
	{
		y_increment = 1;
	}
	else
	{
		y_increment = -1;
	}
	int prev_direction =
		4 + (m.src_neuron->core->id % ARCH_MAX_CORES_PER_TILE);
	for (int x = m.src_x; x != m.dest_x; x += x_increment)
	{
		int direction;
		if (x_increment > 0)
		{
			direction = EAST;
		}
		else
		{
			direction = WEST;
		}
		if (x == m.src_x)
		{
			int link = 4 + (m.src_neuron->core->id %
				ARCH_MAX_CORES_PER_TILE);
			noc.noc_messages_in[x][m.src_y][link] += adjust;
		}
		else
		{
			noc.noc_messages_in[x][m.src_y][direction] += adjust;
		}
		prev_direction = direction;
	}
	for (int y = m.src_y; y != m.dest_y; y += y_increment)
	{
		int direction;
		if (y_increment > 0)
		{
			direction = NORTH;
		}
		else
		{
			direction = SOUTH;
		}
		if ((m.src_x == m.dest_x) && (y == m.src_y))
		{
			int link = 4 + (m.src_neuron->core->id %
				ARCH_MAX_CORES_PER_TILE);
			noc.noc_messages_in[m.dest_x][y][link] += adjust;
		}
		else
		{
			noc.noc_messages_in[m.dest_x][y][prev_direction] +=
				adjust;
		}

		prev_direction = direction;
	}

	if ((m.src_x == m.dest_x) && (m.src_y == m.dest_y))
	{
		int link = 4 +
			(m.src_neuron->core->id % ARCH_MAX_CORES_PER_TILE);
		noc.noc_messages_in[m.dest_x][m.dest_y][link] += adjust;
	}
	else
	{
		noc.noc_messages_in[m.dest_x][m.dest_y][prev_direction] += adjust;
	}

	// Update rolling averages and message counters
	if (message_in)
	{
		// Message entering NoC
		noc.mean_in_flight_receive_delay +=
			(m.receive_delay - noc.mean_in_flight_receive_delay) /
				(noc.messages_in_noc + 1);
		noc.messages_in_noc++;
	}
	else
	{
		// Message leaving the NoC
		if ((noc.messages_in_noc) > 1)
		{
			noc.mean_in_flight_receive_delay +=
			(noc.mean_in_flight_receive_delay - m.receive_delay) /
				(noc.messages_in_noc - 1);
		}
		else
		{
			noc.mean_in_flight_receive_delay = 0.0;
		}

		noc.messages_in_noc--;
	}

	return;
}

double sim_calculate_messages_along_route(Message &m, NocInfo &noc)
{
	double flow_density;
	int x_increment, y_increment;
	int direction;

	flow_density = 0.0;
	if (m.src_x < m.dest_x)
	{
		x_increment = 1;
	}
	else
	{
		x_increment = -1;
	}
	if (m.src_y < m.dest_y)
	{
		y_increment = 1;
	}
	else
	{
		y_increment = -1;
	}
	int prev_direction = 4 + (m.src_neuron->core->id % ARCH_MAX_CORES_PER_TILE);
	for (int x = m.src_x; x != m.dest_x; x += x_increment)
	{
		int direction = 0;
		if (x_increment > 0)
		{
			direction = EAST;
		}
		else
		{
			direction = WEST;
		}
		if (x == m.src_x)
		{
			flow_density += noc.noc_messages_in[x][m.src_y][
				4+m.src_neuron->core->id % ARCH_MAX_CORES_PER_TILE];
		}
		else
		{
			flow_density +=
				noc.noc_messages_in[x][m.src_y][direction];
		}
		prev_direction = direction;
	}

	for (int y = m.src_y; y != m.dest_y; y += y_increment)
	{
		if (y_increment > 0)
		{
			direction = NORTH;
		}
		else
		{
			direction = SOUTH;
		}
		if (m.src_x == m.dest_x && y == m.src_y)
		{
			flow_density +=
				noc.noc_messages_in[m.dest_x][y][
				4+m.src_neuron->core->id % ARCH_MAX_CORES_PER_TILE];
		}
		else
		{
			flow_density +=
				noc.noc_messages_in[m.dest_x][y][prev_direction];
		}
		prev_direction = direction;
	}
	// Handle the last tile
	if ((m.src_x == m.dest_x) && (m.src_y == m.dest_y))
	{
		flow_density +=
			noc.noc_messages_in[m.dest_x][m.dest_y][
				4+m.src_neuron->core->id % ARCH_MAX_CORES_PER_TILE];
	}
	else
	{
		flow_density +=
			noc.noc_messages_in[m.dest_x][m.dest_y][prev_direction];
	}

	assert(flow_density >= -0.1);
	return flow_density;
}

void sim_update_noc(const double t, NocInfo &noc)
{
	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		struct MessageFifo *q;
		Message *m;

		q = &(noc.messages_received[i]);
		m = q->tail;
		while (m != NULL)
		{
			if (m->in_noc && (t >= m->received_timestamp))
			{
				// Mark the message as not in the NoC, moving it
				//  from the network to the receiving core
				m->in_noc = 0;
				sim_message_fifo_pop(q);
				// Go along the message path and decrement tile
				//  message counters
				sim_update_noc_message_counts(*m, noc, 0);
			}
			m = m->next;
		}
	}

	return;
}

double sim_schedule_messages(std::vector<MessageFifo> &messages_sent,
				const Scheduler &scheduler)
{
	struct MessageFifo *priority_queue;
	NocInfo noc;
	double last_timestamp;

	noc.noc_width = scheduler.noc_width;
	noc.noc_height = scheduler.noc_height;

	noc.messages_in_noc = 0;
	noc.mean_in_flight_receive_delay = 0.0;
	for (int x = 0; x < noc.noc_width; x++)
	{
		for (int y = 0; y < noc.noc_height; y++)
		{
			for (int k = 0; k < (4 + ARCH_MAX_CORES_PER_TILE); k++)
			{
				noc.noc_messages_in[x][y][k] = 0.0;
			}
		}
	}
	for (int i = 0; i < ARCH_MAX_CORES; i++)
	{
		sim_init_fifo(noc.messages_received[i]);
		noc.core_finished_receiving[i] = 0.0;
	}

	priority_queue = sim_init_timing_priority(messages_sent);
	last_timestamp = 0.0;
	TRACE1("Scheduling global order of messages.\n");

	// Each core has a queue of received messages. A structure tracks how
	//  many in-flight messages are in the NoC and occupy each tile. We
	//  track the number of messages passing through each tile at
	//  the point of sending, and the average processing delay of
	//  all of those messages. When a message is added or removed from the
	//  NoC we update the average counts.
	while (priority_queue != NULL)
	{
		// Queue isn't empty
		struct MessageFifo *q;
		Message *next_message, *m;

		// Get the core's queue with the earliest simulation time
		q = sim_pop_priority_queue(&priority_queue);
		m = sim_message_fifo_pop(q);
		last_timestamp = fmax(last_timestamp, m->sent_timestamp);

		// Update the Network-on-Chip state
		sim_update_noc(m->sent_timestamp, noc);

		// Messages without a destination (neuron) are dummy messages.
		//  Dummy messages account for processing time that does not
		//  result in any spike messages. Otherwise, messages are sent
		//  from a src neuron to a dest neuron
		if (!m->dummy_message)
		{
			const int dest_core = m->dest_core_offset;
			// Figure out if we are able to send a message into the
			//  network i.e., is the route to the dest core
			//  saturated and likely to block? Sum along the route
			//  and see the density of messages along all links.
			double messages_along_route =
				sim_calculate_messages_along_route(*m, noc);

			int path_capacity = (m->hops+1) * scheduler.buffer_size;
			if (messages_along_route > path_capacity)
			{
				m->sent_timestamp +=
					(messages_along_route - path_capacity) *
					noc.mean_in_flight_receive_delay;
			}

			// Now, push the message into the right receiving queue
			//  Calculate the network delay and when the message
			//  is received
			m->in_noc = 1;
			sim_message_fifo_push(
				noc.messages_received[dest_core], *m);

			// Update the rolling average for message
			//  receiving times in-flight in the network
			sim_update_noc_message_counts(*m, noc, 1);

			double network_delay = messages_along_route *
				noc.mean_in_flight_receive_delay / (m->hops+1.0);
			TRACE1("Path capacity:%d messages:%lf delay:%e\n",
				path_capacity, messages_along_route, network_delay);

			double earliest_received_time = m->sent_timestamp + fmax(
				m->network_delay, network_delay);
			m->received_timestamp = fmax(
				noc.core_finished_receiving[dest_core],
				earliest_received_time);
			noc.core_finished_receiving[dest_core] = fmax(
				(noc.core_finished_receiving[dest_core] + m->receive_delay),
				(earliest_received_time + m->receive_delay));
			m->processed_timestamp =
				noc.core_finished_receiving[dest_core];
			last_timestamp = fmax(
				last_timestamp, m->processed_timestamp);
		}

		// Get the next message for this core
		next_message = q->tail;
		if (next_message != NULL)
		{
			// If applicable, schedule this next message immediately
			//  after the current message finishes sending
			next_message->sent_timestamp = m->sent_timestamp +
				next_message->generation_delay;
			last_timestamp = fmax(last_timestamp,
				next_message->sent_timestamp);
			sim_insert_priority_queue(&priority_queue, *q);
		}
		else
		{
			TRACE2("\t(cid:%d.%d) finished simulating\n", c->t->id,
				c->id);
		}

		if (priority_queue != NULL)
		{
			TRACE2("\t(cid:%d.%d) time:%e\n",
				(*priority_queue)->t->id,
				(*priority_queue)->id, (*priority_queue)->time);
		}
	}
	TRACE1("Neurons fired: %ld\n", ts.total_neurons_fired);

	return last_timestamp;
}

void sim_process_neuron(Timestep &ts, Architecture &arch, Neuron &n)
{
	if (!n.is_init)
	{
		return;
	}

	Core &c = *(n.core);
	n.processing_latency = 0.0;

	TRACE1("Processing neuron: %d.%d\n", n.id, n.parent_group_id);

	if (c.buffer_pos == BUFFER_SYNAPSE)
	{
		INFO("Error: Not implemented\n");
		exit(1);
	}
	else if (c.buffer_pos == BUFFER_DENDRITE)
	{
		// Go through all synapses connected to this neuron and update
		//  all synaptic currents into the dendrite
		INFO("Error: Not implemented\n");
		exit(1);
	}
	else if (c.buffer_pos == BUFFER_SOMA)
	{
		n.processing_latency = sim_update_soma(ts, arch, n, n.charge);
	}
	else if (c.buffer_pos == BUFFER_AXON_OUT)
	{
		if (n.fired)
		{
			n.processing_latency = n.soma_hw->latency_spiking;
			sim_neuron_send_spike_message(ts, arch, n);
		}
	}
	TRACE1("Updating neuron %d.%d.\n", n.parent_group_id, n.id);

	c.next_message.generation_delay += n.processing_latency;
	n.update_needed = 0;
	n.spike_count = 0;
}

double sim_pipeline_receive(Timestep &ts, Architecture &arch,
	Core &c, Message &m)
{
	// We receive a spike and process up to the time-step buffer
	double message_processing_latency = 0.0;

	TRACE1("Receiving messages for cid:%d\n", c.id);
	if (c.buffer_pos >= BUFFER_SYNAPSE)
	{
		assert(m.dest_axon_id >= 0);
		assert(static_cast<size_t>(m.dest_axon_id) < c.axons_in.size());
		AxonInModel &a = c.axons_in[m.dest_axon_id];
		for (int s: a.synapse_addresses)
		{
			const int synaptic_lookup = 1;
			message_processing_latency = sim_update_synapse(
				ts, arch, c, s, synaptic_lookup);
		}
	}

	return message_processing_latency;
}

struct MessageFifo *sim_init_timing_priority(
	std::vector<MessageFifo> &message_queues)
{
	struct MessageFifo *priority_queue;
	priority_queue = NULL;

	TRACE1("Initializing priority queue.\n");
	for (auto &q: message_queues)
	{
		if (q.count > 0) // messages
		{
			Message *m = q.tail;
			assert(m != NULL);
			m->sent_timestamp = m->generation_delay;
			sim_insert_priority_queue(&priority_queue, q);
		}
		else
		{
			TRACE1("No messages for core %d\n", i);
		}
	}

#ifdef DEBUG2
	int i = 0;
	for (Core *curr = priority_queue; curr != NULL; curr = curr->next)
	{
		// TODO
	}
#endif

	return priority_queue;
}

MessageFifo *sim_pop_priority_queue(MessageFifo **priority_queue)
{
	struct MessageFifo *curr;

	// Pop the first element from the priority queue
	curr = *priority_queue;
	*priority_queue = (*priority_queue)->next;

	// For safety, remove current element from queue and unlink
	curr->next = NULL;
	return curr;
}

void sim_insert_priority_queue(MessageFifo **priority_queue,
	MessageFifo &core_message_fifo)
{
	MessageFifo *next;

	// TODO: implement heap-based priority queue rather than list-based.
	//  Will achieve O(lg N) insertion time rather than O(N)

	assert(priority_queue != NULL);

	//INFO("Inserting into priority queue.\n");
	if (core_message_fifo.tail == NULL)
	{
		INFO("error?\n");
	}
	if (((*priority_queue) == NULL) ||
		(core_message_fifo.tail->sent_timestamp <=
		(*priority_queue)->tail->sent_timestamp))
	{
		// Queue is empty or this is the earliest time (highest
		//  priority), make this core the head of the queue
		core_message_fifo.next = (*priority_queue);
		*priority_queue = &core_message_fifo;
	}
	else
	{
		MessageFifo *curr = *priority_queue;
		next = curr->next;

		// Reinsert core into the correct place in the priority list
		while (next != NULL)
		{
			if (core_message_fifo.tail->sent_timestamp <
				next->tail->sent_timestamp)
			{
				break;
			}
			curr = next;
			next = curr->next;
		}
		curr->next = &core_message_fifo;
		core_message_fifo.next = next;
	}

#ifdef DEBUG
	TRACE3("*** Priority queue ***\n");
	for (Core *tmp = *priority_queue; tmp != NULL;
		tmp = tmp->next_timing)
	{
		// TRACE3
		TRACE3("tmp->time:%e (id:%d)\n", tmp->time, tmp->id);
		assert((tmp->next_timing == NULL) ||
			tmp->time <= tmp->next_timing->time);
	}
#endif

	return;
}

double sim_update_synapse(Timestep &ts, Architecture &arch, Core &c,
	const int synapse_address, const bool synaptic_lookup)
{
	// Update all synapses to different neurons in one core. If a synaptic
	//  lookup, read and accumulate the synaptic weights. Otherwise, just
	//  update filtered current and any other connection properties
	double latency, min_synaptic_resolution;
	latency = 0.0;
	Connection &con = *(c.synapses[synapse_address]);
	Neuron &post_neuron = *(con.post_neuron);

	TRACE1("Updating synapses for (cid:%d)\n", axon->pre_neuron->id);
	while (con.last_updated <= ts.timestep)
	{
		TRACE1("Updating synaptic current (last updated:%ld, ts:%ld)\n",
			con.last_updated, ts.timestep);
		con.current *= con.synaptic_current_decay;

		// "Turn off" synapses that have basically no
		//  synaptic current left to decay (based on
		//  the weight resolution)
		min_synaptic_resolution = (1.0 /
			con.synapse_hw->weight_bits);
		if (fabs(con.current) < min_synaptic_resolution)
		{
			con.current = 0.0;
		}

		TRACE2("(nid:%d->nid:%d) con->current:%lf\n",
			con.pre_neuron->id,
			con.post_neuron->id, con->current);
		con.last_updated++;

	}

	if (synaptic_lookup)
	{
		con.current += con.weight;
		post_neuron.update_needed = 1;
		post_neuron.spike_count++;

		assert(con.synapse_hw != NULL);
		con.synapse_hw->spikes_processed++;
		TRACE2("Sending spike to nid:%d, current:%lf\n",
			post_neuron->id, con.current);
		latency += con.synapse_hw->latency_spike_op;
	}
	if (c.buffer_pos != BUFFER_DENDRITE)
	{
		latency += sim_update_dendrite(
			ts, arch, post_neuron, con.current);
	}

	return latency;
}

double sim_update_dendrite(Timestep &ts, Architecture &arch, Neuron &n,
	const double charge)
{
	// TODO: Support dendritic operations, combining the current in
	//  different neurons in some way, and writing the result to an output
	double dendritic_current, latency;
	latency = 0.0;

	dendritic_current = 0.0;
	while (n.dendrite_last_updated <= ts.timestep)
	{
		TRACE3("Updating dendritic current (last_updated:%d, ts:%ld)\n",
			n.dendrite_last_updated, sim->timesteps);
		n.charge *= n.dendritic_current_decay;
		n.dendrite_last_updated++;
		dendritic_current = n.charge;
		TRACE2("nid:%d charge:%lf\n", n.id, n.charge);
	}

	// Update dendritic tap currents
	// TODO: implement multi-tap models
	TRACE2("Charge:%lf\n", charge);
	dendritic_current += charge;
	n.charge += charge;

	// Finally, send dendritic current to the soma
	TRACE2("nid:%d updating dendrite, charge:%lf\n", n.id, n.charge);
	if (n.core->buffer_pos != BUFFER_SOMA)
	{
		latency += sim_update_soma(ts, arch, n, dendritic_current);
	}

	return latency;
}

double sim_update_soma(Timestep &ts, Architecture &arch, Neuron &n,
	const double current_in)
{
	struct SomaUnit *soma = n.soma_hw;

	TRACE1("nid:%d updating, current_in:%lf\n", n.id, current_in);
	while (n.soma_last_updated <= ts.timestep)
	{
		n.neuron_status = n.model->update(current_in);
		n.soma_last_updated++;
	}

	double latency = 0.0;

	TRACE1("neuron status:%d\n", n.neuron_status);
	if (n.forced_spikes > 0)
	{
		n.neuron_status = sanafe::FIRED;
		n.forced_spikes--;
	}

	// Check for spiking
	if (n.neuron_status == sanafe::FIRED)
	{
		TRACE1("Neuron %d.%d fired\n", n.parent_group_id, n.id);
		ts.total_neurons_fired++;
		sim_neuron_send_spike_message(ts, arch, n);
	}

	// Update soma, if there are any received spikes, there is a non-zero
	//  bias or we force the neuron to update every time-step
	if ((n.neuron_status == sanafe::UPDATED) ||
		(n.neuron_status == sanafe::FIRED) ||
		(n.force_update))
	{
		latency += n.soma_hw->latency_update_neuron;
		soma->neuron_updates++;
	}

	latency += n.soma_hw->latency_access_neuron;

	return latency;
}

void sim_neuron_send_spike_message(Timestep &ts, Architecture &arch, Neuron &n)
{
	TRACE1("nid:%d.%d sending spike message to %lu axons out\n",
		n.parent_group_id, n.id, n.axon_out_addresses.size());
	for (int address: n.axon_out_addresses)
	{
		Core &src_core = *(n.core);
		const int core_id = src_core.id;
		const Tile &src_tile = arch.tiles[src_core.parent_tile_id];
		AxonOutModel &src_axon = src_core.axons_out[address];
		const int dest_address = src_axon.dest_axon_id;

		// Generate a spike message
		Tile &dest_tile = arch.tiles[src_axon.dest_tile_id];
		Core &dest_core = dest_tile.cores[src_axon.dest_core_offset];
		AxonInModel &dest_axon = dest_core.axons_in[dest_address];

		// TODO: figure some constructor for the message i.e., required
		//  fields?
		Message m;
		m.timestep = ts.timestep;
		m.src_neuron = &n;
		m.dest_axon_id = src_axon.dest_axon_id;
		m.spikes = dest_axon.synapse_addresses.size();
		m.dummy_message = false;
		m.dest_tile_id = dest_tile.id;
		m.dest_core_offset = dest_core.offset;
		m.src_x = src_tile.x;
		m.dest_x = dest_tile.x;
		m.src_y = src_tile.y;
		m.dest_y = dest_tile.y;

		// Add axon access cost to message latency and energy
		AxonOutUnit &axon_out_hw = *(n.axon_out_hw);
		m.generation_delay =
			src_core.next_message.generation_delay +
			axon_out_hw.latency_access;
		axon_out_hw.packets_out++;

		ts.messages[core_id].push_back(m);
		sim_message_fifo_push(ts.message_queues[core_id],
			ts.messages[core_id].back());

		// Reset the next message in this core
		src_core.next_message = Message();
	}

	return;
}

double sim_generate_noise(Neuron *n)
{
	assert(n != NULL);
	struct SomaUnit &soma_hw = *(n->soma_hw);
	int noise_val = 0;
	int ret;

	if (soma_hw.noise_type == NOISE_FILE_STREAM)
	{
		// With a noise stream, we have a file containing a series of
		//  random values. This is useful if we want to exactly
		//  replicate h/w without knowing how the stream is generated.
		//  We can record the random sequence and replicate it here
		char noise_str[MAX_NOISE_FILE_ENTRY];
		// If we get to the end of the stream, by default reset it.
		//  However, it is unlikely the stream will be correct at this
		//  point
		if (feof(soma_hw.noise_stream))
		{
			INFO("Warning: At the end of the noise stream. "
			     "Random values are unlikely to be correct.\n");
			fseek(soma_hw.noise_stream, 0, SEEK_SET);
		}
		fgets(noise_str, MAX_NOISE_FILE_ENTRY, soma_hw.noise_stream);
		ret = sscanf(noise_str, "%d", &noise_val);
		TRACE2("noise val:%d\n", noise_val);
		if (ret < 1)
		{
			INFO("Error: invalid noise stream entry.\n");
		}
	}

	// Get the number of noise bits required TODO: generalize
	int sign_bit = noise_val & 0x100;
	noise_val &= 0x7f; // TODO: hack, fixed for 8 bits
	if (sign_bit)
	{
		// Sign extend
		noise_val |= ~(0x7f);
	}

	return (double) noise_val;
}

double sim_calculate_energy(const Architecture &arch)
{
	// Returns the total energy across the design, for this timestep
	double network_energy, synapse_energy, soma_energy, axon_out_energy;
	double axon_in_energy, total_energy;

	network_energy = 0.0;
	axon_in_energy = 0.0;
	synapse_energy = 0.0;
	soma_energy = 0.0;
	axon_out_energy = 0.0;

	for (auto &t: arch.tiles)
	{
		double total_hop_energy = t.east_hops * t.energy_east_hop;

		total_hop_energy += t.west_hops * t.energy_west_hop;
		total_hop_energy += t.south_hops * t.energy_south_hop;
		total_hop_energy +=  t.north_hops * t.energy_north_hop;
		network_energy += total_hop_energy;

		for (auto &c: t.cores)
		{
			for (std::vector<AxonInUnit>::size_type k = 0;
				k < c.axon_in_hw.size(); k++)
			{
				axon_in_energy +=
					c.axon_in_hw[k].spike_messages_in *
					c.axon_in_hw[k].energy_spike_message;
			}
			for (std::vector<SynapseUnit>::size_type k = 0;
				k < c.synapse.size(); k++)
			{
				synapse_energy +=
					c.synapse[k].spikes_processed *
					c.synapse[k].energy_spike_op;
			}
			for (std::vector<SomaUnit>::size_type k = 0;
				k < c.soma.size(); k++)
			{

				soma_energy += c.soma[k].neuron_count *
					c.soma[k].energy_access_neuron;
				soma_energy += c.soma[k].neuron_updates *
					c.soma[k].energy_update_neuron;
				soma_energy += c.soma[k].neurons_fired *
					c.soma[k].energy_spiking;
			}
			for (std::vector<AxonOutUnit>::size_type k = 0;
				k < c.axon_out_hw.size(); k++)
			{
				axon_out_energy +=
					c.axon_out_hw[k].packets_out *
					c.axon_out_hw[k].energy_access;
			}
		}
	}

	total_energy = axon_in_energy + synapse_energy + soma_energy +
		axon_out_energy + network_energy;

	return total_energy;
}

void sim_reset_measurements(Network &net, Architecture &arch)
{
	for (auto &group: net.groups)
	{
		for (auto &n: group.neurons)
		{
			// Neurons can be manually forced to update, for example
			//  if they have a constant input bias
			n.update_needed |=
				(n.force_update || (n.neuron_status >= 1));
			n.processing_latency = 0.0;
			n.fired = 0;
		}
	}

	// Reset any energy, time latency or other measurements of network
	//  hardware
	for (auto &t: arch.tiles)
	{
		// Reset tile
		t.energy = 0.0;

		t.hops = 0;
		t.east_hops = 0;
		t.west_hops = 0;
		t.south_hops = 0;
		t.north_hops = 0;
		t.messages_received = 0;
		for (auto &c: t.cores)
		{
			// Reset core
			c.energy = 0.0;
			c.next_message = Message();

			for (std::vector<AxonInUnit>::size_type k = 0;
				k < c.axon_in_hw.size(); k++)
			{
				c.axon_in_hw[k].spike_messages_in = 0L;
				c.axon_in_hw[k].energy = 0.0;
				c.axon_in_hw[k].time = 0;
			}

			for (std::vector<DendriteUnit>::size_type k = 0;
				k < c.dendrite.size(); k++)
			{
				c.dendrite[k].energy = 0.0;
				c.dendrite[k].time = 0.0;
			}

			for (std::vector<SynapseUnit>::size_type k = 0;
				k < c.synapse.size(); k++)
			{
				c.synapse[k].energy = 0.0;
				c.synapse[k].time = 0.0;
				c.synapse[k].spikes_processed = 0;
			}

			for (std::vector<SomaUnit>::size_type k = 0;
				k < c.soma.size(); k++)
			{
				c.soma[k].energy = 0.0;
				c.soma[k].time = 0.0;
				c.soma[k].neuron_updates = 0L;
				c.soma[k].neurons_fired = 0L;
			}

			for (std::vector<AxonOutUnit>::size_type k = 0;
				k < c.axon_out_hw.size(); k++)
			{
				c.axon_out_hw[k].energy = 0.0;
				c.axon_out_hw[k].time = 0.0;
				c.axon_out_hw[k].packets_out = 0;
			}

			// Reset the message buffer
			c.messages_in = std::vector<Message *>();
		}
	}
}

void sim_perf_write_header(FILE *fp)
{
	fprintf(fp, "timestep,");
	fprintf(fp, "fired,");
	fprintf(fp, "packets,");
	fprintf(fp, "hops,");
	fprintf(fp, "sim_time,");
	fprintf(fp, "total_energy,");
	fprintf(fp, "\n");
}

void sim_perf_log_timestep(const Timestep &ts, FILE *fp)
{
	fprintf(fp, "%ld,", ts.timestep);
	fprintf(fp, "%ld,", ts.total_neurons_fired);
	fprintf(fp, "%ld,", ts.packets_sent);
	fprintf(fp, "%ld,", ts.total_hops);
	fprintf(fp, "%le,", ts.sim_time);
	fprintf(fp, "%le,", ts.energy);
	fprintf(fp, "\n");
}

void sim_write_summary(FILE *fp, const Simulation &sim)
{
	// Write the simulation summary to file
	fprintf(fp, "git_version: %s\n", GIT_COMMIT);
	fprintf(fp, "energy: %e\n", sim.total_energy);
	fprintf(fp, "sim_time: %e\n", sim.total_sim_time);
	fprintf(fp, "total_spikes: %ld\n", sim.total_spikes);
	fprintf(fp, "total_packets: %ld\n", sim.total_messages_sent);
	fprintf(fp, "total_neurons_fired: %ld\n", sim.total_neurons_fired);
	fprintf(fp, "wall_time: %lf\n", sim.wall_time);
	fprintf(fp, "timesteps: %ld\n", sim.timesteps);
}

void sim_spike_trace_write_header(const Simulation &sim)
{
	assert(sim.spike_trace_fp != NULL);
	fprintf(sim.spike_trace_fp, "neuron,timestep\n");

	return;
}

void sim_potential_trace_write_header(const Simulation &sim, const Network &net)
{
	// Write csv header for probe outputs - record which neurons have been
	//  probed
	fprintf(sim.potential_trace_fp, "timestep,");
	for (auto &group: net.groups)
	{
		for (auto n: group.neurons)
		{
			if (sim.potential_trace_fp && sim.log_potential &&
				n.log_potential)
			{
				fprintf(sim.potential_trace_fp, "neuron %d.%d,",
					group.id, n.id);
			}
		}
	}

	if (sim.potential_trace_fp)
	{
		fputc('\n', sim.potential_trace_fp);
	}

	return;
}

void sim_message_trace_write_header(const Simulation &sim)
{
	assert(sim.message_trace_fp);
	fprintf(sim.message_trace_fp, "timestep,src_neuron,");
	fprintf(sim.message_trace_fp, "src_hw,dest_hw,hops,spikes,");
	fprintf(sim.message_trace_fp, "generation_latency,network_latency,");
	fprintf(sim.message_trace_fp, "processing_latency,blocking_latency\n");
}

void sim_trace_record_spikes(
	const Simulation &sim, const Network &net)
{
	// A trace of all spikes that are generated
	assert(sim.spike_trace_fp != NULL);

	for (auto &group: net.groups)
	{
		for (auto &n: group.neurons)
		{
			if (n.log_spikes && n.fired)
			{
				fprintf(sim.spike_trace_fp, "%d.%d,%ld\n",
					n.parent_group_id, n.id,
					sim.timesteps+1);
			}
		}
	}

	return;
}

// TODO: should potential be a required value in the soma?
void sim_trace_record_potentials(
	const Simulation &sim, const Network &net)
{
	// Each line of this csv file is the potential of all probed neurons for
	//  one time-step
	int potential_probe_count = 0;

	fprintf(sim.potential_trace_fp, "%ld,", sim.timesteps+1);
	for (auto &group: net.groups)
	{
		for (auto &n: group.neurons)
		{
			if (sim.potential_trace_fp && n.log_potential)
			{
				//fprintf(sim.potential_trace_fp, "%lf,",
				//	n.model->get_potential());
				potential_probe_count++;
			}
		}
	}

	// Each timestep takes up a line in the respective csv file
	if (sim.potential_trace_fp && (potential_probe_count > 0))
	{
		fputc('\n', sim.potential_trace_fp);
	}

	return;
}

void sim_trace_record_message(const Simulation &sim, const Message &m)
{
	fprintf(sim.message_trace_fp, "%ld,", m.timestep);
	assert(m.src_neuron != nullptr);
	fprintf(sim.message_trace_fp, "%d.%d,", m.src_neuron->parent_group_id,
		m.src_neuron->id);
	assert(m.src_neuron->core != nullptr);
	fprintf(sim.message_trace_fp, "%d.%d,",
		m.src_neuron->core->parent_tile_id,
		m.src_neuron->core->id);

	fprintf(sim.message_trace_fp, "%d.%d,", m.dest_tile_id,
		m.dest_core_offset);
	fprintf(sim.message_trace_fp, "%d,", m.hops);
	fprintf(sim.message_trace_fp, "%d,", m.spikes);
	fprintf(sim.message_trace_fp, "%le,", m.generation_delay);
	fprintf(sim.message_trace_fp, "%le,", m.network_delay);
	fprintf(sim.message_trace_fp, "%le,", m.receive_delay);
	fprintf(sim.message_trace_fp, "%le,", m.blocked_latency);
	fprintf(sim.message_trace_fp, "%le,", m.sent_timestamp);
	fprintf(sim.message_trace_fp, "%le\n", m.processed_timestamp);

	return;
}

int sim_poisson_input(const double firing_probability)
{
	// Simulate a single external input (as one neuron) for a timestep
	//  Return 1 if the input fires, 0 otherwise
	double rand_uniform;
	int input_fired;

	rand_uniform = (double) rand() / RAND_MAX;
	input_fired = rand_uniform < firing_probability;

	return input_fired;
}

int sim_rate_input(const double firing_rate, double *current)
{
	int input_fired;

	// Note: rate-based input (without randomization) is equivalent to a
	//  neuron with a fixed bias.
	TRACE2("rate input:%lf\n", firing_rate);
	*current += firing_rate;
	if (*current > 255.0)
	{
		*current = 0;
		input_fired = 1;
	}
	else
	{
		input_fired = 0;
	}

	TRACE2("input fired value:%d\n", input_fired);
	return input_fired;
}

/*
PYBIND11_MODULE(sanafe, m)
{
	m.doc() = R"pbdoc(
        SANA-FE CPP Module with Pybind11
        --------------------------------

        .. currentmodule:: simcpp

        .. autosummary::
           :toctree: _generate

		   SANA_FE
	)pbdoc";
	pybind11::class_<SanaFe>(m, "SanaFe")
		.def(pybind11::init())
		.def("init", &SanaFe::init)
        .def("update_neuron", &SanaFe::update_neuron)
        .def("run_timesteps", &SanaFe::run_timesteps, pybind11::arg("timesteps")=1)
		.def("set_perf_trace", &SanaFe::set_perf_trace, pybind11::arg("enable")=true)
		.def("set_spike_trace", &SanaFe::set_spike_trace, pybind11::arg("enable")=true)
		.def("set_potential_trace", &SanaFe::set_potential_trace, pybind11::arg("enable")=true)
		.def("set_message_trace", &SanaFe::set_message_trace, pybind11::arg("enable")=true)
		.def("set_out_dir", &SanaFe::set_out_dir)
		.def("set_gui_flag", &SanaFe::set_gui_flag, pybind11::arg("flag")=true)
		.def("set_arch", &SanaFe::set_arch)
		.def("set_net", &SanaFe::set_net)
		.def("get_power", &SanaFe::get_power)
		.def("get_status", &SanaFe::get_status)
		.def("sim_summary", &SanaFe::sim_summary)
		.def("run_summary", &SanaFe::run_summary)
		.def("clean_up", &SanaFe::clean_up, pybind11::arg("ret") = 0);
}
*/

/*
double sim_calculate_time(const Architecture *const arch)
{
	Core *c = n->core;
	TRACE1("nid:%d sending spike(s).\n", n->id);
	int core_id = n->core->id;


	for (int k = 0; k < n->maps_out_count; k++)
	{
		Axon *dest_axon;
		Message *m;
		int message_index;

		dest_axon = n->maps_out[k];
		message_index = ts->message_queues[core_id].count;

		// Generate a spike message
		m = &(ts->messages[core_id][message_index]);
		arch_init_message(m);
		m->timestep = ts->timestep;
		m->src_neuron = n;
		m->spikes = dest_axon->connection_count;
		m->dest_neuron = dest_axon->connections[0]->post_neuron;
		// Add axon access cost to message latency and energy
		m->generation_latency =
			c->next_message.generation_delay +
			c->axon_out.latency_access;
		sim_message_fifo_push(&(ts->message_queues[core_id]), m);

		c->axon_out.packets_out++;

		// Record a spike message at all the connected cores (axons)
		dest_axon->spikes_received++;
		dest_axon->message = m;

		// Reset the next message in this core
		arch_init_message(&(c->next_message));
	}

	return 0.0;
}
*/