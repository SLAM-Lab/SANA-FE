// Copyright (c) 2024 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
//  sim.c
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <vector>
#include <list>
#include <memory>
#include <sstream>
#include <filesystem>

#include "print.hpp"
#include "sim.hpp"
#include "network.hpp"
#include "arch.hpp"

using namespace sanafe;

Simulation::Simulation(
	Architecture &a,
	Network &n,
	const std::string &output_dir=".",
	const bool record_spikes=false, const bool record_potentials=false,
	const bool record_perf=false, const bool record_messages=false):
	arch(a), net(n), out_dir(output_dir)
{
	INFO("Initializing simulation.\n");
	total_energy = 0.0;   // Joules
	total_sim_time = 0.0; // Seconds
	wall_time = 0.0;      // Seconds
	total_timesteps = 0;
	total_spikes = 0;
	total_messages_sent = 0;
	total_neurons_fired = 0;

	// All logging disabled by default
	spike_trace_enabled = record_spikes;
	potential_trace_enabled = record_potentials;
	perf_trace_enabled = record_perf;
	message_trace_enabled = record_messages;
}

Simulation::~Simulation()
{
	// Close any open trace files
	spike_trace.close();
	potential_trace.close();
	perf_trace.close();
	message_trace.close();
}

/*
int Simulation::update_neuron(
	const std::vector<NeuronGroup>::size_type group_id,
	const std::vector<Neuron>::size_type n_id,
	const std::vector<std::string> kwargs, const int count)
{
	for (std::string item: kwargs)
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

	std::vector<Attribute> attr;
	for (auto &s: kwargs)
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
*/

RunData::RunData(const long int start, const long int steps)
{
	timestep_start = start;
	timesteps_executed = steps;
	energy = 0.0;
	sim_time = 0.0;
	spikes = 0L;
	packets_sent = 0L;
	neurons_fired = 0L;
	return;
}

RunData Simulation::run(const long int timesteps, const long int heartbeat)
{
	RunData rd((total_timesteps+1), timesteps);
	if (total_timesteps <= 0)
	{
		// If no timesteps have been simulated, open the trace files
		//  and simulate.
		// TODO: consider moving this back and initializing the
		//  simulation with a reference to an arch and network
		if (spike_trace_enabled)
		{
			spike_trace = sim_trace_open_spike_trace(out_dir);
		}
		if (potential_trace_enabled)
		{
			potential_trace = sim_trace_open_potential_trace(
				out_dir, net);
		}
		if (perf_trace_enabled)
		{
			perf_trace = sim_trace_open_perf_trace(out_dir);
		}
		if (message_trace_enabled)
		{
			message_trace = sim_trace_open_message_trace(out_dir);
		}
	}

	for (long int timestep = 1; timestep <= timesteps; timestep++)
	{
		if ((timestep % heartbeat) == 0)
		{
			// Print heart-beat every hundred timesteps
			INFO("*** Time-step %ld ***\n", timestep);
		}
		const Timestep ts = step();
		rd.energy += ts.energy;
		rd.sim_time += ts.sim_time;
		rd.spikes += ts.spikes;
		rd.packets_sent += ts.packets_sent;
		rd.neurons_fired += ts.neurons_fired;
		rd.wall_time = wall_time;
	}

	return rd;
}

Timestep Simulation::step()
{
	// TODO: remove the need to pass the network struct, only the arch
	//  should be needed (since it links back to the net anyway)
	// Run neuromorphic hardware simulation for one timestep
	//  Measure the CPU time it takes and accumulate the stats
	total_timesteps++;
	struct Timestep ts = Timestep(total_timesteps, arch.get_core_count());
	struct timespec ts_start, ts_end, ts_elapsed;

	// Measure the wall-clock time taken to run the simulation
	//  on the host machine
	clock_gettime(CLOCK_MONOTONIC, &ts_start);

	sim_timestep(ts, arch, net);

	// Calculate elapsed time
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	ts_elapsed = calculate_elapsed_time(ts_start, ts_end);

	total_energy += ts.energy;
	total_sim_time += ts.sim_time;
	total_spikes += ts.spike_count;
	total_neurons_fired += ts.neurons_fired;
	total_messages_sent += ts.packets_sent;
	if (spike_trace_enabled)
	{
		sim_trace_record_spikes(spike_trace, total_timesteps, net);
	}
	if (potential_trace_enabled)
	{
		sim_trace_record_potentials(
			potential_trace, total_timesteps, net);
	}
	if (perf_trace_enabled)
	{
		sim_trace_perf_log_timestep(perf_trace, ts);
	}
	if (message_trace_enabled)
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
					sim_trace_record_message(
						message_trace, m);
				}
			}
		}
	}
	wall_time += (double) ts_elapsed.tv_sec + (ts_elapsed.tv_nsec / 1.0e9);
	TRACE1("Time-step took: %fs.\n", static_cast<double>(
			ts_elapsed.tv_sec + (ts_elapsed.tv_nsec / 1.0e9)));

	return ts;
}

double Simulation::get_power()
{
	double power; // Watts
	if (total_sim_time > 0.0)
	{
		power = total_energy / total_sim_time;
	}
	else
	{
		// Avoid divide by 0
		power = 0.0;
	}

	return power;
}

RunData Simulation::get_run_summary()
{
	// Store the summary data in a string to string mapping
	RunData run_data(0, total_timesteps);

	run_data.energy = total_energy;
	run_data.sim_time = total_sim_time;
	run_data.spikes = total_spikes;
	run_data.packets_sent = total_messages_sent;
	run_data.wall_time = wall_time;
	run_data.neurons_fired = total_neurons_fired;

	return run_data;
}

void sanafe::sim_output_run_summary(
	const std::filesystem::path &out_dir, const RunData &run_data)
{
	// Summarize and output the run data using a YAML format to the console
	sim_format_run_summary(std::cout, run_data);

	// Output the same YAML-formatted summary to the given output file
	const std::filesystem::path summary_filename("run_summary.yaml");
	const std::filesystem::path summary_path = out_dir / summary_filename;
	std::ofstream summary_file(summary_path);
	if (summary_file.is_open())
	{
		sim_format_run_summary(summary_file, run_data);
	}
	else
	{
		INFO("Summary file %s couldn't open.\n", summary_path.c_str());
	}
}

void sanafe::sim_format_run_summary(std::ostream &out,
	const RunData &run_data)
{
	out << "build_git_version: '" << GIT_COMMIT << "'" << std::endl;
	out << "energy: " << std::scientific << run_data.energy << std::endl;
	out << "sim_time: " << std::scientific << run_data.sim_time;
	out << std::endl;
	out << "total_spikes: " << run_data.spikes << std::endl;
	out << "total_messages_sent: " << run_data.packets_sent << std::endl;
	out << "wall_time: " << std::fixed << run_data.wall_time << std::endl;
	out << "total_neurons_fired: " << run_data.neurons_fired << std::endl;

	return;
}

std::ofstream sanafe::sim_trace_open_spike_trace(
	const std::filesystem::path &out_dir)
{
	// TODO: warning, this is specific to Linux
	// To be more portable, consider using the filesystem library in C++17
	const std::filesystem::path spike_path = out_dir / "spikes.csv";
	std::ofstream spike_file(spike_path);

	if (!spike_file.is_open())
	{
		throw std::runtime_error(
			"Error: Couldn't open trace file for writing.");
	}
	sim_trace_write_spike_header(spike_file);
	return spike_file;
}

std::ofstream sanafe::sim_trace_open_potential_trace(
	const std::filesystem::path &out_dir, const Network &net)
{
	// TODO: warning, this is specific to Linux
	// To be more portable, consider using the filesystem library in C++17
	const std::filesystem::path potential_path = out_dir / "potential.csv";
	std::ofstream potential_file(potential_path);

	if (!potential_file.is_open())
	{
		throw std::runtime_error(
			"Error: Couldn't open trace file for writing.");
	}
	sim_trace_write_potential_header(potential_file, net);
	return potential_file;
}

std::ofstream sanafe::sim_trace_open_perf_trace(
	const std::filesystem::path &out_dir)
{
	// TODO: warning, this is specific to Linux
	// To be more portable, consider using the filesystem library in C++17
	const std::filesystem::path perf_path = out_dir / "perf.csv";
	std::ofstream perf_file(perf_path);
	if (!perf_file.is_open())
	{
		throw std::runtime_error(
			"Error: Couldn't open trace file for writing.");
	}
	sim_trace_write_perf_header(perf_file);

	return perf_file;
}

std::ofstream sanafe::sim_trace_open_message_trace(
	const std::filesystem::path &out_dir)
{
	const std::filesystem::path message_path = out_dir / "/messages.csv";
	std::ofstream message_file(message_path);
	if (!message_file.is_open())
	{
		throw std::runtime_error(
			"Error: Couldn't open trace file for writing.");
	}
	sim_trace_write_message_header(message_file);
	return message_file;
}

void sanafe::sim_init_fifo(MessageFifo &f)
{
	f.count = 0;
	f.head = NULL;
	f.tail = NULL;
	f.next = NULL;
}

void sanafe::sim_timestep(Timestep &ts, Architecture &arch, Network &net)
{
	Scheduler s;

	// Start the next time-step
	ts = Timestep(ts.timestep, arch.get_core_count());
	sim_reset_measurements(net, arch);

	sim_process_neurons(ts, net, arch);
	sim_receive_messages(ts, arch);

	s.noc_width = arch.noc_width;
	s.noc_height = arch.noc_height;
	s.buffer_size = arch.noc_buffer_size;

	ts.sim_time = sim_schedule_messages(ts.message_queues, s);
	// Performance statistics for this time step
	ts.energy = sim_calculate_energy(arch);

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
				ts.neurons_fired += c.soma[k].neurons_fired;
			}
			for (std::vector<AxonOutUnit>::size_type k = 0;
				k < c.axon_out_hw.size(); k++)
			{
				ts.packets_sent += c.axon_out_hw[k].packets_out;
			}
		}
	}

	TRACE1("Spikes sent: %ld\n", ts.spike_count);
	return;
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
	neurons_fired = 0L;
	spikes = 0L;
	total_hops = 0L;
	energy = 0.0;
	sim_time = 0.0;
	packets_sent = 0L;
}

void sanafe::sim_process_neurons(Timestep &ts, Network &net, Architecture &arch)
{
#pragma omp parallel for schedule(dynamic)
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

void sanafe::sim_receive_messages(Timestep &ts, Architecture &arch)
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
#pragma omp parallel for schedule(dynamic)
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

double sanafe::sim_estimate_network_costs(Tile &src, Tile &dest)
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

Message *sanafe::sim_message_fifo_pop(MessageFifo *queue)
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

void sanafe::sim_message_fifo_push(MessageFifo &queue, Message &m)
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

void sanafe::sim_update_noc_message_counts(
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

double sanafe::sim_calculate_messages_along_route(Message &m, NocInfo &noc)
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

void sanafe::sim_update_noc(const double t, NocInfo &noc)
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

double sanafe::sim_schedule_messages(std::vector<MessageFifo> &messages_sent,
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

void sanafe::sim_process_neuron(Timestep &ts, Architecture &arch, Neuron &n)
{
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
		n.charge = 0.0;
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

double sanafe::sim_pipeline_receive(Timestep &ts, Architecture &arch,
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

struct MessageFifo *sanafe::sim_init_timing_priority(
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

MessageFifo *sanafe::sim_pop_priority_queue(MessageFifo **priority_queue)
{
	struct MessageFifo *curr;

	// Pop the first element from the priority queue
	curr = *priority_queue;
	*priority_queue = (*priority_queue)->next;

	// For safety, remove current element from queue and unlink
	curr->next = NULL;
	return curr;
}

void sanafe::sim_insert_priority_queue(MessageFifo **priority_queue,
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

double sanafe::sim_update_synapse(Timestep &ts, Architecture &arch, Core &c,
	const int synapse_address, const bool synaptic_lookup)
{
	// Update all synapses to different neurons in one core. If a synaptic
	//  lookup, read and accumulate the synaptic weights. Otherwise, just
	//  update filtered current and any other connection properties
	double latency, min_synaptic_resolution;
	latency = 0.0;
	Connection &con = *(c.synapses[synapse_address]);
	Neuron &post_neuron = *(con.post_neuron);

	TRACE1("Updating synapses for (cid:%d)\n", c.id);
	while (con.last_updated < ts.timestep)
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

	TRACE1("(nid:%d.%d->nid:%d.%d) con->current:%lf\n",
		con.pre_neuron->parent_group_id, con.pre_neuron->id,
		con.post_neuron->parent_group_id, con.post_neuron->id,
		con.current);

	if (c.buffer_pos != BUFFER_DENDRITE)
	{
		latency += sim_update_dendrite(
			ts, arch, post_neuron, con.current);
	}

	return latency;
}

double sanafe::sim_update_dendrite(Timestep &ts, Architecture &arch, Neuron &n,
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

double sanafe::sim_update_soma(Timestep &ts, Architecture &arch, Neuron &n,
	const double current_in)
{
	SomaUnit *const soma = n.soma_hw;

	TRACE1("nid:%d updating, current_in:%lf\n", n.id, current_in);
	while (n.soma_last_updated < ts.timestep)
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
		soma->neurons_fired++;
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

void sanafe::sim_neuron_send_spike_message(
	Timestep &ts, Architecture &arch, Neuron &n)
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

double sanafe::sim_generate_noise(Neuron *n)
{
	assert(n != NULL);
	struct SomaUnit &soma_hw = *(n->soma_hw);
	int noise_val = 0;

	if (soma_hw.noise_type == NOISE_FILE_STREAM)
	{
		// With a noise stream, we have a file containing a series of
		//  random values. This is useful if we want to exactly
		//  replicate h/w without knowing how the stream is generated.
		//  We can record the random sequence and replicate it here
		char noise_str[MAX_NOISE_FILE_ENTRY];
		char *str = &(noise_str[0]);
		// If we get to the end of the stream, by default reset it.
		//  However, it is unlikely the stream will be correct at this
		//  point
		if (feof(soma_hw.noise_stream))
		{
			INFO("Warning: At the end of the noise stream. "
			     "Random values are unlikely to be correct.\n");
			fseek(soma_hw.noise_stream, 0, SEEK_SET);
		}
		char *result = fgets(
			str, MAX_NOISE_FILE_ENTRY, soma_hw.noise_stream);
		if (result != NULL)
		{
			const int ret = sscanf(noise_str, "%d", &noise_val);
			TRACE2("noise val:%d\n", noise_val);

			if (ret < 1)
			{
				INFO("Error: invalid noise stream entry.\n");
			}
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

double sanafe::sim_calculate_energy(const Architecture &arch)
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

void sanafe::sim_reset_measurements(Network &net, Architecture &arch)
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

void sanafe::sim_trace_write_spike_header(std::ofstream &spike_trace_file)
{
	assert(spike_trace_file.is_open());
	spike_trace_file << "neuron,timestep" << std::endl;

	return;
}

void sanafe::sim_trace_write_potential_header(
	std::ofstream &potential_trace_file, const Network &net)
{
	// Write csv header for probe outputs - record which neurons have been
	//  probed
	assert(potential_trace_file.is_open());
	potential_trace_file << "timestep,";
	for (auto &group: net.groups)
	{
		for (auto &n: group.neurons)
		{
			if (n.log_potential)
			{
				potential_trace_file << "neuron " << group.id;
				potential_trace_file << "." << n.id << ",";
			}
		}
	}
	potential_trace_file << std::endl;

	return;
}

void sanafe::sim_trace_write_perf_header(std::ofstream &perf_trace_file)
{
	assert(perf_trace_file.is_open());
	perf_trace_file << "timestep,";
	perf_trace_file << "fired,";
	perf_trace_file << "packets,";
	perf_trace_file << "hops,";
	perf_trace_file << "sim_time,";
	perf_trace_file << "total_energy,";
	perf_trace_file << std::endl;

	return;
}

void sanafe::sim_trace_write_message_header(std::ofstream &message_trace_file)
{
	assert(message_trace_file.is_open());
	message_trace_file << "timestep,";
	message_trace_file << "src_neuron,";
	message_trace_file << "src_hw,";
	message_trace_file << "dest_hw,";
	message_trace_file << "hops,";
	message_trace_file << "spikes,";
	message_trace_file << "generation_latency,";
	message_trace_file << "network_latency,";
	message_trace_file << "processing_latency,";
	message_trace_file << "blocking_latency";
	message_trace_file << std::endl;

	return;
}

void sanafe::sim_trace_record_spikes(
	std::ofstream &out, const long int timestep, const Network &net)
{
	// A trace of all spikes that are generated
	assert(out.is_open());

	for (auto &group: net.groups)
	{
		for (auto &n: group.neurons)
		{
			if (n.log_spikes && (n.neuron_status == sanafe::FIRED))
			{
				out << n.parent_group_id << "." << n.id << ",";
				out << timestep;
				out << std::endl;
			}
		}
	}

	return;
}

// TODO: should potential be a required value in the soma?
void sanafe::sim_trace_record_potentials(
	std::ofstream &out, const int timestep, const Network &net)
{
	// Each line of this csv file is the potential of all probed neurons for
	//  one time-step
	assert(out.is_open());
	TRACE1("Recording potential for timestep: %d\n", timestep);
	out << timestep << ",";

	long int potential_probe_count = 0;
	for (auto &group: net.groups)
	{
		for (auto &n: group.neurons)
		{
			if (n.log_potential)
			{
				out << n.model->get_potential() << ",";
				potential_probe_count++;
			}
		}
	}

	// Each timestep takes up a line in the respective csv file
	if (potential_probe_count > 0)
	{
		out << std::endl;
	}

	return;
}

void sanafe::sim_trace_perf_log_timestep(std::ofstream &out, const Timestep &ts)
{
	out << ts.timestep << ",";
	out << ts.neurons_fired << ",";
	out << ts.packets_sent << ",";
	out << ts.total_hops << ",";
	out << std::scientific << ts.sim_time << ",";
	out << std::scientific << ts.energy << ",";
	out << std::endl;
}

void sanafe::sim_trace_record_message(std::ofstream &out, const Message &m)
{
	assert(out.is_open());

	out << m.timestep << ",";
	assert(m.src_neuron != nullptr);
	out << m.src_neuron->parent_group_id << ".";
	out << m.src_neuron->id << ",";

	assert(m.src_neuron->core != nullptr);
	out << m.src_neuron->core->parent_tile_id << ".";
	out << m.src_neuron->core->offset << ",";
	out << m.dest_tile_id << "." << m.dest_core_offset << ",";

	out << m.hops << ",";
	out << m.spikes << ",";
	out << m.generation_delay << ",";
	out << m.network_delay << ",";
	out << m.receive_delay << ",";
	out << m.blocked_latency << ",";
	out << m.hops << ",";
	out << m.hops << ",";

	out << std::endl;

	return;
}

timespec sanafe::calculate_elapsed_time(const timespec &ts_start,
	const timespec &ts_end)
{
	// Calculate elapsed wall-clock time between ts_start and ts_end
	timespec ts_elapsed;

	ts_elapsed.tv_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
	ts_elapsed.tv_sec = ts_end.tv_sec - ts_start.tv_sec;
	if (ts_end.tv_nsec < ts_start.tv_nsec)
	{
		ts_elapsed.tv_sec--;
		ts_elapsed.tv_nsec += 1000000000UL;
	}

	return ts_elapsed;
}

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


/*
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
*/