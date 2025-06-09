// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// chip.hpp - Neuromorphic simulator kernel
//
// Time-step based simulation, based on loop:
// 1) seed any input spikes
// 2) route spikes
// 3) update neurons and check firing
#ifndef CHIP_HEADER_INCLUDED_
#define CHIP_HEADER_INCLUDED_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

#include "arch.hpp"
#include "core.hpp"
#include "fwd.hpp"
#include "pipeline.hpp"
#include "tile.hpp"

namespace sanafe
{

constexpr long int default_heartbeat_timesteps = 100L;

enum TimingModel : uint8_t
{
    timing_model_simple, // analytical model
    timing_model_detailed, // semi-analytical model
    timing_model_cycle_accurate, // Booksim2 simulator
};

class SpikingChip
{
public:
    std::vector<Tile> tiles;
    // Keep a reference to the different neuron groups mapped to the H/W
    std::map<std::string, std::vector<std::reference_wrapper<MappedNeuron>>> mapped_neuron_groups;

    SpikingChip(const Architecture &arch);
    ~SpikingChip();
    // Do not allow copying
    SpikingChip(const SpikingChip &copy) = delete;
    SpikingChip(SpikingChip &&other) = delete;
    SpikingChip &operator=(const SpikingChip &copy) = delete;
    SpikingChip &operator=(SpikingChip &&other) = delete;
    RunData sim(long int timesteps = 1, TimingModel timing_model = timing_model_detailed, int scheduler_thread_count = 1, bool record_spikes = false, bool record_potentials = false, bool record_perf = false, bool record_messages = false, std::string output_dir = "");
    void step(Scheduler &scheduler);
    void load(const SpikingNetwork &net);
    void reset();
    void flush_timestep_data(RunData &rd, Scheduler &scheduler);
    void retire_timestep(const Timestep &ts);

    void set_heartbeat(long int timesteps) { heartbeat = timesteps; }
    double get_power() const noexcept;
    std::vector<NeuronAddress> get_spikes() const;
    std::vector<double> get_potentials() const;
    long int get_total_timesteps() const noexcept { return total_timesteps; }
    static void update_run_data(RunData &rd, const Timestep &ts);
    std::map<std::string, double> sim_trace_get_optional_traces();

    static void sim_trace_write_spike_header(std::ostream &spike_trace_file);
    void sim_trace_write_potential_header(std::ostream &potential_trace_file);
    void sim_trace_write_perf_header(std::ostream &perf_trace_file);
    static void sim_trace_write_message_header(std::ostream &message_trace_file);

    void sim_trace_record_spikes(std::ostream &spike_trace_file, long int timesteps);
    void sim_trace_record_potentials(std::ostream &potential_trace_file, long int timestep);
    void sim_trace_record_perf(std::ostream &perf_trace_file, const Timestep &ts);
    void sim_output_run_summary(const std::filesystem::path &output_dir, const RunData &run_data) const;

    std::vector<std::reference_wrapper<Core>> cores();
    LookupTable<double> ts_sync_delay_table{};

    size_t core_count{0UL};
    size_t mapped_cores{0UL};
    size_t mapped_tiles{0UL};
    size_t max_cores_per_tile{0UL};
    size_t noc_width_in_tiles{1UL};
    size_t noc_height_in_tiles{1UL};
    size_t noc_buffer_size{1UL};

private:
    std::unique_ptr<BookSimConfig> booksim_config;
    std::atomic<long int> total_messages_sent{0L};
    size_t total_neurons_mapped{0UL};
    long int total_neurons_updated{0L};
    long int total_neurons_fired{0L};
    long int total_timesteps{0L};
    long int total_spikes{0L};
    long int heartbeat{default_heartbeat_timesteps};
    double total_sim_time{0.0};
    double total_energy{0.0};
    double synapse_energy{0.0};
    double dendrite_energy{0.0};
    double soma_energy{0.0};
    double network_energy{0.0};

    // Performance and other misc tracking
    static std::atomic<int> chip_count;
    double neuron_processing_wall{0.0};
    double message_processing_wall{0.0};
    double energy_stats_wall{0.0};
    double scheduler_wall{0.0};

    // Filestreams
    std::ofstream spike_trace;
    std::ofstream potential_trace;
    std::ofstream message_trace;
    std::ofstream perf_trace;

    void map_neurons(const SpikingNetwork &net);
    void map_connections(const SpikingNetwork &net);
    void forward_connection_attributes(const SpikingNetwork &net);
    MappedConnection &map_connection(const Connection &con);
    void map_axons();
    void track_mapped_neurons();
    void track_mapped_tiles_and_cores() noexcept;

    void sim_reset_measurements();
    TimestepHandle sim_hw_timestep(long int timestep, Scheduler &scheduler);
    void sim_timestep_sync(Scheduler &scheduler) const;
    void sim_update_ts_counters(Timestep &ts);
    static double sim_estimate_network_costs(const Tile &src, Tile &dest);
    void sim_calculate_ts_energy(Timestep &ts);
    static double sim_calculate_tile_energy(Timestep &ts, Tile &tile);
    static double sim_calculate_core_energy(Timestep &ts, Core &core);
    void sim_update_total_energy_and_counts(const Timestep &ts);

    void sim_format_run_summary(std::ostream &out, const RunData &run_data) const;
    void sim_print_axon_summary() const noexcept;

    static void sim_create_neuron_axons(MappedNeuron &pre_neuron);
    static void sim_allocate_axon(MappedNeuron &pre_neuron, Core &post_core);
    static void sim_add_connection_to_axon(MappedConnection &con, Core &post_core);

    void process_neurons(Timestep &ts);
    void process_messages(Timestep &ts);
    void forced_updates(const Timestep &ts);

    void process_neuron(Timestep &ts, MappedNeuron &n);
    void receive_message(Message &m);
    static double process_message(Timestep &ts, Core &c, Message &m);
    static PipelineResult execute_pipeline(const std::vector<PipelineUnit *> &pipeline, Timestep &ts, MappedNeuron &n, std::optional<MappedConnection *> con, const PipelineResult &input);

    static double pipeline_process_axon_in(Core &core, const Message &m);
    PipelineResult pipeline_process_axon_out(Timestep &ts, MappedNeuron &n);

    std::ofstream sim_trace_open_perf_trace(const std::filesystem::path &out_dir);
    static std::ofstream sim_trace_open_spike_trace(const std::filesystem::path &out_dir);
    std::ofstream sim_trace_open_potential_trace(const std::filesystem::path &out_dir);
    static std::ofstream sim_trace_open_message_trace(const std::filesystem::path &out_dir);

    static void check_booksim_compatibility(const Scheduler &scheduler, int sim_count);
};

void sim_trace_record_message(std::ostream &message_trace_file, const Message &m);
TimingModel parse_timing_model(const std::string_view &timing_model_str);

struct RunData
{
    double wall_time{0.0};
    long int spikes{0L};
    long int packets_sent{0L};
    long int neurons_updated{0L};
    long int neurons_fired{0L};
    long int timestep_start;
    long int timesteps_executed{0L};

    double total_energy{0.0};
    double synapse_energy{0.0};
    double dendrite_energy{0.0};
    double soma_energy{0.0};
    double network_energy{0.0};
    double sim_time{0.0};

    RunData(long int start);
};
}

#endif
