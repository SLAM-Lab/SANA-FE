// Tests for splitting the hardware pipeline into stages based on the
//  time-step buffer position, including buffers in or before the synapse
#include <algorithm>
#include <filesystem>
#include <memory>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "arch.hpp"
#include "chip.hpp"
#include "network.hpp"
#include "pipeline.hpp"

namespace
{
std::string read_file(const std::filesystem::path &path)
{
    std::ifstream file(path);
    std::stringstream contents;
    contents << file.rdbuf();
    return contents.str();
}

std::filesystem::path write_temp(
        const std::string &name, const std::string &contents)
{
    const std::filesystem::path path =
            std::filesystem::temp_directory_path() / name;
    std::ofstream(path) << contents;
    return path;
}

// Example chip, with the time-step buffer moved to the given position
std::filesystem::path write_arch(
        const std::string &buffer_position, const bool inside = false)
{
    const std::filesystem::path root(SANAFE_ROOT_PATH);
    std::string arch = read_file(root / "arch" / "example_chip.yaml");
    const std::string key = "buffer_position: soma";
    arch.replace(arch.find(key), key.size(),
            "buffer_position: " + buffer_position +
                    "\n            buffer_inside_unit: " +
                    (inside ? "true" : "false"));
    return write_temp("sanafe_stages_arch_" + buffer_position +
                    (inside ? "_inside" : "") + ".yaml",
            arch);
}

// Current-based synapse that manages the time-step buffer itself. Spikes
//  (read=true) are latched, and delivered when the synapse is next updated
//  without a spike (read=false) in the following time-step
class DoubleBufferedSynapse : public sanafe::SynapseUnit
{
public:
    std::vector<double> weights, next, current;
    long int timestep{0L};
    long int spike_reads{0L};
    long int deliveries{0L};

    sanafe::PipelineResult update(const size_t address, const bool read,
            const long int t) override
    {
        if (t > timestep)  // New time-step, swap the double buffer
        {
            std::swap(current, next);
            std::fill(next.begin(), next.end(), 0.0);
            timestep = t;
        }
        sanafe::PipelineResult out;
        if (read)
        {
            ++spike_reads;
            next.at(address) += weights.at(address);
            out.current = 0.0;
        }
        else
        {
            ++deliveries;
            out.current = current.at(address);
            current.at(address) = 0.0;
        }
        return out;
    }
    void set_attribute_hw(const std::string & /*name*/,
            const sanafe::ModelAttribute & /*param*/) override {}
    void set_attribute_edge(const size_t address, const std::string &name,
            const sanafe::ModelAttribute &param) override
    {
        if (address >= weights.size())
        {
            weights.resize(address + 1UL, 0.0);
            next.resize(address + 1UL, 0.0);
            current.resize(address + 1UL, 0.0);
        }
        if ((name == "weight") || (name == "w"))
        {
            weights.at(address) = static_cast<double>(param);
        }
    }
    void reset() override
    {
        std::fill(next.begin(), next.end(), 0.0);
        std::fill(current.begin(), current.end(), 0.0);
        timestep = 0L;
    }
};

// Replace every core's "demo_synapse" with the double-buffered synapse,
//  keeping its default costs. Must happen before the network is loaded
std::vector<std::shared_ptr<DoubleBufferedSynapse>> use_buffered_synapses(
        sanafe::SpikingChip &hw)
{
    std::vector<std::shared_ptr<DoubleBufferedSynapse>> units;
    for (auto &core_ref : hw.cores())
    {
        for (auto &hw_unit : core_ref.get().pipeline_hw)
        {
            if (hw_unit->name != "demo_synapse")
            {
                continue;
            }
            auto unit = std::make_shared<DoubleBufferedSynapse>();
            unit->name = hw_unit->name;
            unit->register_attributes(std::set<std::string>{"weight", "w"});
            unit->default_energy_process_spike =
                    hw_unit->default_energy_process_spike;
            unit->default_latency_process_spike =
                    hw_unit->default_latency_process_spike;
            hw_unit = unit;
            units.push_back(unit);
        }
    }
    return units;
}

// Small multi-core network with fan-in, recurrence and inhibition
constexpr const char *test_network = R"(
network:
  name: stage_test
  groups:
    - name: in
      attributes: [log_spikes: true]
      neurons:
        - 0: [spikes: [1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1]]
        - 1: [spikes: [0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0]]
        - 2: [spikes: [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1]]
    - name: out
      attributes: [log_spikes: true, log_potential: true]
      neurons:
        - 0..5: [soma: {threshold: 3, leak_decay: 0.9}]
  edges:
    - in -> out: [type: dense, weight: [2, 1, -1, 3, 1, 2, 1, 2, 2, -1, 1, 3, 1, 1, 2, 2, 1, -1]]
    - out.0 -> out.3: [weight: 2]
    - out.1 -> out.4: [weight: 2]
    - out.2 -> out.5: [weight: -2]
    - out.3 -> out.0: [weight: 1]
mappings:
  - in.0: [core: 0.0, soma: "demo_input[0]"]
  - in.1: [core: 0.1, soma: "demo_input[0]"]
  - in.2: [core: 1.0, soma: "demo_input[0]"]
  - out.0: [core: 0.0]
  - out.1: [core: 0.1]
  - out.2: [core: 1.0]
  - out.3: [core: 1.1]
  - out.4: [core: 0.0]
  - out.5: [core: 1.1]
)";

struct StepTrace
{
    std::vector<std::vector<double>> potentials;
    std::vector<size_t> spike_counts;
    double synapse_energy{0.0};
    double soma_energy{0.0};
};

StepTrace run_steps(const std::string &buffer_position, const long int steps,
        const bool inside = false, const bool buffered_synapses = false)
{
    sanafe::Architecture arch =
            sanafe::load_arch(write_arch(buffer_position, inside));
    const auto net_path = write_temp("sanafe_stages_net.yaml", test_network);
    const sanafe::SpikingNetwork net = sanafe::load_net(net_path, arch);
    sanafe::SpikingChip hw(arch);
    if (buffered_synapses)
    {
        use_buffered_synapses(hw);
    }
    hw.load(net);

    StepTrace trace;
    for (long int t = 0; t < steps; ++t)
    {
        const sanafe::RunData rd = hw.sim(1);
        trace.potentials.push_back(hw.get_potentials());
        trace.spike_counts.push_back(hw.get_spikes().size());
        trace.synapse_energy += rd.synapse_energy;
        trace.soma_energy += rd.soma_energy;
    }
    return trace;
}
}

TEST(PipelineStagesTest, ParseBufferBeforeSynapse)
{
    EXPECT_EQ(sanafe::pipeline_parse_buffer_position_str("synapse", false),
            sanafe::buffer_before_synapse_unit);
}

TEST(PipelineStagesTest, ParseBufferInsideSynapse)
{
    EXPECT_EQ(sanafe::pipeline_parse_buffer_position_str("synapse", true),
            sanafe::buffer_inside_synapse_unit);
}

// A synapse that buffers internally must give the same behavior as SANA-FE
//  buffering the spikes before the synapse
TEST(PipelineStagesTest, InsideSynapseMatchesBeforeSynapseFunctionally)
{
    constexpr long int steps = 20;
    const StepTrace before_synapse = run_steps("synapse", steps);
    const StepTrace inside_synapse = run_steps("synapse", steps, true, true);

    for (long int t = 0; t < steps; ++t)
    {
        EXPECT_EQ(before_synapse.potentials.at(t),
                inside_synapse.potentials.at(t))
                << "potentials differ at step " << t + 1;
        EXPECT_EQ(before_synapse.spike_counts.at(t),
                inside_synapse.spike_counts.at(t))
                << "spike counts differ at step " << t + 1;
    }
    size_t total_spikes = 0;
    for (const size_t count : inside_synapse.spike_counts)
    {
        total_spikes += count;
    }
    EXPECT_GT(total_spikes, 0UL);

    // Default per-spike costs are charged when spikes are read, not when
    //  buffered inputs are delivered. Inside the synapse, spikes are read when
    //  received, so the final step's spikes are also charged
    EXPECT_GE(inside_synapse.synapse_energy, before_synapse.synapse_energy);
}

TEST(PipelineStagesTest, InsideSynapseDeliversEachSpikeInTheNextUpdate)
{
    sanafe::Architecture arch = sanafe::load_arch(write_arch("synapse", true));
    const auto net_path = write_temp("sanafe_stages_net.yaml", test_network);
    const sanafe::SpikingNetwork net = sanafe::load_net(net_path, arch);
    sanafe::SpikingChip hw(arch);
    const auto units = use_buffered_synapses(hw);
    hw.load(net);

    constexpr long int steps = 10;
    const sanafe::RunData rd = hw.sim(steps);

    long int spike_reads = 0L;
    long int deliveries = 0L;
    long int connections = 0L;
    for (const auto &unit : units)
    {
        spike_reads += unit->spike_reads;
        deliveries += unit->deliveries;
        connections += unit->connection_count;
    }
    // Spikes are read once per synaptic event, and the reported spike count
    //  doesn't include the deliveries of buffered inputs
    EXPECT_GT(spike_reads, 0L);
    EXPECT_EQ(rd.spikes, spike_reads);
    // Only synapses that read a spike deliver it, in the next update. Spikes
    //  read in the final time-step are never delivered
    EXPECT_GT(connections, 0L);
    EXPECT_GT(deliveries, 0L);
    EXPECT_LT(deliveries, spike_reads);
    EXPECT_LT(deliveries, connections * steps);
}

// Moving the buffer from before the soma to before the synapse changes when
//  synapses are processed (and how their cost is accounted), but must not
//  change the network's behavior when synapses are static
TEST(PipelineStagesTest, BeforeSynapseMatchesBeforeSomaFunctionally)
{
    constexpr long int steps = 20;
    const StepTrace before_soma = run_steps("soma", steps);
    const StepTrace before_synapse = run_steps("synapse", steps);

    ASSERT_EQ(before_soma.potentials.size(), before_synapse.potentials.size());
    for (long int t = 0; t < steps; ++t)
    {
        EXPECT_EQ(before_soma.potentials.at(t), before_synapse.potentials.at(t))
                << "potentials differ at step " << t + 1;
        EXPECT_EQ(before_soma.spike_counts.at(t),
                before_synapse.spike_counts.at(t))
                << "spike counts differ at step " << t + 1;
    }

    // The network must actually do something for this test to mean anything
    size_t total_spikes = 0;
    for (const size_t count : before_soma.spike_counts)
    {
        total_spikes += count;
    }
    EXPECT_GT(total_spikes, 0UL);
    EXPECT_FALSE(before_soma.potentials.front().empty());

    // Synapse events from the final step are latched for the next (never
    //  simulated) step, so before-synapse processes no more than before-soma
    EXPECT_GT(before_synapse.synapse_energy, 0.0);
    EXPECT_LE(before_synapse.synapse_energy, before_soma.synapse_energy);
}

TEST(PipelineStagesTest, ResetClearsBufferedSpikes)
{
    sanafe::Architecture arch = sanafe::load_arch(write_arch("synapse"));
    const auto net_path = write_temp("sanafe_stages_net.yaml", test_network);
    const sanafe::SpikingNetwork net = sanafe::load_net(net_path, arch);
    sanafe::SpikingChip hw(arch);
    hw.load(net);

    hw.sim(5);
    bool spikes_buffered = false;
    for (auto &core_ref : hw.cores())
    {
        spikes_buffered |= !core_ref.get().message_buffer.empty();
    }
    EXPECT_TRUE(spikes_buffered);

    hw.reset();
    for (auto &core_ref : hw.cores())
    {
        EXPECT_TRUE(core_ref.get().message_buffer.empty());
    }
}
