import sanafe

# Included examples of architectures and applications
example_arch, example_snn = sanafe.load_example()
loihi_arch = sanafe.load_loihi()  # Intel's Loihi 1 chip
truenorth_arch = sanafe.load_truenorth()  # IBM's TrueNorth chip

# Introspecting architectural features
print(example_arch)
print(example_arch.tiles)
print(example_arch.cores())
print(example_arch.tiles[0].cores)
print(example_arch.tiles[1].cores[0].pipeline_hw)

# Introspecting SNN features
print(example_snn)
print(example_snn.groups)
print(example_snn.groups["out"].neurons)

# Creating a SNN programmatically
python_snn = sanafe.Network()
python_snn.create_neuron_group("foo", 2)
python_snn.create_neuron_group("bar", 3, model_attributes={"bias": 0})
src, dst = python_snn.groups["foo"][0], python_snn.groups["bar"][1]
src.connect_to_neuron(dst, {"weight": 2.5})
# Attributes can be set, using either a small subset of simulator built-in
#  attributes (e.g. log_spikes), or arbitrary model-defined attributes.
python_snn.groups["bar"].neurons[0].set_attributes(
    log_spikes=True, model_attributes={"bias": 0.5})

# Mapping a SNN to h/w
for n in python_snn.groups["foo"]:
    n.map_to_core(example_arch.tiles[0].cores[0])
for n in python_snn.groups["bar"]:
    n.map_to_core(example_arch.tiles[0].cores[1])
# Neurons and edges can be mapped to specific h/w pipeline units. Unit names
#  much match exactly with the corresponding pipeline unit in the
#  architecture. If not provided, SANA-FE defaults by mapping to the first
#  defined synapse/dendrite/soma unit.
python_snn.groups["foo"].neurons[0].set_attributes(soma_hw_name="demo_soma_alt")
# Replicated h/w units use an indexed notation, e.g., u[0..2] ->u[0],u[1],u[2]
python_snn.groups["foo"].neurons[1].set_attributes(soma_hw_name="demo_input[1]")

# Creating and initializing a spiking chip
chip = sanafe.SpikingChip(example_arch)
chip.load(example_snn)
chip.load(python_snn)  # Default behaviour is to reprogram with a new SNN
import json
# Get a complete view and specification of the chip hardware, including all
#  its pipeline h/w and supported attributes
print(json.dumps(chip.describe(), indent=1))

# Running and dynamically controlling simulations
chip.sim(2)
# Mapped neurons can also be reconfigured mid-sim e.g., for dynamic inputs
mapped_input_neurons = chip.mapped_neuron_groups["foo"]
mapped_input_neurons[0].set_attributes(model_attributes={"rate": 0.5})
# Simulations pick up where they left off, unless reset between sim() calls
chip.sim(1)
chip.reset()
# SANA-FE supports multiple ways of providing trace output, including as a
#  Python object, default generated CSV files, or via file handles
with open("message_file.csv", "w") as message_handle:
    final_results = chip.sim(3, spike_trace=True,
            potential_trace="potential_file.csv",
            message_trace=message_handle,
            processing_threads=1,
            scheduler_threads=0)
print(final_results)
