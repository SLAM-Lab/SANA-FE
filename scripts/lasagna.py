import csv
import torch
import os
import sys
import dill
import numpy as np
import pandas as pd
# PyTorch and Benchmark data
from torchvision import datasets
import tonic
import tonic.transforms as transforms
# Plotting
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# Try importing the installed sanafe library. If not installed, require a
#  fall-back to a local build of the sanafe Python library
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir)))
    print(f"Project dir: {PROJECT_DIR}")
    sys.path.insert(0, PROJECT_DIR)
    import sanafe

#dataset = "mnist"
dataset = "shd"
analog_neurons = False

if dataset == "mnist":
    timesteps = 100
    num_inputs = 10
    #num_inputs = 10000
elif dataset == "shd":
    timesteps = None
    num_inputs = 1
    #num_inputs = 10
    #num_inputs = 100
    #num_inputs = 2264

print(f"Loading models for {dataset}")
weights = {}
if dataset == "mnist":
    # Load the MNIST network
    mnist_model = torch.load(
            os.path.join(PROJECT_DIR, "etc", "mnist.pt"),
            pickle_module=dill,
            map_location=torch.device("cpu"))
    for attribute_name, param in mnist_model.named_parameters():
        weights[attribute_name] = param.detach().numpy()

    # Load the MNIST test inputs which were stored during training in a CSV file
    inputs = np.loadtxt(os.path.join(PROJECT_DIR, "etc", "mnist.csv"),
                          delimiter=',')
    # I didn't store the labels during training, so load these from the original
    #  dataset
    # TODO: scrap the csv and just use this for inputs too..
    test_dataset = datasets.MNIST(root="./runs/lasagna/data", train=False,
                                  download=True)
    labels = test_dataset.targets
elif dataset == "shd":
    # Load the Spiking Digits network
    spiking_digits_model = torch.load(
            os.path.join(PROJECT_DIR, "etc", "spiking_digits.pt"),
            pickle_module=dill,
            map_location=torch.device("cpu"))
    weights = {}
    for attribute_name, param in spiking_digits_model.named_parameters():
        weights[attribute_name] = param.detach().numpy()

    # Load the spiking digits test inputs from the raw dataset, applying the
    #  same transformations as the training scripts
    frame_transform = transforms.Compose([
        transforms.Downsample(time_factor=0.1, spatial_factor=0.1),
        transforms.ToFrame(sensor_size=(70, 1, 1), time_window=1000)
    ])
    testset = tonic.datasets.SHD(save_to=os.path.join(PROJECT_DIR, "runs", "lasagna", "data"),
                                 transform=frame_transform, train=False)
    dataloader = torch.utils.data.DataLoader(testset, batch_size=1)

    inputs = []
    labels = []
    count = 0
    for input_frame, label in dataloader:
        input_frame = input_frame.cpu()
        label = label.cpu()
        inputs.append(input_frame.numpy())
        labels.append(int(label))
        count += 1
        if count >= num_inputs:
            break

# PyTorch stores weights in an array with dims (num out x num in)
in_neurons = weights["fc1.weight"].shape[1]
hidden_neurons = weights["fc1.weight"].shape[0]
out_neurons = weights["fc2.weight"].shape[0]

print(f"in:{in_neurons}, hidden:{hidden_neurons}, out:{out_neurons}")
#"""
# Load the LASAGNA architecture with analog neurons
arch = sanafe.load_arch("/home/james/code/lasagna/lasagna.yaml")

print("Creating network in SANA-FE")
network = sanafe.Network()
print(f"in:{in_neurons} hidden:{hidden_neurons} out:{out_neurons}")
in_layer = network.create_neuron_group("in", in_neurons)
hidden_layer = network.create_neuron_group("hidden", hidden_neurons)
out_layer = network.create_neuron_group("out", out_neurons)

print("Creating input layer")
for id, neuron in enumerate(in_layer):
    neuron.configure(soma_hw_name=f"input[{id}]",
                            log_spikes=True)
    neuron.map_to_core(arch.tiles[0].cores[0])

print("Creating hidden layer")
for id, neuron in enumerate(hidden_layer):
    hidden_parameters = {
        "threshold": 1.0,
        "leak_decay": 0.85,
        "reset_mode": "hard",
    }
    if dataset == "shd":
        hidden_parameters["input_decay"] = weights["lif1.alpha"][id]

    if analog_neurons:
        neuron.configure(soma_hw_name=f"analog_lif[{id + hidden_neurons}]",
                         log_spikes=True, log_potential=True,
                         model_attributes=hidden_parameters)
    else:
        if dataset == "shd":
            # Clip anything below 0
            hidden_parameters["reverse_threshold"] = 0.0
            hidden_parameters["reverse_reset_mode"] = "saturate"
        neuron.configure(soma_hw_name=f"loihi",
                         log_spikes=True,
                         model_attributes=hidden_parameters)
    neuron.map_to_core(arch.tiles[0].cores[0])

print("Creating output layer")
for id, neuron in enumerate(out_layer):
    if dataset == "shd":
        neuron.configure(soma_hw_name=f"loihi",
                         log_potential=True,
                         model_attributes={
                            "leak_decay": 0.85,
                            "reset_mode": "none",
                            "input_decay": weights["lif2.alpha"][id]
                        })
    elif dataset == "mnist":
        if analog_neurons:
            neuron.configure(soma_hw_name=f"analog_lif[{id}]", log_spikes=True)
        else:
            neuron.configure(soma_hw_name=f"loihi", log_spikes=True,
                             model_attributes={
                                 "threshold": 1.0,
                                 "leak_decay": 0.85,
                                 "reset_mode": "hard"
                            })
    else:
        raise Exception(f"Dataset: {dataset} not recognized!")
    neuron.map_to_core(arch.tiles[0].cores[0])

# Connect neurons in both layers
print("Connecting neurons")
min_weight = 1.0e-10
print("Adding first layer connections")
for src in range(in_neurons):
    for dst in range(hidden_neurons):
        weight = weights["fc1.weight"][dst, src]
        if abs(weight) > min_weight:
            network["in"][src].connect_to_neuron(
                network["hidden"][dst], {"weight": weight})

if dataset == "shd":  # Add recurrent connections only for spiking digits
    print("Adding recurrent connections")
    for src in range(hidden_neurons):
        for dst in range(hidden_neurons):
            weight = weights["fcr.weight"][dst, src]
            if abs(weight) > min_weight:
                network["hidden"][src].connect_to_neuron(
                    network["hidden"][dst], {"weight": weight}
                )

print("Adding second layer connections")
for src in range(hidden_neurons):
    for dst in range(out_neurons):
        weight = weights["fc2.weight"][dst, src]
        if abs(weight) > min_weight:
            network["hidden"][src].connect_to_neuron(
                network["out"][dst], {"weight": weight})

# Run a simulation
print("Building h/w")
hw = sanafe.SpikingChip(arch, record_spikes=True,
                        record_potentials=True, record_perf=True)
print("Loading SNN")
hw.load(network)
print(f"Running simulation for {timesteps} timesteps")

timesteps_per_input = []
mapped_inputs = hw.mapped_neuron_groups["in"]
print("Setting inputs")
for input in range(num_inputs):
    print(f"Simulating input: {input}")

    if dataset == "mnist":
        mnist_input = inputs[input, :]
        for id, mapped_neuron in enumerate(mapped_inputs):
            mapped_neuron.set_model_attributes(
                model_attributes={"poisson": mnist_input[id]})
    elif dataset == "shd":
        spiking_digit_input = inputs[input].squeeze()
        for id, mapped_neuron in enumerate(mapped_inputs):
            spiketrain = list(spiking_digit_input[:, id])
            mapped_neuron.set_model_attributes(
                model_attributes={"spikes": spiketrain})
        timesteps = len(spiketrain)
    else:
        raise Exception(f"Dataset: {dataset} not recognized!")

    print(f"Simulating for {timesteps} timesteps")
    results = hw.sim(timesteps)
    timesteps_per_input.append(timesteps)
    print(results)
    hw.reset()
#"""

# Create an array tracking the input corresponding to each timestep
timestep_inputs = []
for t, steps in enumerate(timesteps_per_input):
    timestep_inputs.extend([t,]*steps)
timestep_inputs = np.array(timestep_inputs)

if dataset == "shd":
    fig = plt.figure(figsize=(12.0, 8.0))
    gs = GridSpec(3, 1, height_ratios=[4, 4, 4], hspace=0.1)
    columns_to_extract = ["neuron out.0", "neuron out.1", "neuron out.2", "neuron out.3",
                    "neuron out.4", "neuron out.5", "neuron out.6", "neuron out.7",
                    "neuron out.8", "neuron out.9", "neuron out.10", "neuron out.11",
                    "neuron out.12", "neuron out.13", "neuron out.14", "neuron out.15",
                    "neuron out.16", "neuron out.17", "neuron out.18", "neuron out.19"]
    potentials_df = pd.read_csv("potentials.csv")
    potentials = potentials_df[columns_to_extract].to_numpy()

    #potentials_df.plot(x="timestep", y=columns_to_extract, ax=ax)
    #potentials_df.plot(x="timestep", ax=ax)
    #plt.savefig("runs/lasagna/shd.png")

    correct = 0
    start_timestep = 0
    for i, timesteps in enumerate(timesteps_per_input):
        timestep_potentials = potentials[start_timestep:start_timestep+timesteps, :]
        max_potentials = np.max(timestep_potentials, axis=0)
        category = np.argmax(max_potentials, axis=0)

        print(f"Potentials per class for inference:\n{max_potentials}\n"
              f"out:{category} actual:{labels[i]}")

        if category == labels[i]:
            correct += 1
        start_timestep += timesteps

    with open("spikes.csv") as spike_csv:
        spike_data = csv.DictReader(spike_csv)

        in_spikes = [[] for _ in range(0, in_neurons)]
        hidden_spikes = [[] for _ in range(0, hidden_neurons)]
        out_spikes = [[] for _ in range(0, out_neurons)]
        for spike in spike_data:
            # Spike entry has the format <group_name.neuron_id,timestep>
            timestep, neuron = int(spike["timestep"]), spike["neuron"]
            group_name, neuron_id = neuron.split(".")
            neuron_id = int(neuron_id)
            # Track spikes for all three layers
            if group_name == "in":
                in_spikes[neuron_id].append(timestep)
            elif group_name == "hidden":
                hidden_spikes[neuron_id].append(timestep)
            elif group_name == "out":
                out_spikes[neuron_id].append(timestep)
            else:
                print(f"Warning: Group {group_name} not recognized!")

    ax_spikes = fig.add_subplot(gs[0])
    ax_spikes.set_title('Analog Neurons Classifying Spiking Digits')
    ax_spikes.set_xlim((0, timesteps_per_input[0]))
    ax_spikes.set_ylim((0, in_neurons + hidden_neurons))

    total_neurons = 0
    for neuron_id in range(0, in_neurons):
        ax_spikes.scatter(in_spikes[neuron_id],
                    [total_neurons]*len(in_spikes[neuron_id]), c='r', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1

    for neuron_id in range(0, hidden_neurons):
        ax_spikes.scatter(hidden_spikes[neuron_id],
                    [total_neurons]*len(hidden_spikes[neuron_id]), c='b', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1
    ax_spikes.set_ylabel("Neuron")
    ax_spikes.set_xticks([])

    # Plot the output neuron spike counts
    ax_potentials = fig.add_subplot(gs[1])
    ax_potentials.set_xlim((0, timesteps_per_input[0]))
    ax_potentials.set_ylabel("Output Potentials")
    ax_potentials.set_xticks([])
    ax_potentials.set_xlabel("")  # Remove label
    ax_potentials.plot(timestep_potentials)
    #ax_potentials.legend(("0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
    #                      "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), ncol=10)

    ax_perf = fig.add_subplot(gs[2])
    ax_perf.set_xlim((0, timesteps_per_input[0]))
    perf_df = pd.read_csv("perf.csv")
    perf_df["total_energy_pj"] = perf_df["total_energy"] * 1.0e12
    perf_df.plot(x="timestep", y=["total_energy_pj"], ax=ax_perf)
    ax_perf.set_ylabel("Simulated Soma Energy (pJ)")
    ax_perf.get_legend().remove()

    ax_perf.set_xlabel("Time-step")
    plt.savefig("runs/lasagna/shd_raster.png")


elif dataset == "mnist":
    # Read in simulation results
    with open("spikes.csv") as spike_csv:
        spike_data = csv.DictReader(spike_csv)

        in_spikes = [[] for _ in range(0, in_neurons)]
        hidden_spikes = [[] for _ in range(0, hidden_neurons)]
        out_spikes = [[] for _ in range(0, out_neurons)]
        for spike in spike_data:
            # Spike entry has the format <group_name.neuron_id,timestep>
            timestep, neuron = int(spike["timestep"]), spike["neuron"]
            group_name, neuron_id = neuron.split(".")
            neuron_id = int(neuron_id)
            # Track spikes for all three layers
            if group_name == "in":
                in_spikes[neuron_id].append(timestep)
            elif group_name == "hidden":
                hidden_spikes[neuron_id].append(timestep)
            elif group_name == "out":
                out_spikes[neuron_id].append(timestep)
            else:
                print(f"Warning: Group {group_name} not recognized!")


    # 1) Count the total output spikes
    counts = np.zeros((out_neurons, num_inputs), dtype=int)
    for digit, spikes in enumerate(out_spikes):
        for spike_timestep in spikes:
            input_idx = (spike_timestep - 1) // timesteps
            ##input_idx = timestep_inputs[spike_timestep - 1]
            assert(input_idx < num_inputs)
            counts[digit, input_idx] += 1


    def calculate_cumulative_spikes(out_spikes, out_neurons, total_timesteps, timesteps_per_input):
        # Initialize array to store cumulative counts for each output neuron at each timestep
        cumulative_counts = np.zeros((out_neurons, total_timesteps))

        # For each output neuron
        for neuron_idx, spikes in enumerate(out_spikes):
            # Sort spikes by timestep for cumulative counting
            sorted_spikes = sorted([s for s in spikes if s <= total_timesteps])
            # Process each timestep
            for t in range(total_timesteps):
                # Calculate which input period we're in
                input_image = t // timesteps_per_input
                period_start = input_image * timesteps_per_input

                # Count spikes from start of current input period up to current timestep
                count = sum(1 for spike in sorted_spikes
                        if spike <= t + 1 and spike > period_start)

                cumulative_counts[neuron_idx, t] = count

        return cumulative_counts

    cumulative_counts = calculate_cumulative_spikes(
        out_spikes, out_neurons, timesteps * num_inputs, timesteps)

    correct = 0
    for i in range(0, num_inputs):
        print(f"Spike counts per class for inference of digit:{counts[:, i]} "
            f"out:{np.argmax(counts[:, i])} actual:{labels[i]}")
        if np.argmax(counts[:, i]) == labels[i]:
            correct += 1

    # 2) Create a raster plot of all spikes
    fig = plt.figure(figsize=(12.0, 8.0))
    gs = GridSpec(4, 1, height_ratios=[2, 4, 4, 4], hspace=0.1)

    # Top subplot for MNIST digits
    ax_digits = fig.add_subplot(gs[0])
    ax_digits.set_xticks([])
    ax_digits.set_yticks([])

    # Calculate positions for each digit
    digit_width = 28  # MNIST digits are 28x28
    time_per_digit = timesteps  # Each digit is presented for 50 timesteps
    num_digits = len(inputs)

    # Create a blank canvas for all digits
    total_timesteps = timesteps * num_inputs
    display_height = 20  # Adjust this factor to change digit height

    # Place each digit at its corresponding time position
    for i, digit in enumerate(inputs[0:num_inputs, :]):
        start_time = i * time_per_digit
        # Calculate extent for each digit: [left, right, bottom, top]
        # Width of each digit display is set to match the time window
        digit_extent = [start_time, start_time + time_per_digit, 0, display_height]
        ax_digits.imshow(digit.reshape(digit_width, digit_width),
                        cmap='gray', aspect='auto', extent=digit_extent)

    # Display the digits
    ax_digits.set_xlim(0, total_timesteps)
    ax_digits.set_ylim(0, display_height)
    ax_digits.set_title('Analog Neurons Classifying MNIST')

    ax_spikes = fig.add_subplot(gs[1])
    ax_spikes.set_xlim((0, total_timesteps))
    ax_spikes.set_ylim((0, in_neurons + hidden_neurons + out_neurons))

    total_neurons = 0
    for neuron_id in range(0, in_neurons):
        ax_spikes.scatter(in_spikes[neuron_id],
                    [total_neurons]*len(in_spikes[neuron_id]), c='r', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1

    for neuron_id in range(0, hidden_neurons):
        ax_spikes.scatter(hidden_spikes[neuron_id],
                    [total_neurons]*len(hidden_spikes[neuron_id]), c='b', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1

    for neuron_id in range(0, out_neurons):
        ax_spikes.scatter(out_spikes[neuron_id],
                    [total_neurons]*len(out_spikes[neuron_id]), c='k', s=2,
                    marker='.', linewidths=0.5)
        total_neurons += 1

    # Add vertical lines to show digit presentation boundaries
    for i in range(num_inputs + 1):
        ax_spikes.axvline(x=i * time_per_digit, color='gray', linestyle='--', alpha=0.3)
    ax_spikes.set_ylabel("Neuron")
    ax_spikes.set_xticks([])

    # Plot the output neuron spike counts
    ax_potentials = fig.add_subplot(gs[2])
    ax_potentials.set_xlim((0, total_timesteps))
    ax_potentials.plot(cumulative_counts.transpose())
    # Add vertical lines to show digit presentation boundaries
    for i in range(num_inputs + 1):
        ax_potentials.axvline(x=i*time_per_digit + 1, color='gray', linestyle='--', alpha=0.3)
    ax_potentials.set_ylabel("Spike Counts")
    ax_potentials.set_xticks([])
    ax_potentials.set_xlabel("")

    ax_perf = fig.add_subplot(gs[3])
    ax_perf.set_xlim((0, total_timesteps))
    perf_df = pd.read_csv("perf.csv")
    perf_df["total_energy_pj"] = perf_df["total_energy"] * 1.0e12
    perf_df.plot(x="timestep", y=["total_energy_pj"], ax=ax_perf)
    for i in range(num_inputs + 1):
        ax_perf.axvline(x=i*time_per_digit + 1, color='gray', linestyle='--', alpha=0.3)
    ax_perf.set_ylabel("Simulated Soma Energy (pJ)")
    ax_perf.get_legend().remove()

    ax_perf.set_xlabel("Time-step")
    plt.savefig("runs/lasagna/mnist_raster.png")

accuracy = (correct / num_inputs) * 100
print(f"Accuracy: {accuracy}%")

print("Finished.")
plt.show()
