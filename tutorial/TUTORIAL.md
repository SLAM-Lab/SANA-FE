# SANA-FE Tutorial #

Presentation slides (and additional files) can also be found at:
https://drive.google.com/drive/folders/1GzjXAFouakm3b6GcFIHsw67H8t6l3BtY?usp=sharing
As of Oct 2024.

Note: This tutorial was presented at NICE Apr 2024, ICONS Aug 2024 and ESWEEK
Sep 2024. It used an older version of SANA-FE (v1), and is based on a Docker
image (described below). There is now a newer version of SANA-FE (v2), which
has more features and an updated tutorial (which is similar but uses a
Jupyter-notebook). See the updated tutorial here.

## Docker Setup ##

For this tutorial, we recommend using our Docker image that includes SANA-FE
and the architecture files, network files and scripts needed for this tutorial.
The first step is to get the Docker image and launch the SANA-FE container. You
can also optionally mount a directory inside the image to get access to the
network, architecture and script files on your local desktop (e.g. if you want
to edit files using a text editor outside of Linux).

Now launch a Docker container using the SANA-FE image. The easiest way to do
this is using the Docker Desktop GUI.

![Docker Desktop](https://github.com/SLAM-Lab/SANA-FE/blob/main/tutorial/docker_1.png)

For Docker Desktop, do this by going to "Search for images, containers, volumes"
in the top toolbar. Type in the image name here, which is:

    jamesaboyle/sana-fe

Click on the first image for SANA-FE, now you should see the screen below.

![SANA-FE image](https://github.com/SLAM-Lab/SANA-FE/blob/main/tutorial/docker_2.png)


Click on the Pull button and the download should start (around 300 MB). This may
take some time. Once the download has finished, in the same place click on the
"Run" button to open another dialog box. Now you should see a screen like below:

![Docker Desktop](https://github.com/SLAM-Lab/SANA-FE/blob/main/tutorial/docker_3.png)

Before launching, Docker gives you a few optional settings. Leave the container
name blank. If you want to edit files locally, create a new folder on your host machine
which will be accessible inside the Docker container. Note that this is where
Docker will write a number of files to. Then, in Docker expand the
"Optional Settings" and go to "Volumes". Here, set the "Host path" to the newly
created folder. Manually set the "Container path" to `/tutorial`. Docker will
link these two folders so that changes from Docker appear in your Host folder
and vice versa. Now we can run this container. If a mounted folder was setup,
you should see new files appear. These files can now be edited in either the Docker
container or on the host machine. The Docker container will be left running,
allowing us to connect to it via the Docker GUI exec shell.

## Trying out SANA-FE ##

After running the SANA-FE container using the steps above, you should now be able
to run SANA-FE, edit input files and view simulator outputs. Inside Docker GUI,
go to the "Exec" tab inside the "Container" section, which will launch a Linux
shell. The shell is where we interact with the container's Linux environment and
run SANA-FE.

Inside this docker container, the SANA-FE simulator is installed at the default
(root) directory. Additional tutorial files and scripts are saved inside
`/tutorial`.

To run SANA-FE for the first time in the Docker shell, run the simulation script
using the command below. Note that Linux is case-sensitive so reproduce these
commands exactly (e.g., copy and paste from here):

    python3 sim.py tutorial/arch.yaml tutorial/snn.net 1000

SANA-FE should run and print the following to the console.

    Command: /home/james/code/sana-fe/sim /home/james/code/sana-fe/arch.yaml.parsed tutorial/snn.net 1000
    [main.c:46:main()] Initializing simulation.
    [main.c:200:main()] Reading network from file.
    [network.c:194:network_create_neuron_group()] Created neuron group gid:0 count:2 threshold:1.000000 neg threshold:-1.000000 pos reset:0.000000 reset mode:2 neg reset:0.000000 neg reset mode:0
    [network.c:194:network_create_neuron_group()] Created neuron group gid:1 count:1 threshold:2.000000 neg threshold:-1.000000 pos reset:0.000000 reset mode:2 neg reset:0.000000 neg reset mode:0
    [arch.c:552:arch_print_connection_map_summary()] ** Mapping summary **
    [arch.c:581:arch_print_connection_map_summary()] Total cores: 1
    [arch.c:582:arch_print_connection_map_summary()] Average in map count: 1.000000
    [arch.c:583:arch_print_connection_map_summary()] Average out map count: 1.000000
    [main.c:210:main()] Creating probe and perf data files.
    [main.c:228:main()] Running simulation.
    [main.c:235:main()] *** Time-step 100 ***
    [main.c:235:main()] *** Time-step 200 ***
    [main.c:235:main()] *** Time-step 300 ***
    [main.c:235:main()] *** Time-step 400 ***
    [main.c:235:main()] *** Time-step 500 ***
    [main.c:235:main()] *** Time-step 600 ***
    [main.c:235:main()] *** Time-step 700 ***
    [main.c:235:main()] *** Time-step 800 ***
    [main.c:235:main()] *** Time-step 900 ***
    [main.c:235:main()] *** Time-step 1000 ***
    [main.c:240:main()] ***** Run Summary *****
    git_version:
    energy: 3.498000e-09
    time: 3.498000e-06
    total_spikes: 166
    total_packets: 166
    total_neurons_fired: 166
    wall_time: 0.041942
    timesteps: 1000
    [main.c:250:main()] Average power consumption: 0.001000 W.
    [main.c:259:main()] Run finished.
    sim finished

SANA-FE dumps information about the initialized SNN and design, followed by a
heart-beat after executing every 100 time-steps, and finally summarizes run
statistics before exiting.

To check your run executed correctly, check the functional output of the
simulator against correct reference values by running:

    diff -I wall run_summary.yaml run_results

If you ran correctly, then the files should match i.e., you should see no output
from diff. Note, we ignore the wall_time entry, since this is measuring
execution time and will differ from run to run i.e., is not deterministic.

Next, we will look at an example of an SNN-based architecture in SANA-FE.
Inside `tutorial/arch.yaml` there is a minimal architecture, based on the
diagram below. Follow the three exercises given at the top of the `arch.yaml`
file to extend the design and add the dashed orange elements (shown below).
Run the previous command again:

![Example architecture](https://github.com/SLAM-Lab/SANA-FE/blob/main/tutorial/example_arch_small.png)


    python3 sim.py tutorial/arch.yaml tutorial/snn.net 1000

Similar to before, results can be checked by running the command and confirming
there are no differences between your results and a provided solution:

    diff -I wall run_summary.yaml arch_results

Next, we will finish a partially defined SNN and map it to our new architecture
The SNN is shown in the diagram below, with pre-defined elements drawn in black.
Complete the four exercises at the top of `tutorial/snn.net` to add the ornage
elements shown and complete this small SNN.

![Example SNN](https://github.com/SLAM-Lab/SANA-FE/blob/main/tutorial/example_snn_small.png)

To check your solutions, run the previous command again and check your results:

    python3 sim.py tutorial/arch.yaml tutorial/snn.net 1000
    diff -I wall run_summary.yaml snn_results

Note that this part requires both the architecture and SNN to be correct.

(If you get stuck, correct run logs can be found under `arch.log` and
`snn.log`. Also, solutions for both the architecture and net files can also be
found at `completed_arch.yaml` and `completed_snn.net`)

Now that both an architecture and SNN have been defined, we can look at the
outputs SANA-FE can generate for these files. We use -o to specify an output
directory. First, we will just look at the run summary that SANA-FE writes
to a file, saving all output to the `tutorial` directory.

    python3 sim.py -o tutorial tutorial/arch.yaml tutorial/snn.net 10
    cat tutorial/run_summary.yaml

The file `run_summary.yaml` matches the printed summary at the end of the run.

Next, run another simulation but with the spike and potential trace flags set.
To enable spike traces, prepend (`-s`) to your run command. If the spike trace
flag is set, the simulator will record spikes for neurons with `log_spikes` set
to 1 in the SNN description. Similarly, the potential trace is set by
prepending (`-v`) to the run command. If enabled, SANA-FE logs the neuron
potentials of any voltage probes. Voltage probes are set by defining
`log_potential` to 1 for neurons in the SNN description.

In this example, two different flags are set simultaneously, i.e., (`-s -v`).
Therefore, two traces will be generated: `spikes.csv` and `potential.csv`.

    python3 sim.py -s -v -o tutorial tutorial/arch.yaml tutorial/snn.net 10
    cat tutorial/spikes.csv
    cat tutorial/potentials.csv

Every line in the spike trace (`tutorial/spikes.csv`) has the format:

    <neuron>,<timestep spiked>

The potential trace (`tutorial/potentials.csv`) has a column per probe and
one line per time-step i.e.:

    <neuron 0, timestep 0>,<neuron 1, timestep 0>
    <neuron 0, timestep 1>,<neuron 1, timestep 1>
    ...

Try visualizing `potentials.csv` by plotting the dynamics of neurons 0.0 and
0.1 for 10 simulated timesteps, using your preferred plotting application.

Next, run the same simulation but with message and performance traces
enabled instead. Message tracing is set by prepending `-m` to the command.
This will cause the simulator to record information about spike messages sent
by hardware. Performance tracing is set by prepending `-p` to the command.
Performance traces contain more detailed per-timestep information about
activity on the simulated design. Similar to the last example, we can combine
multiple traces, using `-m -p`.

    python3 sim.py -m -p -o tutorial tutorial/arch.yaml tutorial/snn.net 10
    cat tutorial/perf.csv
    cat tutorial/messages.csv

Two different traces will be generated: `perf.csv` and `messages.csv`.
Both files have similar formats, with one line per time-step and different
columns for various statistics. Try visualizing the `fired` and `total_energy`
fields over time.

Finally, run a script that launches a simulation a larger application
(DVS gesture categorization) running on a real-world architecture (Loihi).

    python3 tutorial/dvs_challenge.py

The script loads the convolutional kernel weights and generates the SNN.
If running the script in Docker these weights are already included. If running
the script on your own SANA-FE install, you will need to get the weights
manually (see ```dvs_challenge.py``` for more information).
The script maps each layer of the SNN across one or more cores on Loihi,
saves this mapped SNN in SANA-FE's format, and launches a simulation
from the Python script with performance traces. Finally, the script checks
results from the run and makes sure the SNN executed correctly.

## Challenge ##

Now that we have showcased a real-world example, there is an open-ended
challenge using the gesture categorization application used above:

Optimize the mapping of the DVS gesture SNN to Loihi cores to get the lowest
Energy-Delay product (total energy * total run-time) for a valid mapping.

One quick way to try new mappings is by changing the number of cores each layer
is mapped to (note that each core can only have 1024 neurons mapped to it).
This code is highlighted in the challenge script. Notice how different mappings
of neurons to cores produces different energy-delay product values.
Note that the SNN application must be executed successfully, the simulation must
finish without simulator errors, and there must be no script errors from the
run validation checks.

The end.

## (Building From Source) ##

If you want to build SANA-FE from source, you can get the code from the online
github repository:

    git clone https://github.com/SLAM-Lab/SANA-FE sana-fe
    cd sana-fe

Building requires make and a compiler that supports the C99 standard. To build
simply run the command:

    make

To run SANA-FE in your own Linux environment, you need Python 3.6 or later and
the pyyaml library. The challenge script also requires NumPy to load the
convolutional kernel weights. To create a venv and install these dependencies,
run:

    python -m venv ./venv && source ./venv/bin/activate
    pip install --upgrade pip && pip install pyyaml numpy

## (Docker from the Command-line) ##

If you wish to run Docker from its CLI interface instead of Docker desktop use
the following commands (assuming Linux):

    docker pull jamesaboyle/sana-fe
    docker run --mount type=bind,source=<src dir here>,target="/tutorial" --entrypoint "/bin/bash" -it jamesaboyle/sana-fe:latest

Note that where it says "src dir here", you will have to replace this with
directory you want Docker to use for file sharing.

For example, to choose the current directory in a bash environment, you can run:

    docker run --mount type=bind,source="$(pwd)",target="/tutorial" --entrypoint "/bin/bash" -it jamesaboyle/sana-fe:latest
