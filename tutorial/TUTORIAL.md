# SANA-FE Tutorial #

Presentation slides (and additional files) can also be found at:
https://drive.google.com/drive/folders/1GzjXAFouakm3b6GcFIHsw67H8t6l3BtY?usp=sharing
As of Apr 2024.

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

    jamesaboyle/sana-fe:latest


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
allowing us to connect to it via the Docker GUI again.

## Trying out SANA-FE ##

After running the SANA-FE container using the steps above, you should now be able
to run SANA-FE, edit input files and view simulator outputs. Inside Docker GUI,
go to the "Exec" tab inside the "Container" section, which will launch a Linux
shell. The shell is where we interact with the container's Linux environment and
run SANA-FE.

Inside this docker container, the SANA-FE simulator is installed at the default
(root) directory. Additional tutorial files and scripts are saved inside
`/tutorial`.

To run SANA-FE for the first time in the Docker shell, go to the tutorial
directory and run the simulation script with the commands below. Note that
Linux is case-sensitive so reproduce these commands exactly (e.g., copy and
paste from here):

    python3 sim.py tutorial/arch.yaml tutorial/snn.net 1000

SANA-FE should run and print to the console. SANA-FE dumps information about the
initialized SNN and design, followed by a heart-beat after
every 100 time-steps, and finally summarizes run statistics before exiting.

First, we will look at an example of an SNN-based architecture in SANA-FE.
Inside `tutorial/arch.yaml` there is a minimal architecture. Follow the
three exercises given at the top of the `arch.yaml` file to extend the design.
The extended design will be used in the next step.

For the next step, we will finish an SNN to map to the architecture (shown in
the diagram below). Complete the three exercises at the top of
`tutorial/snn.net` to add the elements shown in blue and finish defining a
small SNN.

To check your solutions, run the previous command again:

    python3 sim.py tutorial/arch.yaml tutorial/snn.net 1000

If all exercises were correctly completed for both the architecture and SNN,
your new architecture should run for 1000 time-steps and print the following
statistics after running.

    energy: 1.331600e-08
    time: 1.048500e-05
    total_spikes: 665
    total_packets: 665
    total_neurons_fired: 666

(Solutions for both arch and net can also be found at tutorial/solutions/.)

Now that both an architecture and SNN have been defined, we can look at the
outputs SANA-FE can generate for these files. Use -o to specify an output
directory. First, we will just look at the run summary that SANA-FE saves
to file.

    python3 sim.py -o tutorial tutorial/arch.yaml tutorial/snn.net 10
    cat tutorial/run_summary.yaml

The file `run_summary.yaml` matches the printed summay at the end of the run.

Next, run another simulation but with the spike and potential trace flags set.
To enable spike traces, prepend (`-s`) to your run command. If the spike trace
flag is set, the simulator will record spikes for neurons with log_spikes set
to 1 in the SNN description. Similarly, the potential trace is set by
prepending (`-v`) to the run command. If set, SANA-FE logs the neuron potentials
of any voltage probes. Voltage probes are set by setting log_potential of
neurons to 1 in the SNN description.

In this example, two different flags are being set, so (`-s -v`).
Two traces will be generated: `spikes.csv` and `potential.csv`.

    python3 sim.py -s -v -o tutorial tutorial/arch.yaml tutorial/snn.net 10
    cat tutorial/spikes.csv
    cat tutorial/potential.csv

Every line in the spike trace (`tutorial/spikes.csv`) has the format:

    <neuron>,<timestep spiked>

The potential trace (`tutorial/potential.csv`) has a column per probe and
one line per time-step i.e.:

    <neuron 0, timestep 0>,<neuron 1, timestep 0>
    <neuron 0, timestep 1>,<neuron 1, timestep 1>
    ...

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
columns for various statistics.

Finally, run a script that launches a simulation a larger application
(DVS gesture categorization) running on a real-world architecture (Loihi).

    python3 tutorial/dvs_challenge.py

The script loads the convolutional kernel weights and generates the SNN.
If running the script in Docker these weights are already included. If running
the script on your own SANA-FE install, you will need to get the weights
manually (see ```dvs_challenge.py``` for more information).
The script maps each layer of the SNN across one or more cores on Loihi,
saves this mapped SNN in SANA-FE's format, and launches a simulation
from the Python script with performance traces.

## Challenge ##

Now that we have showcased a real-world example, there is an open-ended
challenge using this application:

Optimize the mapping of the DVS gesture SNN to Loihi cores to get the lowest
Energy-Delay product (total energy * run-time)

One quick way to try new mappings is by changing the number of cores each layer
is mapped to (note that each core can only have 1024 neurons mapped to it).
Notice how different mappings of neurons to cores produces different numbers.

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
the pyyaml library. The challenge script also requires NumPy, to load the
convolutional kernel weights.

// docker pull jamesaboyle/sana-fe:latest
// docker run --mount type=bind,source="$(pwd)",target="/home/tutorial/sana-fe" --entrypoint "/bin/bash" -it jamesaboyle/sana-fe:latest