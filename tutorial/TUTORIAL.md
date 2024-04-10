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

Either using the Docker desktop GUI or command line, pull the docker image:

    jamesaboyle/sana-fe:latest

Now launch a Docker container using the SANA-FE image. The easiest way to do
this is using the Docker Desktop GUI. If you want to edit files locally rather
than in a Linux terminal, follow the steps below to mount a directory.

To edit files locally, create a new folder on your host machine which will be
accessible inside the Docker container. Note that this is where Docker will
write a number of files to. After clicking run on the image, expand
"Optional Settings" to create a new volume. Set the "Host path" to the newly
created folder. Manually set the "Container path" to `/tutorial`.
Run this container. You should see files appear in this new folder.
These files can now be edited in either the Docker container or on the host
machine.

## Trying out SANA-FE ##

After running a Docker container using the steps above, you should now be able
to run SANA-FE, edit input files and view simulator outputs. Inside Docker GUI,
go to the "Exec" tab inside the "Container" section, which will launch a Linux
shell. The shell is where we interact with the Linux environment and run
SANA-FE.

To run SANA-FE for the first time in the Docker shell, go to the tutorial
directory and run the simulation script with the commands below. Note that
Linux is case-sensitive so reproduce these commands exactly (e.g., copy and
paste):

    python3 sim.py tutorial/arch.yaml tutorial/snn.net 1000

SANA-FE should run and print to the console. SANA-FE dumps information about the
initialized SNN and design, followed by a heart-beat after
every 100 time-steps, and finally summarizes run statistics before exiting.

Next, we will finish writing an architecture description. Inside
`tutorial/tutorial.yaml` there is an incomplete architecture. Follow the three
exercises to complete the design.

Then, we will finish a corresponding mapped SNN. Again, complete the three
exercises to define the SNN in `tutorial/tutorial.net`

Once both architecture and SNN have been completed, we will look at the
outputs SANA-FE can generate for these files. Start by running and looking at
the summary file. This file matches the summary printed by SANA-FE after
it finishes a simulation.

    python3 sim.py tutorial/arch.yaml tutorial/snn.net 10
    cat run_summary.yaml

Next, run another simulation but with the spike and potential trace flags set.
Two traces will be generated: `spikes.csv` and `potential.csv`.

    python3 sim.py -s -v tutorial/arch.yaml tutorial/snn.net 10
    cat spikes.csv
    cat potential.csv

Then, run the same simulation but with message and performance traces
enabled instead. Two different traces will be generated: `perf.csv` and
`messages.csv`.

    python3 sim.py -m -p tutorial/arch.yaml tutorial/snn.net 10
    cat perf.csv
    cat messages.csv

Finally, run a script that launches a simulation a larger application
(DVS gesture categorization) running on a real-world architecture (Loihi).

    python3 tutorial/dvs_challenge.py

The script loads the convolutional kernel weights and generates the SNN.
Then, the script maps each layer of the SNN across one or more cores on Loihi.
The script saves this mapped SNN in SANA-FE's format, and launches a simulation
from the Python script with performance tracing.

## Challenge ##

After running the DVS gesture challenge, there is an open-ended challenge:

Optimize the mapping of the DVS gesture SNN to Loihi cores to get the lowest
power usage.

One quick way to try new mappings is by changing the number of cores each layer
is mapped to (note that each core can only have 1024 neurons mapped to it).
Notice how different mappings of neurons to cores produces different numbers.

Otherwise, it is also possible to map neurons individually. By using a mix of
profiling (using spike traces), it might be possible to optimize even further.

SANA-FE can be used to rapidly explore design-spaces and mapping spaces.

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