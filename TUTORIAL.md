# SANA-FE Tutorial #

## Docker Setup ##

For this tutorial, we recommend using our Docker image which comes with SANA-FE
and the architectures, networks and scripts needed for this tutorial.
The first step is to get the Docker image and launch the SANA-FE container. You
can also optionally mount a directory inside the image to get access to the
network, architecture and script files on your local desktop (e.g. if you want
to edit files using a text editor outside of Linux).

Either using the Docker desktop GUI or command line, pull the docker image:

    jamesaboyle/sana-fe:latest

Now launch a docker container using this image. The easiest way to do this is
using the Docker GUI. When launching a container with the GUI, it is possible
to export the /tutorial/files directory to anywhere on your local machine.
This will give you access to all the tutorial files outside the docker
environment.

## Trying out SANA-FE ##

To run SANA-FE for the first time, load the Python dependencies and run
the simulation script:

    cd /tutorial
    source venv/bin/activate
    python3 sim.py files/example.yaml files/example.net 1000

Next, we will look at an incomplete architecture. Inside this file are three
exercises to complete, to finish defining the architecture.

    cat files/tutorial.yaml

Then, we look at an incomplete SNN description. Again, there are three exercises
to define the SNN.

    cat files/tutorial.net

Now that we have finished defining an architecture and SNN, we will look at the
outputs SANA-FE can generate for these files.

    python3 sim.py files/tutorial.yaml files/tutorial.net
    cat run_summary.yaml

Next we enable spike and potential traces and run another simulation. Two traces
will be generated.

    python3 sim.py –s -v files/tutorial.yaml files/tutorial.net 10
    cat spikes.csv
    cat potential.csv

Then, we do the same but for message and performance traces.

    python3 sim.py –m -p files/tutorial.yaml files/tutorial.net 10
    cat perf.csv
    cat messages.csv

Finally, we go onto a larger application running on a real-world architecture
(Loihi). As part of this, we will look at a challenge inside the run script,
where we can optimize mappings to get the lowest power usage.

    python3 files/tutorial.py

The end.

## Building From Source ##

If you want to build SANA-FE from source, you can get the code from the online
github repository:

    git clone https://github.com/SLAM-Lab/SANA-FE sana-fe
    cd sana-fe

Building requires make and a compiler that supports the C99 standard. To build
simply run the command:

    make

To run SANA-FE in your own Linux environment, you need Python 3.6 or later and
the pyyaml library.

// docker pull jamesaboyle/sana-fe:latest
// docker run --mount type=bind,source="$(pwd)",target="/home/tutorial/sana-fe" --entrypoint "/bin/bash" -it jamesaboyle/sana-fe:latest