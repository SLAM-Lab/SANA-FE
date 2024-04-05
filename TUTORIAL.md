# SANA-FE Tutorial #

To get the code from the online repository:

    git clone https://github.com/SLAM-Lab/SANA-FE sana-fe
    cd sana-fe

For this tutorial we will use Docker, which comes with dependencies installed.
Next we get the Docker image and launch a container. Here we also mount the
current directory in Docker, so you can access the SANA-FE code in both your own
environment and the Docker virtual environment!

(It is also possible to setup SANA-FE on your own environment. SANA-FE
requires GCC, make, Python 3.6 or later and the pyyaml library to build and run.
Setup should be straightforward on any standard Linux distro.)

    docker pull jamesaboyle/sana-fe:v1
    docker run --mount type=bind,source="$(pwd)",target="/home/tutorial/sana-fe" --entrypoint "/bin/bash" -it jamesaboyle/sana-fe:v1

To run SANA-FE for the first time:

    make
    source tutorial/venv/bin/activate
    python3 sim.py snn/example.net arch/example.yaml 1000

Next, we will look at an incomplete architecture.

    cat arch/tutorial.yaml

Then, we look at an incomplete SNN description.

    cat snn/tutorial.net

Now that we have finished the architecture and SNN, we will look at the outputs
SANA-FE can generate, using these as examples:
    python3 sim.py arch/tutorial.yaml snn/tutorial.net
    cat run_summary.yaml

Now we look at the per-neuron statistics, looking at probed spikes and
membrane potentials.

    python3 sim.py –s -v arch/tutorial.yaml snn/tutorial.net
    cat spikes.csv
    cat potential.csv

Next, we look at the per-timestep statistics

    python3 sim.py –m -p arch/tutorial.yaml snn/tutorial.net
    cat perf.csv
    cat messages.csv

Finally, we go onto a larger application running on a real-world architecture
(Loihi). As part of this, we will look at a challenge inside the run script,
where we can optimize mappings to get the lowest power usage.

    cat tutorial.py
    python3 tutorial.py

The end.

