# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""example_neurofem, modified to actually simulate on the SANA-FE backend.

Note the notebook only builds the network and never runs it. So there's not
really much to compare. Here I input some random stimulus and run so that
something actually happens within SANA-FE.
"""
import itertools
import sys
import warnings

import numpy as np
import meshpy.triangle as tri
from matplotlib.tri import Triangulation

import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
sys.path.insert(0, PROJECT_DIR)

from sango import Network, NodeGroup, EdgeGroup, NodeList
from sango.model import IN, LIF, pLIF, PSP
from sango.backend.brian.brian import SimBrian
from sanafe.sango.sango import SimSANAFE

ARCH = f"{PROJECT_DIR}/arch/sango.yaml"
MAX_V = 0.04
N_BDRY = 64
NPM = 12
TIMESTEPS = 60

## build mesh, exactly as the notebook builds it ##
bdry_idxs = np.linspace(0, 2 * np.pi, N_BDRY, endpoint=False)
bdry_pts = [(np.cos(x), np.sin(x)) for x in bdry_idxs]
bdry_fcs = [sorted([x, (x + 1) % N_BDRY]) for x in range(N_BDRY)]

mesh_info = tri.MeshInfo()
mesh_info.set_points(bdry_pts)
mesh_info.set_facets(bdry_fcs)
mesh = tri.build(mesh_info, verbose=False, min_angle=30,
                 volume_constraints=True, max_volume=MAX_V,
                 allow_boundary_steiner=False)
ts = Triangulation([x[0] for x in mesh.points],
                   [x[1] for x in mesh.points], mesh.elements)

mesh_edges = []
for edge in ts.edges:
    if edge[0] < N_BDRY and edge[1] < N_BDRY:
        pass
    elif edge[0] < N_BDRY:
        mesh_edges.append([edge[0], edge[1]])
    elif edge[1] < N_BDRY:
        mesh_edges.append([edge[1], edge[0]])
    else:
        mesh_edges.append([edge[0], edge[1]])
        mesh_edges.append([edge[1], edge[0]])

mesh_points = list(mesh.points)
interior = [i for i in range(len(mesh_points)) if i >= N_BDRY]
print(f"mesh: {len(mesh_points)} points ({len(interior)} interior), "
      f"{len(mesh_edges)} directed edges")


## build the network, as the notebook defines it ##
class MeshPoint(Network):
    def __init__(self, num_neuron=8, neuron=None):
        super().__init__()
        self.num_neuron = num_neuron
        self.half_neuron = num_neuron // 2
        self.neuron = neuron if neuron is not None else pLIF()

    def build(self):
        self.pos = NodeGroup(self.neuron, self.half_neuron)
        self.neg = NodeGroup(self.neuron, self.half_neuron)
        self.all = NodeList(self.pos + self.neg)

        pairs = list(itertools.product(range(self.half_neuron),
                                       range(self.half_neuron)))
        self.p_p = EdgeGroup(self.pos, self.pos, PSP(weight=1.0),
                             edges=[(i, j) for i, j in pairs if i != j])
        self.p_n = EdgeGroup(self.pos, self.neg, PSP(weight=1.0), edges=pairs)
        self.n_p = EdgeGroup(self.neg, self.pos, PSP(weight=-1.0), edges=pairs)
        self.n_n = EdgeGroup(self.neg, self.neg, PSP(weight=-1.0),
                             edges=[(i, j) for i, j in pairs if i != j])


class Mesh(Network):
    def __init__(self, mesh_points, mesh_edges, npm=8, neuron=None,
                 driven=(), stim_times=()):
        super().__init__()
        self.mesh_points = mesh_points
        self.mesh_edges = mesh_edges
        self.npm = npm
        self.hnpm = npm // 2
        self.neuron = neuron
        self.driven = list(driven)
        self.stim_times = list(stim_times)

    def build(self):
        self.mp = [MeshPoint(self.npm, self.neuron)
                   for _ in range(len(self.mesh_points))]
        self.readout = NodeGroup(LIF(threshold=1.0, leak=0.3),
                                 len(self.mesh_points))

        self.p_r = [EdgeGroup(self.mp[i].pos, self.readout, PSP(weight=1.0),
                              edges=[(j, i) for j in range(self.hnpm)])
                    for i in range(len(self.mesh_points))]
        self.n_r = [EdgeGroup(self.mp[i].neg, self.readout, PSP(weight=-1.0),
                              edges=[(j, i) for j in range(self.hnpm)])
                    for i in range(len(self.mesh_points))]

        cross = list(itertools.product(range(self.hnpm), range(self.npm)))
        self.p_m, self.n_m = [], []
        for edge in self.mesh_edges:
            self.p_m.append(EdgeGroup(self.mp[edge[0]].pos,
                                      self.mp[edge[1]].all,
                                      PSP(weight=1.0), edges=cross))
            self.n_m.append(EdgeGroup(self.mp[edge[0]].neg,
                                      self.mp[edge[1]].all,
                                      PSP(weight=-1.0), edges=cross))

        # Stimulus. The notebook never runs the network, so this is invented:
        #  a regular drive on a few interior mesh points
        self.stim = [NodeGroup(IN(), self.hnpm,
                               times=[self.stim_times[k]] * self.hnpm)
                     for k in range(len(self.driven))]
        self.drive = [EdgeGroup(self.stim[k], self.mp[point].pos,
                                PSP(weight=1.5),
                                edges=[(j, j) for j in range(self.hnpm)])
                      for k, point in enumerate(self.driven)]


def build_net():
    # A regular drive on a few interior points, to give the mesh something
    #  to propagate
    driven = interior[:4]
    times = [[t for t in range(0, TIMESTEPS, 3)] for _ in driven]

    net = Network()
    net.m = Mesh(mesh_points, mesh_edges, NPM,
                 neuron=pLIF(threshold=1.0, leak=0.4, prob=0.85),
                 driven=driven, stim_times=times)
    net.build()
    return net


net = build_net()
print(net)

# Run simulations using Brian 2 and SANA-FE and compare. Note that they will
#  absolutely not match 1-1 because they use random neurons (with no easy
#  way to match seeds etc.). The behaviour is quite different run to run, so
#  it's hard to really assess the match until we actually have the NeuroFEM
#  application here
results = {}
for label, Sim in (("brian", SimBrian), ("sanafe", SimSANAFE)):
    sim = Sim(build_net())
    sim.compile(**({"arch": ARCH} if Sim is SimSANAFE else {}))
    sim.run(TIMESTEPS)
    spikes = sim.get_spikes()
    results[label] = spikes
    print(f"{label:7s} total spikes: {sum(len(s) for s in spikes)}, "
          f"active neurons: {sum(1 for s in spikes if s)}/{len(spikes)}")
    if label == "sanafe":
        print(f"        energy: {sim.energy['total']:.3e} J, "
              f"latency: {sim.latency:.3e} s")

b, s = results["brian"], results["sanafe"]
nb, ns = sum(len(x) for x in b), sum(len(x) for x in s)
print(f"\nspike count ratio sanafe/brian: {ns / max(nb, 1):.3f}")
agree = sum(1 for x, y in zip(b, s) if bool(x) == bool(y))
print(f"neurons agreeing on active/silent: {agree}/{len(b)}")