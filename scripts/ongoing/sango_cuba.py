# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""example_brian.ipynb modified to run on the SANA-FE backend.

This notebook simply creates some interesting/complex connectivity and then
allows neurons to desynchronize and spike randomly. I don't think it actually
does anything useful.

Add a synapse channel attribute to distinguish synapse types in the dendrite.
This notebook requires custom dendrite and soma models that are application
specific.

TODO: separate those models into plugins, move the architecture YAML file into
 the script dir to keep it separate from SANA-FE's supported architectures.
"""
import sys
import warnings
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
sys.path.insert(0, PROJECT_DIR)

from sango import Network, NodeGroup, EdgeGroup, NodeList
from sango.model import Neuron, PSP, shared
from sango.backend.brian.brian import SimBrian
from sanafe.sango.sango import SimSANAFE
from brian2 import ms

ARCH = f"{PROJECT_DIR}/arch/sango_cuba.yaml"
N = 4000
P_CONNECT = 0.02
N_EXC = (N * 4) // 5
DT = 0.1  # ms per time-step
TIMESTEPS = 10000  # 1000 ms


# ## models, as the notebook defines them ##
@dataclass
class CUBA(Neuron):
    model: str = 'CUBA'
    v: float = 0.0
    ge: float = 0.0
    gi: float = 0.0
    Vt: float = shared(10.0)
    Vr: float = shared(0.0)
    El: float = shared(5.0)
    taum: float = shared(20.0)
    taue: float = shared(5.0)
    taui: float = shared(10.0)
    tref: float = shared(5.0)
    dt: float = shared(DT)          # added: SANA-FE needs the tick duration


@dataclass
class PSPe(PSP):
    model: str = 'PSPe'


@dataclass
class PSPi(PSP):
    model: str = 'PSPi'


brian_registry = {
    'CUBA': {'graph_type': 'neuron',
             'model_eqs': '''
                          dv/dt  = (ge+gi-(v-El))/taum : 1 (unless refractory)
                          dge/dt = -ge/taue : 1
                          dgi/dt = -gi/taui : 1
                          ''',
             'method': 'exact',
             'threshold': 'v>Vt',
             'reset': 'v=Vr',
             'refractory': 'tref',
             'events': {},
             'state': {'v':  {'dsl': 'v',  'default': -60.0},
                       'ge': {'dsl': 'ge', 'default': 0.0},
                       'gi': {'dsl': 'gi', 'default': 0.0}},
             'param': {'Vt':   {'dsl': 'Vt',   'default': -50.0},
                       'Vr':   {'dsl': 'Vr',   'default': -60.0},
                       'El':   {'dsl': 'El',   'default': -49.0},
                       'taum': {'dsl': 'taum', 'default': 20.0, 'unit': 'ms'},
                       'taue': {'dsl': 'taue', 'default': 5.0,  'unit': 'ms'},
                       'taui': {'dsl': 'taui', 'default': 10.0, 'unit': 'ms'},
                       'tref': {'dsl': 'tref', 'default': 5.0,  'unit': 'ms'}}},
    'PSPe': {'graph_type': 'synapse',
             'model_eqs': 'weight : 1',
             'on_pre': 'ge+=weight',
             'state': {'delay':  {'dsl': 'delay',  'default': 1.0, 'unit': 'ms'},
                       'weight': {'dsl': 'weight', 'default': 1.0}}},
    'PSPi': {'graph_type': 'synapse',
             'model_eqs': 'weight : 1',
             'on_pre': 'gi+=weight',
             'state': {'delay':  {'dsl': 'delay',  'default': 1.0, 'unit': 'ms'},
                       'weight': {'dsl': 'weight', 'default': 1.0}}},
}

# The SANA-FE equivalent. The neuron becomes a dendrite/soma pair, and the two
#  synapse types become one synapse model plus a per-edge channel flag
sanafe_registry = {
    'CUBA': {'graph_type': 'neuron',
             'hw_type': 'soma',
             'model_type': 'sango_cuba_soma',
             # CUBA spans two stages, so name both units
             'hw_name': {'soma': 'cuba_somas', 'dendrite': 'cuba_dendrites'},
             # taue/taui/dt are read by the dendrite, Vt/Vr/El/taum/tref/dt by
             #  the soma. Group attributes reach both units
             'param': {'Vt':   {'dsl': 'Vt',   'default': 10.0},
                       'Vr':   {'dsl': 'Vr',   'default': 0.0},
                       'El':   {'dsl': 'El',   'default': 5.0},
                       'taum': {'dsl': 'taum', 'default': 20.0},
                       'taue': {'dsl': 'taue', 'default': 5.0},
                       'taui': {'dsl': 'taui', 'default': 10.0},
                       'tref': {'dsl': 'tref', 'default': 5.0},
                       'dt':   {'dsl': 'dt',   'default': DT}},
             'state': {'v': {'dsl': 'v', 'default': 0.0}}},
    'PSPe': {'graph_type': 'synapse',
             'hw_type': 'synapse',
             'model_type': 'current_based',
             'param': {},
             'state': {'weight':  {'dsl': 'weight', 'default': 1.0},
                       'delay':   {'dsl': 'delay',  'default': 1.0,
                                   'rep': 'tick'},
                       'channel': {'dsl': None,     'default': 0}}},
    'PSPi': {'graph_type': 'synapse',
             'hw_type': 'synapse',
             'model_type': 'current_based',
             'param': {},
             'state': {'weight':  {'dsl': 'weight', 'default': 1.0},
                       'delay':   {'dsl': 'delay',  'default': 1.0,
                                   'rep': 'tick'},
                       'channel': {'dsl': None,     'default': 1}}},
}


def rand_edges(s1, s2, prob, rng):
    index = rng.choice(range(s1 * s2), size=int(s1 * s2 * prob), replace=False)
    return [(int(idx) // s2, int(idx) % s2) for idx in index]


Vt, Vr, El = -50.0, -60.0, -49.0
we = 60 * 0.27 / 10
wi = -20 * 4.5 / 10


def build_net(seed=1):
    rng = np.random.default_rng(seed)
    net = Network()
    net.p = NodeGroup(CUBA(Vt=Vt, Vr=Vr, El=El), N)
    net.p.v = Vr + rng.uniform(size=(N,)) * (Vt - Vr)
    net.pe = NodeList(net.p[:N_EXC])
    net.pi = NodeList(net.p[N_EXC:])
    net.pe_p = EdgeGroup(net.pe, net.p, PSPe(weight=we),
                         edges=rand_edges(N_EXC, N, P_CONNECT, rng))
    net.pi_p = EdgeGroup(net.pi, net.p, PSPi(weight=wi),
                         edges=rand_edges(N - N_EXC, N, P_CONNECT, rng))
    net.build()
    return net


print(f"CUBA: {N} neurons ({N_EXC} exc), p={P_CONNECT}, "
      f"{TIMESTEPS} steps of {DT} ms")

results = {}
for label, Sim, registry in (("brian", SimBrian, brian_registry),
                             ("sanafe", SimSANAFE, sanafe_registry)):
    # Run each backend in a loop
    print(f"Running sim using {label} backend")
    sim = Sim(build_net())
    sim.model_registry.update(registry)

    if Sim is SimBrian:
        sim.tstep = DT * ms  # Scale the spike time back to continuous (s)
        sim.compile()
    else:
        # TODO: maybe timestep/DT timing should be stored in the backend and applied to the parsed spike data there?
        sim.compile(arch=ARCH, dt=DT)
    sim.run(TIMESTEPS)

    spike_list = sim.get_spikes()
    results[label] = spike_list
    total = sum(len(s) for s in spike_list)
    rate = total / (N * TIMESTEPS * DT * 1e-3)
    print(f"{label:7s} spikes: {total:7d}  mean rate: {rate:6.2f} Hz  "
          f"active: {sum(1 for s in spike_list if s)}/{N}")

    if Sim is SimSANAFE:
        print(f"        energy: {sim.energy['total']:.3e} J  "
              f"latency: {sim.latency:.3e} s")

    # Plot spike list as raster
    spike_table = [[],[]] # t, n
    for n, spk in enumerate(spike_list):
        for t in spk:
            spike_table[0].append(t)
            spike_table[1].append(n)
    plt.figure()
    plt.plot(spike_table[0], spike_table[1], ',k')
    plt.xlabel('Time (ms)')
    plt.ylabel('Neuron index')
    print(len(spike_table[0]))
    plt.savefig(f"{PROJECT_DIR}/runs/sango/{label}.png")

b, s = results["brian"], results["sanafe"]
nb, ns = sum(len(x) for x in b), sum(len(x) for x in s)
print(f"\nspike count ratio comparing sanafe/brian: {ns / max(nb, 1):.3f}")
