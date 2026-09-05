# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
import sys
import itertools
import random
import warnings
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
sys.path.insert(0, PROJECT_DIR)
from sango.network import Network, NodeGroup, EdgeGroup
from sango.model import LIF, PSP, IN
from sango.backend.brian.brian import SimBrian

import sanafe
from sanafe.sango.sango import SimSANAFE

class Input(Network):
    def __init__(self, t):
        super().__init__()
        self.t = t

    def build(self):
        self.spikegen = NodeGroup(IN(), len(self.t), times=self.t)

random.seed(11)
T, fails = 20, 0
for trial in range(15):
    nin, nout = random.randint(1,3), random.randint(2,5)
    times = [sorted(random.sample(range(T), random.randint(1,5))) for _ in range(nin)]
    edges = [(i,j) for i,j in itertools.product(range(nin), range(nout))
             if random.random() < 0.7]

    if not edges:
        continue

    weights = [round(random.uniform(-1.5, 2.0), 3) for _ in edges]
    delays = [random.randint(1,6) for _ in edges]
    lif = dict(threshold=round(random.uniform(0.3, 2.5),3),
               leak=round(random.uniform(0.0, 1.0),3),
               reset=round(random.uniform(-1.0, 0.5),3),
               bias=round(random.uniform(-0.2, 0.4),3),
               voltage=round(random.uniform(-0.5, 0.5),3))

    def build():
        n = Network(); n.inp = Input(times)
        n.layer = NodeGroup(LIF(**lif), nout)
        n.dense = EdgeGroup(n.inp.spikegen, n.layer, PSP(),
                            edges=edges, weight=list(weights), delay=list(delays))
        n.build()
        return n

    b = SimBrian(build())
    b.compile()
    b.run(T)

    s = SimSANAFE(build())
    s.compile(arch=f"{PROJECT_DIR}/arch/sango.yaml")
    s.run(T)

    ok = all(sorted(x) == sorted(y) for x,y in zip(b.get_spikes(), s.get_spikes()))
    if not ok:
        fails += 1
        print("FAIL", lif, edges, weights, delays)
        print("  brian ", b.get_spikes()); print("  sanafe", s.get_spikes())

print(f"\n{15-fails}/15 random networks (with random per-edge delays) match Brian exactly")
