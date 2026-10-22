"""Split per-core step time into waiting (dependencies) and implied latency.

Usage: python3 decompose.py <data_dir> [tadv]

Applies the NeuroScale rule T[c,t] = L[c,t] + max(T[c,t-1], T[fanin,t-1],
T[fanout,t-tadv]) to both actsim and SANA-FE timings, and writes
decomposed.csv with, per (core, step):
  delta  = T[c,t] - T[c,t-1]       (observed step time)
  wait   = start - T[c,t-1]        (time held back by neighbours)
  L      = T[c,t] - start          (implied per-core latency)
"""
import json
import os
import sys

import numpy as np
import pandas as pd

d = sys.argv[1]
tadv = int(sys.argv[2]) if len(sys.argv) > 2 else 6
act = pd.read_csv(os.path.join(d, "activity.csv"))
fan_in = {int(k): v for k, v in json.load(open(os.path.join(d, "fan_in.json"))).items()}
fan_out = {int(k): v for k, v in json.load(open(os.path.join(d, "fan_out.json"))).items()}

# Time matrices T[core, step], step 0 = start (time 0)
aT = pd.read_csv(os.path.join(d, "actsim_T.csv"), index_col=0)
aT.columns = aT.columns.astype(int)
A = aT.sort_index().to_numpy()                       # actsim: steps 0..99
sT = act.pivot_table(index="core", columns="step", values="sanafe_T")
S = np.zeros((64, sT.shape[1] + 1))                  # SANA-FE: steps 0..100
S[:, 1:] = sT.sort_index().to_numpy()


def decompose(T, label):
    rows = []
    n_steps = T.shape[1]
    for t in range(1, n_steps):
        for c in range(64):
            prev = T[c, t - 1]
            fi = max((T[f, t - 1] for f in fan_in.get(c, [])), default=0.0)
            fo = max((T[g, t - tadv] for g in fan_out.get(c, []) if t - tadv >= 0),
                     default=0.0)
            start = max(prev, fi, fo)
            binding = ("self" if start == prev else
                       "fan_in" if start == fi else "fan_out")
            rows.append(dict(source=label, core=c, step=t, T=T[c, t],
                             delta=T[c, t] - prev, wait=start - prev,
                             L=T[c, t] - start, binding=binding))
    return pd.DataFrame(rows)


dec = pd.concat([decompose(A, "actsim"), decompose(S, "sanafe")], ignore_index=True)
dec.to_csv(os.path.join(d, "decomposed.csv"), index=False)

# Sanity check: SANA-FE's implied L should equal its modelled N + M + sync
s = dec[dec.source == "sanafe"].merge(act[["core", "step", "sanafe_gen", "sanafe_proc"]],
                                      on=["core", "step"])
model_L = s["sanafe_gen"] + s["sanafe_proc"] + 14.0e-6
print(f"tadv={tadv}: SANA-FE implied L vs modelled N+M+sync, max abs error "
      f"{(s['L'] - model_L).abs().max() * 1e9:.3f} ns")

for src in ("actsim", "sanafe"):
    x = dec[dec.source == src]
    print(f"{src:7s} mean delta {x.delta.mean()*1e6:6.2f} us | mean L {x.L.mean()*1e6:6.2f} us "
          f"| mean wait {x.wait.mean()*1e6:5.2f} us | steps waiting {(x.wait > 1e-12).mean()*100:5.1f}% "
          f"| binding: " + ", ".join(f"{k} {v*100:.0f}%" for k, v in
                                     x.binding.value_counts(normalize=True).items()))
