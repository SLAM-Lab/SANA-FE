"""Characterize large unexplained stalls in actsim per-core latency.

Usage: python3 stalls.py <data_dir> [threshold_us]
"""
import json
import os
import sys

import numpy as np
import pandas as pd

d = sys.argv[1]
thr = float(sys.argv[2]) * 1e-6 if len(sys.argv) > 2 else 3.0e-6
X = pd.read_csv(os.path.join(d, "features.csv"))
fan_in = {int(k): v for k, v in json.load(open(os.path.join(d, "fan_in.json"))).items()}
fan_out = {int(k): v for k, v in json.load(open(os.path.join(d, "fan_out.json"))).items()}
X["pop"] = X["core"] // 4 + 1
stall = X["resid"] > thr
S = X[stall]
print(f"Stalls (residual > {thr*1e6:.1f} us): {len(S)} of {len(X)} core-steps "
      f"({len(S)/len(X)*100:.1f}%), carrying {S.resid.sum()/X.resid[X.resid > 0].sum()*100:.0f}% "
      f"of all positive residual time")
print(f"  stall size: median {S.resid.median()*1e6:.2f} us, max {S.resid.max()*1e6:.2f} us")

# Where: by population and by step
print("  stalls per population: " + " ".join(
    f"P{p}:{n}" for p, n in S.groupby("pop").size().reindex(range(1, 17), fill_value=0).items()))
early = (S.step <= 10).sum()
print(f"  stalls in steps 1-10: {early}, steps 11-99: {len(S) - early}")

# Clustering in space/time
key = set(zip(S.core, S.step))


def frac(pred):
    return np.mean([pred(c, t) for c, t in key]) * 100


same_tile = frac(lambda c, t: any((o, t) in key for o in range(c // 4 * 4, c // 4 * 4 + 4) if o != c))
upstream = frac(lambda c, t: any((f, t) in key or (f, t - 1) in key for f in fan_in.get(c, []) if f // 4 != c // 4))
downstream = frac(lambda c, t: any((g, t) in key or (g, t + 1) in key for g in fan_out.get(c, []) if g // 4 != c // 4))
own_next = frac(lambda c, t: (c, t + 1) in key)
# Baseline: same tests with random positions at the same stall rate
rng = np.random.default_rng(0)
rand = set(zip(rng.integers(0, 64, len(key)), rng.integers(1, 100, len(key))))
base_tile = np.mean([any((o, t) in rand for o in range(c // 4 * 4, c // 4 * 4 + 4) if o != c)
                     for c, t in rand]) * 100
print(f"  another core on the same tile also stalls that step: {same_tile:.0f}% "
      f"(random baseline {base_tile:.0f}%)")
print(f"  an upstream (other-tile fan-in) core stalls at t or t-1: {upstream:.0f}%")
print(f"  a downstream (other-tile fan-out) core stalls at t or t+1: {downstream:.0f}%")
print(f"  the same core stalls again next step: {own_next:.0f}%")

# What the chip is doing during stalls vs normal steps
cols = ["events_in@cur", "msgs_in@cur", "msgs_out@cur", "remote_msgs_in@cur",
        "remote_msgs_out@cur", "tile_remote_msgs_in@cur", "tile_remote_msgs_out@cur",
        "dest_cores_out@cur", "src_cores_in@cur", "wait"]
cmp = pd.DataFrame({"stall": S[cols].mean(), "normal": X[~stall][cols].mean()})
cmp["ratio"] = cmp["stall"] / cmp["normal"]
print("\nMean activity, stall vs normal steps:")
print(cmp.round(3).to_string())

print("\nLargest 12 stalls:")
show = ["core", "pop", "step", "L", "resid", "msgs_in@cur", "events_in@cur", "msgs_out@cur",
        "tile_remote_msgs_out@cur", "tile_remote_msgs_in@cur"]
t = S.sort_values("resid", ascending=False)[show].head(12).copy()
t["L"] *= 1e6
t["resid"] *= 1e6
print(t.round(2).to_string(index=False))
