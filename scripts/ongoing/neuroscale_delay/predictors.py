"""Search for predictors of actsim's per-core step latency.

Usage: python3 predictors.py <data_dir>

Alignment: actsim step k <-> SANA-FE step k+1. Features are taken from the
SANA-FE step matching the actsim step ("@cur") and the one before it ("@prev").
"""
import os
import sys

import numpy as np
import pandas as pd

d = sys.argv[1]
act = pd.read_csv(os.path.join(d, "activity.csv"))
dec = pd.read_csv(os.path.join(d, "decomposed.csv"))
a = dec[dec.source == "actsim"][["core", "step", "L", "wait", "delta"]].copy()
a["sf_step"] = a["step"] + 1                      # matching SANA-FE step

# Derived activity features (per core, per SANA-FE step)
act["local_msgs_out"] = act["msgs_out"] - act["remote_msgs_out"]
act["local_msgs_in"] = act["msgs_in"] - act["remote_msgs_in"]
act["tile"] = act["core"] // 4
tile = act.groupby(["tile", "step"])[["msgs_in", "msgs_out", "events_in", "remote_msgs_in",
                                      "remote_msgs_out", "fired"]].sum().add_prefix("tile_")
act = act.join(tile, on=["tile", "step"])
# Traffic of the other three cores on the same tile (shared router / NoC port)
for f in ("msgs_in", "msgs_out", "remote_msgs_in", "remote_msgs_out"):
    act[f"others_{f}"] = act[f"tile_{f}"] - act[f]

base = ["msgs_in", "msgs_out", "events_in", "fired", "local_msgs_in", "remote_msgs_in",
        "local_msgs_out", "remote_msgs_out", "src_cores_in", "dest_cores_out",
        "tile_msgs_in", "tile_msgs_out", "tile_remote_msgs_in", "tile_remote_msgs_out",
        "others_msgs_in", "others_msgs_out", "others_remote_msgs_in", "others_remote_msgs_out"]
cur = act[["core", "step"] + base].rename(columns={f: f"{f}@cur" for f in base})
prev = act[["core", "step"] + base].copy()
prev["step"] += 1
prev = prev.rename(columns={f: f"{f}@prev" for f in base})
X = a.merge(cur, left_on=["core", "sf_step"], right_on=["core", "step"], suffixes=("", "_c"))
X = X.merge(prev, left_on=["core", "sf_step"], right_on=["core", "step"], suffixes=("", "_p"))
feats = [c for c in X.columns if "@" in c]
X.to_csv(os.path.join(d, "features.csv"), index=False)


def ols(cols, target="L"):
    M = np.column_stack([np.ones(len(X))] + [X[c].to_numpy(float) for c in cols])
    y = X[target].to_numpy(float)
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    res = y - M @ coef
    return coef, 1 - (res ** 2).sum() / ((y - y.mean()) ** 2).sum(), res


print("Single-feature R^2 for actsim L (top 12):")
single = sorted(((ols([f])[1], f) for f in feats), reverse=True)[:12]
for r2, f in single:
    print(f"  {f:30s} R^2={r2:.3f}")

# Forward stepwise selection
chosen, best = [], 0.0
print("\nForward selection:")
for _ in range(8):
    cand = [(ols(chosen + [f])[1], f) for f in feats if f not in chosen]
    r2, f = max(cand)
    if r2 - best < 0.005:
        break
    chosen.append(f)
    best = r2
    print(f"  + {f:30s} R^2={r2:.3f}")
coef, r2, res = ols(chosen)
print(f"\nSelected model, R^2={r2:.3f}, residual sd {res.std()*1e9:.0f} ns:")
print(f"  intercept {coef[0]*1e6:.3f} us")
for c, v in zip(chosen, coef[1:]):
    print(f"  {c:30s} {v*1e9:8.2f} ns/unit")
X["resid"] = res
X.to_csv(os.path.join(d, "features.csv"), index=False)

# What's left: structure of the residual
print("\nResidual structure:")
lag1 = X.sort_values(["core", "step"]).groupby("core")["resid"].apply(
    lambda r: r.autocorr(1)).mean()
print(f"  mean lag-1 autocorrelation within a core: {lag1:+.3f}")
by_core = X.groupby("core")["resid"].mean()
print(f"  per-core mean residual: sd {by_core.std()*1e9:.0f} ns "
      f"(vs residual sd {res.std()*1e9:.0f} ns)")
by_step = X.groupby("step")["resid"].mean()
print(f"  per-step mean residual (common to all cores): sd {by_step.std()*1e9:.0f} ns")
X["pos"] = X["core"] % 4
print("  mean residual by core position in tile: " +
      ", ".join(f"{p}:{v*1e9:+.0f} ns" for p, v in X.groupby("pos")["resid"].mean().items()))
q = X["resid"].quantile([0.5, 0.9, 0.99])
print(f"  residual quantiles: median {q[0.5]*1e9:+.0f}, p90 {q[0.9]*1e9:+.0f}, p99 {q[0.99]*1e9:+.0f} ns")
