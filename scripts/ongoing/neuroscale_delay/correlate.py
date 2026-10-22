"""Correlate actsim's implied per-core latency with chip activity.

Usage: python3 correlate.py <data_dir>
"""
import os
import sys

import numpy as np
import pandas as pd

d = sys.argv[1]
act = pd.read_csv(os.path.join(d, "activity.csv"))
dec = pd.read_csv(os.path.join(d, "decomposed.csv"))
metrics = ["msgs_in", "events_in", "remote_msgs_in", "src_cores_in",
           "msgs_out", "events_out", "remote_msgs_out", "dest_cores_out", "fired"]

a = dec[dec.source == "actsim"][["core", "step", "delta", "wait", "L"]]
s = dec[dec.source == "sanafe"][["core", "step", "L"]].rename(columns={"L": "L_sanafe"})


def joined(lag, frame=a):
    """actsim step k joined with activity of step k + lag."""
    x = act.copy()
    x["step"] = x["step"] - lag
    return frame.merge(x, on=["core", "step"])


# 1. Alignment: which activity step best explains actsim step k?
print("Lag search (Pearson r of actsim L vs activity at step k+lag):")
for lag in (-2, -1, 0, 1, 2):
    j = joined(lag)
    print(f"  lag {lag:+d}: " + "  ".join(
        f"{m}={j['L'].corr(j[m]):+.3f}" for m in ("events_in", "msgs_in", "msgs_out", "fired")))

print("\nSame profile for SANA-FE's own L (known structure, step k = SANA-FE step k):")
for lag in (-1, 0, 1):
    j = joined(lag, s.rename(columns={"L_sanafe": "L"}))
    print(f"  lag {lag:+d}: " + "  ".join(
        f"{m}={j['L'].corr(j[m]):+.3f}" for m in ("events_in", "msgs_in", "msgs_out", "fired")))


def design(frame, lag_feats):
    """Activity features at the given lags, joined onto frame (core, step)."""
    out = frame.copy()
    for lag, feats in lag_feats.items():
        x = act[["core", "step"] + feats].copy()
        x["step"] = x["step"] - lag
        x = x.rename(columns={f: f"{f}@{lag:+d}" for f in feats})
        out = out.merge(x, on=["core", "step"], how="inner")
    return out


def fit(frame, target, cols, label):
    X = np.column_stack([np.ones(len(frame))] + [frame[c].to_numpy(float) for c in cols])
    y = frame[target].to_numpy(float)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    rmse = np.sqrt(((y - pred) ** 2).mean())
    print(f"\n{label}: R^2={r2:.3f}  RMSE={rmse*1e9:.0f} ns  (target sd {y.std()*1e9:.0f} ns)")
    print(f"    intercept            {coef[0]*1e6:9.3f} us")
    for c, v in zip(cols, coef[1:]):
        print(f"    {c:22s} {v*1e9:9.2f} ns/unit")
    return coef, pred


# SANA-FE (before synapse): axon_in on messages received this step, synapses
#  on events received last step, axon_out on messages sent this step
sf = design(s, {0: ["msgs_in", "msgs_out"], -1: ["events_in"]})
fit(sf, "L_sanafe", ["msgs_in@+0", "events_in@-1", "msgs_out@+0"],
    "SANA-FE L ~ known cost terms (YAML: axon_in 11.2, synapse 4.7, axon_out 3.0 ns)")

# actsim, if its step k corresponds to SANA-FE step k+1
af = design(a, {1: ["msgs_in", "msgs_out"], 0: ["events_in"]})
fit(af, "L", ["msgs_in@+1", "events_in@+0", "msgs_out@+1"],
    "actsim L ~ same terms, actsim step k = SANA-FE step k+1")
af2 = design(a, {0: ["msgs_in", "msgs_out"], -1: ["events_in"]})
fit(af2, "L", ["msgs_in@+0", "events_in@-1", "msgs_out@+0"],
    "actsim L ~ same terms, actsim step k = SANA-FE step k")
