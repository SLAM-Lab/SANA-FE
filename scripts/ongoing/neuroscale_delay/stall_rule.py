"""Fit a conditional stall term on top of the linear per-core latency model.

Usage: python3 stall_rule.py <data_dir> <fig_dir>

L = linear(core activity) + stall * 1[tile remote msgs out >= q and sender ahead
of its destination tile by >= a]. The condition uses only state a NeuroScale
scheduler already tracks (per-core timing history and per-step traffic).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d, fig_dir = sys.argv[1], sys.argv[2]
X = pd.read_csv(os.path.join(d, "features.csv"))
J = pd.read_csv(os.path.join(d, "tile_dest.csv"))[["tile", "step", "ahead_us"]]
X["tile"] = X["core"] // 4
X = X.merge(J, on=["tile", "step"], how="left")
X["ahead_us"] = X["ahead_us"].fillna(-1e9)        # last tile has no remote destination
y = X["L"].to_numpy()
core_cols = ["msgs_in@cur", "events_in@prev", "msgs_out@cur"]
burst = X["tile_remote_msgs_out@cur"].to_numpy(float)
ahead = X["ahead_us"].to_numpy(float)


def fit(extra=None):
    cols = [np.ones(len(X))] + [X[c].to_numpy(float) for c in core_cols]
    if extra is not None:
        cols.append(extra)
    M = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    pred = M @ coef
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum(), coef, pred


r2_lin, coef_lin, pred_lin = fit()
best = (r2_lin, None, None, coef_lin, pred_lin)
for q in range(50, 301, 10):
    for a in np.arange(0.0, 12.1, 0.5):
        r2, coef, pred = fit(((burst >= q) & (ahead >= a)).astype(float))
        if r2 > best[0]:
            best = (r2, q, a, coef, pred)
r2, q, a, coef, pred = best
print(f"Linear core model:             R^2={r2_lin:.3f}")
print(f"+ conditional stall term:      R^2={r2:.3f}   (burst >= {q} remote msgs, sender >= {a:.1f} us ahead)")
print(f"  intercept {coef[0]*1e6:.2f} us; " + ", ".join(
    f"{c} {v*1e9:.1f} ns" for c, v in zip(core_cols, coef[1:4])) + f"; stall {coef[4]*1e6:.2f} us")
res = y - pred
print(f"  residual sd {res.std()*1e9:.0f} ns (linear {(y-pred_lin).std()*1e9:.0f} ns), "
      f"p99 {np.quantile(res, .99)*1e9:+.0f} ns (linear {np.quantile(y-pred_lin, .99)*1e9:+.0f} ns)")
flag = (burst >= q) & (ahead >= a)
big = (y - pred_lin) > 3e-6
print(f"  rule fires on {flag.mean()*100:.1f}% of core-steps; catches {flag[big].mean()*100:.0f}% of "
      f"stalls > 3 us; precision {big[flag].mean()*100:.0f}%")

# Per-core totals: does it recover the finish-time pattern?
X["pred_lin"], X["pred_rule"] = pred_lin, pred
tot = X.groupby("core")[["L", "pred_lin", "pred_rule"]].sum() * 1e3
print("\nSum of per-core latency over the run (ms), by population (mean of 4 cores):")
pop = tot.groupby(tot.index // 4).mean()
pop.index = pop.index + 1
pop.columns = ["actsim", "linear", "linear+stall"]
print(pop.round(3).T.to_string())
print(f"spread across populations: actsim {np.ptp(pop.actsim)*1e3:.0f} us, "
      f"linear {np.ptp(pop.linear)*1e3:.0f} us, linear+stall {np.ptp(pop['linear+stall'])*1e3:.0f} us")

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(pop.index, pop.actsim, "o-", color="#222222", label="actsim")
ax.plot(pop.index, pop.linear, "s-", color="#3b7dd8", label="linear in activity")
ax.plot(pop.index, pop["linear+stall"], "^-", color="#e07b39", label="linear + conditional stall")
ax.set_xlabel("Population (bias 1-16)"); ax.set_ylabel("Total step latency over run (ms)")
ax.set_xticks(range(1, 17)); ax.grid(alpha=0.3); ax.legend(frameon=False)
ax.set_title("Per-population latency: actsim vs fitted models")
fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "stall_rule_per_population.png"), dpi=200)
