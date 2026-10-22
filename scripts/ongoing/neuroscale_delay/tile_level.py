"""Tile-level view: does stall size follow traffic leaving the tile?

Usage: python3 tile_level.py <data_dir> <fig_dir>
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d, fig_dir = sys.argv[1], sys.argv[2]
os.makedirs(fig_dir, exist_ok=True)
X = pd.read_csv(os.path.join(d, "features.csv"))
X["tile"] = X["core"] // 4
cols = {"L": "mean", "resid": "mean", "tile_remote_msgs_out@cur": "first",
        "tile_remote_msgs_in@cur": "first", "tile_msgs_out@cur": "first",
        "tile_msgs_in@cur": "first", "events_in@cur": "sum", "msgs_out@cur": "sum"}
T = X.groupby(["tile", "step"]).agg(cols).reset_index()
T["tile_local_msgs_out"] = T["tile_msgs_out@cur"] - T["tile_remote_msgs_out@cur"]
T["stall"] = T["resid"] > 3e-6

print("Correlation of tile-mean residual with tile traffic (Pearson / Spearman):")
for f in ("tile_remote_msgs_out@cur", "tile_local_msgs_out", "tile_remote_msgs_in@cur",
          "tile_msgs_in@cur", "events_in@cur"):
    print(f"  {f:28s} {T['resid'].corr(T[f]):+.3f} / {T['resid'].corr(T[f], method='spearman'):+.3f}")

# Binned: stall probability and mean residual vs remote messages out
bins = [0, 25, 50, 100, 150, 200, 250, 300, 400, 800]
T["bin"] = pd.cut(T["tile_remote_msgs_out@cur"], bins, right=False)
g = T.groupby("bin", observed=True).agg(n=("resid", "size"), p_stall=("stall", "mean"),
                                        mean_resid_us=("resid", lambda r: r.mean() * 1e6),
                                        p90_resid_us=("resid", lambda r: r.quantile(0.9) * 1e6))
print("\nBy remote messages leaving the tile this step:")
print(g.round(2).to_string())

# Does size scale with the burst once stalled?
S = T[T.stall]
slope = np.polyfit(S["tile_remote_msgs_out@cur"], S["resid"], 1)
print(f"\nAmong stalled tile-steps: residual ~ {slope[0]*1e9:.1f} ns per outgoing remote msg "
      f"+ {slope[1]*1e6:.2f} us, r={S['resid'].corr(S['tile_remote_msgs_out@cur']):+.3f}")

# Figure: residual vs outgoing remote traffic, coloured by stall
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
ax.scatter(T["tile_remote_msgs_out@cur"], T["resid"] * 1e6, s=8, alpha=0.5,
           c=np.where(T.stall, "#d1495b", "#3b7dd8"))
ax.axhline(3, color="0.5", lw=0.8, ls="--")
ax.set_xlabel("Remote messages leaving the tile in the step")
ax.set_ylabel("Unexplained step latency (us)")
ax.set_title("Tile-mean residual vs outgoing remote traffic")
ax = axes[1]
ax.bar(range(len(g)), g["p_stall"] * 100, color="#d1495b")
ax.set_xticks(range(len(g)))
ax.set_xticklabels([f"{int(b.left)}-{int(b.right)}" for b in g.index], rotation=45)
ax.set_xlabel("Remote messages leaving the tile")
ax.set_ylabel("Steps with a stall > 3 us (%)")
ax.set_title("Stall probability by outgoing traffic")
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "stalls_vs_outgoing_traffic.png"), dpi=200)
T.drop(columns="bin").to_csv(os.path.join(d, "tile_steps.csv"), index=False)
