"""Do tile stalls depend on the receiving tile (backpressure / handshake)?

Usage: python3 receiver.py <data_dir> <fig_dir>

Each population's tile sends its remote spikes to the next population's tile.
For every tile-step, compare stalls against the destination tile's state.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d, fig_dir = sys.argv[1], sys.argv[2]
T = pd.read_csv(os.path.join(d, "tile_steps.csv"))          # per (tile, step)
aT = pd.read_csv(os.path.join(d, "actsim_T.csv"), index_col=0)
aT.columns = aT.columns.astype(int)
tileT = aT.groupby(aT.index // 4).max()                     # tile finishes when its last core does

T = T.sort_values(["tile", "step"]).reset_index(drop=True)
T["pop"] = T["tile"] + 1
dest = T[["tile", "step", "resid", "stall", "L", "tile_msgs_in@cur", "tile_msgs_out@cur",
          "tile_remote_msgs_out@cur"]].copy()
dest["tile"] -= 1                                           # dest of tile p is tile p+1
dest = dest.add_prefix("dest_").rename(columns={"dest_tile": "tile", "dest_step": "step"})
J = T.merge(dest, on=["tile", "step"], how="inner")         # drops the last tile (no dest)


def lead(row):
    """Sender's finish time of step k-1 minus the destination's (negative = sender earlier)."""
    t, k = int(row.tile), int(row.step)
    return tileT.loc[t, k - 1] - tileT.loc[t + 1, k - 1]


J["sender_lead_us"] = J.apply(lead, axis=1) * 1e6
# Positive: the sender finished step k-1 this much earlier than its destination
J["ahead_us"] = -J["sender_lead_us"]
busy = J["tile_remote_msgs_out@cur"] >= 150
print(f"High-traffic tile-steps (>=150 remote msgs out): {busy.sum()}, "
      f"stall rate {J.stall[busy].mean()*100:.0f}%")

H = J[busy]
for f, label in (("dest_tile_msgs_out@cur", "destination's own outgoing burst"),
                 ("dest_tile_msgs_in@cur", "destination's incoming load"),
                 ("dest_L", "destination's step latency"),
                 ("ahead_us", "sender ahead of destination (us)")):
    med = H[f].median()
    hi, lo = H[H[f] > med], H[H[f] <= med]
    print(f"  {label:38s} above median: {hi.stall.mean()*100:4.0f}%  "
          f"at/below: {lo.stall.mean()*100:4.0f}%  (median {med:.3g})")

print(f"\n  destination also stalls in the same step: "
      f"{H[H.stall].dest_stall.mean()*100:.0f}% of sender stalls "
      f"(vs {H[~H.stall].dest_stall.mean()*100:.0f}% when the sender doesn't)")

# All tile-steps: stall rate vs sender lead
J["lead_bin"] = pd.cut(J["ahead_us"], [-50, -5, -2, -0.5, 0.5, 2, 5, 50])
g = J.groupby("lead_bin", observed=True).agg(n=("stall", "size"), stall_pct=("stall", "mean"),
                                             busy_n=("tile_remote_msgs_out@cur", lambda x: (x >= 150).sum()))
g["stall_pct"] *= 100
print("\nStall rate by how far the sender was ahead of its destination at the previous step:")
print(g.round(1).to_string())

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.scatter(J["ahead_us"], J["resid"] * 1e6, s=8, alpha=0.5,
           c=np.where(J["tile_remote_msgs_out@cur"] >= 150, "#d1495b", "#3b7dd8"))
ax.axvline(0, color="0.5", lw=0.8)
ax.set_xlabel("Sender ahead of destination tile at previous step (us)")  # positive = sender earlier
ax.set_ylabel("Unexplained step latency (us)")
ax.set_title("Stalls vs sender/receiver timing (red: >=150 remote msgs out)")
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "stalls_vs_receiver_timing.png"), dpi=200)
J.to_csv(os.path.join(d, "tile_dest.csv"), index=False)
