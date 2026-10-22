"""Unconsumed-message backlog as the stall trigger.

Usage: python3 backlog.py <data_dir> <fig_dir>

A destination consumes a sender's step-j spikes during its own step j+1. When
the sender starts step k, everything it sent for steps the destination hasn't
consumed yet is still queued ahead of the new burst:
    backlog[k] = sum(burst[j] for j in dest_done .. k-1) + burst[k]
where dest_done is the last step the destination had completed.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d, fig_dir = sys.argv[1], sys.argv[2]
aT = pd.read_csv(os.path.join(d, "actsim_T.csv"), index_col=0)
aT.columns = aT.columns.astype(int)
T = aT.groupby(aT.index // 4).max().to_numpy()
J = pd.read_csv(os.path.join(d, "tile_dest.csv"))
burst = J.pivot_table(index="tile", columns="step", values="tile_remote_msgs_out@cur").reindex(
    columns=range(0, 100)).fillna(0.0).to_numpy()


def dest_done(tile, t):
    done = np.nonzero(T[tile + 1] <= t + 1e-15)[0]
    return int(done.max()) if len(done) else -1


backlog, queued = [], []
for t, k in zip(J.tile, J.step):
    dd = max(dest_done(t, T[t, k - 1]), 0)
    prior = burst[t, dd:k].sum()                 # sent, not yet consumed
    queued.append(prior)
    backlog.append(prior + burst[t, k])
J["queued_before"] = queued
J["backlog"] = backlog


def auc(score, label):
    """Probability a random stall scores higher than a random non-stall."""
    r = pd.Series(score).rank().to_numpy()
    pos = label.astype(bool)
    n1, n0 = pos.sum(), (~pos).sum()
    return (r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


y = J["stall"].to_numpy()
print("How well each quantity separates stall from non-stall tile-steps (AUC, 0.5 = chance):")
for name, s in (("burst this step", J["tile_remote_msgs_out@cur"]),
                ("steps ahead", J["steps_ahead"]),
                ("time ahead (us)", J["ahead_us"]),
                ("queued before burst", J["queued_before"]),
                ("backlog (queued + burst)", J["backlog"]),
                ("burst x steps ahead", J["tile_remote_msgs_out@cur"] * J["steps_ahead"])):
    print(f"  {name:26s} AUC={auc(s.to_numpy(), y):.3f}")

bins = [0, 100, 200, 300, 400, 500, 600, 800, 1200, 3000]
J["backlog_bin"] = pd.cut(J["backlog"], bins, right=False)
g = J.groupby("backlog_bin", observed=True).agg(
    n=("stall", "size"), stall_pct=("stall", lambda s: s.mean() * 100),
    mean_resid_us=("resid", lambda r: r.mean() * 1e6),
    stall_size_us=("resid", lambda r: r[r > 3e-6].median() * 1e6 if (r > 3e-6).any() else np.nan))
print("\nStall rate and size by backlog (messages):")
print(g.round(2).to_string())

S = J[J.stall]
print(f"\nAmong stalls: size vs backlog r={S.resid.corr(S.backlog):+.3f}, "
      f"vs queued-before r={S.resid.corr(S.queued_before):+.3f}, vs burst "
      f"r={S.resid.corr(S['tile_remote_msgs_out@cur']):+.3f}")
J.to_csv(os.path.join(d, "tile_dest.csv"), index=False)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
ax = axes[0]
ax.scatter(J.backlog[~J.stall], J.resid[~J.stall] * 1e6, s=6, alpha=0.4, color="#3b7dd8", label="no stall")
ax.scatter(J.backlog[J.stall], J.resid[J.stall] * 1e6, s=12, alpha=0.8, color="#d1495b", label="stall > 3 us")
ax.set_xscale("symlog", linthresh=50)
ax.set_xlabel("Unconsumed remote messages queued for the destination (backlog)")
ax.set_ylabel("Unexplained step latency (us)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Stalls vs backlog")
ax = axes[1]
ax.bar(range(len(g)), g["stall_pct"], color="#d1495b")
ax.set_xticks(range(len(g)))
ax.set_xticklabels([f"{int(b.left)}-{int(b.right)}" for b in g.index], rotation=45)
ax.set_xlabel("Backlog (messages)")
ax.set_ylabel("Stall rate (%)")
ax.set_title("Stall probability by backlog")
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "stalls_vs_backlog.png"), dpi=200)
