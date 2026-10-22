"""What releases a stalled core: run-ahead limit (steps) or drained buffers?

Usage: python3 release.py <data_dir> <fig_dir>

Uses actsim times only (tile = max over its 4 cores). Each tile's remote
destination is the next population's tile.
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
T = aT.groupby(aT.index // 4).max().to_numpy()      # T[tile, step], steps 0..99
J = pd.read_csv(os.path.join(d, "tile_dest.csv"))    # per (tile, step) with stall, resid
J = J[J.tile < 15].copy()
n_steps = T.shape[1]


def steps_completed(tile, t):
    """Last step the tile had finished by wall-clock time t."""
    done = np.nonzero(T[tile] <= t + 1e-15)[0]
    return int(done.max()) if len(done) else -1


# 1. Steps ahead of the destination when the sender starts step k
J["steps_ahead"] = [(k - 1) - steps_completed(t + 1, T[t, k - 1]) for t, k in zip(J.tile, J.step)]
J["busy"] = J["tile_remote_msgs_out@cur"] >= 150
print("Steps the sender is ahead of its destination when it starts the step:")
g = J.groupby("steps_ahead").agg(n=("stall", "size"), stall_pct=("stall", lambda s: s.mean() * 100),
                                 busy_n=("busy", "sum"),
                                 busy_stall_pct=("stall", lambda s: s[J.loc[s.index, "busy"]].mean() * 100))
print(g.round(1).to_string())

# 2. Release event: for stalled steps, how long after each destination step
#    completion does the sender finish? A consistent small gap = release event.
S = J[J.stall].copy()
print(f"\nStalled tile-steps: {len(S)}")
rows = []
for off in range(-3, 2):
    gaps = []
    for t, k in zip(S.tile, S.step):
        j = k + off
        if 0 <= j < n_steps:
            gaps.append(T[t, k] - T[t + 1, j])
    gaps = np.array(gaps) * 1e6
    rows.append(dict(dest_step=f"k{off:+d}", median_gap_us=np.median(gaps),
                     iqr_us=np.subtract(*np.percentile(gaps, [75, 25])),
                     share_positive=(gaps > 0).mean() * 100))
print("Sender finish time minus destination finishing step j (stalled steps):")
print(pd.DataFrame(rows).round(2).to_string(index=False))

# Control: the same statistic on non-stalled busy steps
C = J[~J.stall & J.busy]
rows = []
for off in range(-3, 2):
    gaps = np.array([T[t, k] - T[t + 1, k + off] for t, k in zip(C.tile, C.step)
                     if 0 <= k + off < n_steps]) * 1e6
    rows.append(dict(dest_step=f"k{off:+d}", median_gap_us=np.median(gaps),
                     iqr_us=np.subtract(*np.percentile(gaps, [75, 25]))))
print("Same, busy steps without a stall (control):")
print(pd.DataFrame(rows).round(2).to_string(index=False))

# 3. Within the same steps-ahead, does burst size still matter?
print("\nStall rate by burst size, among steps 0 and 1 ahead:")
J["burst_bin"] = pd.cut(J["tile_remote_msgs_out@cur"], [0, 50, 150, 200, 300, 800], right=False)
for sa in (0, 1):
    sub = J[J.steps_ahead == sa]
    r = sub.groupby("burst_bin", observed=True)["stall"].agg(["size", "mean"])
    print(f"  {sa} ahead: " + ", ".join(f"{b}: {m*100:.0f}% (n={n})" for b, (n, m) in r.iterrows()))
J.to_csv(os.path.join(d, "tile_dest.csv"), index=False)

# Figure: stall end vs destination completing step k-1
fig, ax = plt.subplots(figsize=(6.5, 4))
gap = np.array([T[t, k] - T[t + 1, k - 1] for t, k in zip(J.tile, J.step)]) * 1e6
ax.scatter(gap[~J.stall], J.resid[~J.stall] * 1e6, s=6, alpha=0.4, color="#3b7dd8", label="no stall")
ax.scatter(gap[J.stall], J.resid[J.stall] * 1e6, s=12, alpha=0.8, color="#d1495b", label="stall > 3 us")
ax.axvline(0, color="0.5", lw=0.8)
ax.set_xlabel("Sender finishes step k minus destination finishes step k-1 (us)")
ax.set_ylabel("Unexplained step latency (us)")
ax.set_title("When do stalled tiles finish, relative to their destination?")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, "stall_release_timing.png"), dpi=200)
