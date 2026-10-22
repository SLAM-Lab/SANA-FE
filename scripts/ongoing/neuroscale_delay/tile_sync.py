"""Do a tile's cores behave as one timing unit in actsim?

Usage: python3 tile_sync.py <data_dir>

1. Spread of the 4 cores' finish times per tile-step (stalled, busy, quiet).
2. Tile step latency fitted on the busiest core's work (cores in parallel)
   vs the sum of all 4 cores' work (cores serialized).
"""
import os
import sys

import numpy as np
import pandas as pd

d = sys.argv[1]
aT = pd.read_csv(os.path.join(d, "actsim_T.csv"), index_col=0)
aT.columns = aT.columns.astype(int)
A = aT.sort_index().to_numpy()
F = pd.read_csv(os.path.join(d, "features.csv"))
F["tile"] = F.core // 4
cols = ["msgs_in@cur", "events_in@prev", "msgs_out@cur"]

g = F.groupby(["tile", "step"])
st = g.agg(stall=("resid", lambda r: (r > 3e-6).any()),
           work_diff=("msgs_out@cur", lambda x: x.max() - x.min()),
           burst=("tile_remote_msgs_out@cur", "first")).reset_index()
st["finish_spread_ns"] = [(A[t*4:t*4+4, k].max() - A[t*4:t*4+4, k].min()) * 1e9
                          for t, k in zip(st.tile, st.step)]
busy = st.burst >= 150
print("Spread of the 4 cores' finish times within a tile (actsim):")
for label, m in (("stall steps", st.stall), ("busy, no stall", ~st.stall & busy),
                 ("quiet steps", ~st.stall & ~busy)):
    s = st[m]
    print(f"  {label:15s} n={m.sum():4d}  median {s.finish_spread_ns.median():7.1f} ns, "
          f"max {s.finish_spread_ns.max():9.1f} ns; cores' msgs_out differ by median {s.work_diff.median():.0f}")

y = g["L"].mean().to_numpy()
train = g["step"].first().to_numpy() < 50


def fit(Xf, label):
    X = np.column_stack([np.ones(len(Xf))] + [Xf[c].to_numpy(float) for c in cols])
    c, *_ = np.linalg.lstsq(X[train], y[train], rcond=None)
    p = X @ c

    def r2(m):
        return 1 - ((y[m] - p[m]) ** 2).sum() / ((y[m] - y[m].mean()) ** 2).sum()
    print(f"  {label:32s} train R^2 {r2(train):.3f}, test R^2 {r2(~train):.3f} | fixed {c[0]*1e6:.2f} us, "
          + ", ".join(f"{k.split('@')[0]} {v*1e9:.1f} ns" for k, v in zip(cols, c[1:])))


print("\nTile step latency, fitted on steps 1-49, tested on 50-99:")
fit(g[cols].max(), "busiest core's work (parallel)")
fit(g[cols].sum(), "sum over 4 cores (serialized)")
