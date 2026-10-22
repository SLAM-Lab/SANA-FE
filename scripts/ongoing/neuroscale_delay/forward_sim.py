"""Forward-simulate NeuroScale timing with candidate stall models; score vs actsim.

Usage: python3 forward_sim.py <data_dir> <fig_dir>

T[c,k] = L[c,k] + max(T[c,k-1], T[fanin,k-1], T[fanout,k-tadv])
'ahead' for the stall rule is computed from the simulated times, as SANA-FE's
analytical scheduler would. Activity for actsim step k is SANA-FE step k+1.
Parameters are fitted on steps 1-49 and scored on 50-99 (and overall).
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

d, fig_dir = sys.argv[1], sys.argv[2]
TADV, NC, NS = 6, 64, 100
F = pd.read_csv(os.path.join(d, "features.csv"))
fan_in = {int(k): v for k, v in json.load(open(os.path.join(d, "fan_in.json"))).items()}
fan_out = {int(k): v for k, v in json.load(open(os.path.join(d, "fan_out.json"))).items()}
aT = pd.read_csv(os.path.join(d, "actsim_T.csv"), index_col=0)
aT.columns = aT.columns.astype(int)
A = aT.sort_index().to_numpy()


def grid(col):
    g = np.zeros((NC, NS))
    p = F.pivot_table(index="core", columns="step", values=col).reindex(index=range(NC), columns=range(NS))
    g[:, :] = p.fillna(0.0).to_numpy()
    return g


mi, ev, mo = grid("msgs_in@cur"), grid("events_in@prev"), grid("msgs_out@cur")
burst = grid("tile_remote_msgs_out@cur")


def simulate(coef, stall=None):
    """coef = (fixed, per msg in, per event, per msg out); stall = fn(c, k, T) -> seconds."""
    T = np.zeros((NC, NS))
    for k in range(1, NS):
        for c in range(NC):
            start = max(T[c, k - 1],
                        max((T[f, k - 1] for f in fan_in.get(c, [])), default=0.0),
                        max((T[g, k - TADV] for g in fan_out.get(c, []) if k - TADV >= 0), default=0.0))
            L = coef[0] + coef[1] * mi[c, k] + coef[2] * ev[c, k] + coef[3] * mo[c, k]
            if stall is not None:
                L += stall(c, k, T)
            T[c, k] = start + L
    return T


def ahead_us(c, k, T):
    """Sim time by which this core's tile finished step k-1 before its destination tile."""
    tile = c // 4
    if tile == 15:
        return -1e9
    src = T[tile * 4:tile * 4 + 4, k - 1].max()
    dst = T[(tile + 1) * 4:(tile + 1) * 4 + 4, k - 1].max()
    return (dst - src) * 1e6


def score(T, steps):
    dA = np.diff(A, axis=1)[:, steps - 1]
    dT = np.diff(T, axis=1)[:, steps - 1]
    mape = np.mean(np.abs(dT - dA) / np.abs(dA)) * 100
    fin = T[:, steps[-1]]
    pop_fin = fin.reshape(16, 4).max(axis=1)
    pop_act = A[:, steps[-1]].reshape(16, 4).max(axis=1)
    return mape, (pop_fin.max() - pop_fin.min()) * 1e6, np.mean(np.abs(pop_fin - pop_act) / pop_act) * 100


train = F.step < 50
X = np.column_stack([np.ones(len(F)), F["msgs_in@cur"], F["events_in@prev"], F["msgs_out@cur"]])
lin, *_ = np.linalg.lstsq(X[train], F.L[train], rcond=None)
calm = train & (F.resid < 3e-6)
lin_calm, *_ = np.linalg.lstsq(X[calm], F.L[calm], rcond=None)
yaml = np.array([14.38e-6, 11.2e-9, 4.7e-9, 3.0e-9])


def rule(q, a, size):
    return lambda c, k, T: size if (burst[c, k] >= q and ahead_us(c, k, T) >= a) else 0.0


def expected(q, a, size, p):
    return lambda c, k, T: p * size if (burst[c, k] >= q and ahead_us(c, k, T) >= a) else 0.0


# Fit the stall rule's thresholds on the training steps (forward-simulated)
train_steps, test_steps, all_steps = np.arange(1, 50), np.arange(50, 100), np.arange(1, 100)
best = None
for q in (150, 180, 200, 220, 250):
    for a in (2.0, 3.0, 4.0, 4.5, 5.0, 6.0, 8.0):
        for size in (4e-6, 5e-6, 6e-6, 8e-6):
            T = simulate(lin_calm, rule(q, a, size))
            m = score(T, train_steps)[0]
            if best is None or m < best[0]:
                best = (m, q, a, size)
_, q, a, size = best
print(f"Fitted stall rule (train steps): burst >= {q}, ahead >= {a} us, stall {size*1e6:.0f} us")
p_hit = 0.35   # observed stall probability when the rule's conditions hold

models = {
    "SANA-FE YAML costs (current)": simulate(yaml),
    "refit costs, no stall": simulate(lin),
    "refit (calm) + stall rule": simulate(lin_calm, rule(q, a, size)),
    "refit (calm) + expected stall": simulate(lin_calm, expected(q, a, size, p_hit)),
}
rows = []
for name, T in models.items():
    m_all, spread, fin_err = score(T, all_steps)
    m_test = score(T, test_steps)[0]
    rows.append(dict(model=name, mape_all=m_all, mape_test=m_test, pop_finish_spread_us=spread,
                     pop_finish_err_pct=fin_err))
pa = A[:, 99].reshape(16, 4).max(axis=1)
rows.append(dict(model="actsim", mape_all=0, mape_test=0,
                 pop_finish_spread_us=(pa.max() - pa.min()) * 1e6, pop_finish_err_pct=0))
R = pd.DataFrame(rows)
print(R.round(2).to_string(index=False))
R.to_csv(os.path.join(d, "forward_sim_scores.csv"), index=False)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(range(1, 17), pa * 1e3, "o-", color="#222222", label="actsim")
for (name, T), col in zip(models.items(), ("#888888", "#3b7dd8", "#e07b39", "#2a9d8f")):
    ax.plot(range(1, 17), T[:, 99].reshape(16, 4).max(axis=1) * 1e3, ".-", color=col, label=name)
ax.set_xlabel("Population (bias 1-16)"); ax.set_ylabel("Finish time, step 99 (ms)")
ax.set_xticks(range(1, 17)); ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=8)
ax.set_title("Forward-simulated NeuroScale timing vs actsim")
fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "forward_sim_finish.png"), dpi=200)
