"""Build a per-(core, step) dataset joining actsim timing with chip activity.

Usage: python3 build_dataset.py <sanafe_project> <arch_yaml> <actsim_csv> <out_dir>

Activity comes from a SANA-FE run of the same network (1.json) using the
NeuroScale protocol. Core IDs are converted to the paper (snake) order used by
actsim. All times are in seconds.
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

project, arch_path, actsim_csv, out_dir = sys.argv[1:5]
os.makedirs(out_dir, exist_ok=True)
os.environ.setdefault("NEUROSCALE_LOG_DIR", out_dir)
spec = importlib.util.spec_from_file_location(
    "ns", os.path.join(project, "scripts/ongoing/neuroscale.py"))
ns = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ns)

# --- Run SANA-FE with spike + message traces -------------------------------
arch = ns.make_arch_with_protocol(ns.load_arch_dict(arch_path), "neuroscale")
tiles, core_id_map = ns.compute_tile_order(arch)
net, _ = ns.build_network_from_json(arch, ns.USE_JSON_NETWORK)
chip = ns.sanafe.SpikingChip(arch)
chip.load(net)
msg_path = os.path.join(out_dir, "messages.csv")
res = chip.sim(ns.TIMESTEPS, timing_model="simple", message_trace=msg_path,
               spike_trace=True)


def paper_core(hw):
    tile, core = str(hw).split(".")
    return core_id_map[int(tile) * 4 + int(core)]


# --- Messages -> per (core, step) activity ---------------------------------
m = pd.read_csv(msg_path)
m["src"] = m["src_hw"].map(paper_core)
real = m[m["mid"] >= 0].copy()
real["dst"] = real["dest_hw"].map(paper_core)
real["remote"] = (real["src"] // 4) != (real["dst"] // 4)   # crosses tiles
real["step"] = real["timestep"]

out_g = real.groupby(["src", "step"])
act_out = pd.DataFrame({
    "msgs_out": out_g.size(),
    "events_out": out_g["spikes"].sum(),
    "hops_out": out_g["hops"].sum(),
    "remote_msgs_out": out_g["remote"].sum(),
    "dest_cores_out": out_g["dst"].nunique(),
})
in_g = real.groupby(["dst", "step"])
act_in = pd.DataFrame({
    "msgs_in": in_g.size(),
    "events_in": in_g["spikes"].sum(),
    "hops_in": in_g["hops"].sum(),
    "remote_msgs_in": in_g["remote"].sum(),
    "src_cores_in": in_g["src"].nunique(),
})
act_out.index.names = act_in.index.names = ["core", "step"]

# SANA-FE's own per-core latency terms (NeuroScale analytical model)
gen = m.groupby(["src", "timestep"])["generation_delay"].sum()
proc = (real.groupby(["dst", "timestep"])["processing_delay"].sum()
        + real.groupby(["dst", "timestep"])["network_delay"].sum())
gen.index.names = proc.index.names = ["core", "step"]
sanafe_T = m.groupby(["src", "timestep"])["sent_timestamp"].first()
sanafe_T.index.names = ["core", "step"]

# Neurons fired per core (pop p -> paper tile p-1, 50 neurons per core)
fired = {}
for step, spiking in enumerate(res["spike_trace"], start=1):
    for s in spiking:
        pop = int(s.group_name[1:])
        core = (pop - 1) * 4 + s.neuron_offset // 50
        fired[(core, step)] = fired.get((core, step), 0) + 1
fired = pd.Series(fired, name="fired")
fired.index.names = ["core", "step"]

# --- Static connectivity: core-level fan-in / fan-out ----------------------
pairs = real[["src", "dst"]].drop_duplicates()
fan_in = pairs.groupby("dst")["src"].apply(set).to_dict()
fan_out = pairs.groupby("src")["dst"].apply(set).to_dict()

# --- Assemble the grid -------------------------------------------------------
cores = range(64)
steps = range(1, ns.TIMESTEPS + 1)
idx = pd.MultiIndex.from_product([cores, steps], names=["core", "step"])
df = pd.DataFrame(index=idx)
for part in (act_out, act_in, fired.to_frame()):
    df = df.join(part)
df = df.fillna(0.0)
df["sanafe_gen"] = gen.reindex(idx).fillna(0.0)
df["sanafe_proc"] = proc.reindex(idx).fillna(0.0)
df["sanafe_T"] = sanafe_T.reindex(idx)
df = df.reset_index()
df["pop"] = df["core"] // 4 + 1
df["fan_in_cores"] = df["core"].map(lambda c: len(fan_in.get(c, ())))
df["fan_out_cores"] = df["core"].map(lambda c: len(fan_out.get(c, ())))

# --- actsim timing -------------------------------------------------------------
a = pd.read_csv(actsim_csv)
aT = a.pivot_table(index="core_id", columns="sim_step", values="wall_clock")
df.to_csv(os.path.join(out_dir, "activity.csv"), index=False)
aT.to_csv(os.path.join(out_dir, "actsim_T.csv"))
pd.Series({k: sorted(v) for k, v in fan_in.items()}).to_json(
    os.path.join(out_dir, "fan_in.json"))
pd.Series({k: sorted(v) for k, v in fan_out.items()}).to_json(
    os.path.join(out_dir, "fan_out.json"))
print(f"rows={len(df)} msgs={len(real)} fired={int(df.fired.sum())} "
      f"actsim steps={aT.columns.min()}..{aT.columns.max()}")
