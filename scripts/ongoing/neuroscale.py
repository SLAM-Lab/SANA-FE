"""
neuroscale.py
==================

Single-file recreation of Figure 4 from
  Li, Imam & Manohar, "A deterministic neuromorphic architecture with scalable
  time synchronization" (Nature Communications, 2025)
using the SANA-FE neuromorphic simulator.

The experiment runs the SAME recurrent SNN under three time-synchronization
protocols and compares how simulation time advances across cores in wall-clock
time:
  * "clock"      -> TrueNorth-style global fixed-interval barrier.
  * "barrier"    -> Loihi-style two-phase mesh barrier.
  * "neuroscale" -> distributed local synchronization (this paper).

Network (paper Methods, "Network configurations" + Fig. 4):
  * 16 recurrently-connected populations P1..P16, 200 neurons each.
  * Populations 1..16 have biases 1..16; fixed threshold 100. Higher bias =>
    denser spiking => that core runs "slower".
  * Intra-population connections: uniform random, prob 0.1.
  * Inter-population (P_i -> P_{i+1}) connections: prob 0.05.
  * All synaptic weights = 1.
  * Each population maps onto one tile (a 2x2 block of cores == 1 tile here).

Per-core timing extraction (how we recover Fig. 4b/c dots):
  * neuroscale: read the message trace. We group timings by core to get async advance.
  * barrier:    read the perf trace `sim_time` per timestep (all cores advance
    in lockstep, so one global timeline suffices).
  * clock:      analytic -- fixed period, no need to simulate for timing. We
    still run it once to confirm SANA-FE supports the mode.

Outputs:
  * results/fig4_core_timing.csv  : tidy (protocol, core_id, sim_step, wall_clock)
  * results/fig4_summary.csv      : per-protocol aggregate sim_time / spikes
  * results/fig4b_progression.png : full-range progression (Fig. 4b style)
  * results/fig4c_zoom.png        : magnified sync behavior (Fig. 4c style)

Usage:
  python neuroscal.py
"""

import copy
import io
import os
import sys
import random

import matplotlib.pyplot as plt
import pandas as pd
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
# Setup paths
try:
    import sanafe
except ImportError:
    # Not installed, fall-back to local build
    sys.path.insert(0, PROJECT_DIR)
    import sanafe

# ===========================================================================
# Configuration
# ===========================================================================
ARCH_PATH = "arch/neuroscale.yaml"
PROTOCOLS = ["clock", "barrier", "neuroscale"]
PROTOCOL_LABELS = {
    "clock": "TrueNorth (clock)",
    "barrier": "Loihi (barrier)",
    "neuroscale": "NeuroScale",
}
LOG_DIR = "runs/neuroscale"

# Network constants (from the paper's Methods section)
NUM_POPULATIONS = 16
NEURONS_PER_POP = 200
THRESHOLD = 100.0
INTRA_POP_PROB = 0.10
INTER_POP_PROB = 0.05
SYNAPSE_WEIGHT = 1.0
# SYNAPSE_WEIGHT = 0.0
# SYNAPSE_WEIGHT = 100.0

TIMESTEPS = 100
SEED = 41

CLOCK_PERIOD_FALLBACK = None  # seconds

# ===========================================================================
# Architecture: swap sync_protocol without touching the base file
# ===========================================================================
def load_arch_dict(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def make_arch_with_protocol(base_arch_dict, protocol):
    """Write a temp arch YAML with sync_protocol overridden, then load it."""
    arch_dict = copy.deepcopy(base_arch_dict)
    # Protocol lives under architecture.attributes.sync_protocol.
    arch_dict["architecture"]["attributes"]["sync_protocol"] = protocol

    os.makedirs(LOG_DIR, exist_ok=True)
    tmp_path = os.path.join(LOG_DIR, f"_arch_{protocol}.yaml")
    with open(tmp_path, "w") as f:
        yaml.safe_dump(arch_dict, f, sort_keys=False)
    return sanafe.load_arch(tmp_path)


# ===========================================================================
# Network construction (1 population -> 1 tile)
# ===========================================================================
def compute_tile_order(arch):
    """
    Reproduce the paper tile order, which snakes from the top of the chip to bottom.

    Returns:
        tiles:        list of Tile objects in snake order
                      (index == population index == paper tile index)
        core_id_map:  dict {raw_hw_core_id: paper_core_id}
    """
    # First get the order of tiles to assign
    tiles_by_xy = {(t.x, t.y): t for t in arch.tiles}
    width = max(x for x, _ in tiles_by_xy) + 1
    height = max(y for _, y in tiles_by_xy) + 1

    tiles = []
    for row in range(height):
        cols = range(width) if row % 2 == 0 else range(width - 1, -1, -1)
        for col in cols:
            key = (col, row)
            if key in tiles_by_xy:
                tiles.append(tiles_by_xy[key])

    # Now build a map of SANA-FE core id to paper core ID
    cores_per_tile = len(tiles[0].cores)
    core_id_map = {}
    for paper_tile_idx, tile in enumerate(tiles):
        for offset in range(len(tile.cores)):
            raw_id = (tile.id * cores_per_tile) + offset
            paper_id = (paper_tile_idx * cores_per_tile) + offset
            core_id_map[raw_id] = paper_id

    return tiles, core_id_map


def build_network(arch, tiles, seed=SEED):
    """
    Build the Fig. 4 recurrent SNN and map each population onto one tile.

    Returns (net, pop_groups) where pop_groups[i] is P(i+1).
    """
    rng = random.Random(seed)
    net = sanafe.Network()

    if len(tiles) < NUM_POPULATIONS:
        raise RuntimeError(
            f"Architecture has {len(tiles)} tiles but the experiment needs "
            f"{NUM_POPULATIONS} (one per population)."
        )

    # --- Create the 16 populations, each with its own bias ----------------
    pop_groups = []
    for p in range(NUM_POPULATIONS):
        bias = p + 1  # P1 -> 1, ..., P16 -> 16
        group = net.create_neuron_group(
            f"P{p + 1}",
            NEURONS_PER_POP,
            model_attributes={
                "bias": float(bias),
                "threshold": THRESHOLD,
            },
            log_spikes=True,
        )
        pop_groups.append(group)

    # --- Intra-population recurrent connections (prob 0.1) ----------------
    for group in pop_groups:
        neurons = group.neurons
        for src in neurons:
            for dst in neurons:
                if rng.random() < INTRA_POP_PROB:
                    src.connect_to_neuron(dst, {"weight": SYNAPSE_WEIGHT})

    # --- Inter-population connections P_i -> P_{i+1} (prob 0.05) -----------
    for p in range(NUM_POPULATIONS - 1):
    # for p in range(NUM_POPULATIONS):
        nxt = (p + 1) % NUM_POPULATIONS
        for src in pop_groups[p].neurons:
            for dst in pop_groups[nxt].neurons:
                if rng.random() < INTER_POP_PROB:
                    src.connect_to_neuron(dst, {"weight": SYNAPSE_WEIGHT})

    # --- Map each population onto one tile ---------------------------------
    # 1 population == 1 tile. Spread the 200 neurons across the cores of that
    # tile as evenly as possible.
    for p, group in enumerate(pop_groups):
        tile_cores = tiles[p].cores
        neurons = group.neurons
        per_core = -(-len(neurons) // len(tile_cores))  # ceil division
        for i, neuron in enumerate(neurons):
            core = tile_cores[min(i // per_core, len(tile_cores) - 1)]
            neuron.map_to_core(core)

    return net, pop_groups


def build_network_from_json(arch, json_path):
    import json
    """
    Reconstruct the exact connectivity recorded in a SANA-FE hardware trace
    (per-core axon/fanout tables), instead of re-sampling the paper's random
    connection probabilities. Guarantees bit-identical topology to whatever
    run produced `json_path`.
    """
    with open(json_path) as f:
        data = json.load(f)

    cores = data["cores"]
    core_by_xyid = {(c["x"], c["y"], c["id"]): c for c in cores}

    # bias is written straight into the trace per tile -> tells us which
    # P_i this tile is without needing to re-derive the snake order.
    tile_bias = {
        (c["x"], c["y"]): int(c["neuron_models"][0]["bias"])
        for c in cores if c["id"] == 0
    }
    pop_of_tile = {xy: bias - 1 for xy, bias in tile_bias.items()}
    n_pops = len(pop_of_tile)

    # per-tile core ids sorted, with cumulative neuron-offset base per core
    # (robust to uneven core/neuron splits, not just 4x50)
    core_base_offset = {}
    for (x, y) in pop_of_tile:
        tile_cores = sorted(
            (c for c in cores if c["x"] == x and c["y"] == y),
            key=lambda c: c["id"],
        )
        offset = 0
        for c in tile_cores:
            core_base_offset[(x, y, c["id"])] = offset
            offset += len(c["neurons"])

    net = sanafe.Network()
    pop_groups = [None] * n_pops
    for (x, y), p in pop_of_tile.items():
        pop_groups[p] = net.create_neuron_group(
            f"P{p + 1}",
            NEURONS_PER_POP,
            model_attributes={"bias": float(tile_bias[(x, y)]), "threshold": THRESHOLD},
            log_spikes=True,
        )

    # Decode every source neuron's fanout via the destination core's axon
    # table and issue the equivalent connect_to_neuron() call.
    for c in cores:
        src_pop = pop_of_tile[(c["x"], c["y"])]
        base = core_base_offset[(c["x"], c["y"], c["id"])]
        for n_idx, neuron in enumerate(c["neurons"]):
            src_neuron = pop_groups[src_pop].neurons[base + n_idx]
            for fo in neuron["fanout"]:
                remote = c["remote"][fo["remote"]]
                dst_core = core_by_xyid[(remote["x"], remote["y"], remote["id"])]
                dst_pop = pop_of_tile[(remote["x"], remote["y"])]
                dst_base = core_base_offset[(remote["x"], remote["y"], remote["id"])]
                for target in dst_core["input_axon"][fo["axon"]]:
                    dst_neuron = pop_groups[dst_pop].neurons[dst_base + target["neuron"]]
                    src_neuron.connect_to_neuron(dst_neuron, {"weight": target["fixed-wt"]})

    # Map populations onto tiles exactly as recorded in the trace.
    for (x, y), p in pop_of_tile.items():
        tile = next(t for t in arch.tiles if (t.x, t.y) == (x, y))
        for c in sorted((c for c in cores if c["x"] == x and c["y"] == y), key=lambda c: c["id"]):
            core = tile.cores[c["id"]]
            base = core_base_offset[(x, y, c["id"])]
            for i in range(len(c["neurons"])):
                pop_groups[p].neurons[base + i].map_to_core(core)

    return net, pop_groups

# ===========================================================================
# Per-core timing extraction
# ===========================================================================
def extract_neuroscale_timing(message_df, core_id_map):
    """
    NeuroScale: messages carry each timestep's
    wall-clock time in `sent_timestamp`, per core. Group by core to get the async
    per-core advance timeline.

    Returns list of (core_id, sim_step, wall_clock).
    """
    if message_df is None or message_df.empty:
        return []

    rows = []
    grouped = message_df.groupby(["src_hw", "timestep"], as_index=False).first()
    for _, r in grouped.iterrows():
        core_address = str(r["src_hw"])
        tile_id, core_offset = core_address.split(".")
        raw_core_id = (int(tile_id) * 4) + int(core_offset)
        core_id = core_id_map[raw_core_id]     # <-- translate to paper order
        step = int(r["timestep"])
        rows.append((core_id, step, float(r["sent_timestamp"])))
    return rows


def extract_barrier_timing(perf_df, n_cores):
    """
    Barrier (Loihi): all cores advance in lockstep. perf trace `sim_time` gives
    the latency for each timestep; replicate across cores.

    Returns list of (core_id, sim_step, wall_clock).
    """
    global CLOCK_PERIOD_FALLBACK # TODO: hack

    if perf_df is None or perf_df.empty:
        return []

    timestep_latencies = perf_df["sim_time"].tolist()
    print(f"Max ts latency: {max(timestep_latencies)}")
    CLOCK_PERIOD_FALLBACK = max(timestep_latencies)
    rows = []
    cumulative_time = 0.0
    for step, latency in enumerate(timestep_latencies, start=1):
        cumulative_time += float(latency)
        for core_id in range(n_cores):
            rows.append((core_id, step, cumulative_time))
    return rows


def extract_clock_timing(n_cores, n_steps):
    """
    Clock (TrueNorth): strict fixed period. Purely analytic -- every core
    advances every `period` of wall-clock time, in lockstep.

    Returns list of (core_id, sim_step, wall_clock).
    """
    rows = []
    for step in range(1, n_steps + 1):
        wc = step * CLOCK_PERIOD_FALLBACK
        for core_id in range(n_cores):
            rows.append((core_id, step, wc))
    return rows


# ===========================================================================
# Run one protocol
# ===========================================================================
USE_JSON_NETWORK = "1.json"  # set to None to fall back to random generation\
USE_JSON_NETWORK = None
def run_protocol(protocol, base_arch_dict):
    """
    Run one protocol and return (timing_rows, summary_dict).

    timing_rows: list of (core_id, sim_step, wall_clock).
    """
    print(f"\n=== Running protocol: {protocol} ===")
    arch = make_arch_with_protocol(base_arch_dict, protocol)
    tiles, core_id_map = compute_tile_order(arch)

    if USE_JSON_NETWORK:
        net, _pops = build_network_from_json(arch, USE_JSON_NETWORK)
    else:
        net, _pops = build_network(arch, tiles=tiles, seed=SEED)
    net.save(os.path.join(LOG_DIR, "snn.yaml"))

    chip = sanafe.SpikingChip(arch)
    chip.load(net)

    n_cores = len(core_id_map)

    if protocol == "clock":
        # We still run it to confirm SANA-FE supports the mode, but timing is
        # analytic (strict fixed period). Prefer the barrier-derived period so
        # panels are comparable; fall back to a constant otherwise.
        results = chip.sim(TIMESTEPS, timing_model="simple")
        timing = extract_clock_timing(n_cores, TIMESTEPS)

    elif protocol == "barrier":
        perf_path = os.path.join(LOG_DIR, "perf_barrier.csv")
        results = chip.sim(
            TIMESTEPS,
            timing_model="simple",
            # scheduler_threads=16,
            # processing_threads=1,
            perf_trace=perf_path,
            spike_trace=True
        )
        perf_df = pd.read_csv(perf_path)
        timing = extract_barrier_timing(perf_df, n_cores)
        plot_raster(results, os.path.join(LOG_DIR, "raster.png"))

    else:  # neuroscale
        message_path = os.path.join(LOG_DIR, "messages_neuroscale.csv")
        results = chip.sim(
            TIMESTEPS,
            timing_model="simple",
            message_trace=message_path,  # placeholder mid==-1 rows carry sent_timestamp
            spike_trace=True,
        )
        message_df = pd.read_csv(message_path)
        timing = extract_neuroscale_timing(message_df, core_id_map)

    summary = {
        "protocol": protocol,
        "sim_time": results.get("sim_time"),
        "spikes": results.get("spikes"),
        "timesteps_executed": results.get("timesteps_executed"),
    }
    print(f"  spikes={summary['spikes']}")
    return timing, summary


# ===========================================================================
# Plotting (Fig. 4b / 4c style)
# ===========================================================================
def select_representative_cores(df_proto, n=NUM_POPULATIONS):
    """
    One representative core per population for a clean raster (Fig. 4b/c show
    16 cores, one per population). With 1 tile == 1 population and 4 cores per
    tile, take the first core of every tile => every 4th core id.
    """
    core_ids = sorted(df_proto["core_id"].unique())
    return core_ids[::4][:n]


def plot_progression(df, global_max, out_path):
    """Full-range progression (Fig. 4b): 3 stacked panels."""
    fig, axes = plt.subplots(len(PROTOCOLS), 1, figsize=(7, 8), sharex=True)
    for ax, protocol in zip(axes, PROTOCOLS):
        dfp = df[df["protocol"] == protocol]
        reps = select_representative_cores(dfp)
        for y_index, core_id in enumerate(reps):
            pts = dfp[dfp["core_id"] == core_id]
            xs = pts["wall_clock"] / global_max
            # xs = pts["wall_clock"]
            ax.scatter(xs, [y_index] * len(xs), s=3, marker="|", color="#1f77b4")
        ax.set_ylabel("Core ID")
        ax.set_title(PROTOCOL_LABELS[protocol], loc="left", fontsize=10)
        ax.set_ylim(-1, len(reps))
    axes[-1].set_xlabel("Wall-clock time (normalized)")
    fig.suptitle("Fig. 4b: simulation-time progression across cores")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


def plot_zoom(df, global_max, out_path, x_lo=0.1, x_hi=0.15):
    """Magnified view (Fig. 4c): 3 side-by-side panels over a narrow window."""
    fig, axes = plt.subplots(1, len(PROTOCOLS), figsize=(11, 4), sharey=True)
    for ax, protocol in zip(axes, PROTOCOLS):
        dfp = df[df["protocol"] == protocol]
        reps = select_representative_cores(dfp)
        for y_index, core_id in enumerate(reps):
            pts = dfp[dfp["core_id"] == core_id]
            xs = pts["wall_clock"] / global_max
            ax.scatter(xs, [y_index] * len(xs), s=3, marker="|", color="#1f77b4")
        ax.set_xlim(x_lo, x_hi)
        ax.set_title(PROTOCOL_LABELS[protocol], fontsize=10)
        ax.set_xlabel("Wall-clock time (normalized)")
    axes[0].set_ylabel("Core ID")
    fig.suptitle("Fig. 4c: magnified synchronization behavior")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")

def plot_raster(results, out_path):
    """Raster of spikes over time. y = global neuron index (grouped by pop)."""
    spike_trace = results.get("spike_trace")
    if not spike_trace:
        print("no spike_trace; enable spike_trace=True")
        return

    xs, ys = [], []
    for t, spiking in enumerate(spike_trace):
        for mn in spiking:
            group = mn.group_name
            nid = mn.neuron_offset
            group_number = int(group.split("P")[1])
            xs.append(t)
            ys.append((group_number * NEURONS_PER_POP) + nid)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(xs, ys, s=1, marker=".", linewidths=0)
    # population band boundaries
    for p in range(1, NUM_POPULATIONS):
        ax.axhline(p * NEURONS_PER_POP - 0.5, color="0.85", lw=0.5)
    ax.set_xlabel("Simulation timestep")
    ax.set_ylabel("Neuron (grouped by population P1..P16)")
    ax.set_title("Spike Raster")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")

# ===========================================================================
# Main
# ===========================================================================
def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    base_arch_dict = load_arch_dict(ARCH_PATH)

    # Run barrier first so we can derive the clock's fixed period (paper sets
    # the TrueNorth tick to the longest Loihi barrier interval).
    order = ["barrier", "neuroscale", "clock"]
    # order = ["barrier", "neuroscale",]

    timing_rows = []       # (protocol, core_id, sim_step, wall_clock)
    summaries = []

    for protocol in order:
        rows, summary = run_protocol(protocol, base_arch_dict)
        for core_id, sim_step, wall_clock in rows:
            timing_rows.append((protocol, core_id, sim_step, wall_clock))
        summaries.append(summary)

    # --- Serialize to CSV ------------------------------
    timing_df = pd.DataFrame(
        timing_rows, columns=["protocol", "core_id", "sim_step", "wall_clock"]
    )
    timing_csv = os.path.join(LOG_DIR, "fig4_core_timing.csv")
    timing_df.to_csv(timing_csv, index=False)
    print(f"\nSaved {timing_csv}")

    summary_df = pd.DataFrame(summaries)
    summary_csv = os.path.join(LOG_DIR, "fig4_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Saved {summary_csv}")

    # --- Plots -------------------------------------------------------------
    global_max = timing_df["wall_clock"].max()
    plot_progression(timing_df, global_max,
                     os.path.join(LOG_DIR, "fig4b_progression.pdf"))
    plot_zoom(timing_df, global_max,
              os.path.join(LOG_DIR, "fig4c_zoom.pdf"))


if __name__ == "__main__":
    main()
    plt.show()
