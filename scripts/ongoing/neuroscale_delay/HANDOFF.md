# NeuroScale timing: SANA-FE vs actsim — handoff

## Goal

Find where SANA-FE's NeuroScale timing model diverges from actsim (the
reference NeuroScale simulator), by reading actsim's implementation. In
particular: handshakes and time-synchronization messages, network flow control
and buffering, and how cores within a tile share resources. The output should
be a list of concrete mechanisms that SANA-FE does not model, each with the
actsim code that implements it and a proposal for how SANA-FE could model it
(an equation or rule, plus its parameters).

Analysis so far is statistical: it shows *where* and *when* actsim is slower,
but actsim only exports per-step completion times, so it cannot see *why*.
Reading the actsim code, and ideally adding instrumentation, is the next step.

## Setup

SANA-FE, `neuroscale` branch, with these patches applied (`git am *.patch`):

1. Slot-by-slot pipeline execution; buffers in or before the synapse
2. NeuroScale arch fixes (`arch/neuroscale3.yaml`), `neuroscale.py` path options
3. `neuroscale3.yaml`: 16 tiles (was 32, which broke the core-ID mapping)
4. This analysis (`scripts/ongoing/neuroscale_delay/`)

Build the Python module and run:

```bash
mkdir build_py && cd build_py
cmake .. -DSTANDALONE_BUILD_ENABLED=OFF -DPYTHON_BUILD_ENABLED=ON -DCMAKE_BUILD_TYPE=Release
make -j sanafecpp && cd ..
export PYTHONPATH=$PWD/build_py

# Inputs (not in git): network and actsim reference timing
#   runs/neuroscale/1.json                 network/mapping exported from actsim's setup
#   runs/neuroscale/ticks_core_timing.csv  actsim per-core tick times (from fig4.py)

# Fig. 4 reproduction + actsim comparison
NEUROSCALE_ACTSIM_CSV=runs/neuroscale/ticks_core_timing.csv python3 scripts/ongoing/neuroscale.py

# Full delay analysis (about 5 min); writes data/ and figs/
scripts/ongoing/neuroscale_delay/run_all.sh runs/neuroscale/ticks_core_timing.csv
```

## The experiment

Li, Imam & Manohar (Nature Communications 2025), Fig. 4: 16 recurrently
connected populations of 200 neurons, biases 1–16 (higher bias = more
spiking), threshold 100, intra-population p=0.1, population i → i+1 p=0.05,
weights 1. One population per tile, 4 cores per tile, 50 neurons per core,
4×4 mesh, tiles in snake order. Remote traffic of population i therefore goes
over exactly one link, to the adjacent tile of population i+1 (1 hop).
100 time steps.

## What SANA-FE models (and doesn't)

**Per-core pipeline** (`arch/neuroscale3.yaml`, buffer before the synapse):
each step, a core

- receives spike messages: `latency_message_in` 11.2 ns per message (axon in);
- processes the *previous* step's latched synaptic events: 4.7 ns per event
  (`latency_process_spike`);
- updates all its neurons: 50 × 7.6 ns (`latency_access_neuron`);
- sends one message per (spiking neuron, destination core): 3 ns per message
  (`latency_message_out`).

**Timing recurrence** (`src/schedule.cpp`,
`schedule_messages_timestep_analytical_neuroscale`), per core c and step t:

```
L[c,t] = N[c,t] + M[c,t] + latency_sync           # send-side + receive-side work + 14 us
T[c,t] = L[c,t] + max( T[c,t-1],                  # own previous step
                       max over fan-in  f: T[f,t-1],
                       max over fan-out g: T[g,t-Tadv] )   # Tadv = max_steps_async = 6
```

Fan-in/fan-out are core-to-core connectivity sets (`Core::fan_in/fan_out`).

**Not modelled:**

- network contention, blocking, or link buffers: the analytical model ignores
  `link_buffer_size`, and hop latency is 0 in the YAML;
- any explicit handshake, token, or acknowledgement message for time
  synchronization; the dependency rule above is the whole protocol model;
- latency of synchronization messages (the rule is instantaneous);
- sharing between cores of a tile: each core is fully independent;
- message multicast: SANA-FE sends one message per destination core.

## Alignment facts (needed to compare with actsim)

- actsim `sim_step` k corresponds to SANA-FE step k+1 (0- vs 1-based). The
  activity lag profiles confirm this.
- actsim `core_id` is the paper core ID: tile in snake order × 4 + core.
  `neuroscale.py:compute_tile_order` maps SANA-FE cores to it (this needs the
  16-tile arch).
- Neuron spikes from `1.json` are assumed identical between the simulators.
  This has not been checked against actsim spike output and would be worth
  confirming.

## Findings so far

All from `run_all.sh`. "L" is a core's own step latency, after removing time
spent waiting on neighbours with the recurrence above.

1. **Waiting on neighbours explains almost nothing.** Cores wait on fan-in in
   77% of steps but only 0.06 µs on average (SANA-FE: the same). The Tadv
   fan-out limit never binds: windows of 3 or less make the fit collapse
   (R² 0.61 → ≤ 0.15). The whole gap is in L: 17.5 µs mean in actsim vs
   15.6 µs in SANA-FE.
2. **Activity explains about 60% of L.** Fitted per-unit costs vs YAML:
   fixed 15.3–15.6 µs vs 14.38 µs, 8.5–12.8 ns per message in vs 11.2,
   6.1–6.9 ns per event vs 4.7, **5.7–17 ns per message out vs 3.0** (the
   high value absorbs stalls; 5.7 ns once stalls are modelled).
3. **The remaining error is rare, large, tile-wide stalls.** About 3.6% of
   core-steps hold 46% of the excess time (median 4.9 µs, max 18 µs). In 99%
   of stalls all four cores of the tile stall together. They do not propagate
   to neighbouring tiles or repeat on the next step.
4. **Stalls need a big outgoing burst *and* the sender running ahead of its
   destination.**
   - ≥ about 200 remote messages leaving the tile in the step, and
   - the tile finished the previous step ≥ about 4.5 µs before the destination
     tile.

   Both together: 34.5% stall. Either alone: about 1%. Stall size is bounded
   by the time ahead. Stalls happen 1–3 steps ahead, never near Tadv.
5. **It is not a buffer filling across steps.** Messages queued from earlier,
   unconsumed steps have no predictive power (AUC 0.52). The current burst is
   the best single signal (AUC 0.77). There is no capacity threshold.
6. **No destination event cleanly releases a stalled sender.** Its finish time
   relative to the destination's step completions looks like that of
   non-stalled busy steps.
7. **A tile's four cores always finish within about 9.5 ns of each other**, at
   every step (stalled, busy or quiet), even when their own workloads differ by
   about 30 messages (about 0.5 µs at the fitted costs). Either actsim reports
   tile-level timing, or a tile's cores wait for each other or share a
   serialized resource. Fitting with the sum of the four cores' work vs the
   busiest core's work can't separate the two, since cores in a tile are
   similarly loaded. The fitted fixed cost (15.3 µs) is close to
   14 µs + 200 × 7.6 ns = 15.5 µs, i.e. the whole tile's neuron updates done
   one after another.
8. **Calibration alone recovers most of the accuracy.** Forward-simulating the
   recurrence with refitted costs gives 5.0% MAPE on held-out steps (YAML:
   10.7%) and population finish times within 1.4%. An explicit stall rule adds
   little, because identical step-level conditions stall only about 35% of the
   time. The remaining ~5% is per-step jitter, including about 2 µs on the
   quietest population.

## Questions for actsim, in priority order

1. **Tile vs core timing (finding 7).** How is `ticks_core_timing.csv`
   produced? Are the per-core ticks really per core, or per tile/cluster? In
   the hardware model, do the 4 cores of a tile share a controller, router
   port, neuron-update engine, or a local barrier? Is anything serialized
   across them (neuron updates, spike output, synaptic memory access)?
2. **What blocks a sender during a burst (findings 4–6).** Router/NI buffer
   depths, flow control (credit-based, handshake per flit or packet), packet
   format (one packet per spike? per destination core? multicast?), and how
   fast the receiving tile drains its input while it is busy with its own
   step. Does the receiver accept spikes for future steps at full rate, or
   only at certain phases? Is there a per-link or per-tile output queue
   between cores and the router?
3. **Time-synchronization handshakes.** What messages implement NeuroScale's
   time advance (tokens/acks per step, per fan-in/fan-out pair, per tile)?
   Their size, latency, and whether a core blocks on them. Do they share
   links/queues with spike traffic (so a burst could delay a token)? Confirm
   TIME_ADVANCE / TIME_HISTORY (the YAML comment says 7/8; SANA-FE uses
   max_steps_async=6).
4. **Fixed per-step cost.** What makes up the ~15.3 µs baseline step (vs 14 µs
   `latency_sync` + neuron updates), and why the quietest populations run
   faster (population 1: 14.9 µs mean, with ±2 µs jitter even though it gets
   no remote input).
5. **Per-message send cost** (fitted 5.7 ns, YAML 3 ns) and receive cost
   (8.5–12.8 ns vs 11.2 ns): what operations they correspond to.

## Instrumentation that would settle it

Per core (or tile) and per step, from actsim:

- time blocked sending (back-pressure), time waiting on sync tokens, and time
  computing;
- maximum input/output buffer occupancy, and per-link stall cycles;
- spike messages sent and received (to confirm activity matches SANA-FE's);
- optionally, per-message inject, arrive and consume times for one busy step.

With these, `build_dataset.py` can join actsim's breakdown onto the existing
per-(core, step) table directly.

## Files

| File | Purpose |
|---|---|
| `build_dataset.py` | Runs SANA-FE on `1.json`; per-(core, step) activity in paper core IDs |
| `decompose.py` | Splits step time into waiting vs own latency L; validates on SANA-FE |
| `correlate.py` | Step-alignment lag search; per-unit cost fits |
| `predictors.py` | Single-feature and stepwise predictors of L |
| `stalls.py` | Stall clustering (tile-wide, propagation, activity during stalls) |
| `tile_level.py` | Stall rate vs outgoing remote traffic |
| `receiver.py` | Stalls vs destination state and sender run-ahead |
| `release.py` | Steps ahead at stalls; which destination event ends a stall |
| `backlog.py` | Cross-step backlog vs current burst as the trigger |
| `stall_rule.py` | Conditional stall term on top of the linear model |
| `tile_sync.py` | Within-tile finish spread; parallel vs serialized tile model |
| `forward_sim.py` | Forward-simulates the recurrence with candidate models vs actsim |
| `run_all.sh` | Runs everything in order |
