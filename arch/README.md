## Architecture Description

The architecture description format is also based on the YAML file format.

Different architectures are defined using a hierarchical description.
This tool models neuromorphic designs with several assumptions, in order to
simplify the tool.

1) The chip is time-step based. A time-step is a small discrete amount of time.
    This is as opposed to a purely event driven simulation e.g. ROSS.
2) The neural cores adhere to some common design patterns

At the top level, the description begins with the "architecture" keyword. Any
other top-level sections will be ignored. This defines anything at the chip
level, including the NoC interconnect.

A chip contains one or more network tiles, representing some shared network
resources e.g., a router. Each `tile` contains one or more cores, where a core
performs computation. Each neuromorphic `core` contains a fixed spike processing
hardware pipeline. It is assumed that tiles and cores are all parallel processing
elements.

Each core is assumed to have a neuromorphic pipeline which processes the updates
for one or more neurons. The pipeline is a fixed sequence of niche hardware
units. Those hardware units could contain digital logic, analog circuits or
even novel devices.

The pipeline contains the following units:

* The input axon unit receive spike packets from the network and generate
  synaptic addresses for memory lookups.

* The synaptic unit looks up connectivity for incoming spikes and updates any
  relevant synaptic currents.

* The dendritic unit combines currents based some internal structure and a set
  of operations.

* The soma unit updates membrane potentials based on the dendritic current and
  neuron model. If the firing criteria is met, it generates a spike for that
  neuron.

* The output axon unit send spikes from the soma out to the network to go to
  other cores' pipelines.

For an example, see `arch/loihi.yaml`. There are a nested series of
keywords, where keywords define required hardware units. Each block must be
contain a `name` keyword, which may optionally specify the number of instances.
Units are duplicated the number specified in the range, for example:

    # Define 8 cores, 0 through 7
    -name: neuromorphic_core[0..7]

Units must also have both an `attributes` section and the next hardware units
in the hierarchy. The attributes section will generate one or more parameters
that are passed to be parsed by the simulator and the relevant hardware models
implemented either internally (`models.cpp`) or externally (plugins).

Timestep baselines and synchronization activity costs use separate attributes:

```yaml
energy_timestep_fixed: 6.874e-07     # Chip-wide joules per timestep
latency_timestep_fixed: 1.334e-06    # Seconds per timestep
energy_sync:
  per_sending_core: 3.550e-10
  per_receiving_core: 3.550e-10
  per_message_sent: 0.0
  per_message_received: 1.709e-10
latency_sync:
  per_sending_core: 1.731e-08
  per_receiving_core: 1.731e-08
  per_message_sent: 7.317e-09
  per_message_received: 7.317e-09
```

Both fixed fields accept a scalar, list, or map indexed by mapped tile count.
Lookup uses the largest key at or below that count, or the first entry when the
count is below all keys. For example, Loihi's existing calibration becomes:

```yaml
latency_timestep_fixed: {1: 0.6e-6, 2: 1.0e-6, 4: 1.4e-6, 29: 1.8e-6}
```

No `sync_model` selector is needed. Fixed energy is charged once chip-wide per
timestep; fixed latency is added once on the critical timing path. Lists and
maps select that chip-wide value, rather than multiplying it by the tile count.
For compatibility, legacy scalar/table `latency_sync`, scalar `energy_sync`,
`energy_sync.fixed`, `energy_sync.per_sync_message`, and the previous
`latency_sync_processing` map remain readable. Explicit new fixed fields take
precedence over legacy baselines. Energy reports retain the `sync` category,
which includes both the fixed timestep energy and activity-dependent sync energy.

The `energy_sync` and `latency_sync` dictionaries have identical endpoint fields,
in joules and seconds respectively. A core pays the per-core cost for each
nonempty direction, plus the per-message cost multiplied by its logical activity
count. Missing costs default to zero; all costs must be finite and nonnegative.
The chip computes activity directly, including source-only cores without an
active synapse model. Energy sums endpoint costs; latency combines local work
through the scheduler's dependency graph.

For analytical NeuroScale timing, let `core_latency = max(neuron_work,
message_work)`. Only the sync work not hidden by this existing processing
contributes to `sync_delay`:

```text
sync_delay = clamp(message_work + sync_cost - core_latency, 0, sync_cost)
completion = latency_timestep_fixed + core_latency
           + max(previous_time, fanin_time + sync_delay, fanout_time + sync_delay)
```

The async history lag is unchanged, and zero sync costs reproduce existing
behavior. This is an aggregate approximation: control packets, arbitration,
hash collisions and shared-resource stalls are not explicitly simulated.
Detailed barrier timing treats sync work as initial receiver occupancy; external
cycle-level timing rejects nonzero sync costs.

The example coefficients fit two effective costs using silent recurrent groups.
Equal sender/receiver latency values assume symmetry; they are not four
independently identified parameters. NeuroScale's logical activity counts
represent directed core links rather than individual physical done/advance
packets. These local processing costs are independent of energy coefficients
and packet propagation delays.
