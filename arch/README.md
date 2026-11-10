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

Synchronization processing latency can be specified separately from the fixed
or table-based `latency_sync` setting:

```yaml
latency_sync_processing:
  per_sending_core: 1.731e-08
  per_receiving_core: 1.731e-08
  per_message_sent: 7.317e-09
  per_message_received: 7.317e-09
```

These optional, finite, nonnegative costs default to zero. They follow the chip
sync protocol's logical send/receive activity counts. A core pays the first two
costs when its corresponding direction is nonempty, plus the latter two costs
multiplied by its respective message counts. Unlike energy, these local endpoint
costs are not summed over the chip to obtain elapsed time. Source-only cores
require no active synapse model for this accounting.

Analytical NeuroScale timing overlaps this budget with neuron processing on
fan-in/fan-out dependency paths. Self progress retains the existing neuron/spike
processing cost. Zero costs reproduce the existing recurrence, and the async
history lag is unchanged. This is an aggregate approximation: control packets,
queue arbitration, shared-resource stalls and hash collisions are not explicitly
simulated. Detailed barrier timing treats the budget as initial receiver
occupancy; external cycle-level network timing rejects nonzero budgets.

The example coefficients fit two effective costs using silent, recurrent
intra-tile groups. Equal sender/receiver values assume symmetry; they are not
four independently identified coefficients. In the NeuroScale activity model,
the logical counts represent directed core links rather than separate physical
done/advance packets. Do not interpret these costs as packet propagation delays
or derive them from the energy coefficients.
