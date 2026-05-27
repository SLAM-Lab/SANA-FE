# YAML SNN Description Format

Different mapped SNNs can be defined hierarchically using our custom YAML-based
SNN description file format. An SNN file describes a given `network`, and its
`mapping` to hardware in corresponding sections.

Under a `network`, is two sections: the neuron `groups` and `edges` section,
describing how those neurons connect. Within the neuron `groups`, we define one
or more `neurons`, e.g., that shares some common `attributes`. This is similar
to how other frameworks may define populations, or layers of similar neurons.
Attributes are simply any named parameters associated with the neuron, a user
can tag any relevant data that can be represented in YAML: strings, scalar
values, lists and mappings of these types. Attributes can be complex, e.g.,
nested structures, but ultimately limited to the types mentioned.

In each `neurons` subsection, we list all sets of neurons belonging to the
group. For conciseness we support specifying multiple neurons using a range
(..) notation. Following each neuron, give an ordered or unordered list of
attributes e.g.,

```yaml
- 0..2: [attribute1: value1]
- 3: {attribute1: value1}
```

Within the `edges` subsection, we define neuron to neuron connections or group
to group hyper-edges, including any edge attributes. The edge format uses a
notation similar to the graph DOT format e.g.,

```yaml
- layer1.0 -> layer2.1: [weight: 1]
- layer1 -> layer2: [weight: 1]
```

Note that by default, attributes are forwarded onto all relevant units.
Neurons forward their attributes onto synapse, dendrite, and soma units, whereas
edges forward their attributes on synapse and dendrite units. However, if you
want to forward an attribute onto one specific unit, you can define `synapse`,
`dendrite`, or `soma` subsections.

E.g.
```yaml
- 0..2: [soma: {threshold: 2, bias: 0.5}]
```
Or,
```
- layer1.0 -> layer2.1: [synapse: {weight: 2}]
```

In the `mappings` top-level section, we map neurons to hardware cores.
Under the section heading is a list of mappings, with an example of one mapping
as follows:

```yaml
- layer1.0..1: [core: 0.0]
```

Similar to before, neurons may provide a range. This example maps two neurons to
tile 0, core 0 (the first core in the chip).

Note as long as valid YAML syntax is used, SANA-FE does not distinguish between
different YAML styles (block style, flow style, or a mix of the two). For a
complete example of a small SNN defined using our format, see
`snn/example_snn.yaml`.

## Legacy SNN Description (Netlist) Format

Version 1 of SANA-FE (written in C) defined a simpler, less capable SNN
description format (compared to the current YAML-based format). For
back-compatability, the netlist-style format is still supported. To use this
format, use the command-line flag (-n).

In the netlist format each line defines a new entry, which may either be a
neuron group (g), neuron (n), edge (e), or mapping (&). Each line starts with
the type of entry followed by one required field and then any number of named
attributes. Fields are separated by one or more spaces.

Attributes are defined using the syntax: `<attribute>=<value>`. Note, there
is no space before or after the equals. The attribute `soma_hw_name` is required
to be set for every neuron or neuron group.

A neuron group helps reduce the number of repeated, shared parameters for a
population of neurons e.g., for a layer of neurons in a deep SNN.

`g <number of neurons> <common attributes>`

Neurons are addressed (using the group number followed by neuron number), and
then all attributes are specified. Note the group must be defined first.

`n group_id.neuron_id <unique attributes>`

An edge connects one source neuron (presynaptic) to one destination neuron
(postsynaptic). The edge may also have attributes such as synaptic weight.

`e src_group_id.src_neuron_id->dest_group_id.dest_neuron_id <edge attributes>`

Finally, mappings place predefined neurons on a hardware core. Here we specify
the neuron and the core.

`& group_id.neuron_id@tile_id.core_id`

An example of how to use the netlist format is given in `snn/example.net`
