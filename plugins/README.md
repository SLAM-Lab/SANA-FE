# Plugins

As part of SANA-FE, the user can implement different hardware models using
custom plugins. Models for synapses, dendrites and somas are all supported.
SANA-FE supports a base hardware model base class, with which it implements
all of its synaptic, dendritic and somatic hardware models. Using SANA-FE's
`PipelineUnit` base class, you can implement your own models as hardware
plugins.

## Using Plugins

There is one example already provided in the `/plugins` folder implementing a
Hodgkin-Huxley neuron model (`hodgkin_huxley.cpp`). There are a few steps
required to use plugins in SANA-FE:

1. Specify the plugin path in the architecture yaml file, in the corresponding
`synapse`, `dendrite` or `soma` hardware section. Specify the plugin path using
the attribute `plugin: <pathname>`.
2. Specify the model name using the attribute `model: <name>`.
3. Map neurons to the hardware unit as usual with the attribute: `soma_hw_name`.

For example, for the Hodgkin-Huxley example provided with SANA-FE, you could use
it as follows:

    # Rest of arch description
    ...
    soma:
    - name: plugin_example_soma
      attributes:
        plugin: plugins/hodgkin_huxley.cpp
        model: HodgkinHuxley
    ...

## Creating a New Plugin

SANA-FE can run any models provided as user plugins. The plugin must be compiled
as a shared library containing one or more hardware models.
Models can execute arbitrary code, but interfaces must be derived either from
the  general `PipelineUnit` class, or one of the specialized `SynapseUnit`,
`DendriteUnit` or `SomaUnit` base classes.

SANA-FE's plugin mechanism makes it easy to integrate plugins with your
architectural simulations. However, a few steps are needed to get plugins
running:

1. You must make sure your plugin has been built as a shared library (`.so`),
either by updating the plugin CMake file or providing your own build scripts.
2. Your new plugin must implement a hardware model class with the hardware
functionality you want. The model class you implement must be derived from
`PipelineUnit` in `chip.hpp`, which defines the required interfaces. These are
enforced by pure virtual methods, including attribute parsing methods update
methods. For examples of different hardware models, see either `models.cpp` or
the `plugins` folder.
3. Finally, provide a class factory function that returns a new instance of
your model class. This has to be in the format `create_<modelname>`. For
example, for a `HodgkinHuxley` model, we would specify the following code in
the plugin C++ file:

```
extern "C" sanafe::PipelineUnit *create_HodgkinHuxley()
{
    return (sanafe::PipelineUnit *) new HodgkinHuxley();
}
```

It is recommended new users look through the rest of the `hodgkin_huxley.cpp`
file to see what an example plugin looks like.
