===========================================
SANA-FE Trace Data Conversion (sanafe.data)
===========================================

The ``sanafe.data`` module converts SANA-FE trace outputs into other useful
formats, e.g., Pandas.

Conversion functions accept inputs using many different types:

* a path or :class:`~pathlib.Path` to a CSV produced by ``chip.sim()``
* the dict returned by ``chip.sim()``
* the raw in-memory values
* a :class:`~pandas.DataFrame` (returned as-is)

Quickstart
----------

Build a small network, run a simulation with all traces enabled, and convert
each trace into a DataFrame:

.. code-block:: python

   import sanafe
   import sanafe.data
   import pandas

   arch, _ = sanafe.load_example()
   snn = sanafe.Network()
   group = snn.create_neuron_group(
       "in", 2, {"bias": 0.5, "threshold": 1.0, "reset": 0.0},
       log_spikes=True, log_potential=True)
   for neuron in group:
       neuron.map_to_core(arch.tiles[0].cores[0])

   chip = sanafe.SpikingChip(arch)
   chip.load(snn)
   results_dict = chip.sim(5, spike_trace=True, potential_trace=True,
                           perf_trace=True, message_trace=True,
                           neuron_trace=True)

   sanafe.data.spikes_to_raster(results_dict, n_timesteps=5)
   sanafe.data.spikes_to_dataframe(results_dict)
   sanafe.data.potentials_to_dataframe(results_dict)
   sanafe.data.performance_to_dataframe(results_dict)
   sanafe.data.messages_to_dataframe(results_dict)
   sanafe.data.neuron_traces_to_dataframe(results_dict)

.. autofunction:: spikes_to_raster

.. autofunction:: spikes_to_dataframe

.. autofunction:: potentials_to_dataframe

.. autofunction:: performance_to_dataframe

.. autofunction:: messages_to_dataframe
