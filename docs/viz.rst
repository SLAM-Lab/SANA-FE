========================================
SANA-FE Visualization (sanafe.viz)
========================================

The ``sanafe.viz`` module provides plotting helpers for SANA-FE trace outputs,
built on top of matplotlib. Each plotting function accepts the same input
forms as ``sanafe.data``:

* a path or :class:`~pathlib.Path` to a CSV produced by ``chip.sim()``
* the dict returned by ``chip.sim()``
* the raw in-memory values
* a :class:`~pandas.DataFrame` (returned as-is)

Plots share a common styling system (see :ref:`styling`) so that figures look
consistent across spike, potential, and performance views.

Quickstart
----------

Run a simulation and generate the standard set of plots:

.. code-block:: python

    import sanafe
    import sanafe.viz
    import matplotlib.pyplot as plt

    arch, snn = sanafe.load_example()
    chip = sanafe.SpikingChip(arch)
    chip.load(snn)
    results_dict = chip.sim(5, spike_trace=True, potential_trace=True,
                            perf_trace=True, message_trace=True)
    sanafe.viz.plot_raster(results_dict)
    sanafe.viz.plot_potential(results_dict)
    sanafe.viz.plot_energy(results_dict)
    sanafe.viz.plot_throughput(results_dict)
    sanafe.viz.plot_message_latency(results_dict)
    plt.show()


Spike Plots
-----------

.. autofunction:: plot_raster

.. autofunction:: raster_plot_matrix


Membrane Potential Plots
------------------------

The default rendering is a heatmap, which is consistent with the raster view
for spikes. ``plot_potential_lines`` is available for traditional per-neuron
line plots, which is more useful when only a handful of neurons are being
inspected.

.. autofunction:: plot_potential

.. autofunction:: plot_potential_lines


Performance Plots
-----------------

Hardware performance plots cover energy breakdowns, throughput, and message
latency distributions. Energy and time values are automatically rescaled to
sensible units (e.g., pJ, nJ, µs, ms) based on the magnitude of the data.

.. autofunction:: plot_energy

.. autofunction:: plot_throughput

.. autofunction:: plot_message_latency


.. _styling:

Styling
-------

All plots share a common styling system controlled by :class:`SANAFEStyle`.
A few preset styles are provided for common use cases:

* ``PUBLICATION_STYLE`` — serif fonts, 300 DPI, compact figure sizes
* ``PRESENTATION_STYLE`` — large fonts and thick lines for slides
* ``NOTEBOOK_STYLE`` — larger figures with a light grid for interactive use

Apply a preset globally with :func:`set_default_style`, or pass a
:class:`SANAFEStyle` instance to any individual plotting function via its
``style`` argument.

.. code-block:: python

   import sanafe.viz
   from sanafe.viz.styles import PUBLICATION_STYLE

   sanafe.viz.set_default_style(PUBLICATION_STYLE)
   sanafe.viz.apply_style()

   sanafe.viz.plot_raster(results_dict)

.. autoclass:: SANAFEStyle
   :members:

.. autofunction:: get_default_style

.. autofunction:: set_default_style

.. autofunction:: apply_style

.. autofunction:: get_group_colors

.. autofunction:: get_colormap

.. autofunction:: create_figure

.. autofunction:: style_axis
