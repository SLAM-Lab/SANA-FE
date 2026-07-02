======================
Models (sanafe.models)
======================

This page documents every pipeline model available in SANA-FE. It begins with
the framework attributes shared by all models, followed by a section for each
model listing its model-specific attributes.


Framework attributes
--------------------

These attributes are accepted by every pipeline model in SANA-FE. They control
shared behavior such as energy/latency reporting, hardware unit selection, and
update scheduling, and apply equally to built-in and plugin models.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``connections_out``
     - (int) Connections outgoing from a neuron (deprecated)
   * - ``dendrite_hw_name``
     - (str) Unique name of the dendrite H/W unit.
   * - ``energy_access_neuron``
     - (float) Energy cost for a soma to access a neuron (J).
   * - ``energy_message_in``
     - (float) Energy cost of receiving a spike message (J).
   * - ``energy_message_out``
     - (float) Energy cost of sending a spike message (J)
   * - ``energy_process_spike``
     - (float) Energy cost for one synapse look-up/access (J).
   * - ``energy_spike_out``
     - (float) Energy cost for a soma to spike (J).
   * - ``energy_update``
     - (float) Energy cost of updating a dendrite (s)
   * - ``energy_update_neuron``
     - (float) Energy cost for a soma to update (J).
   * - ``force_update``
     - (bool) Force updates every time-step.
   * - ``latency_access_neuron``
     - (float) Latency cost for a soma to access a neuron (s).
   * - ``latency_message_in``
     - (float) Latency cost of receiving a spike message (s).
   * - ``latency_message_out``
     - (float) Latency cost of sending a spike message (s)
   * - ``latency_process_spike``
     - (float) Latency cost for one synapse look-up/access (s).
   * - ``latency_spike_out``
     - (float) Latency cost for a soma to spike (s).
   * - ``latency_update``
     - (float) Latency cost of updating a dendrite (s)
   * - ``latency_update_neuron``
     - (float) Energy cost for a soma to update (s).
   * - ``model``
     - (str) Unique model name, either built-in or plugin.
   * - ``plugin``
     - (str) Plug-in library path.
   * - ``soma_hw_name``
     - (str) Unique name of the soma H/W unit.
   * - ``synapse_hw_name``
     - (str) Unique name of the synapse H/W unit.


accumulator
-----------

*(none)*


accumulator_with_delay
----------------------

*(none)*


current_based
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``d``
     - (float) Time-steps that a spike is delayed.
   * - ``delay``
     - (float) Time-steps that a spike is delayed.
   * - ``w``
     - (float) Synaptic weight associated with connection.
   * - ``weight``
     - (float) Synaptic weight associated with connection.



input
-----

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``poisson``
     - (float) Randomized Poisson input encoding, i.e., random > poisson: spike, else no spike.
   * - ``rate``
     - (float) Rate-based input encoding.
   * - ``spikes``
     - (list[bool]) A per-time-step spike-train.



leaky_integrate_fire
------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``bias``
     - Bias current applied every step: v[t+1] = bias + v[t]*leak_decay + u[t]
   * - ``force_potential``
     - *(no description)*
   * - ``force_update``
     - (bool) Force soma to update every step, regardless of inputs.
   * - ``force_update_every_timestep``
     - (bool) Force soma to update every step, regardless of inputs.
   * - ``leak_decay``
     - (float) Decay term applied every step: v[t+1] = bias + v[t]*leak_decay + u[t]
   * - ``log_u``
     - (bool) Record input current (u) for Loihi soma.
   * - ``noise``
     - (str) Noise source. Only file-based noise stream supported
   * - ``noise_bits``
     - (int) The number of noise bits simulated.
   * - ``refractory_delay``
     - (int) The number of refractory steps after a spike, default=0.
   * - ``reset``
     - (float) The potential to reset to after a spike. Default=0.0
   * - ``reset_mode``
     - (str) The type of reset to apply on spikes [none/soft/hard/saturate]. Default=hard
   * - ``reverse_reset``
     - (float) The potential to reset to after a reverse spike.
   * - ``reverse_reset_mode``
     - (str) The type of reset to apply on negative/reverse spikes [none/soft/hard/saturate]. Default=None
   * - ``reverse_threshold``
     - (float) Fires when V < reverse_threshold, then V <- f(reverse_reset_mode, reverse_reset).
   * - ``threshold``
     - (float) Fires when V > threshold, then V <- f(reset_mode, reset).



taps
----

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``space_constants``
     - (list[float]) List of RC line space constants, N-1 in total.
   * - ``taps``
     - (int) Number of N dendritic 'taps' or compartments.
   * - ``time_constants``
     - (list[float]) List of RC line time constants, one per tap.



truenorth
---------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``bias``
     - Additive bias current applied every step.
   * - ``leak``
     - (float) Subtractive leak term applied every step
   * - ``leak_towards_zero``
     - (bool) Leak towards zero if true, leak away otherwise.
   * - ``random_mask``
     - (int) Positive mask to apply to randomly generated noise.
   * - ``reset``
     - (float) The potential to reset to after a spike. Default=0.0
   * - ``reset_mode``
     - (str) The type of reset to apply on spikes [none/soft/hard/saturate]. Default=hard
   * - ``reverse_reset``
     - (float) The potential to reset to after a reverse spike.
   * - ``reverse_reset_mode``
     - (str) The type of reset to apply on negative/reverse spikes [none/soft/hard/saturate]. Default=None
   * - ``reverse_threshold``
     - (float) The potential at which a reverse spike is triggered.
   * - ``threshold``
     - (float) The potential at which a spike is triggered.


