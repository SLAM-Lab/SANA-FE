Framework attributes
====================

These attributes are accepted by every pipeline model in SANA-FE. They control
shared behavior such as energy/latency reporting, hardware unit selection, and
update scheduling, and apply equally to built-in and plugin models.

For attributes specific to a particular model, see that model's documentation
page.

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

