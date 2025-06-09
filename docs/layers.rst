=====================================
SANA-FE Layers Module (sanafe.layers)
=====================================

The ``sanafe.layers`` module provides higher-level machine learning abstractions
that simplify the construction of deep spiking neural network applications in
SANA-FE. These layers wrap the low-level C++ kernel with PyTorch/Keras-style
interfaces.

Quick Start
===========

.. code-block:: python

    import sanafe
    from sanafe import layers
    import numpy as np

    # Create network
    net = sanafecpp.Network()

    # Build a simple CNN
    input_layer = layers.Input2D(net, 28, 28, 1)

    # Add convolutional layer
    conv_weights = np.random.random((3, 3, 1, 32))
    conv1 = layers.Conv2D(net, input_layer, conv_weights, threshold=1.0)

    # Add dense layer
    dense_weights = np.random.random((conv1.width * conv1.height * conv1.channels, 10))
    output = layers.Dense(net, conv1, 10, dense_weights)

For a more in-depth real-world example, see our tutorial implementing
`DVS gesture categorization <https://github.com/SLAM-Lab/SANA-FE/blob/main/tutorial/tutorial_5_dvs.ipynb>`_

Layer Base Class
================

.. autoclass:: sanafe.layers.Layer
    :members:
    :undoc-members:
    :show-inheritance:

Input Layers
============

.. autoclass:: sanafe.layers.Input2D
    :members:
    :undoc-members:
    :show-inheritance:

Convolutional Layers
====================

.. autoclass:: sanafe.layers.Conv2D
    :members:
    :undoc-members:
    :show-inheritance:

Dense Layers
============

.. autoclass:: sanafe.layers.Dense
    :members:
    :undoc-members:
    :show-inheritance:

Application Examples
=====================

MNIST CNN
---------

.. code-block:: python

    import sanafecpp
    from sanafe import layers
    import numpy as np

    def create_mnist_cnn(net):
        """Create a simple CNN for MNIST classification."""

        # Input layer: 28x28 grayscale images
        input_layer = layers.Input2D(net, 28, 28, 1, threshold=1.0)

        # First conv layer: 3x3 kernels, 32 filters
        conv1_weights = np.random.normal(0, 0.1, (3, 3, 1, 32))
        conv1 = layers.Conv2D(net, input_layer, conv1_weights,
                             stride_width=1, stride_height=1,
                             threshold=1.2, leak=0.05)

        # Second conv layer: 3x3 kernels, 64 filters
        conv2_weights = np.random.normal(0, 0.1, (3, 3, 32, 64))
        conv2 = layers.Conv2D(net, conv1, conv2_weights,
                             stride_width=2, stride_height=2,
                             threshold=1.1, leak=0.05)

        # Flatten and classify
        flatten_size = conv2.width * conv2.height * conv2.channels
        dense_weights = np.random.normal(0, 0.1, (flatten_size, 10))
        output = layers.Dense(net, conv2, 10, dense_weights, threshold=2.0)

        return input_layer, conv1, conv2, output

Multi-Layer Perceptron
----------------------

.. code-block:: python

    def create_mlp(net, input_size, hidden_sizes, output_size):
        """Create a multi-layer perceptron."""

        # Input layer
        input_layer = layers.Input2D(net, input_size, 1, 1, threshold=0.5)

        prev_layer = input_layer
        layers_list = [input_layer]

        # Hidden layers
        for i, hidden_size in enumerate(hidden_sizes):
            weights = np.random.normal(0, 0.1, (len(prev_layer), hidden_size))
            hidden = layers.Dense(net, prev_layer, hidden_size, weights,
                                 threshold=1.0, leak=0.1)
            layers_list.append(hidden)
            prev_layer = hidden

        # Output layer
        output_weights = np.random.normal(0, 0.1, (len(prev_layer), output_size))
        output = layers.Dense(net, prev_layer, output_size, output_weights,
                             threshold=1.5)
        layers_list.append(output)

        return layers_list

Error Handling
==============

The layers module includes error checking:

- **Dimension validation**: Ensures weight matrices match layer sizes
- **Parameter validation**: Checks for positive dimensions and valid strides
- **Channel compatibility**: Verifies input/output channel consistency
- **Output size validation**: Prevents zero or negative output dimensions

Common error scenarios:

.. code-block:: python

    # This will raise ValueError: weight matrix size mismatch
    try:
        wrong_weights = np.random.random((100, 50))  # Wrong input size
        layer = layers.Dense(net, prev_layer, 50, wrong_weights)
    except ValueError as e:
        print(f"Error: {e}")

    # This will raise ValueError: invalid stride
    try:
        layer = layers.Conv2D(net, input_layer, weights, stride_width=0)
    except ValueError as e:
        print(f"Error: {e}")

See Also
========

* :doc:`api` - API reference