"""
SANA-FE Neural Network Layers Module

High-level layer abstractions for building spiking neural networks with
automatic neuron group creation and connectivity patterns.

This module provides PyTorch/Keras-style layer interfaces that wrap the
low-level SANA-FE C++ kernel for easier network construction.

Example:
    >>> import sanafe
    >>> from sanafe import layers
    >>>
    >>> # Create network and architecture
    >>> net = sanafecpp.Network()
    >>>
    >>> # Build CNN with high-level layers
    >>> input_layer = layers.Input2D(net, 28, 28, 1)
    >>> conv1 = layers.Conv2D(net, input_layer, conv_weights, stride_width=1)
    >>> dense1 = layers.Dense(net, conv1, 128, dense_weights)
"""

class Layer:
    """
    Base class for all neural network layers.

    Provides common functionality for layer objects including indexing,
    iteration, and length operations that delegate to the underlying
    neuron group.

    Attributes:
        group (NeuronGroup): Underlying SANA-FE neuron group

    Example:
        >>> # Access individual neurons
        >>> neuron = layer[0]
        >>>
        >>> # Iterate over neurons
        >>> for neuron in layer:
        ...     neuron.set_attributes(log_spikes=True)
        >>>
        >>> # Get layer size
        >>> print(f"Layer has {len(layer)} neurons")
    """

    def __init__(self):
        """Initialize base layer with empty neuron group."""
        self.group = None
        return

    def __getitem__(self, key):
        """
        Access neuron(s) by index or slice.

        Args:
            key (int or slice): Index or slice for neuron access

        Returns:
            Neuron or list: Single neuron or list of neurons
        """
        return self.group[key]

    def __len__(self):
        """
        Get number of neurons in this layer.

        Returns:
            int: Number of neurons in the layer
        """
        return len(self.group)

    def __iter__(self):
        """
        Iterate over all neurons in the layer.

        Yields:
            Neuron: Individual neurons in the layer
        """
        i = 0
        while i < len(self.group):
            yield self.group[i]
            i += 1
        return


class Input2D(Layer):
    """
    2D input layer for image-like data with configurable dimensions.

    Creates a flattened neuron group representing a 2D input space with
    optional channel dimension. Neurons are arranged in row-major order
    (channels-last format).

    Attributes:
        width (int): Input width in pixels
        height (int): Input height in pixels
        channels (int): Number of input channels
        group (NeuronGroup): Underlying neuron group

    Example:
        >>> # Create 28x28 grayscale input layer
        >>> input_layer = Input2D(net, 28, 28, 1, threshold=1.0)
        >>>
        >>> # Create 32x32 RGB input layer
        >>> rgb_input = Input2D(net, 32, 32, 3, leak=0.01)
        >>>
        >>> # Access specific neuron at position (y=10, x=15, channel=0)
        >>> neuron_idx = 10 * input_layer.width + 15
        >>> neuron = input_layer[neuron_idx]
    """

    _count = 0

    def __init__(self, snn, width, height, channels=1, **kwargs):
        """
        Create 2D input layer with specified dimensions.

        Args:
            snn (Network): SANA-FE network to add this layer to
            width (int): Input width in pixels
            height (int): Input height in pixels
            channels (int, optional): Number of input channels. Defaults to 1.
            **kwargs: Additional neuron model parameters (threshold, leak, etc.)

        Raises:
            ValueError: If width, height, or channels are non-positive
        """
        super().__init__()

        if width <= 0 or height <= 0 or channels <= 0:
            raise ValueError("Width, height, and channels must be positive")

        neuron_count = width * height * channels
        self.width = width
        self.height = height
        self.channels = channels

        self.group = snn.create_neuron_group(
            f"input_{Input2D._count}",
            neuron_count,
            model_attributes=kwargs
        )

        Input2D._count += 1


class Conv2D(Layer):
    """
    2D convolutional layer implementing CNN-style feature extraction.

    Applies learnable filters across the input to produce feature maps.
    Automatically calculates output dimensions based on input size, kernel
    size, stride, and padding parameters.

    Attributes:
        width (int): Output feature map width
        height (int): Output feature map height
        channels (int): Number of output channels (filters)
        group (NeuronGroup): Underlying neuron group

    Example:
        >>> # Create 3x3 convolution with 32 filters (3*3*1*32 = 288 weights)
        >>> conv_weights = [0.1 if i % 2 == 0 else -0.1 for i in range(288)]
        >>> conv_layer = Conv2D(net, input_layer, conv_weights,
        ...                     stride_width=1, stride_height=1)
        >>>
        >>> # Create 5x5 convolution with stride 2 (5*5*32*64 = 51200 weights)
        >>> conv_weights = [0.05] * 51200
        >>> conv2 = Conv2D(net, conv_layer, conv_weights,
        ...                 stride_width=2, stride_height=2,
        ...                 threshold=1.2, leak=0.05)
    """

    _count = 0

    def __init__(self, snn, prev_layer, weights, stride_width=1,
                 stride_height=1, pad_width=0, pad_height=0, **kwargs):
        """
        Create 2D convolutional layer.

        Args:
            snn (Network): SANA-FE network to add this layer to
            prev_layer (Layer): Previous layer to connect from
            weights (np.ndarray): 4D weight tensor with shape (W, H, C_in, C_out)
                where W=kernel width, H=kernel height, C_in=input channels,
                C_out=output channels
            stride_width (int, optional): Horizontal stride. Defaults to 1.
            stride_height (int, optional): Vertical stride. Defaults to 1.
            pad_width (int, optional): Horizontal padding. Defaults to 0.
            pad_height (int, optional): Vertical padding. Defaults to 0.
            **kwargs: Additional neuron model parameters

        Raises:
            ValueError: If weights don't have expected 4D shape
            ValueError: If stride or padding values are invalid
        """
        super().__init__()

        if weights.ndim != 4:
            raise ValueError(
                "Expected weights kernel with 4 dimensions in the "
                "order 'WHCN' (Width, Height, Channels_in, Channels_out)"
            )

        if stride_width <= 0 or stride_height <= 0:
            raise ValueError("Stride values must be positive")

        if pad_width < 0 or pad_height < 0:
            raise ValueError("Padding values cannot be negative")

        (kernel_width, kernel_height, filter_channels, filter_count) = weights.shape

        # Validate input channel compatibility
        if hasattr(prev_layer, 'channels') and prev_layer.channels != filter_channels:
            raise ValueError(
                f"Input channels mismatch: prev_layer has {prev_layer.channels} "
                f"channels but weights expect {filter_channels}"
            )

        weights_flat = weights.flatten()

        # Calculate output dimensions using standard convolution formula
        self.width = 1 + ((prev_layer.width + (2 * pad_width) - kernel_width) //
                          stride_width)
        self.height = 1 + ((prev_layer.height + (2 * pad_height) - kernel_height) //
                           stride_height)
        self.channels = filter_count

        if self.width <= 0 or self.height <= 0:
            raise ValueError(
                f"Invalid output dimensions ({self.width}x{self.height}). "
                "Check kernel size, stride, and padding parameters."
            )

        neuron_count = self.width * self.height * self.channels

        # Create neuron group with user-specified parameters
        self.group = snn.create_neuron_group(
            f"conv2d_{Conv2D._count}",
            neuron_count,
            model_attributes=kwargs
        )

        # Connect to previous layer with convolutional connectivity
        attributes = {"w": weights_flat}
        prev_layer.group.connect_neurons_conv2d(
            self.group,
            attributes,
            prev_layer.width,
            prev_layer.height,
            prev_layer.channels,
            kernel_width,
            kernel_height,
            filter_count,
            stride_width,
            stride_height
        )

        Conv2D._count += 1


class Dense(Layer):
    """
    Fully-connected (dense) layer with all-to-all connectivity.

    Connects every neuron in the previous layer to every neuron in this
    layer, implementing a standard neural network dense layer.

    Attributes:
        group (NeuronGroup): Underlying neuron group

    Example:
        >>> # Create dense layer with 128 neurons
        >>> prev_size = len(prev_layer)
        >>> weights = [0.1 if i % 3 == 0 else 0.05 for i in range(prev_size * 128)]
        >>> dense_layer = Dense(net, prev_layer, 128, weights,
        ...                     threshold=1.0, leak=0.1)
        >>>
        >>> # Create output layer with 10 classes
        >>> output_weights = [0.08] * (128 * 10)
        >>> output_layer = Dense(net, dense_layer, 10, output_weights)
    """

    _count = 0

    def __init__(self, snn, prev_layer, neuron_count, weights, **kwargs):
        """
        Create fully-connected layer.

        Args:
            snn (Network): SANA-FE network to add this layer to
            prev_layer (Layer): Previous layer to connect from
            neuron_count (int): Number of neurons in this layer
            weights (np.ndarray): 2D weight matrix with shape
                (prev_layer_size, neuron_count)
            **kwargs: Additional neuron model parameters

        Raises:
            ValueError: If weight matrix dimensions don't match layer sizes
            ValueError: If neuron_count is non-positive
        """
        super().__init__()

        if neuron_count <= 0:
            raise ValueError("Neuron count must be positive")

        expected_shape = (len(prev_layer), neuron_count)
        if weights.shape != expected_shape:
            raise ValueError(
                f"Weight matrix shape {weights.shape} doesn't match expected "
                f"shape {expected_shape} for connection from {len(prev_layer)} "
                f"to {neuron_count} neurons"
            )

        self.group = snn.create_neuron_group(
            f"dense_{Dense._count}",
            neuron_count,
            model_attributes=kwargs
        )

        # Connect with dense connectivity pattern
        attributes = {"w": weights.flatten()}
        prev_layer.group.connect_neurons_dense(self.group, attributes)

        Dense._count += 1