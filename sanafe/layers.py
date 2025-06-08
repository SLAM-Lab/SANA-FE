
class Layer:
    def __init__(self):
        self.group = None
        return

    def __getitem__(self, key):
        return self.group[key]

    def __len__(self):
        return len(self.group)

    def __iter__(self):
        i = 0
        while i < len(self.group):
            yield self.group[i]
            i += 1
        return


class Input2D(Layer):
    _count = 0
    def __init__(self, snn, width, height, channels=1, **kwargs):
        neuron_count = width * height
        self.width = width
        self.height = height
        self.channels = channels
        self.group = snn.create_neuron_group(f"input_{Input2D._count}",
                                             neuron_count,
                                             model_attributes=kwargs)


class Conv2D(Layer):
    _count = 0
    def __init__(self, snn, prev_layer, weights, stride_width=1,
                 stride_height=1, pad_width=0, pad_height=0, **kwargs):
        if weights.ndim != 4:
            print("Error: Expected weights kernel with 4-dimensions in the "
                  "order 'WHCN' (Width, Height, Channels, Number)")
        (kernel_width, kernel_height, filter_channels, filter_count) = \
            weights.shape
        weights = weights.flatten()

        self.width = 1 + ((prev_layer.width + (2*pad_width) - kernel_width) //
                          stride_width)
        self.height = 1 + ((prev_layer.height + (2*pad_height) - kernel_height) //
                           stride_height)
        self.channels = filter_count
        neuron_count = self.width * self.height * self.channels

        # Pass all other keyword arguments as default neuron parameters
        self.group = snn.create_neuron_group(f"conv2d_{Conv2D._count}",
                                             neuron_count,
                                             model_attributes=kwargs)

        attributes = {"w": weights}
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
    _count = 0
    def __init__(self, snn, prev_layer, neuron_count, weights, **kwargs):
        self.group = snn.create_neuron_group(f"dense_{Dense._count}",
                                             neuron_count,
                                             model_attributes=kwargs)
        attributes = {"w": weights.flatten()}
        prev_layer.group.connect_neurons_dense(self.group, attributes)

        Dense._count += 1
        return
