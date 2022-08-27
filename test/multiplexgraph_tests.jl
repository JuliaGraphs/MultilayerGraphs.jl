# Specify layers
layers_u = [
    Layer(:layer_1, get_SimpleGraph(); U=Float64),
    Layer(:layer_2, get_SimpleWeightedGraph(); U=Float64),
    Layer(:layer_3, get_SimpleWeightedGraph(); U=Float64),
]

# Test instantiation. This also tests add_layer! and specify_interlayer!
multilayergraph = MultilayerGraph(layers_u)
