using Revise
using MultilayerGraphs

# The objects that `MultilayerVertex`s represent are `Node`s
const n_nodes = 100
const nodes = [Node("node_$i") for i in 1:n_nodes]

# We next define some empty layers (i.e. layers without edges). These layers will be the building blocks The nodes that each layer represents are randomly sampled for each layer
const simple_layer = layer_simplegraph(
                                        :simple_layer,                # Name of the layer
                                        rand(nodes, rand(1:n_nodes)), # The random sample of nodes that this layer will represent
                                        0                             # Number of (randomly sampled) edges. Here it is set to 0 to get an empty layer.
)

