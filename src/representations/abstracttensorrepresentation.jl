"""
    AbstractTensorRepresentation{U}

An abstract type encoding a generic tensorial representation of the links and metadata of a multilayer graph. 

Concrete subtypes must have an `array` field (a 4-dimensional tensor of eltype U, indexes as [source_node_idx, destination_node_idx, source_layer_idx, destination_layer_idx ]).

# PARAMETRIC TYPES

- `U`: the weight type of the multilayer graph.
"""
abstract type AbstractTensorRepresentation end

"""
    getindex(atr::AbstractTensorRepresentation, src_mv::MultilayerVertex, dst_mv::MultilayerVertex)
"""
function Base.getindex(
    atr::AbstractTensorRepresentation, src_mv::MultilayerVertex, dst_mv::MultilayerVertex
)
    # Find the index of the source node in the array of node
    src_n_idx = atr.idx_N_associations(src_mv.node)
    # Find the index of the destination node in the array of node
    dst_n_idx = atr.idx_N_associations(dst_mv.node)
    # Find the index of the source layer in the array of layers
    src_layer_idx = findfirst(name -> name == src_mv.layer, atr.layers_names)
    # Find the index of the destination layer in the array of layers
    dst_layer_idx = findfirst(name -> name == dst_mv.layer, atr.layers_names)
    # Return the value in the array at the position of the source node, destination node, source layer, and destination layer
    return atr.array[src_n_idx, dst_n_idx, src_layer_idx, dst_layer_idx]
end

"""
    getindex(atr::AbstractTensorRepresentation, src_tup::Tuple{String, Symbol}, dst_tup::Tuple{String, Symbol})
"""
function Base.getindex(
    atr::AbstractTensorRepresentation,
    src_tup::Tuple{String,Symbol},
    dst_tup::Tuple{String,Symbol},
)
    # Find the index of the source node in the adjacency tensor
    src_n_idx = atr.idx_N_associations(Node(src_tup[1]))
    # Find the index of the destination node in the adjacency tensor
    dst_n_idx = atr.idx_N_associations(Node(dst_tup[1]))
    # Find the index of the source layer in the adjacency tensor
    src_layer_idx = findfirst(name -> name == src_tup[2], atr.layers_names)
    # Find the index of the destination layer in the adjacency tensor
    dst_layer_idx = findfirst(name -> name == dst_tup[2], atr.layers_names)
    # Return the value of the adjacency tensor at the given indices
    return atr.array[src_n_idx, dst_n_idx, src_layer_idx, dst_layer_idx]
end

"""
    array(atr::AbstractTensorRepresentation) 

Return the array representation of `atr`.
"""
array(atr::AbstractTensorRepresentation) = atr.array
