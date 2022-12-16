"""
    AbstractMultilayerVertex{S} <: AbstractVertex

An abstract type representing an abstract MultilayerGraph vertex.
"""
abstract type AbstractMultilayerVertex{S} <: AbstractVertex end

"""
    MultilayerVertex{N <: Integer} <: AbstractMultilayerVertex{N}

A struct representing a vertex of a MultilayerGraph.

# FIELDS
- `node::Node`: the `Node` that the `MultilayerVertex` represents;
- `layer::Union{Nothing, Symbol}`: the name of the `Layer` the `MultilayerVertex` lies in;
- `metadata::Union{<: NamedTuple, <: Tuple}`: the metadata associated to this `MultilayerVertex`.

# CONSTRUCTORS

    MultilayerVertex(node::Node, layer::Union{Nothing, Symbol},  metadata::Union{<: NamedTuple, <: Tuple})

Constructs a `MultilayerVertex` representing [`Node`](@ref) `node` in [`Layer`](@ref) with metadata `metadata`.
"""
struct MultilayerVertex{S} <: AbstractMultilayerVertex{S}
    node::Node
    layer::Union{Nothing,Symbol}
    metadata::Union{<:NamedTuple,<:Tuple}

    function MultilayerVertex(
        node::Node, layer::Union{Nothing,Symbol}, metadata::Union{<:NamedTuple,<:Tuple}
    )
        return new{layer}(node, layer, metadata)
    end
end

"""
    MultilayerVertex(node::Node, layer::Symbol)

Return `MultilayerVertex(node, layer, NamedTuple())`
"""
MultilayerVertex(node::Node, layer::Symbol) = MultilayerVertex(node, layer, NamedTuple())

"""
    MultilayerVertex(node::Node, not::Nothing)

Return `MultilayerVertex(node, nothing, NamedTuple())`.
"""
MultilayerVertex(node::Node, not::Nothing) = MultilayerVertex(node, nothing, NamedTuple())

"""
    MultilayerVertex(node::Node, metadata::Union{Tuple,NamedTuple})

Return `MultilayerVertex(node, nothing, metadata)`.
"""
function MultilayerVertex(node::Node, metadata::Union{Tuple,NamedTuple})
    return MultilayerVertex(node, nothing, metadata)
end

"""
    MultilayerVertex(node::Node)

Return `MultilayerVertex(node, nothing, NamedTuple())`.
"""
MultilayerVertex(node::Node) = MultilayerVertex(node, nothing, NamedTuple())

# Shorter alias for MultilayerVertex
"""
    MV

Alias for `MultilayerVertex`
"""
const MV = MultilayerVertex

# Base overrides
function Base.:(==)(lhs::MultilayerVertex, rhs::MultilayerVertex)
    return (lhs.node == rhs.node) &&
           (lhs.layer == rhs.layer) &&
           (lhs.metadata == rhs.metadata)
end
Base.isequal(lhs::MultilayerVertex, rhs::MultilayerVertex) = lhs == rhs

"""
    get_bare_mv(mv::MultilayerVertex)
"""
get_bare_mv(mv::MultilayerVertex) = MV(mv.node, mv.layer)

"""
    node(mv::MultilayerVertex)

Returns the Node represented by `mv`.
"""
node(mv::MultilayerVertex) = mv.node

"""
    layer(mv::MultilayerVertex)

Return the name of the layer which the `MultilayerVertex` belongs to.
"""
layer(mv::MultilayerVertex) = mv.layer

"""
    metadata(mv::MultilayerVertex)

Return the metadata associated to `mv`.
"""
metadata(mv::MultilayerVertex) = mv.metadata

# Compare multilayer vertices 
function compare_multilayervertices(
    V1::MultilayerVertex, V2::MultilayerVertex; check_metadata=false
)
    # Check if the two nodes are the same
    V1.node == V2.node || return false
    # Check if the two layers are the same
    V1.layer == V2.layer || return false
    # Check if the two metadata are the same
    if check_metadata
        V1.metadata == V2.metadata || return false
    end
    # Return true if all tests passed
    return true
end

# Console print utilities
function to_string(x::M) where {S,M<:MultilayerVertex{S}}
    return "MV($(x.node), :$(x.layer), $(x.metadata))"
end
Base.show(io::IO, x::MultilayerVertex) = print(io, to_string(x))
