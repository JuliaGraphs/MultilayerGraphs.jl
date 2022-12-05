"""
    AbstractMultilayerVertex{S}

An abstract type representing an abstract MultilayerGraph vertex.
"""
abstract type AbstractMultilayerVertex{S} <: AbstractVertex end #<: AbstractVertex{S} end

"""
    MultilayerVertex{N <: Integer} <: AbstractMultilayerVertex{N}

A struct representing a vertex of a MultilayerGraph.

# FIELDS
- `node::Node`: the `Node` that the `MultilayerVertex` represents;
- `layer::Union{Nothing, Symbol}`: the name of the `Layer` the `MultilayerVertex` lies in;
- `layer::Symbol`: the layer the `MultilayerVertex` belongs to.

# CONSTRUCTORS

    MultilayerVertex(node::Node, layer::Union{Nothing, Symbol},  metadata::Union{<: NamedTuple, <: Tuple})

Constructs a `MultilayerVertex` representing [`Node`](@ref) `node` in [`Layer`](@ref) with metadata `metadata`.
"""
struct MultilayerVertex{S} <: AbstractMultilayerVertex{S}
    node::Node
    layer::Union{Nothing, Symbol}
    metadata::Union{<: NamedTuple, <: Tuple}

    MultilayerVertex(node::Node, layer::Union{Nothing, Symbol},  metadata::Union{<: NamedTuple, <: Tuple}) = new{layer}(node, layer, metadata)
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
MultilayerVertex(node::Node, metadata::Union{Tuple,NamedTuple}) = MultilayerVertex(node, nothing, metadata)

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
Base.:(==)(lhs::MultilayerVertex, rhs::MultilayerVertex) = (lhs.node == rhs.node) && (lhs.layer == rhs.layer) && (lhs.metadata == rhs.metadata)
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

function compare_multilayervertices(V1::MultilayerVertex, V2::MultilayerVertex; check_metadata = false)
    V1.node == V2.node || return false
    V1.layer == V2.layer || return false

    if check_metadata
        V1.metadata == V2.metadata || return false
    end
    
    return true
end

# Console print utilities
to_string( x::M) where {S  , M <: MultilayerVertex{S}} = "MV($(x.node), :$(x.layer), $(x.metadata))" #{:$S}
Base.show(io::IO, x::MultilayerVertex) = print(io, to_string(x))