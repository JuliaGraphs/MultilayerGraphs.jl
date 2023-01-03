"""
    abstract type AbstractHalfEdge{T ,U} end

An abstract type representing an HelfEdge (i.e. the edge structure to be used in fadjlists/badjlists).
"""
abstract type AbstractHalfEdge{T,U} end

"""
    mutable struct HalfEdge{T <: AbstractMultilayerVertex  , U<:Union{<:Real,Nothing}} <: AbstractHalfEdge{T,U}

# FIELDS
- `vertex::T`: the `MultilayerVertex` which is source/destination (depending on the backward/forward nature of the adjacency list) of the `HalfEdge`.
- `weight::U`: the weight of the edge.
- `metadata::Union{Tuple, NamedTuple}`: metadata associated to the edge.
"""
mutable struct HalfEdge{T<:AbstractMultilayerVertex,U<:Union{<:Real,Nothing}} <:
               AbstractHalfEdge{T,U}
    vertex::T
    weight::U
    metadata::Union{Tuple,NamedTuple}
end

"""
    vertex(he::HalfEdge)

Return the `MultilayerVertex` which is source/destination (depending on the backward/forward nature of the adjacency list) of `he`.
"""
vertex(he::HalfEdge) = he.vertex

"""
    weight(he::HalfEdge) 

Return the weight of the edge.
"""
SimpleWeightedGraphs.weight(he::HalfEdge) = he.weight
"""
    metadata(he::HalfEdge) 

Return the metadata associated to the edge.
"""
metadata(he::HalfEdge) = he.metadata

# Console print overrides
to_string(x::HalfEdge) = "HalfEdge($(to_string(vertex(x))), $(weight(x)), $(metadata(x)))"
Base.show(io::IO, x::HalfEdge) = print(io, to_string(x))
