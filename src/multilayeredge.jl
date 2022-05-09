"""
    AbstractMultilayerEdge{T} <: AbstractEdge{T}

An abstract type representing a `MultilayerGraph` edge.

It must have fields: `src`, `dst`, `weight`.
"""
abstract type  AbstractMultilayerEdge{T} <: AbstractEdge{T} end

"""
    struct MultilayerEdge{ T <: MultilayerVertex, U <: Union{ <: Real, Nothing}} <: AbstractMultilayerEdge{T}

Default concrete subtype of AbstractMultilayerEdge.

# FIELDS

- `src::T`: the source vertex of the edge ;
- `dst::T`: the destination vertex of the edge;
- `weight::U`: the edge weight.

# CONSTRUCTORS

    MultilayerEdge(src::T, dst::T, weight::U) where { T <: MultilayerVertex, U <: Union{ <: Real, Nothing}}

Default constructor.

    MultilayerEdge(src::T, dst::T) where {T <: MultilayerVertex}

Unweighted edge.
"""
struct MultilayerEdge{ T <: MultilayerVertex, U <: Union{ <: Real, Nothing}} <: AbstractMultilayerEdge{T}
    src::T
    dst::T
    weight::U
end


MultilayerEdge(src::T, dst::T) where {T <: MultilayerVertex}  = MultilayerEdge{T,Nothing}(src, dst, nothing)

"""
    src(e::AbstractMultilayerEdge)
"""
src(e::AbstractMultilayerEdge)    = e.src

"""
    dst(e::AbstractMultilayerEdge)
"""
Graphs.dst(e::AbstractMultilayerEdge)    = e.dst

"""
    weight(e::AbstractMultilayerEdge)
"""
weight(e::AbstractMultilayerEdge) = e.weight

