"""
    AbstractMultilayerEdge{T} <: AbstractEdge{T}

An abstract type representing a `MultilayerGraph` edge.

It must have fields: `src`, `dst`, `weight`.
"""
abstract type AbstractMultilayerEdge{U} <: AbstractEdge{AbstractMultilayerVertex} end

"""
    struct MultilayerEdge{ T <: MultilayerVertex, U <: Union{ <: Real, Nothing}} <: AbstractMultilayerEdge{T}

Default concrete subtype of AbstractMultilayerEdge.

# FIELDS

- `src::T`: the source vertex of the edge;
- `dst::T`: the destination vertex of the edge;
- `weight::U`: the edge weight.

# CONSTRUCTORS

    MultilayerEdge(src::T, dst::T, weight::U) where { T <: MultilayerVertex, U <: Union{ <: Real, Nothing}}

Default constructor.
"""
struct MultilayerEdge{U<:Union{<:Real,Nothing}} <: AbstractMultilayerEdge{U}
    src::AbstractMultilayerVertex
    dst::AbstractMultilayerVertex
    weight::U
    metadata::Union{NamedTuple,Tuple}
end

"""
    MultilayerEdge(src::AbstractMultilayerVertex, dst::AbstractMultilayerVertex) 

Convert to `MultilayerEdge{Nothing}(src, dst, nothing, NamedTuple())`.
"""
function MultilayerEdge(src::AbstractMultilayerVertex, dst::AbstractMultilayerVertex)
    return MultilayerEdge{Nothing}(src, dst, nothing, NamedTuple())
end

"""
    MultilayerEdge(src::AbstractMultilayerVertex, dst::AbstractMultilayerVertex, weight::U) where {U <: Real}

Convert to `MultilayerEdge{U}(src, dst, weight, NamedTuple())`.
"""
function MultilayerEdge(
    src::AbstractMultilayerVertex, dst::AbstractMultilayerVertex, weight::U
) where {U<:Real}
    return MultilayerEdge{U}(src, dst, weight, NamedTuple())
end

"""
    MultilayerEdge(src::AbstractMultilayerVertex, dst::AbstractMultilayerVertex, metadata::NamedTuple)

Convert to `MultilayerEdge{Nothing}(src, dst, nothing, metadata)`.
"""
function MultilayerEdge(
    src::AbstractMultilayerVertex, dst::AbstractMultilayerVertex, metadata::Union{Tuple,NamedTuple}
)
    return MultilayerEdge{Nothing}(src, dst, nothing, metadata)
end

"""
    ME

Shorter alias for MultilayerEdge.
"""
const ME = MultilayerEdge

"""
    src(e::AbstractMultilayerEdge)

Return the source `MultilayerVertex` of `e`.
"""
Graphs.src(e::AbstractMultilayerEdge) = e.src

"""
    dst(e::AbstractMultilayerEdge)

Return the destination `MultilayerVertex` of `e`.
"""
Graphs.dst(e::AbstractMultilayerEdge) = e.dst

"""
    weight(e::AbstractMultilayerEdge)

Return the weight of `e`.
"""
SimpleWeightedGraphs.weight(e::AbstractMultilayerEdge) = e.weight

"""
    metadata(e::AbstractMultilayerEdge)

Return the metadata of `e`.
"""
metadata(e::AbstractMultilayerEdge) = e.metadata

"""
    reverse(e::MultilayerEdge)

Return and edge between `dst(e)` and `src(e)` with same `weight(e)` and `metadata(e)`.
"""
Base.reverse(e::MultilayerEdge) = MultilayerEdge(dst(e), src(e), weight(e), metadata(e))

# Compare multilayer edges
function compare_multilayeredges(
    lhs::MultilayerEdge,
    rhs::MultilayerEdge;
    check_weight::Bool=false,
    check_metadata::Bool=false,
)
    # Check source
    lhs.src != rhs.src && return false
    # check destination
    lhs.dst != rhs.dst && return false
    # Check weight
    _check_weight = false
    if check_weight
        _check_weight = lhs.weight == rhs.weight ? true : return false
    end
    # Check metadata
    _check_metadata = false
    if check_metadata
        _check_metadata = lhs.metadata == rhs.metadata ? true : return false
    else
        _check_metadata = true
    end
    return true
end

# Base overrides
function Base.:(==)(lhs::MultilayerEdge, rhs::MultilayerEdge)
    return (lhs.src == rhs.src) &&
           (lhs.dst == rhs.dst) &&
           (lhs.weight == rhs.weight) &&
           (lhs.metadata == rhs.metadata)
end # 

function Base.isequal(lhs::E1, rhs::E2) where {E1<:MultilayerEdge,E2<:MultilayerEdge}
    # println("here")
    return lhs == rhs
end

# Console print utilities
function to_string(x::MultilayerEdge)
    return "ME($(to_string(src(x))) --> $(to_string(dst(x))),\tweight = $(weight(x)),\tmetadata = $(metadata(x)))"
end
Base.show(io::IO, x::MultilayerEdge) = print(io, to_string(x))
