function __add_vertex!(
    g::Graphs.SimpleGraphs.AbstractSimpleGraph{T};
    metadata::Union{Tuple,NamedTuple}=NamedTuple(),
) where {T<:Integer}
    !isempty(metadata) && println(
        "Trying to add a vertex with metadata to a graph of type $(typeof(g)). Metadata $(metadata) will be ignored.",
    )
    return add_vertex!(g)
end

function _get_vertex_metadata(
    g::Graphs.SimpleGraphs.AbstractSimpleGraph{T}, vertex::T
) where {T}
    return NamedTuple()
end

function _add_edge!(
    g::Graphs.SimpleGraphs.AbstractSimpleGraph{T},
    src::T,
    dst::T;
    weight::W=nothing,
    metadata::Union{Tuple,NamedTuple}=NamedTuple(),
) where {T<:Integer,W<:Union{<:Real,Nothing}}
    (isnothing(weight) || weight == 1) ||
        @warn "Trying to add a weighted edge to an unweighted graph of type $(typeof(g)). Weight $(weight) will be ignored."
    !isempty(metadata) && @warn (
        "Trying to add an edge with metadata to a graph of type $(typeof(g)). Metadata $(metadata) will be ignored."
    )
    return add_edge!(g, src, dst)
end

function _get_edge_weight(
    g::Graphs.SimpleGraphs.AbstractSimpleGraph{T}, src::T, dst::T, weighttype::Type{U}
) where {T,U<:Real}
    return one(U)
end

# `_get_edge_metadata` must return NamedTuple() (or maybe `nothing` would be better?) instead of throwing an exception in order for `edges` to consistently work on layers and interlayers 
function _get_edge_metadata(
    g::Graphs.SimpleGraphs.AbstractSimpleGraph{T}, src::T, dst::T
) where {T}
    return NamedTuple()
end
