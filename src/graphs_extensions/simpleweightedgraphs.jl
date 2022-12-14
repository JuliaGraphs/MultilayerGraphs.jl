Graphs.weights(g::G) where {T, G <: AbstractSimpleWeightedGraph{T}} = Graphs.weights(g)


function __add_vertex!(g::AbstractSimpleWeightedGraph{T}; metadata::Union{Tuple,NamedTuple} = NamedTuple()) where {T <: Integer}
    !isempty(metadata) && println("Trying to add a vertex with metadata to a graph of type $(typeof(g)). Metadata $(metadata) will be ignored.")
    add_vertex!(g)
end

 _get_vertex_metadata(g::AbstractSimpleWeightedGraph{T}, vertex::T) where T = NamedTuple()

function _add_edge!(g::AbstractSimpleWeightedGraph{T}, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple,NamedTuple} = NamedTuple()) where {T <: Integer, W<: Union{<: Real, Nothing}}
    !isempty(metadata) && println("Trying to add an edge with metadata to a graph of type $(typeof(g)). Metadata $metadata will be ignored.")
    isnothing(weight) ? add_edge!(g, src, dst) : add_edge!(g, src, dst, weight)
end

_get_edge_weight(g::AbstractSimpleWeightedGraph{T}, src::T, dst::T, weighttype::Type{U} ) where {T, U <: Real} = U(weights(g)[src,dst])

_get_edge_metadata(g::AbstractSimpleWeightedGraph{T}, src::T, dst::T ) where T = NamedTuple()


_set_weight!!(g::AbstractSimpleWeightedGraph{T}, src::T, dst::T, weight::U) where {T,U} = add_edge!(g, src, dst, weight)

