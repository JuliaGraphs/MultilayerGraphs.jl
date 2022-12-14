#= # Some custom extensions of Graphs.jl and SimpleWeightedGraphs.jl
"""
    SimpleWeightedGraph(n_vertices::Integer, n_edges::Integer; T::Type = Int64, U::Type = Float64)

Random  SimpleWeightedGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedGraph(
    n_vertices::Integer, n_edges::Integer; T::Type=Int64, U::Type=Float64
)
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        adjm[cart_idx] = U(rand())
    end

    adjm .= (adjm .+ adjm') ./ U(2)

    return SimpleWeightedGraph{T,U}(adjm)
end

"""
    SimpleWeightedGraph{T}(n_vertices::Integer, n_edges::Integer; U::Type = Float64) where { T }

Random  SimpleWeightedGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedGraph{T}(
    n_vertices::Integer, n_edges::Integer; U::Type=Float64
) where {T}
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        adjm[cart_idx] = U(rand())
    end

    adjm .= (adjm .+ adjm') ./ U(2)

    return SimpleWeightedGraph{T,U}(adjm)
end

"""
    SimpleWeightedGraph{T,U}(n_vertices::Integer, n_edges::Integer) where { T, U }

Random  SimpleWeightedGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedGraph{T,U}(
    n_vertices::Integer, n_edges::Integer; weight_range::Tuple{U,U} = (zero(U),one(U)),
) where {T,U}
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        if U <: Integer
            adjm[cart_idx] = U(rand(weight_range[1]:weight_range[2]))
        elseif U <: Real
            adjm[cart_idx] = U(rand(Uniform(weight_range[1],weight_range[2])))
        else
            throw(ErrorException("The rand function suited for U = $U has not yet been implemented. Please file an issue."))
        end
    end

    adjm .= (adjm .+ adjm') ./ U(2)

    return SimpleWeightedGraph{T,U}(adjm)
end

"""
    SimpleWeightedDiGraph(n_vertices::Integer, n_edges::Integer; T::Type = Int64, U::Type = Float64)

Random  SimpleWeightedDiGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedDiGraph(
    n_vertices::Integer, n_edges::Integer; T::Type=Int64, U::Type=Float64
)
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        adjm[cart_idx] = U(rand())
    end

    return SimpleWeightedDiGraph{T,U}(adjm)
end

"""
    SimpleWeightedGraphs.SimpleWeightedDiGraph{T}(n_vertices::Integer, n_edges::Integer; U::Type = Float64) where {T}

Random  SimpleWeightedDiGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedDiGraph{T}(
    n_vertices::Integer, n_edges::Integer; U::Type=Float64
) where {T}
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        adjm[cart_idx] = U(rand())
    end

    return SimpleWeightedDiGraph{T,U}(adjm)
end

"""
    SimpleWeightedDiGraph{T,U}(n_vertices::Integer, n_edges::Integer)

Random  SimpleWeightedDiGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedDiGraph{T,U}(
    n_vertices::Integer, n_edges::Integer, weight_range::Tuple{U,U} = (zero(U),one(U))
) where {T,U}
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        if U <: Integer
            adjm[cart_idx] = U(rand(weight_range[1]:weight_range[2]))
        elseif U <: Real
            adjm[cart_idx] = U(rand(Uniform(weight_range[1],weight_range[2])))
        else
            throw(ErrorException("The rand function suited for U = $U has not yet been implemented. Please file an issue."))
        end
    end
    return SimpleWeightedDiGraph{T,U}(adjm)
end =#

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

