# Some custom extensions of Graphs.jl and SimpleWeightedGraphs.jl

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
        adjm[cart_idx] = rand()
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
        adjm[cart_idx] = rand()
    end

    adjm .= (adjm .+ adjm') ./ U(2)

    return SimpleWeightedGraph{T,U}(adjm)
end

"""
    SimpleWeightedGraph{T,U}(n_vertices::Integer, n_edges::Integer) where { T, U }

Random  SimpleWeightedGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedGraph{T,U}(
    n_vertices::Integer, n_edges::Integer
) where {T,U}
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        adjm[cart_idx] = rand()
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
        adjm[cart_idx] = rand()
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
        adjm[cart_idx] = rand()
    end

    return SimpleWeightedDiGraph{T,U}(adjm)
end

"""
    SimpleWeightedDiGraph{T,U}(n_vertices::Integer, n_edges::Integer)

Random  SimpleWeightedDiGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. Edge weights are uniformly extracted between 0 and 1. 
"""
function SimpleWeightedGraphs.SimpleWeightedDiGraph{T,U}(
    n_vertices::Integer, n_edges::Integer
) where {T,U}
    adjm = zeros(U, n_vertices, n_vertices)
    rand_nnz_cart_idxs = rand(CartesianIndices(adjm), n_edges)

    for cart_idx in rand_nnz_cart_idxs
        adjm[cart_idx] = rand()
    end

    return SimpleWeightedDiGraph{T,U}(adjm)
end
