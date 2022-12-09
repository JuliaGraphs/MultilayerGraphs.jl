"""
    mutable struct GraphOfGraphs{T} <: AbstractGraphOfGraphs{T}

Default undirected concrete implementation of `AbstractGraphOfGraphs`.

# FIELDS

- `nodes::OrderedDict{Int64,T}`: An `OrderedDict` that associates graphs (values) to Int64 nodes (keys).
- `graph::SimpleWeightedGraph{Int64,Float64}`: The weighted undirected graph representing the graph of graphs. Its Int64 nodes are associated to graphs via the `nodes` field.
"""
mutable struct GraphOfLayers{T,G} <: AbstractGraphOfGraphs{T,G}
    v_L_associations::Bijection{T, <: Layer}
    graph::G
end

"""
    GraphOfGraphs(graphs::Vector{T}, adjacency_matrix::Matrix{Float64})

Outer constructor for GraphOfGraphs.

# ARGUMENTS

- `graphs::Vector{T}`: The vector of nodes of the graphs of graphs. Each of them will be assigned to an integer, from 1 to length(graphs), and their order is assumed to be the same as the rows and columns of the `adjacency_matrix` argument.
- `adjacency_matrix::Matrix{Float64}`: The adjacency matrix of the graph of graphs. Its (i,j)-element is the weight of the link between `graphs[i]` and `graphs[j]`.
"""
function GraphOfGraphs(
    graphs::Vector{T},
    adjacency_matrix::Union{Symmetric{Float64,Matrix{Float64}},Matrix{Float64}},
) where {T<:AbstractGraph}
    size(adjacency_matrix, 1) == length(graphs) || throw(
        ErrorException(
            "The `adjacency_matrix` has more columns than the number of graphs specified in `graphs`. Found $(size(adjacency_matrix,1)) and $(length(graphs)).",
        ),
    )
    nodes = OrderedDict{Int64,T}(i => graph for (i, graph) in enumerate(graphs))
    return GraphOfGraphs{T}(nodes, SimpleWeightedGraph(adjacency_matrix))
end