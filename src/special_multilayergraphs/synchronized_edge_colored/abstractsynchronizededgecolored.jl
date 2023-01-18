"""
"""
abstract type AbstractSynchronizedEdgeColoredGraph{T,U}  <: AbstractMultilayerGraph{T,U} end

# Nodes
"""
    add_node!(mg::AbstractSynchronizedEdgeColoredGraph, n::Node; add_vertex_to_layers::Union{Vector{Symbol}, Symbol} = Symbol[])

Add node `n` to `mg`. Return true if succeeds. Additionally, add a corresponding vertex to all layers.
"""
add_node!(mg::AbstractSynchronizedEdgeColoredGraph, n::Node) = _add_node!(mg, n; add_vertex_to_layers = :all) 


"""
    rem_node!(mg::AbstractSynchronizedEdgeColoredGraph, n::Node)

Remove node `n` to `mg`. Return true if succeeds.
"""
rem_node!(mg::AbstractSynchronizedEdgeColoredGraph, n::Node) = _rem_node!(mg, n)


# Vertices
## add_vertex! and rem_vertex! are not implemented, since oly nodes can be added/removed to/from synchronized graphs





