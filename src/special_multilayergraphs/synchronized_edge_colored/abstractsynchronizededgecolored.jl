"""
    AbstractSynchronizedEdgeColoredGraph{T,U}  <: AbstractMultilayerGraph{T,U} 

An abstract type representing an edge-colored and synchronized (i.e. every node is represented in each layer) graph. As such:

- `add_node!` will always add the corresponding vertex in all layers;
- `add_vertex!` and `rem_vertex!` are not available for this type;
-  All Interlayers automatically added by `add_layer!` are empty simple graphs.
- `specify_interlayer!` is not available.
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


# Edges
"""
    add_edge!(mg::M, me::E) where {T,U, M <: AbstractSynchronizedEdgeColoredGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add a MultilayerEdge between `src` and `dst` with weight `weight` and metadata `metadata`. Return true if succeeds, false otherwise.
"""
function Graphs.add_edge!(mg::M, me::E) where {T,U, M <: AbstractSynchronizedEdgeColoredGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}} 
    layer(src(me)) == layer(dst(me)) && throw(ArgumentErorr("Cannot add an interlayer edge to an edge colored graph"))
    _add_edge!(mg, me)

"""
    rem_edge!(mg::AbstractSynchronizedEdgeColoredGraph, me::MultilayerEdge)

Remove edge from `src(me)` to `dst(me)` from `mg`. Return true if succeeds, false otherwise.
"""
Graphs.rem_edge!(mg::AbstractSynchronizedEdgeColoredGraph, src::MultilayerVertex, dst::MultilayerVertex) = _rem_edge!(mg, src, dst)

# Layers and Interlayers
"""
    add_layer!( mg::M,
        new_layer::L;
    ) where {T,U,G<:AbstractGraph{T},L<:Layer{T,U,G}, H <: AbstractGraph{T}, M<:AbstractSynchronizedEdgeColoredGraph{T,U}; !IsDirected{M}}

Add layer `layer` to the synchronized edge-colored graph `mg`.

# ARGUMENTS

- `mg::M`: the `AbstractSynchronizedEdgeColoredGraph` which the new layer will be added to;
- `new_layer::L`: the new `Layer` to add to `mg`;
"""
function add_layer!(
    mg::M,
    new_layer::L
) where {
    T,
    U,
    G<:AbstractGraph{T},
    L<:Layer{T,U,G},
    M<:AbstractSynchronizedEdgeColoredGraph{T,U}
}

    return add_layer_directedness!(
        mg,
        new_layer;
        default_interlayers_structure="empty"
    )

end



