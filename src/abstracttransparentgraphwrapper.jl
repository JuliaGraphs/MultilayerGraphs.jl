"""
    abstract type AbstractTrasparentGraphWrapper{T <: Integer,G <: AbstractGraph{T}} end

An abstract type representing a graph with vertices that are not Integer and for which it may be inappropriate to use a graph with vertex-level metadata. Its concrete subtypes MUST have the following fields:

- `graph::AbstractGraph{T}`: the underlyng graph used to store the relations (edges) between the vertices;
- `v_V_associations::Bijection`: A bijection from Bijections.jl that keepps track of the associations between the internal Integer representation of each vertex in the `graph` and the `AbstractVertex`s..
"""
abstract type AbstractTrasparentGraphWrapper{T <: Integer} <: AbstractGraph{T} end


# Vertices
"""
    Base.eltype(subgraph::S) where { S <: AbstractSubGraph}

Return the vertex type of `subgraph`.
"""
Base.eltype(subgraph::S) where {S<:AbstractSubGraph} = typeof(subgraph).parameters[1] 


"""
    has_vertex(subgraph::S, v::MultilayerVertex) where { T,U,G, S <: AbstractSubGraph{T,U,G}}

Return `true` if `v` is a vertex of `subgraph`.
"""
Graphs.has_vertex( subgraph::S, v::T ) where {T,U,G,S<:AbstractSubGraph{T,U,G}} =  has_vertex(subgraph.graph, v)


"""
    nv(subgraph::S) where { S <: AbstractSubGraph}

Return the number of vertices in `subgraph`.
"""
Graphs.nv(subgraph::S) where {S<:AbstractSubGraph} = nv(subgraph.graph)


"""
    vertices(subgraph::S) where {S <: AbstractSubGraph{ <: Integer, <: AbstractSimpleGraph}}

Return the collection of the vertices of `subgraph`.
"""
Graphs.vertices(subgraph::S) where {S<:AbstractSubGraph} = vertices(subgraph.graph)


"""
    inneighbors(subgraph::S, v::T) where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of inneighbors of `v` within `subgraph`.
"""
Graphs.inneighbors(subgraph::S, v::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = inneighbors(subgraph.graph, v)



"""
    outneighbors(subgraph::S, v::T) where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of outneighbors of `v` within `subgraph`.
"""
Graphs.outneighbors(subgraph::S, v::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = outneighbors(subgraph.graph, v)



"""
    neighbors(subgraph::S, v::T) where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of neighbors of `v` within `subgraph`.
"""
Graphs.neighbors(subgraph::S, v::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = neighbors(subgraph.graph, v)

"""
    get_V

The methods of this function must have a signature like `get_V(::AbstractTransparentGraphWrapper, ::Integer)` and return the `AbstractVertex` associated to `v` via the internal field `v_V_associations`. 
"""
function get_V end


# Edges
"""
    has_edge(subgraph::S, s::MultilayerVertex, d::MultilayerVertex) where { T,U,G, S <: AbstractSubGraph{T,U,G}}

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
Graphs.has_edge( subgraph::S, s::T, d::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = has_edge( subgraph.graph, s, d )


"""
    ne(subgraph::S) where { S <: AbstractSubGraph}

Return the number of edges in `subgraph`.
"""
Graphs.ne(subgraph::S) where {S<:AbstractSubGraph} = ne(subgraph.graph)


"""
    add_edge!(subgraph::S, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple()) 

Add edge from `src` to `dst` with weight `weight` and metadata `metadata` to `subgraph`.
"""
Graphs.add_edge!(subgraph::S, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U<: Real, W<:Union{ U, Nothing},G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} = _add_edge!(subgraph.graph, src, dst, weight = weight, metadata = metadata)


"""
    rem_edge!(subgraph::S, src::V, dst::V) where { T, U, S <: AbstractSubGraph{T,U}, V <: MultilayerVertex} 

Remove edge from `src` to `dst` in a directed `subgraph`.
"""
function Graphs.rem_edge!(
    subgraph::S, src::T, dst::T
) where {T,U, G, S<:AbstractSubGraph{T,U,G}}
    has_vertex(subgraph, src) && has_vertex(subgraph, dst) || throw(
        ErrorException(
            "One of the two vertices ($(subgraph.v_V_associations[src]) , $(subgraph.v_V_associations[dst])) (or both) does not belong to the subgraph.",
        ),
    )

    success = false
    if !has_edge(
        subgraph, src, dst
    )
        @warn "Edge from vertex $(subgraph.v_V_associations[src]) to vertex $(subgraph.v_V_associations[dst]) already doesn't exists in subgraph $(subgraph.name)."
        
        return false
    else
        success = rem_edge!(
            subgraph.graph,
            src,
            dst
        )
    end
end


"""
    get_edge_metadata(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex)

Return the metadata of the edge between the source vertex `src` and the destination vertex `dst` in `subgraph`. 
"""
get_metadata(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex) where {T,U,G, S <: AbstractSubGraph{T,U,G}} =  _get_edge_metadata(subgraph.graph, get_v(subgraph, src), get_v(subgraph, dst)) # subgraph.edge_dict[subgraph.v_V_associations(src), subgraph.v_V_associations(dst)].metadata

"""
    get_weight(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex)

Return the weight of the edge between the source vertex `src` and the destination vertex `dst` in `subgraph`. 
"""
SimpleWeightedGraphs.get_weight(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex) where {T,U,G, S <: AbstractSubGraph{T,U,G}} =  _get_edge_weight(subgraph.graph, get_v(subgraph,src), get_v(subgraph,dst), U) 

# Graph


"""
    weights(subgraph::S) where { T,U, G <: AbstractGraph{T}, S <:AbstractSubGraph{T,U,G}} 

Return the weights of `subgraph.graph`, with the eltype converted to `U`.
"""
weights(subgraph::S) where {T,U,G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} = U.(weights(subgraph.graph))




