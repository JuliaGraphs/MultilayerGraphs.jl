

# Nodes
"""
    nodes(subgraph::AbstractSubGraph)

Return the collection of the nodes of `subgraph`.
"""
nodes(subgraph::AbstractSubGraph) = unique([mv.node for mv in mv_vertices(subgraph)])

# Vertices
"""
    Base.eltype(subgraph::AbstractSubGraph)

Return the vertex type of `subgraph`.
"""
Base.eltype(subgraph::AbstractSubGraph) = typeof(subgraph).parameters[1]

"""
    has_vertex( subgraph::S, v::T ) where {T,S<:AbstractSubGraph{T}}

Return `true` if `v` is a vertex of `subgraph`.
"""
function Graphs.has_vertex(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
    return has_vertex(subgraph.graph, v)
end

"""
    nv(subgraph::AbstractSubGraph)

Return the number of vertices in `subgraph`.
"""
Graphs.nv(subgraph::AbstractSubGraph) = nv(subgraph.graph)

"""
    vertices(subgraph::AbstractSubGraph)

Return the collection of the vertices of `subgraph`.
"""
Graphs.vertices(subgraph::AbstractSubGraph) = vertices(subgraph.graph)

"""
    mv_vertices(subgraph::AbstractSubGraph)

Return the collection of the `MultilayerVertex`s of `subgraph`.
"""
mv_vertices(subgraph::AbstractSubGraph) = get_rich_mv.(Ref(subgraph), vertices(subgraph))

"""
    inneighbors(subgraph::S, v::T) where {T, S <: AbstractSubGraph{T}} 

Return the list of inneighbors of `v` within `subgraph`.
"""
function Graphs.inneighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
    return inneighbors(subgraph.graph, v)
end

"""
    inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)

Return the list of inneighbors of `mv` within `subgraph`.
"""
function Graphs.inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
    return inneighbors(subgraph, get_v(subgraph, mv))
end

"""
    mv_inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex) 

Return the `MultilayerVertex`s inneighbors of `mv` within `subgraph`.
"""
function mv_inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
    return get_rich_mv.(Ref(subgraph), inneighbors(subgraph, mv))
end

"""
    outneighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}} 

Return the list of outneighbors of `v` within `subgraph`.
"""
function Graphs.outneighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
    return outneighbors(subgraph.graph, v)
end

"""
    outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
 
Return the list of outneighbors of `mv` within `subgraph`.
"""
function Graphs.outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
    return outneighbors(subgraph, subgraph.v_V_associations(get_bare_mv(mv)))
end

"""
    mv_outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)

Return the `MultilayerVertex`s outneighbors of `mv` within `subgraph`.
"""
function mv_outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
    return get_rich_mv.(Ref(subgraph), outneighbors(subgraph, mv))
end

"""
    neighbors(subgraph::S, v::T) where {T, S <: AbstractSubGraph{T}} 

Return the list of neighbors of `v` within `subgraph`.
"""
function Graphs.neighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
    return neighbors(subgraph.graph, v)
end

"""
    neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)

Return the list of neighbors of `mv` within `subgraph`.
"""
function Graphs.neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
    return outneighbors(subgraph, get_v(subgraph, mv))
end

"""
    mv_neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)  

Defaults to [`mv_outneighbors`](@ref).
"""
function mv_neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
    return mv_outneighbors(subgraph, mv)
end

"""
    get_v(subgraph::AbstractSubGraph, V::MultilayerVertex)

Return `v` associated with `V`. 
"""
function get_v(subgraph::AbstractSubGraph, V::MultilayerVertex)
    # Convert V to a bare vertex
    bare_V = get_bare_mv(V)
    # Check if subgraph has this vertex
    has_vertex(subgraph, bare_V) || return nothing
    # Get the list of edges
    return subgraph.v_V_associations(bare_V)
end

"""
    get_V(subgraph::S, v::T; perform_checks::Bool = false) where {T,U, S <: AbstractSubGraph{T,U}}

Return the `MultilayerVertex` whose internal representation is `v`.
"""
function get_V(
    subgraph::S, v::T; perform_checks::Bool=false
) where {T,U,S<:AbstractSubGraph{T,U}}
    # Check if v is a vertex label in subgraph
    if perform_checks
        haskey(subgraph.v_V_associations, v) || throw(
            ErrorException("$v is not an integer label of any vertex in the subgraph")
        )
    end
    # Return the vertex object
    return subgraph.v_V_associations[v]
end

"""
    get_rich_mv(mg::M, i::T) where {T,U, M <: AbstractSubGraph{T,U}}

Return `V` together with its metadata.
"""
function get_rich_mv(
    subgraph::S, i::T; perform_checks::Bool=false
) where {T,U,S<:AbstractSubGraph{T,U}}
    # Make sure that i is a valid vertex in the subgraph
    if perform_checks
        haskey(subgraph.v_V_associations, i) || throw(
            ErrorException("$i is not an integer label of any vertex in the subgraph")
        )
    end
    # Get the bare vertex corresponding to i
    bare_V = subgraph.v_V_associations[i]
    # Return the vertex as a multilayer vertex
    return MV(bare_V.node, bare_V.layer, get_metadata(subgraph, bare_V))
end

"""
    get_rich_mv(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex) 

Return `V` together with its metadata.
"""
function get_rich_mv(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex)
    has_vertex(subgraph, bare_mv) ||
        throw(ErrorException("$bare_mv is not part of $(subgraph.name)"))
    return MV(bare_mv.node, bare_mv.layer, get_metadata(subgraph, bare_mv))
end

"""
    get_metadata(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex)

Return the metadata of the vertex `bare_mv` in `subgraph` (metadata assigned to `bare_mv` will be discarded). 
"""
function get_metadata(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex)
    return _get_vertex_metadata(subgraph.graph, get_v(subgraph, bare_mv))
end

# Edges
"""
    edgetype(::S) where {T,U,S<:AbstractSubGraph{T,U}}

Return the edge type for `subgraph`.
"""
Graphs.edgetype(::S) where {T,U,S<:AbstractSubGraph{T,U}} = MultilayerEdge{U}

"""
    has_edge(subgraph::AbstractSubGraph,me::MultilayerEdge)

Return `true` if there is an edge from `src(me)` to `dst(me)` within subgraph, `false` otherwise.
"""
function Graphs.has_edge(subgraph::AbstractSubGraph, me::MultilayerEdge)
    return has_edge(subgraph, src(me), dst(me))
end

"""
    has_edge( subgraph::AbstractSubGraph, s::MultilayerVertex, d::MultilayerVertex)

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
function Graphs.has_edge(
    subgraph::AbstractSubGraph, s::MultilayerVertex, d::MultilayerVertex
)
    return has_edge(subgraph, get_v(subgraph, s), get_v(subgraph, d))
end

"""
    has_edge( subgraph::S, s::T, d::T) where {T,S<:AbstractSubGraph{T}}

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
function Graphs.has_edge(subgraph::S, s::T, d::T) where {T,S<:AbstractSubGraph{T}}
    return has_edge(subgraph.graph, s, d)
end

"""
    ne(subgraph::AbstractSubGraph)

Return the number of edges in `subgraph`.
"""
Graphs.ne(subgraph::AbstractSubGraph) = ne(subgraph.graph)

"""
    edges(subgraph::S) where {T,U,S<:AbstractSubGraph{T,U}} 

Return an iterator over all the edges of `subgraph`.
"""
function Graphs.edges(subgraph::S) where {T,U,S<:AbstractSubGraph{T,U}}
    return (
        MultilayerEdge(
            get_rich_mv(subgraph, src), get_rich_mv(subgraph, dst), weight, metadata
        ) for (src, dst, weight, metadata) in _edges(subgraph.graph, U)
    )
end

"""
    add_edge!( subgraph::S, me::E) where {T,U<:Real,S<:AbstractSubGraph{T,U},E<:MultilayerEdge{ <: Union{U, Nothing}}}

Add unweighted edge `me` to `subgraph`. Its `weight` and `metadata` fields are passed to the uniform interface of [`add_edge!(::Layer, ::MultilayerVertex, ::MultilayerVertex, ::Tuple)`](@ref).
"""
function Graphs.add_edge!(
    subgraph::S, me::E
) where {T,U<:Real,S<:AbstractSubGraph{T,U},E<:MultilayerEdge{<:Union{U,Nothing}}}
    return add_edge!(subgraph, src(me), dst(me); weight=weight(me), metadata=metadata(me))
end

"""
    add_edge_standard!(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple() ) where {T,U <: Real, W<:Union{U, Nothing},S<:AbstractSubGraph{T,U}}

Add edge from `src` to `dst` with weight `weight` to `subgraph`.
"""
function add_edge_standard!(
    subgraph::S,
    src::MultilayerVertex,
    dst::MultilayerVertex;
    weight::W=nothing,
    metadata::Union{Tuple,NamedTuple}=NamedTuple(),
) where {T,U<:Real,W<:Union{U,Nothing},S<:AbstractSubGraph{T,U}}
    return add_edge!(
        subgraph,
        get_v(subgraph, src),
        get_v(subgraph, dst);
        weight=weight,
        metadata=metadata,
    )
end

"""
    add_edge!(subgraph::S, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U<: Real, W<:Union{ U, Nothing},S<:AbstractSubGraph{T,U}}

Add edge from `src` to `dst` with weight `weight` and metadata `metadata` to `subgraph`.
"""
function Graphs.add_edge!(
    subgraph::S,
    src::T,
    dst::T;
    weight::W=nothing,
    metadata::Union{Tuple,NamedTuple}=NamedTuple(),
) where {T,U<:Real,W<:Union{U,Nothing},S<:AbstractSubGraph{T,U}}
    return _add_edge!(subgraph.graph, src, dst; weight=weight, metadata=metadata)
end

"""
    rem_edge!(subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex) 

Remove edge from `src` to `dst` in `subgraph`.
"""
function Graphs.rem_edge!(
    subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex
)
    return rem_edge!(subgraph, get_v(subgraph, src), get_v(subgraph, dst))
end

"""
    rem_edge!(subgraph::AbstractSubGraph, me::MultilayerEdge)

Remove edge from `src(me)` to `dst(me)` in `subgraph`.
"""
function Graphs.rem_edge!(subgraph::AbstractSubGraph, me::MultilayerEdge)
    return rem_edge!(subgraph, src(me), dst(me))
end

"""
    rem_edge!(subgraph::S, src::T, dst::T) where {T, S<:AbstractSubGraph{T}}

Remove edge from `src` to `dst` in a directed `subgraph`.
"""
function Graphs.rem_edge!(subgraph::S, src::T, dst::T) where {T,S<:AbstractSubGraph{T}}
    has_vertex(subgraph, src) && has_vertex(subgraph, dst) || throw(
        ErrorException(
            "One of the two vertices ($(subgraph.v_V_associations[src]) , $(subgraph.v_V_associations[dst])) (or both) does not belong to the subgraph.",
        ),
    )

    success = false
    if !has_edge(subgraph, src, dst)
        @warn "Edge from vertex $(subgraph.v_V_associations[src]) to vertex $(subgraph.v_V_associations[dst]) already doesn't exists in subgraph $(subgraph.name)."

        return false
    else
        success = rem_edge!(subgraph.graph, src, dst)
    end
end

"""
    get_edge_metadata(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex)

Return the metadata of the edge between the source vertex `src` and the destination vertex `dst` in `subgraph`. 
"""
function get_metadata(
    subgraph::S, src::MultilayerVertex, dst::MultilayerVertex
) where {T,U,G,S<:AbstractSubGraph{T,U,G}}
    return _get_edge_metadata(subgraph.graph, get_v(subgraph, src), get_v(subgraph, dst))
end # subgraph.edge_dict[subgraph.v_V_associations(src), subgraph.v_V_associations(dst)].metadata

"""
    get_weight(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex)

Return the weight of the edge between the source vertex `src` and the destination vertex `dst` in `subgraph`. 
"""
function SimpleWeightedGraphs.get_weight(
    subgraph::S, src::MultilayerVertex, dst::MultilayerVertex
) where {T,U,G,S<:AbstractSubGraph{T,U,G}}
    return _get_edge_weight(subgraph.graph, get_v(subgraph, src), get_v(subgraph, dst), U)
end

# Graph
"""
    is_directed(subgraph::AbstractSubGraph)

Return `true` if `subgraph` is directed, `false` otherwise. 
"""
Graphs.is_directed(subgraph::AbstractSubGraph) = is_directed(subgraph.descriptor.null_graph)

"""
    is_directed(::Type{S}) where {T,U,G,S <: AbstractSubGraph{T,U,G}}

Return `true` if instances of `S` are directed, `false` otherwise. 
"""
Graphs.is_directed(::Type{S}) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = is_directed(G)

"""
    adjacency_matrix(subgraph::AbstractSubGraph)

Return the adjacency matrix of `subgraph.graph`.
"""
Graphs.adjacency_matrix(subgraph::AbstractSubGraph) = adjacency_matrix(subgraph.graph)

"""
    weights(subgraph::S) where { T,U, G <: AbstractGraph{T}, S <:AbstractSubGraph{T,U,G}} 

Return the weights of `subgraph.graph`, with the eltype converted to `U`.
"""
function Graphs.weights(subgraph::S) where {T,U,S<:AbstractSubGraph{T,U}}
    return U.(weights(subgraph.graph))
end

# Utilities
"""
    Base.(==)(x::AbstractSubGraph, y::AbstractSubGraph)

Overload equality for `AbstractSubGraph`s.
"""
function Base.:(==)(x::AbstractSubGraph, y::AbstractSubGraph)
    # Check that each field in AbstractSubGraph is equal in x and y
    for field in fieldnames(AbstractSubGraph)
        # If the field is not equal in x and y, return false
        if @eval $x.$field != $y.$field
            return false
        end
    end
    # If all fields are equal in x and y, return true
    return true
end

"""
    name(subgraph::AbstractSubGraph)

Return the name of `subgraph`. 
"""
name(subgraph::AbstractSubGraph) = subgraph.name
