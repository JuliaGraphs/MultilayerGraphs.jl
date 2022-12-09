abstract type AbstractSubGraph{T <: Integer,U <: Real,G <: AbstractGraph{T}} end

# Nodes
"""
    nodes(subgraph::S) where { S <: AbstractSubGraph}

Return the collection of the nodes of `subgraph`.
"""
nodes(subgraph::S) where {S<:AbstractSubGraph} = unique([mv.node for mv in mv_vertices(subgraph)])

# Vertices
"""
    Base.eltype(subgraph::S) where { S <: AbstractSubGraph}

Return the vertex type of `subgraph`.
"""
Base.eltype(subgraph::S) where {S<:AbstractSubGraph} = typeof(subgraph).parameters[1] 

"""
    has_vertex(subgraph::S, v::T ) where {T,U,G,S<:AbstractSubGraph{T,U,G}}

Return `true` if `v` is a vertex of `subgraph`.
"""
Graphs.has_vertex( subgraph::S, v::T ) where {T,U,G,S<:AbstractSubGraph{T,U,G}} =  has_vertex(subgraph.graph, v)

"""
    nv(subgraph::S) where { S <: AbstractSubGraph}

Return the number of vertices in `subgraph`.
"""
Graphs.nv(subgraph::S) where {S<:AbstractSubGraph} = nv(subgraph.graph) #length(vertices(subgraph))

"""
    vertices(subgraph::S) where {S <: AbstractSubGraph{ <: Integer, <: AbstractSimpleGraph}}

Return the collection of the vertices of `subgraph`.
"""
Graphs.vertices(subgraph::S) where {S<:AbstractSubGraph} = vertices(subgraph.graph)

"""
    mv_vertices(subgraph::S) where {S <: AbstractSubGraph{ <: Integer, <: AbstractSimpleGraph}}

Return the collection of the `MultilayerVertex`s of `subgraph`.
"""
mv_vertices(subgraph::S) where {S<:AbstractSubGraph} = get_rich_mv.(Ref(subgraph), vertices(subgraph)) 

"""
    inneighbors(subgraph::S, v::T) where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of inneighbors of `v` within `subgraph`.
"""
Graphs.inneighbors(subgraph::S, v::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = inneighbors(subgraph.graph, v)

"""
    inneighbors(subgraph::S, mv::MultilayerVertex)  where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of inneighbors of `mv` within `subgraph`.
"""
Graphs.inneighbors(subgraph::S, mv::MultilayerVertex) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = inneighbors(subgraph, get_v(subgraph,mv))

"""
    mv_inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex) 

Return the `MultilayerVertex`s inneighbors of `mv` within `subgraph`.
"""
mv_inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex) = get_rich_mv.(Ref(subgraph), inneighbors(subgraph,mv))

"""
    outneighbors(subgraph::S, v::T) where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of outneighbors of `v` within `subgraph`.
"""
Graphs.outneighbors(subgraph::S, v::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = outneighbors(subgraph.graph, v)

"""
    outneighbors(subgraph::S, mv::MultilayerVertex)  where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of outneighbors of `mv` within `subgraph`.
"""
Graphs.outneighbors(subgraph::S, mv::MultilayerVertex) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = outneighbors(subgraph, subgraph.v_V_associations(get_bare_mv(mv)))

"""
    mv_outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)

Return the `MultilayerVertex`s outneighbors of `mv` within `subgraph`.
"""
mv_outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex) = get_rich_mv.(Ref(subgraph), outneighbors(subgraph,mv))


"""
    neighbors(subgraph::S, v::T) where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of neighbors of `v` within `subgraph`.
"""
Graphs.neighbors(subgraph::S, v::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = neighbors(subgraph.graph, v)

"""
    neighbors(subgraph::S, mv::MultilayerVertex)  where {T,U,G, S <: AbstractSubGraph{T,U,G}} 

Return the list of neighbors of `mv` within `subgraph`.
"""
Graphs.neighbors(subgraph::S, mv::MultilayerVertex) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = outneighbors(subgraph, get_v(subgraph,mv))


"""
    mv_neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)  

Defaults to [`mv_outneighbors`](@ref).
"""
mv_neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex) = mv_outneighbors(subgraph, mv)


"""
    get_v(subgraph::AbstractSubGraph, V::MultilayerVertex)

Return `v` associated with `V`. 
"""
function get_v(subgraph::AbstractSubGraph, V::MultilayerVertex)
    bare_V = get_bare_mv(V)
    has_vertex(subgraph, bare_V) || return nothing

    subgraph.v_V_associations(bare_V)
end


"""
    get_V(subgraph::S, v::T; perform_checks::Bool = false) where {T,U, S <: AbstractSubGraph{T,U}}

Return the `MultilayerVertex` whose internal representation is `v`.
"""
function get_V(subgraph::S, v::T; perform_checks::Bool = false) where {T,U, S <: AbstractSubGraph{T,U}}
    if perform_checks
        haskey(subgraph.v_V_associations,v) || throw(ErrorException("$v is not an integer label of any vertex in the subgraph"))
    end

    return subgraph.v_V_associations[v]
end




"""
    get_rich_mv(mg::M, i::T) where {T,U, M <: AbstractSubGraph{T,U}}

Return `V` together with its metadata.
"""
function get_rich_mv(subgraph::S, i::T; perform_checks::Bool = false) where {T,U, S <: AbstractSubGraph{T,U}}
    if perform_checks
        haskey(subgraph.v_V_associations,i) || throw(ErrorException("$i is not an integer label of any vertex in the subgraph"))
    end

    bare_V = subgraph.v_V_associations[i]
    return MV(bare_V.node, bare_V.layer,  get_metadata(subgraph, bare_V))
end

"""
    get_rich_mv(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex) 

Return `V` together with its metadata.
"""
function get_rich_mv(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex) 
    has_vertex(subgraph, bare_mv) || throw(ErrorException("$bare_mv is not part of $(subgraph.name)"))
    MV(bare_mv.node, bare_mv.layer, get_metadata(subgraph, bare_mv))
end

"""
    get_metadata(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex)

Return the metadata of `V`. 
"""
get_metadata(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex)  =_get_vertex_metadata(subgraph.graph, get_v(subgraph, bare_mv))

# Edges
"""
    edgetype(subgraph::S) where {T,U,G,S <: AbstractSubGraph{T,U,G} }

Return the edge type for `subgraph`.
"""
Graphs.edgetype(::S) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = MultilayerEdge{U}

"""
    has_edge(subgraph::S, s::MultilayerVertex, d::MultilayerVertex) where { T,U,G, S <: AbstractSubGraph{T,U,G}}

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
Graphs.has_edge( subgraph::S,me::MultilayerEdge) where {T,U,G,S<:AbstractSubGraph{T,U,G}} =  has_edge(subgraph, src(me), dst(me))

"""
    has_edge(subgraph::S, s::MultilayerVertex, d::MultilayerVertex) where { T,U,G, S <: AbstractSubGraph{T,U,G}}

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
Graphs.has_edge( subgraph::S, s::MultilayerVertex, d::MultilayerVertex) where {T,U,G,S<:AbstractSubGraph{T,U,G}} =  has_edge( subgraph, get_v(subgraph, s), get_v(subgraph, d))

"""
    has_edge(subgraph::S, s::MultilayerVertex, d::MultilayerVertex) where { T,U,G, S <: AbstractSubGraph{T,U,G}}

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
Graphs.has_edge( subgraph::S, s::T, d::T) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = has_edge( subgraph.graph, s, d )

"""
    ne(subgraph::S) where { S <: AbstractSubGraph}

Return the number of edges in `subgraph`.
"""
Graphs.ne(subgraph::S) where {S<:AbstractSubGraph} = ne(subgraph.graph) # length(subgraph.edge_list)



"""
    edges(subgraph::S) where {T,U,G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}}

Return an iterator over all the edges of `subgraph`.
"""
function Graphs.edges(subgraph::S) where {T,U,G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} # = subgraph.edge_list
    (
        MultilayerEdge(
            get_rich_mv(subgraph, src),
            get_rich_mv(subgraph, dst),
            weight,
            metadata
        ) for (src,dst,weight,metadata) in _edges(subgraph.graph,U)
    )
end

"""
    add_edge!(subgraph::S, me::E) where {T,U <: Real,G <: AbstractGraph{T} , S <: AbstractSubGraph{T,U,G}, E <: MultilayerEdge{U}}

Add unweighted edge `me` to `subgraph`. Its `weight` and `metadata` fields are passed to the uniform interface of [`add_edge!(::Layer, ::MultilayerVertex, ::MultilayerVertex, ::Tuple)`](@ref).
"""
Graphs.add_edge!( subgraph::S, me::E) where {T,U<:Real,G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G},E<:MultilayerEdge{ <: Union{U, Nothing}}} = add_edge!(subgraph, src(me), dst(me); weight =  weight(me), metadata = metadata(me)) 

"""
    add_edge!(subgraph::S,src::V, dst::V, weight::U) where { T, U <: Real, G <: AbstractGraph{T}, S <: AbstractSubGraph{T,U,G}}

Add edge from `src` to `dst` with weight `weight` to `subgraph`.
"""
add_edge_standard!(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple() ) where {T,U <: Real, W<:Union{U, Nothing},G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} = add_edge!( subgraph, get_v(subgraph,src),get_v(subgraph,dst), weight =  weight, metadata = metadata)

"""
    add_edge!(subgraph::S, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple()) 

Add edge from `src` to `dst` with weight `weight` and metadata `metadata` to `subgraph`.
"""
Graphs.add_edge!(subgraph::S, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U<: Real, W<:Union{ U, Nothing},G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} = _add_edge!(subgraph.graph, src, dst, weight = weight, metadata = metadata)

"""
    rem_edge!(subgraph::S, src::V, dst::V) where { T, U, S <: AbstractSubGraph{T,U}, V <: MultilayerVertex} 

Remove edge from `src` to `dst` in `subgraph`.
"""
Graphs.rem_edge!(subgraph::S, src::MultilayerVertex, dst::MultilayerVertex) where {T,U,S<:AbstractSubGraph{T,U}} = rem_edge!(subgraph, get_v(subgraph, src), get_v(subgraph, dst))

"""
    rem_edge!(subgraph::S, me::E) where { T, U, S <: AbstractSubGraph{T,U}, E<:MultilayerEdge} 

Remove edge `me` in `subgraph`.
"""
Graphs.rem_edge!(subgraph::S, me::E) where {T,U,S<:AbstractSubGraph{T,U},E<:MultilayerEdge} = rem_edge!(subgraph, src(me), dst(me) )

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
    is_directed(subgraph::S) where { S <: AbstractSubGraph} 

Return `true` if `subgraph` is directed, `false` otherwise. 
"""
Graphs.is_directed(subgraph::S) where {S<:AbstractSubGraph} = is_directed(subgraph.descriptor.null_graph)

"""
    is_directed(::Type{S}) where {T,U,G,S <: AbstractSubGraph{T,U,G}}

Return `true` if instances of `S` are directed, `false` otherwise. 
"""
Graphs.is_directed(::Type{S}) where {T,U,G,S<:AbstractSubGraph{T,U,G}} = is_directed(G)

"""
    adjacency_matrix(subgraph::S) where { T,U, G <: AbstractGraph{T}, S <:AbstractSubGraph{T,U,G}} 

Return the adjacency matrix of `subgraph.graph`.
"""
Graphs.adjacency_matrix(subgraph::S) where {T,U,G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} = adjacency_matrix(subgraph.graph) # weights(subgraph) .> zero(U)

"""
    weights(subgraph::S) where { T,U, G <: AbstractGraph{T}, S <:AbstractSubGraph{T,U,G}} 

Return the weights of `subgraph.graph`, with the eltype converted to `U`.
"""
weights(subgraph::S) where {T,U,G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} = U.(weights(subgraph.graph))

# Utilities
"""
    Base.(==)(x::AbstractSubGraph, y::AbstractSubGraph)

Overload equality for `AbstractSubGraph`s.
"""
function Base.:(==)(x::AbstractSubGraph, y::AbstractSubGraph)
    for field in fieldnames(AbstractSubGraph)
        if @eval $x.$field != $y.$field
            return false
        end
    end
    return true
end

"""
    name(subgraph::AbstractSubGraph)

Return the name of `subgraph`. 
"""
name(subgraph::AbstractSubGraph) = subgraph.name