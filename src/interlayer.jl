# CUSTOM TYPE/CLASS
# This file defines the custom type `Interlayer` and straightforwardly makes it compatible with the Graphs.jl ecosystem. The reason to have this custom type is to have a way not to break other people's code when modification to interlayers are required.

"""
    AbstractInterlayer{T,U,G}

An abstract type representing a generic Interlayer.

# FIELDS

- `T`: the node type;
- `U`: the adjacency matrix/tensor eltype;
- `G`: the underlying graph type.
"""
abstract type AbstractInterlayer{T,U,G} end # <: Interlayer{T,U,G}

"""
    mutable struct Interlayer{G <: AbstractGraph}

Represents an interlayer in a `Multilayer(Di)Graph`. 

# FIELDS

- `nv`::Int64 : the number of nodes;
- `name``::Symbol : name of the Interlayer;
- `layer_1`::Symbol: name of one of the Layers connected by this Interlayer;
- `layer_2`::Symbol: name of one of the Layers connected by this Interlayer;
- `graph_type``::Type{G}: type of the graph underlying the Interlayer;
- `edge_list`::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}: edge list for the Interlayer;
- `forbidden_vertices`::Vector{MultilayerVertex{T}}: nodes of the MultilayerGraph that are not part of this Interlayer (they will be formally present in the Interlayer but it will be checked that they aren't adjacacent to any other node);
- `forbidden_edges`::Vector{NTuple{2, MultilayerVertex{T}}}: edges that are required not to exist in this Interlayer.


# Constructors
    Interlayer(name::Symbol,layer_1::Symbol,layer_2::Symbol, graph::G, forbidden_vertices::Vector{MultilayerVertex{T}}, forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}; U::Union{Type{ <: Real}, Nothing} = nothing  ) where { T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}}
         

Overridden inner constructor. Return an `Interlayer` between `layer_1` and `layer_2`, whose underlying graph is `graph`. All `Layer`s and `Interlayer`s of a `Multilayer(Di)Graph` need to formally have the same vertices, but in real applications it may be that some vertices are excluded from som layers. Such vertices should be specified in `forbiddes_vertices`. Similarly for `forbiddes_edges`. This constructor (to which all the other eventually fall back to) will check that `forbidden_vertices` have no neighbors in `graph`, and that `forbidden_edges` actually correspond to zero entries in the adjacency matrix of `graph`. Also check that it is a proper bipartite graph.

"""
mutable struct Interlayer{T <: Union{ <: Integer, AbstractVertex}, U <: Real, G <: AbstractGraph{T}} <: AbstractInterlayer{T,U,G}
    name::Symbol
    layer_1::Symbol
    layer_2::Symbol
    graph::G
    forbidden_vertices::Vector{MultilayerVertex{T}}
    forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}

    # Override the inner constructor to check that `forbidden_vertices` have no neighbors in `graph`, and that `forbidden_edges` actually correspond to zero entries in the adjacency matrix of `graph`. Also check that it is a proper bipartite graph.
    function Interlayer(name::Symbol,layer_1::Symbol,layer_2::Symbol, graph::G, forbidden_vertices::Vector{MultilayerVertex{T}}, forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}; U::Union{Type{ <: Real}, Nothing} = nothing  ) where { T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}}

        is_bipartite(graph) || throw(ErrorException("The specified underlying `graph` is not bipartite."))

        for forbidden_node in forbidden_vertices
            isempty(all_neighbors(graph, forbidden_node)) || throw(ErrorException("The node $(forbidden_node) has been found to actually have neighbors $(all_neighbors(graph, forbidden_node)).\nThis is ny definition of forbidden nodes not allowed."))
        end

        _adjacency_matrix = Graphs.adjacency_matrix(graph)
        for forbidden_edge in forbidden_edges
            view(_adjacency_matrix, forbidden_edge...) .== 0 || throw(ErrorException("The forbidden edge $(forbidden_edge), actually corresponds to a non-zero entry in the adjacency matrix of the underlying graph."))
        end

        _U = isnothing(U) ? eltype(adjacency_matrix(graph)) :  U

        return new{T,_U,G}(name, layer_1, layer_2, graph, forbidden_vertices, forbidden_edges)
    end
end


# Outer constructors for unweighted and weighted Interlayers
#= """
    Interlayer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}, edge_list::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}; forbidden_vertices::Vector{MultilayerVertex{T}} = MultilayerVertex{T}[], forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}} = NTuple{2, MultilayerVertex{T}}[]) where {T <: Union{ <: Integer, AbstractVertex}, U <: Real, G <: AbstractGraph{T}}

Outer constructor for unweighted `Interlayer`.

## ARGUMENTS

- `nv`::Int64 : the number of nodes;
- `name``::Symbol : name of the Interlayer;
- `layer_1`::Symbol: name of one of the Layers connected by this Interlayer;
- `layer_2`::Symbol: name of one of the Layers connected by this Interlayer;
- `graph_type``::Type{G}: type of the graph underlying the Interlayer;
- `edge_list`::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}: edge list for the Interlayer;
- `forbidden_vertices`::Vector{MultilayerVertex{T}}: nodes of the MultilayerGraph that are not part of this Interlayer (they will be formally present in the Interlayer but it will be checked that they aren't adjacent to any other node);
- `forbidden_edges`::Vector{NTuple{2, MultilayerVertex{T}}}: edges that are required not to exist in this Interlayer.

"""
 =#

"""
    Interlayer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}, edge_list::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}, forbidden_vertices::Vector{MultilayerVertex{T}}, forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}) where {T <: Union{ <: Integer, AbstractVertex}, U <: Real, G <: AbstractGraph{T}; IsWeighted{G}}


Outer constructor for `Interlayer`.  Return an `Interlayer` between `layer_1` and `layer_2`, whose underlying graph type is `graph_type`. Edges are given via an edge list `edge_list`. All `Layer`s and `Interlayer`s of a `Multilayer(Di)Graph` need to formally have the same vertices, but in real applications it may be that some vertices are excluded from som layers. Such vertices should be specified in `forbiddes_vertices`. Similarly for `forbiddes_edges`.

# ARGUMENTS

- `nv`::Int64 : the number of nodes;
- `name``::Symbol : name of the Interlayer;
- `layer_1`::Symbol: name of one of the Layers connected by this Interlayer;
- `layer_2`::Symbol: name of one of the Layers connected by this Interlayer;
- `graph_type``::Type{G}: type of the graph underlying the Interlayer;
- `edge_list`::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}: edge list for the Interlayer;
- `forbidden_vertices`::Vector{MultilayerVertex{T}}: nodes of the MultilayerGraph that are not part of this Interlayer (they will be formally present in the Interlayer but it will be checked that they aren't adjacent to any other node);
- `forbidden_edges`::Vector{NTuple{2, MultilayerVertex{T}}}: edges that are required not to exist in this Interlayer.

"""
function Interlayer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}, edge_list::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}; forbidden_vertices::Vector{MultilayerVertex{T}} = MultilayerVertex{T}[], forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}} = NTuple{2, MultilayerVertex{T}}[]) where {T <: Union{ <: Integer, AbstractVertex}, U <: Real, G <: AbstractGraph{T}}
    
    layer_1 != layer_2 || throw(ErrorException("Argument `layer_1` cannot match argument `layer_2`"))
    layer_1_vertices = [MultilayerVertex(v,layer_1) for v in T.(1:nv)] # Set{MutlilayerVertex{T}}
    layer_2_vertices = [MultilayerVertex(v,layer_2) for v in T.(1:nv)] # Set{MutlilayerVertex{T}}
    _vertices = T.(1:(2*nv))
    associations = Dict{MultilayerVertex{T},T}(V => v for (v,V) in zip(_vertices, vcat(layer_1_vertices,layer_2_vertices)))
  
    adjm = zeros(U,2*nv,2*nv)
    if !istrait(IsDirected{G}) 
        for me in edge_list
            graph_vertices = [associations[me.src],associations[me.dst]]
            adjm[graph_vertices...] = adjm[reverse(graph_vertices)...] =  length(unique(graph_vertices)) == 1 ? me.weight + me.weight : me.weight
        end
        issymmetric(adjm) || throw(ErrorException("The provided `edge_list` does not correspond to a symmetric adjacency matrix. Simmetry is required since the underlying graph type $(G) is undirected."))
    else 
        for me in edge_list
            graph_vertices = [associations[me.src],associations[me.dst]]
            adjm[graph_vertices...] = me.weight
        end
    end

    graph = graph_type(adjm)
    
    Interlayer(name, layer_1,layer_2, graph, forbidden_vertices, forbidden_edges; U = U)

end

"""
    Interlayer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}, ne::Int64, forbidden_vertices::Vector{MultilayerVertex{T}}, forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}) where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}; !IsDirected{G}, IsWeighted{G}}

Random `Interlayer`.

# ARGUMENTS

- `nv::Int64`: number of vertices;
- `name::Symbol`: name of the Interlayer;
- `layer_1::Symbol`: The first Layer it connects;
- `layer_2::Symbol`: The second Layer it connects;
- `graph_type::Type{G}`: the underlying graph type;
- `ne::Int64`: the number of edges;
- `forbidden_vertices::Vector{MultilayerVertex{T}}`: list of vertices that are not considered present in the Interlayer;
- `forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}`: list of edges whose existence is a priori excluded from the Interlayer.
"""
function Interlayer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}, ne::Int64; U::Union{Type{ <: Real},Nothing} = nothing)  where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}} 

    _U = isnothing(U) ?  eltype(adjacency_matrix(graph_type())) : U

    edge_list = nothing
    if istrait(IsDirected{G})
        edge_list = Tuple(MultilayerEdge( MultilayerVertex(rand(1:nv), l1), MultilayerVertex(rand(1:nv), setdiff([layer_1,layer_2], [l1])[1] ), rand(_U) ) for (i,l1) in zip(1:ne,[rand([layer_1,layer_2]) for i in 1:ne ]))
    else
        edge_list = Tuple(MultilayerEdge( MultilayerVertex(rand(1:nv), layer_1), MultilayerVertex(rand(1:nv), layer_2), rand(_U) ) for i in 1:ne)
    end

    return Interlayer(nv,name,layer_1,layer_2,graph_type,edge_list,forbidden_vertices = MultilayerVertex{T}[], forbidden_edges = NTuple{2, MultilayerVertex{T}}[])

end


"""
    multiplex_interlayer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}; forbidden_vertices::Vector{MultilayerVertex{T}}, forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}) where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T};  IsWeighted{G}}

Return an `Interlayer{T,U,G}` that has edges only between vertices that represent the same node.

# ARGUMENTS

- `nv::Int64`: number of vertices;
- `name::Symbol`: name of the Interlayer;
- `layer_1::Symbol`: The first Layer it connects;
- `layer_2::Symbol`: The second Layer it connects;
- `graph_type::Type{G}`: the underlying graph type;
- `U::Type = Float64`: the eltype of the adjacency matrix of the Interlayer. Note: it doesn't have to coincide with the underlying graph's adjacency matrix eltype, since right now there is no guarantee that all Graphs.jl's extension will allow the user to set such eltype on new graph types.
- `forbidden_vertices::Vector{MultilayerVertex{T}}`: list of vertices that are not considered present in the Interlayer;
- `forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}`: list of edges whose existence is a priori excluded from the Interlayer.
"""
function multiplex_interlayer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}; U::Type = Float64 ,forbidden_vertices::Vector{MultilayerVertex{T}} = MultilayerVertex{T}[], forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}} = NTuple{2, MultilayerVertex{T}}[]) where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}}
    
    edge_list = Tuple( MultilayerEdge( MultilayerVertex( i, layer_1), MultilayerVertex( i, layer_2), one(U) ) for i in 1:nv )# Like SimpleWeightedGraphs.jl, we assume that G.parameters[2] is the weight type.

    idxs_tbr = Int64[]
    for forbidden_vertex in forbidden_vertices
        idxs_tbr = [i for (i,edge) in enumerate(edge_list) if !((edge.src == forbidden_vertex || edge.dst == forbidden_vertex) || edge ∈ forbidden_edges ) ]
    end
    deleteat!.(Ref(edge_list),idxs_tbr )

    
    return Interlayer(nv,name,layer_1,layer_2,graph_type,edge_list,forbidden_vertices = forbidden_vertices,forbidden_edges = forbidden_edges) 
end

# Extend all Graphs.jl required methods (https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/)
"""
    Graphs.edges(interlayer::In) where {T,U,G <: AbstractGraph{T}, In <: Interlayer{T,U,G}}

Return an iterator over all the edges of `interlayer`.
"""
function Graphs.edges(interlayer::In) where {T,U,G <: AbstractGraph{T}, In <: Interlayer{T,U,G}}
    adjm = adjacency_matrix(interlayer.graph)
    v_V_associations = interlayer.v_V_associations
    hasweight = hasfield(edgetype(interlayer.graph), :weight)
    return (MultilayerEdge(v_V_associations[edge.src], v_V_associations[edge.dst], hasweight ? edge.weight : adjm[edge.src, edge.dst]) for edge in  edges(interlayer.graph))
end

"""
    Base.eltype(interlayer::In) where { In <: Interlayer}

Return the vertex type of `interlayer`.
"""
Base.eltype(interlayer::In) where { In <: Interlayer}  = Base.eltype(interlayer.graph)

"""
    edgetype(interlayer::In) where {T,U,G,In <: Interlayer{T,U,G} }

Return the edge type for `interlayer`.
"""
Graphs.edgetype(interlayer::In) where {T,U,G,In <: Interlayer{T,U,G} } = MultilayerEdge{MultilayerVertex{T},U}

"""
    has_edge(interlayer::In, s::MultilayerVertex{T}, d::MultilayerVertex{T}) where { T,U,G, In <: Interlayer{T,U,G}}

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
Graphs.has_edge(interlayer::In, s::MultilayerVertex{T}, d::MultilayerVertex{T}) where { T,U,G, In <: Interlayer{T,U,G}} = has_edge(interlayer.graph, interlayer.V_v_associations[s], interlayer.V_v_associations[d])

"""
    has_edge(interlayer::In,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, In <: Interlayer{T,U,G}; IsDirected{G}

Return `true` if there is an edge between `s` and `d` with weight `weight`, `false` otherwise.
"""
@traitfn Graphs.has_edge(interlayer::In,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, In <: Interlayer{T,U,G}; IsDirected{G}} = has_edge(interlayer.graph, interlayer.V_v_associations[s], interlayer.V_v_associations[d]) && adjacency_matrix(interlayer.graph)[interlayer.V_v_associations[s],interlayer.V_v_associations[d]] == weight

"""
    has_edge(interlayer::In,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, In <: Interlayer{T,U,G}; !IsDirected{G}}

    Return `true` if there is an edge between `s` and `d` with weight `weight`, `false` otherwise.
"""
@traitfn Graphs.has_edge(interlayer::In,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, In <: Interlayer{T,U,G}; !IsDirected{G}} = has_edge(interlayer.graph, interlayer.V_v_associations[s], interlayer.V_v_associations[d]) && adjacency_matrix(interlayer.graph)[interlayer.V_v_associations[s],interlayer.V_v_associations[d]] == weight && adjacency_matrix(interlayer.graph)[interlayer.V_v_associations[s],interlayer.V_v_associations[d]] == weight


"""
    has_vertex(interlayer::In, v::MultilayerVertex{T}) where { T,U,G, In <: Interlayer{T,U,G}}

Return `true` if `v` is a vertex of `interlayer`.
"""
Graphs.has_vertex(interlayer::In, v::MultilayerVertex{T}) where { T,U,G, In <: Interlayer{T,U,G}} = !(v in interlayer.forbidden_vertices)

"""
    inneighbors(interlayer::In, mv::MultilayerVertex{T}) where {T,U,G, In <: Interlayer{T,U,G}}

Return the list of inneighbors of `v` within `interlayer`.
"""
Graphs.inneighbors(interlayer::In, mv::MultilayerVertex{T}) where {T,U,G, In <: Interlayer{T,U,G}} = [interlayer.v_V_associations[v] for v in  inneighbors(interlayer.graph, interlayer.V_v_associations[mv])]

"""
    ne(interlayer::In) where { In <: Interlayer}

Return the number of edges in `interlayer`.
"""
Graphs.ne(interlayer::In) where { In <: Interlayer} = ne(interlayer.graph)

"""
    nv(interlayer::In) where { In <: Interlayer}

Return the number of vertices in `interlayer`.
"""
Graphs.nv(interlayer::In) where { In <: Interlayer} = length(vertices(interlayer)) - length(interlayer.forbidden_vertices)

"""
    outneighbors(interlayer::In, v::T) where {In <: Interlayer{T} } where { T <: Integer}

Return the list of outneighbors of `v` within `interlayer`, looping first over all layers (in the order they are given in `interlayer.layers`), then over all interlayers (in the order they are given in `interlayer.interlayers`).
"""
Graphs.outneighbors(interlayer::In, mv::MultilayerVertex{T}) where {T,U,G, In <: Interlayer{T,U,G}} = [interlayer.v_V_associations[v] for v in  outneighbors(interlayer.graph, interlayer.V_v_associations[mv])]

"""
    vertices(interlayer::In) where {In <: Interlayer{ <: Integer, <: AbstractSimpleGraph}}

Return the collection of the vertices of `interlayer`.
"""
Graphs.vertices(interlayer::In) where { In <: Interlayer} = collect(Iterators.flatten((interlayer.layer_1_vertices,interlayer.layer_2_vertices)))

"""
    is_directed(interlayer::In) where { In <: Interlayer} 

Returns `true` if `interlayer` is directed, `false` otherwise. 
"""
Graphs.is_directed(interlayer::In) where { In <: Interlayer} = is_directed(interlayer.graph)

"""
    is_directed(interlayer_type::Type{In}) where {In <: Interlayer}

Returns `true` if `interlayer_type` is directed, `false` otherwise. 
"""
Graphs.is_directed(interlayer_type::Type{In}) where {In <: Interlayer} = is_directed(interlayer_type.parameters[1])

"""
    adjacency_matrix(interlayer::In) where { T,U, G <: AbstractGraph{T}, In <:Interlayer{T,U,G}}

Return the adjacency matrix of `interlayer.graph`, with the eltype converted to `U`.
"""
Graphs.adjacency_matrix(interlayer::In) where { T,U, G <: AbstractGraph{T}, In <:Interlayer{T,U,G}} = U.(adjacency_matrix(interlayer.graph))

"""
    add_edge!(interlayer::In,src::V, dst::V, weight::U) where { T, U <: Real, G <: AbstractGraph{T}, In <: Interlayer{T,U,G}, V <: MultilayerVertex{T}; IsWeighted{G}}

Add edge from `src` to `dst` with weight `weight` to `interlayer`.
"""
@traitfn function Graphs.add_edge!(interlayer::In,src::V, dst::V, weight::U) where { T, U <: Real, G <: AbstractGraph{T}, In <: Interlayer{T,U,G}, V <: MultilayerVertex{T}; IsWeighted{G}}
    has_vertex(interlayer,src) && has_vertex(interlayer,dst)  || throw(ErrorException("One of the two vertices ($(src) , $(dst)) (or both) does not belong to the interlayer."))
    if has_edge(interlayer, src, dst, weight)
        println("Edge from vertex $src to vertex $dst with weight ($weight) already exists in interlayer $(interlayer.name)")
        return false
    else
        add_edge!(interlayer.graph, interlayer.V_v_associations[src], interlayer.V_v_associations[dst], weight)
    end
end

"""
    add_edge!(interlayer::In,src::V, dst::V) where { T, U <: Real, G <: AbstractGraph{T}, In <: Interlayer{T,U,G}, V <: MultilayerVertex{T}; !IsWeighted{G}}

Add edge from `src` to `dst` with weight `weight` to `interlayer`.
"""
@traitfn function Graphs.add_edge!(interlayer::In,src::V, dst::V) where { T, U <: Real, G <: AbstractGraph{T}, In <: Interlayer{T,U,G}, V <: MultilayerVertex{T}; !IsWeighted{G}}

    has_vertex(interlayer,src) && has_vertex(interlayer,dst)  || throw(ErrorException("One of the two vertices ($(src) , $(dst)) (or both) does not belong to the interlayer."))

    if has_edge(interlayer, src, dst)
        println("Edge from vertex $src to vertex $dst already exists in interlayer $(interlayer.name)")
        return false
    else
        add_edge!(interlayer.graph, interlayer.V_v_associations[src], interlayer.V_v_associations[dst])
    end

end

"""
    add_edge!(interlayer::In, me::E) where {T,U <: Real,G, In <: Interlayer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},U}}

Add unweighted edge `me` to `interlayer`.
"""
Graphs.add_edge!(interlayer::In, me::E) where {T,U <: Real,G, In <: Interlayer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},U}} = add_edge!(interlayer, src(me), dst(me), weight(me))

"""
    add_edge!(interlayer::In, me::E) where {T,U,G, In <: Interlayer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},Nothing}}

Add unweighted edge `me` to `interlayer`.
"""
Graphs.add_edge!(interlayer::In, me::E) where {T,U,G, In <: Interlayer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},Nothing}} = add_edge!(interlayer, src(me), dst(me))

"""
    rem_edge!(interlayer::In, src::V, dst::V) where { T, U, In <: Interlayer{T,U}, V <: MultilayerVertex{T}} 


Remove edge from `src` to `dst` in `interlayer`.
"""
function Graphs.rem_edge!(interlayer::In, src::V, dst::V) where { T, U, In <: Interlayer{T,U}, V <: MultilayerVertex{T}} 

    has_vertex(interlayer,src) && has_vertex(interlayer,dst)  || throw(ErrorException("One of the two vertices ($(src) , $(dst)) (or both) does not belong to the interlayer."))

    if !has_edge(interlayer.graph, interlayer.V_v_associations[src], interlayer.V_v_associations[dst])
        println("Edge from vertex $src to vertex $dst already doesn't exists in interlayer $(interlayer.name)")
        return false
    else
        rem_edge!(interlayer.graph, interlayer.V_v_associations[src], interlayer.V_v_associations[dst])
    end
end

# Base overloads
"""
    Base.(==)(x::Interlayer, y::Interlayer)

Overload equality for `Interlayer`s.
"""
function Base.:(==)(x::Interlayer, y::Interlayer)
    for field in fieldnames(Interlayer)
        if @eval $x.$field != $y.$field
            return false
        end
    end
    return true
end


"""
    Base.getproperty(interlayer::In, f::Symbol) where { In <: Interlayer}
"""
function Base.getproperty(interlayer::In, f::Symbol) where { In <: Interlayer }
    if f == :layers
        [Base.getfield(interlayer, :layer_1),Base.getfield(interlayer, :layer_2)]
    elseif f == :layer_1_vertices || f == :U 
        num_vertices_graph = nv(interlayer.graph)
        (MultilayerVertex(v, interlayer.layer_1) for v in vertices(interlayer.graph)[1:(num_vertices_graph ÷ 2)]) 
    elseif f == :layer_2_vertices || f == :V
        num_vertices_graph = nv(interlayer.graph)
        (MultilayerVertex(v, interlayer.layer_2) for v in vertices(interlayer.graph)[1:(num_vertices_graph ÷ 2)]) 
    elseif f == :V_v_associations
        OrderedDict(V => v for (v,V) in zip(vertices(interlayer.graph), Iterators.flatten((interlayer.layer_1_vertices, interlayer.layer_2_vertices) )))
    elseif f == :v_V_associations
        OrderedDict(v => V for (v,V) in zip(vertices(interlayer.graph), Iterators.flatten((interlayer.layer_1_vertices, interlayer.layer_2_vertices) )))
    else
        Base.getfield(interlayer, f)
    end
end

"""
    get_symmetric_interlayer(interlayer::In; symmetric_interlayer_name::String) where{T,U,G, In <: Interlayer{T,U,G}}

Return the `Interlayer` corresponding to `interlayer` where `layer_1` and `layer_2` are swapped. Its name will be `symmetric_interlayer_name` (defaults to `interlayer_(interlayer.layer_2)_(interlayer.layer_1)` ).
"""
function get_symmetric_interlayer(interlayer::In; symmetric_interlayer_name::String = "interlayer_$(interlayer.layer_2)_$(interlayer.layer_1)") where{T,U,G, In <: Interlayer{T,U,G}}

    interlayer_adjm = adjacency_matrix(interlayer.graph)
    symmetric_interlayer_adjm = similar(interlayer_adjm)
    nv = size(interlayer_adjm,1) ÷ 2

    range_1 = 1:nv
    range_2 = (nv+1):2*nv

    @views symmetric_interlayer_adjm[range_1, range_1] = interlayer_adjm[range_2, range_2]
    @views symmetric_interlayer_adjm[range_1, range_2] = interlayer_adjm[range_2, range_1]
    @views symmetric_interlayer_adjm[range_2, range_1] = interlayer_adjm[range_1, range_2]
    @views symmetric_interlayer_adjm[range_2, range_2] = interlayer_adjm[range_1, range_1]

    symmetric_graph = G(symmetric_interlayer_adjm)

    return Interlayer(Symbol(symmetric_interlayer_name),interlayer.layer_2,interlayer.layer_1,symmetric_graph, interlayer.forbidden_vertices, interlayer.forbidden_edges; U = U )
end