# CUSTOM TYPE/CLASS
# This file defines the custom type `Layer` and straightforwardly makes it compatible with the Graphs.jl ecosystem. The reason to have this custom type is to have a way not to break other people's code when modification to layers are required.

"""
    AbstractLayer{T,U,G}

An abstract type representing a generic Layer.

# FIELDS

- `T`: the node type;
- `U`: the adjacency matrix/tensor eltype;
- `G`: the underlying graph type.
"""
abstract type AbstractLayer{T,U,G} end

"""
    mutable struct Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}} <: AbstractLayer{T,U,G}

Represents a layer in a `Multilayer(Di)Graph`. 

# FIELDS

- `name::Symbol`: the name of the layer;
- `graph::G`: underlying graph of the layer;
- `forbidden_vertices::Vector{MultilayerVertex{T}}`: nodes of the MultilayerGraph that are not part of this Layer (they will be formally present in the Layer but it will be checked that they aren't adjacent to any other node);
- `forbidden_edges`::Vector{NTuple{2, MultilayerVertex{T}}}: edges that are required not to exist in this Layer.


# CONSTRUCTORS
    Layer(name::Symbol, graph::G, forbidden_vertices::Tuple{Vararg{T}}, forbidden_edges::Tuple{Vararg{NTuple{2, T}}};  U::Union{Type{ <: Real}, Nothing}  = nothing ) where {T,G <: AbstractGraph{T}}
         
Overridden inner constructor. Return an `Layer` whose underlying graph is `graph`. All `Layer`s and `Layer`s of a `Multilayer(Di)Graph` need to formally have the same nodes, but in real applications it may be that some vertices are excluded from some layers. Such vertices should be specified in `forbiddes_vertices`. Similarly for `forbiddes_edges`. This constructor (to which all the other eventually fall back to) will check that `forbidden_vertices` have no neighbors in `graph`, and that `forbidden_edges` actually correspond to zero entries in the adjacency matrix of `graph`. 

    Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}} <: AbstractLayer{T,U,G}

Incomplete initialization, used to write type-stable functions
"""
mutable struct Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}} <: AbstractLayer{T,U,G}
    name::Symbol
    graph::G
    forbidden_vertices::Vector{MultilayerVertex{T}}
    forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}
    
    # Override the inner constructor to check that `forbidden_vertices` have no neighbors in `graph`, and that `forbidden_edges` actually correspond to zero entries in the adjacency matrix of `graph`
    function Layer(name::Symbol, graph::G, forbidden_vertices::Tuple{Vararg{T}}, forbidden_edges::Tuple{Vararg{NTuple{2, T}}};  U::Union{Type{ <: Real}, Nothing}  = nothing ) where {T,G <: AbstractGraph{T}}

        for forbidden_node in forbidden_vertices
            isempty(all_neighbors(graph, forbidden_node))  || throw(ErrorException("The node $(forbidden_node) has been found to actually have neighbors $(all_neighbors(graph, forbidden_node)).\nThis is ny definition of forbidden nodes not allowed."))
        end

        _adjacency_matrix = adjacency_matrix(graph)
        for forbidden_edge in forbidden_edges
            view(_adjacency_matrix, forbidden_edge...) .== 0 || throw(ErrorException("The forbidden edge $(forbidden_edge), actually corresponds to a non-zero entry in the adjacency matrix of the underlying graph."))
        end

        _U = isnothing(U) ?  eltype(adjacency_matrix(graph)) : U

        _forbidden_vertices = [MultilayerVertex(i, name) for i in forbidden_vertices]
        _forbidden_edges    = [(MultilayerVertex(i, name), MultilayerVertex(j, name) ) for (i,j) in forbidden_edges]

        return new{T,_U,G}(name, graph, _forbidden_vertices,_forbidden_edges)
    end
    # Incomplete inner constructor, used to write type-stable functions
    function Layer(graph_type::Type{G}; U::Union{Type{ <: Real}, Nothing} = nothing) where { G <: AbstractGraph}
        _U = isnothing(U) ? eltype(adjacency_matrix(graph_type())) : U
        return new{graph_type.parameters[1], U, graph_type}()
    end
end

#= """
    Layer(nv::Int64, name::Symbol,layer_1::Symbol, layer_2::Symbol, graph_type::Type{G}, edge_list::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}, forbidden_vertices::Vector{MultilayerVertex{T}}, forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}}) where {T <: Union{ <: Integer, AbstractVertex}, U <: Real, G <: AbstractGraph{T}; IsWeighted{G}}


Outer constructor for `Layer`.  Return an `Layer` named `name` whose underlying graph type is `graph_type`. Edges are given via an edge list `edge_list`. All `Layer`s and `Interlayer`s of a `Multilayer(Di)Graph` need to formally have the same vertices, but in real applications it may be that some vertices are excluded from som layers. Such vertices should be specified in `forbiddes_vertices`. Similarly for `forbiddes_edges`.

# ARGUMENTS

- `nv`::Int64 : the number of nodes;
- `name``::Symbol : name of the Layer;
- `graph_type``::Type{G}: type of the graph underlying the Layer;
- `edge_list`::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}: edge list for the Layer;
- `forbidden_vertices`::Vector{MultilayerVertex{T}}: nodes of the MultilayerGraph that are not part of this Layer (they will be formally present in the Layer but it will be checked that they aren't adjacent to any other node);
- `forbidden_edges`::Vector{NTuple{2, MultilayerVertex{T}}}: edges that are required not to exist in this Layer.
"""
function Layer(nv::Int64, name::Symbol, graph_type::Type{G}, edge_list::Tuple{Vararg{ <: MultilayerEdge{MultilayerVertex{T},U} }}; forbidden_vertices::Vector{MultilayerVertex{T}} = MultilayerVertex{T}[], forbidden_edges::Vector{NTuple{2, MultilayerVertex{T}}} = NTuple{2, MultilayerVertex{T}}[]) where {T <: Union{ <: Integer, AbstractVertex}, U <: Real, G <: AbstractGraph{T}}
    
    @assert layer_1 != layer_2
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
        @assert issymmetric(adjm)
    else 
        for me in edge_list
            graph_vertices = [associations[me.src],associations[me.dst]]
            adjm[graph_vertices...] = me.weight
        end
    end

    graph = graph_type(adjm)
    
    Layer(name, layer_1,layer_2, graph, forbidden_vertices, forbidden_edges; U = U)

end
 =#

"""
    Layer(nv::Int64, name::Symbol, graph_type::Type{G}, ne::Int64; U::Union{Type{ <: Real},Nothing} = nothing)  where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}} 

Random `Layer`.

# ARGUMENTS

- `nv::Int64`: number of vertices;
- `name::Symbol`: name of the layer;
- `graph_type::Type{G}`: the underlying graph type;
- `ne::Int64`: the number of edges;
"""
function Layer(nv::Int64, name::Symbol, graph_type::Type{G}, ne::Int64; U::Union{Type{ <: Real},Nothing} = nothing)  where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}} 

    _U = isnothing(U) ?  eltype(adjacency_matrix(graph_type())) : U

    random_graph = graph_type(nv,ne)

    return Layer(name,random_graph,(), (); U = _U) #NTuple{2, MultilayerVertex{T}}[]

end

"""
    Layer(name::Symbol, graph::G; U::Union{Type{ <: Real}, Nothing} = nothing) where {G <: AbstractGraph} 

Return a Layer with name `name` and graph `graph` with no forbidden nodes or edges.
"""
Layer(name::Symbol, graph::G; U::Union{Type{ <: Real}, Nothing} = nothing) where {G <: AbstractGraph} = Layer(name, graph, (), (); U = U)
    
# Extend all Graphs.jl required methods (https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/)
"""
    Graphs.edges(layer::L) where {T,U,G <: AbstractGraph{T}, L <: Layer{T,U,G}}

Return an iterator over all the edges of `layer`.
"""
function Graphs.edges(layer::L) where {T,U,G <: AbstractGraph{T}, L <: Layer{T,U,G}}
    adjm = adjacency_matrix(layer.graph)

    hasweight = hasfield(edgetype(layer.graph), :weight)
    return (MultilayerEdge(MultilayerVertex(edge.src,  layer.name), MultilayerVertex(edge.dst,  layer.name), hasweight ? edge.weight : adjm[edge.src, edge.dst]) for edge in edges(layer.graph))
end

"""
    Base.eltype(layer::L) where { L <: Layer}

Return the vertex type of `layer`.
"""
Base.eltype(layer::L) where { L <: Layer}  = Base.eltype(layer.graph)

"""
    edgetype(layer::L) where {T,U,G,L <: Layer{T,U,G} }

Return the edge type for `layer`.
"""
Graphs.edgetype(layer::L) where {T,U,G,L <: Layer{T,U,G}} = MultilayerEdge{MultilayerVertex{T},U}

"""
    has_edge(layer::L, s::MultilayerVertex{T}, d::MultilayerVertex{T}) where { T,U,G, L <: Layer{T,U,G}}

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
Graphs.has_edge(layer::L,  s::MultilayerVertex{T}, d::MultilayerVertex{T}) where { T,U,G, L <: Layer{T,U,G}} = has_edge(layer.graph, layer.V_v_associations[s], layer.V_v_associations[d])

"""
    has_edge(layer::L,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, L <: Layer{T,U,G}; IsDirected{G}}

Return `true` if there is an edge between `s` and `d` with weight `weight`, `false` otherwise.
"""
@traitfn Graphs.has_edge(layer::L,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, L <: Layer{T,U,G}; IsDirected{G}} = has_edge(layer.graph, layer.V_v_associations[s], layer.V_v_associations[d]) && adjacency_matrix(layer.graph)[layer.V_v_associations[s],layer.V_v_associations[d]] == weight

"""
    has_edge(layer::L,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, L <: Layer{T,U,G}; IsDirected{G}}

Return `true` if there is an edge between `s` and `d` with weight `weight`, `false` otherwise.
"""
@traitfn Graphs.has_edge(layer::L,  s::MultilayerVertex{T}, d::MultilayerVertex{T}, weight::U) where { T,U,G, L <: Layer{T,U,G}; !IsDirected{G}} = has_edge(layer.graph, layer.V_v_associations[s], layer.V_v_associations[d]) && adjacency_matrix(layer.graph)[layer.V_v_associations[s],layer.V_v_associations[d]] == weight && adjacency_matrix(layer.graph)[layer.V_v_associations[d],layer.V_v_associations[s]] == weight

"""
    has_vertex(layer::L, v::Integer) where { L <: Layer} 

Return `true` if `v` is a vertex of `layer`.
"""
Graphs.has_vertex(layer::L, v::Integer) where { L <: Layer} = has_vertex(layer.graph, v)

"""
    inneighbors(layer::L, mv::MultilayerVertex{T}) where {T,U,G, L <: Layer{T,U,G}} 

Return the list of inneighbors of `v` within `layer`.
"""
Graphs.inneighbors(layer::L, mv::MultilayerVertex{T}) where {T,U,G, L <: Layer{T,U,G}} = [layer.v_V_associations[v] for v in  inneighbors(layer.graph, layer.V_v_associations[mv])] 

"""
    ne(layer::L) where { L <: Layer}

Return the number of edges in `layer`.
"""
Graphs.ne(layer::L) where { L <: Layer} = ne(layer.graph)

"""
    nv(layer::L) where { L <: Layer}

Return the number of vertices in `layer`.
"""
Graphs.nv(layer::L) where { L <: Layer} = nv(layer.graph) - length(layer.forbidden_vertices)


"""
    outneighbors(layer::L, mv::MultilayerVertex{T}) where {T,U,G, L <: Layer{T,U,G}} =

Return the list of outneighbors of `v` within `layer`.
"""
Graphs.outneighbors(layer::L, mv::MultilayerVertex{T}) where {T,U,G, L <: Layer{T,U,G}} = [layer.v_V_associations[v] for v in  outneighbors(layer.graph, layer.V_v_associations[mv])] 

"""
    vertices(layer::L) where { L <: Layer}

Return the collection of the vertices of `layer`.
"""
Graphs.vertices(layer::L) where { L <: Layer} = [MultilayerVertex(v, layer.name) for v in vertices(layer.graph)]

"""
    is_directed(layer::L) where { L <: Layer} 

Returns `true` if `layer` is directed, `false` otherwise. 
"""
Graphs.is_directed(layer::L) where { L <: Layer} = is_directed(layer.graph)

"""
    is_directed(::Type{L}) where {T,U,G,L <: Layer{T,U,G}}

Returns `true` if `layer_type` is directed, `false` otherwise. 
"""
Graphs.is_directed(::Type{L}) where {T,U,G,L <: Layer{T,U,G}} = is_directed(G)

"""
    add_edge!(layer::L,src::V, dst::V, weight::U) where { T, U <: Real, G <: AbstractGraph{T}, L <: Layer{T,U,G}, V <: MultilayerVertex{T}; IsWeighted{G}}
Add edge from `src` to `dst` with weight `weight` to `layer`.
"""
@traitfn function Graphs.add_edge!(layer::L,src::V, dst::V, weight::U) where { T, U <: Real, G <: AbstractGraph{T}, L <: Layer{T,U,G}, V <: MultilayerVertex{T}; IsWeighted{G}}

    has_vertex(layer,layer.V_v_associations[src]) && has_vertex(layer,layer.V_v_associations[dst]) || throw(ErrorException("One of the two vertices ($(src) , $(dst)) (or both) does not belong to the layer."))

    # Mimic Graphs.add_edge! behaviour
    if has_edge(layer, src, dst, weight)
        println("Edge from vertex $src to vertex $dst with the same weight ($weight) already exists in intelrayer $(layer.name)")
        return false
    else
        add_edge!(layer.graph, layer.V_v_associations[src], layer.V_v_associations[dst], weight)
    end
end

"""
    add_edge!(layer::L, src::V, dst::V) where { T, U <: Real, G <: AbstractGraph{T}, L <: Layer{T,U,G}, V <: MultilayerVertex{T}; !IsWeighted{G}}

Add edge from `src` to `dst` with weight `weight` to `layer`.
"""
@traitfn function Graphs.add_edge!(layer::L, src::V, dst::V) where { T, U <: Real, G <: AbstractGraph{T}, L <: Layer{T,U,G}, V <: MultilayerVertex{T}; !IsWeighted{G}}
    has_vertex(layer,layer.V_v_associations[src]) && has_vertex(layer,layer.V_v_associations[dst])  || throw(ErrorException("One of the two vertices ($(src) , $(dst)) (or both) does not belong to the layer."))
    if has_edge(layer, src, dst)
        println("Edge from vertex $src to vertex $dst already exists in intelrayer $(layer.name)")
        return false
    else
        add_edge!(layer.graph, layer.V_v_associations[src], layer.V_v_associations[dst])
    end
end

"""
    add_edge!(layer::L, me::E) where {T,U <: Real,G <: AbstractGraph{T} , L <: Layer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},U}; IsWeighted{G}}

Add unweighted edge `me` to `layer`.
"""
Graphs.add_edge!(layer::L, me::E) where {T,U <: Real,G <: AbstractGraph{T} , L <: Layer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},U}} = add_edge!(layer, src(me), dst(me), weight(me))

"""
    add_edge!(layer::L, me::E) where {T,U,G<: AbstractGraph{T}, L <: Layer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},Nothing}; !IsWeighted{G}}

Add unweighted edge `me` to `layer`.
"""
Graphs.add_edge!(layer::L, me::E) where {T,U,G<: AbstractGraph{T}, L <: Layer{T,U,G}, E <: MultilayerEdge{MultilayerVertex{T},Nothing}} = add_edge!(layer, src(me), dst(me))

"""
    rem_edge!(layer::L, src::V, dst::V) where { T, U, G<: AbstractGraph{T}, L <:Layer{T,U,G}, V <: MultilayerVertex{T}} 


Remove edge from `src` to `dst` in `layer`.
"""
function Graphs.rem_edge!(layer::L, src::V, dst::V) where { T, U, G<: AbstractGraph{T}, L <:Layer{T,U,G}, V <: MultilayerVertex{T}}
    has_vertex(layer,layer.V_v_associations[src]) && has_vertex(layer,layer.V_v_associations[dst])  || throw(ErrorException("One of the two vertices ($(src) , $(dst)) (or both) does not belong to the layer."))
    if !has_edge(layer, src, dst)
        println("The layer doesn't have any edge between $src and $dst")
        return false
    else
        rem_edge!( layer.graph,  layer.V_v_associations[src], layer.V_v_associations[dst] ) 
    end
end

"""
    adjacency_matrix(layer::L) where { T,U, G <: AbstractGraph{T}, L <:Layer{T,U,G}} 

Return the adjacency matrix of `layer.graph`, with the eltype converted to `U`.
"""
Graphs.adjacency_matrix(layer::L) where { T,U, G <: AbstractGraph{T}, L <:Layer{T,U,G}} = U.(Graphs.adjacency_matrix(layer.graph))
    
# Base overloads
"""
    Base.(==)(x::Layer, y::Layer)

Overload equality for `Layer`s.
"""
function Base.:(==)(x::Layer, y::Layer)
    for field in fieldnames(Layer)
        if @eval $x.$field != $y.$field
            return false
        end
    end
    return true
end

"""
    Base.getproperty(layer::L, f::Symbol) where { L <: layer}
"""
function Base.getproperty(layer::L, f::Symbol) where { L <: Layer }
    if f == :V_v_associations
        Dict(V => V.node for V in vertices(layer))
    elseif f == :v_V_associations
        Dict(V.node => V for V in vertices(layer))
    else
        Base.getfield(layer, f)
    end
end