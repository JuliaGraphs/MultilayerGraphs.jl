# CUSTOM TYPE/CLASS
# This file defines the custom type `Layer` and straightforwardly makes it compatible with the Graphs.jl ecosystem. The reason to have this custom type is to have a way not to break other people's code when modification to layers are required.

"""
    AbstractLayer{T,U,G}

An abstract type representing a generic Layer.

# FIELDS

- `T`: the node type;
- `U`: the `MultilayerEdge` weight eltype;
- `G`: the underlying graph type.
"""
abstract type AbstractLayer{T,U,G}  <: AbstractSubGraph{T,U,G} end

"""
    mutable struct Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}} <: AbstractLayer{T,U,G}

Represents a layer in a `Multilayer(Di)Graph`. 

# FIELDS

- `name::Symbol`: the name of the layer;
- `graph::G`: underlying graph of the layer;
- `forbidden_vertices::Vector{MultilayerVertex}`: nodes of the MultilayerGraph that are not part of this Layer (they will be formally present in the Layer but it will be checked that they aren't adjacent to any other node);
- `forbidden_edges`::Vector{NTuple{2, MultilayerVertex}}: edges that are required not to exist in this Layer.

# CONSTRUCTORS

    Layer(name::Symbol, graph::G, forbidden_vertices::Tuple{Vararg{T}}, forbidden_edges::Tuple{Vararg{NTuple{2, T}}};  U::Union{Type{ <: Real}, Nothing}  = nothing ) where {T,G <: AbstractGraph{T}}
         
Overridden inner constructor. Return an `Layer` whose underlying graph is `graph`. All `Layer`s and `Layer`s of a `Multilayer(Di)Graph` need to formally have the same nodes, but in real applications it may be that some vertices are excluded from some layers. Such vertices should be specified in `forbidden_vertices`. Similarly for `forbidden_edges`. This constructor (to which all the other eventually fall back to) will check that `forbidden_vertices` have no neighbors in `graph`, and that `forbidden_edges` actually correspond to zero entries in the adjacency matrix of `graph`. 

    Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}} <: AbstractLayer{T,U,G}

Incomplete initialization, used to write type-stable functions
"""
mutable struct Layer{T<:Integer,U<:Real,G<:AbstractGraph{T}} <: AbstractLayer{T,U,G}
    descriptor::LayerDescriptor{T,U,G}
    graph::G
    v_V_associations::Bijection{T, <: MultilayerVertex}

    # Inner constructor that performs checks on request
    function Layer(descriptor::LayerDescriptor{T,U,G}, graph::G, v_V_associations::Bijection{T, <: MultilayerVertex}; check_consistency = true) where {T,U,G}
        if check_consistency
            typeof(descriptor.null_graph) == typeof(graph) || throw(ErrorException("Graph types between the provided `descriptor` and `graph` cannot differ. Found $(typeof(descriptor.null_graph)) and $(typeof(graph))."))

            all(vertices(graph) .== sort(domain(v_V_associations))) || throw(ErrorException("The graph has a different set of vertices w.r.t. the domain of `v_V_associations`. Found $(vertices(graph)) and $(sort(domain(v_V_associations)))."))
        end

        return new{T,U,G}(descriptor, graph, v_V_associations)
    end
end

# TODO: 
# Need two constructors:
# 1. takes a MultilayerVertex list and a graph as input (for e.g. underlying graphs like metagraphs) and then we match MultilayerVertices and graph's vertices by order (check that the numbers of vertices are the same)
# 2. takes a Multilayer(Weighted)Edge list and a graph type as input
# 3. takes a graph as input. The MultilayerVertices that will go into v_N_associations will have a default name or an uninitialized name
# From these constructors we may later implement the random one

# 1. takes a MultilayerVertex list and a graph as input (for e.g. underlying graphs like metagraphs) and then we match MultilayerVertices and graph's vertices by order (check that the numbers of vertices are the same)
"""
    Layer(name::Symbol, vertices::Vector{<: MultilayerVertex}, edge_list::Vector{ <: MultilayerEdge}, null_graph::G, weighttype::Type{U};  default_vertex_metadata::Function = mv -> NamedTuple(), default_edge_weight::Function = (src, dst) -> one(U), default_edge_metadata::Function = (src, dst) -> NamedTuple()) where {T <: Integer, U <: Real,  G <: AbstractGraph{T}}

...
"""
function Layer(name::Symbol, vertices::Vector{<: MultilayerVertex}, edge_list::Vector{ <: MultilayerEdge}, null_graph::G, weighttype::Type{U};  default_vertex_metadata::Function = mv -> NamedTuple(), default_edge_weight::Function = (src, dst) -> one(U), default_edge_metadata::Function = (src, dst) -> NamedTuple()) where {T <: Integer, U <: Real,  G <: AbstractGraph{T}}
    descriptor = LayerDescriptor(name, null_graph, weighttype,  default_vertex_metadata = default_vertex_metadata, default_edge_weight = default_edge_weight, default_edge_metadata = default_edge_metadata)

    return Layer(descriptor, vertices, edge_list)
end

"""
    Layer(descriptor::LayerDescriptor{T,U,G}, vertices::Vector{<: MultilayerVertex}, edge_list::Vector{<:MultilayerEdge}) where {T <: Integer, U <: Real, G <: AbstractGraph{T}}

...
"""
function Layer(descriptor::LayerDescriptor{T,U,G}, vertices::Vector{<: MultilayerVertex}, edge_list::Vector{<:MultilayerEdge}) where {T <: Integer, U <: Real, G <: AbstractGraph{T}}

    if hasproperty(eltype(vertices), :parameters)
        par = eltype(vertices).parameters[1]
        (isnothing(par) || par == descriptor.name) || throw(ErrorException("`vertices` should be a `Vector{MultilayerVertex{:$(descriptor.name)}}` or a `Vector{MultilayerVertex{nothing}}`. Found $(typeof(vertices))"))
    else
        throw(ErrorException("`vertices` should be a `Vector{MultilayerVertex{:$(descriptor.name)}}` or a `Vector{MultilayerVertex{nothing}}`. Found $(typeof(vertices))"))
    end

    layer = Layer(descriptor, deepcopy(descriptor.null_graph),  Bijection{T, MultilayerVertex{descriptor.name}}(), check_consistency = false)

    for mv in vertices
        add_vertex!(layer, mv)
    end

    for edge in edge_list
        add_edge!(layer, edge)
    end

    return layer
end

"""
    Layer(nv::Int64, name::Symbol, graph_type::Type{G}, ne::Int64; U::Union{Type{ <: Real},Nothing} = nothing)  where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}} 

Return a random `Layer`.

# ARGS

- `name::Symbol`: The name of the Layer
- `vertices::Vector{ <: MultilayerVertex}`: The `MultilayerVertex`s of the Layer
- `ne::Int64`: The number of edges of the Layer
- `null_graph::G`: the Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
- `weighttype::Type{U}`: The type of the `MultilayerEdge` weights (evem when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted);

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `allow_self_loops::Bool`: whether to allow self loops to be geenrated or not. Deafults to `false`.
"""
function Layer(
    name::Symbol,
    vertices::Vector{ <: MultilayerVertex},
    ne::Int64,
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function = mv -> NamedTuple(),
    default_edge_weight::Function = (src, dst) -> nothing,
    default_edge_metadata::Function = (src, dst) -> NamedTuple(),
    allow_self_loops::Bool = false

) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

    descriptor = LayerDescriptor(name, null_graph, weighttype; default_vertex_metadata = default_vertex_metadata, default_edge_weight = default_edge_weight, default_edge_metadata = default_edge_metadata) 
    edge_list = MultilayerEdge[]

    for i in 1:ne
        rand_vertex_1  = rand(vertices)
        rand_vertex_2  = nothing
        if !allow_self_loops
            while isnothing(rand_vertex_2) || rand_vertex_2.node == rand_vertex_1.node
                rand_vertex_2 = rand(vertices)
            end
        else
            rand_vertex_2 = rand(vertices)
        end
        push!(edge_list, MultilayerEdge(MV(rand_vertex_1.node, name), MV(rand_vertex_2.node, name), default_edge_weight(rand_vertex_1,rand_vertex_2), default_edge_metadata(rand_vertex_1,rand_vertex_2 )))
    end

    edge_list = MultilayerEdge[ rand() < 0.5 ? me : reverse(me) for me in  edge_list]
    layer = Layer(descriptor, vertices,  edge_list)

    return layer 
end

"""
    has_node(layer::L, n::Node) where { L <: Layer} 

Return `true` if `n` is a node of `layer`.
"""
has_node(layer::L, n::Node) where {L<:Layer} = MV(n, layer.name) ∈ image(layer.v_V_associations)


"""
    has_vertex(layer::L, v::MultilayerVertex) where { T,U,G, L <: Layer{T,U,G}}

Return `true` if `v` is a vertex of `layer`.
"""
Graphs.has_vertex(layer::L, mv::MultilayerVertex) where {T,U,G,L<:Layer{T,U,G}} = MV(node(mv), name(layer)) ∈ collect(image(layer.v_V_associations))

# TODO:
# Implement a MultilayerVertex constructor that leaves the .layer field unspecified, for ease of use of the following function
"""
    add_vertex!(layer::L, mv::MultilayerVertex) where { L <: Layer} 

Add vertex to layer `layer`. 
"""
function Graphs.add_vertex!(layer::L, mv::MultilayerVertex) where { L <: Layer}  
    (isnothing(mv.layer) || mv.layer == layer.name) || throw(ErrorException("The multilayer vertex $mv cannot belong to layer $(layer.name)."))
    add_vertex!(layer, mv.node; metadata = mv.metadata )
end

"""
    add_vertex!(layer::L, n::Node, args...; kwargs...) where {T, U, G,  L <: Layer{T,U,G}}      

Add vertex associated with node `n` to layer `layer`. This method supports the uniform and transparent interfaces. See the [Vertices](@ref) section of the Tutorial.
"""
function Graphs.add_vertex!(layer::L, n::Node, args...; kwargs...) where {T, U, G,  L <: Layer{T,U,G}} 
    has_node(layer, n) && return false

    success = false
    if isempty(args) && length(kwargs) == 1 && issetequal(Set([:metadata]), Set(keys(kwargs)) )
        success = add_vertex_standard!(layer; metadata = values(kwargs).metadata)
    elseif length(args) == length(kwargs) == 0
        success = add_vertex_standard!(layer, metadata = layer.default_vertex_metadata(MV(n, layer.name)))
    else
        success =  add_vertex!(layer.graph, args...; kwargs... )
    end

    if success
        last_vertex = length(layer.v_V_associations) == 0 ? zero(T) : maximum(domain(layer.v_V_associations))
        layer.v_V_associations[last_vertex + one(T)] = MV(n, layer.name)
        return true
    else
        return false
    end
end

"""
    add_vertex_standard!(layer::Layer; metadata::Union{Tuple, NamedTuple}= NamedTuple())

Add vertex with metadata to layer `layer`.
"""
add_vertex_standard!(layer::Layer; metadata::Union{Tuple, NamedTuple}= NamedTuple()) = _add_vertex!(layer; metadata = metadata)

"""
    _add_vertex!( layer::L; metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U, G, L <: Layer{T,U,G}}

Add vertex with metadata to layer `layer`.
"""
_add_vertex!( layer::L; metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U, G, L <: Layer{T,U,G}} = __add_vertex!(layer.graph; metadata = metadata)
     
"""
    rem_vertex!(layer::Layer, mv::MultilayerVertex) 

Remove vertex `mv` from layer `layer`.
"""
function Graphs.rem_vertex!(layer::Layer, mv::MultilayerVertex) 
    (isnothing(mv.layer) || mv.layer == layer.name) || return false
    rem_vertex!(layer, mv.node)
end


"""
    rem_vertex!(layer::Layer, n::Node)

Remove node `n` from `layer`. Modify `layer.v_N_associations` according to how `rem_vertex!` works in [Graph.jl](https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs/#Graphs.SimpleGraphs.rem_vertex!-Tuple{Graphs.SimpleGraphs.AbstractSimpleGraph,%20Integer}).
"""
function Graphs.rem_vertex!(layer::Layer, n::Node)

    !has_node(layer, n) && return false

    success = rem_vertex!(layer, layer.v_V_associations(MV(n, layer.name)))

    if success
        # Get the key of the node to be removed
        v = layer.v_V_associations(MV(n,layer.name))
        # Get the last node and its key
        last_v = maximum(domain(layer.v_V_associations))

        if v != last_v
            last_V = layer.v_V_associations[last_v]
            # Delete the node to be removed
            delete!(layer.v_V_associations, v)
            # Delete the last node
            delete!(layer.v_V_associations, last_v)
            # Re-add the last node with the key of the node to be removed
            layer.v_V_associations[v] = last_V
            true
        else
            # Delete the node to be removed
            delete!(layer.v_V_associations, v)
            true
        end
    else
        return false
    end
end

"""
    rem_vertex!(layer::L, v::T) where {T,U,G, L <: Layer{T,U,G}}    

Remove vertex `v` from layer `layer`.
"""
Graphs.rem_vertex!(layer::L, v::T) where {T,U,G, L <: Layer{T,U,G}} = rem_vertex!(layer.graph, v)

"""
    add_edge!(layer::L, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...) where {L <: Layer} 

Add edge from vertex `src` to vertex `dst` to layer `layer`. This method supports the uniform and transparent interfaces. See the [Edges](@ref) section of the Tutorial.
"""
function Graphs.add_edge!(layer::L, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...) where {L <: Layer} 
    bare_src = get_bare_mv(src)
    bare_dst = get_bare_mv(dst)

    !has_vertex(layer, src)  && throw( ErrorException( "Vertex $(src) does not belong to the layer."))
    !has_vertex(layer, dst) && throw( ErrorException( "Vertex $(dst) does not belong to the layer."))

    if !has_edge(layer, src, dst)
        if isempty(args) && length(kwargs) == 2 && issetequal(Set([:weight, :metadata]), Set(keys(kwargs)) )
            success = add_edge_standard!(layer, bare_src, bare_dst, weight = values(kwargs).weight, metadata = values(kwargs).metadata)

        elseif isempty(args) && length(kwargs) == 1 && issetequal(Set([:weight]), Set(keys(kwargs)) )
            success = add_edge_standard!(layer, bare_src, bare_dst, weight = values(kwargs).weight, metadata = layer.default_edge_metadata(bare_src, bare_dst))

        elseif isempty(args) && length(kwargs) == 1 && issetequal(Set([:metadata]), Set(keys(kwargs)) )
            success = add_edge_standard!(layer, bare_src, bare_dst, weight = layer.default_edge_weight(bare_src, bare_dst), metadata =  values(kwargs).metadata)

        elseif length(args) == length(kwargs) == 0
            success = add_edge_standard!(layer, bare_src, bare_dst, weight = layer.default_edge_weight(bare_src, bare_dst), metadata = layer.default_edge_metadata(bare_src, bare_dst))
        else
            success = add_edge!(layer.graph, get_v(layer,bare_src), get_v(layer,bare_dst), args...; kwargs... )
        end

        return success
    else
        return false
    end
end

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
    Base.getproperty(layer::L, f::Symbol) where {L<:Layer}

Overload `getproperty` for `Layer`s.
"""
function Base.getproperty(layer::L, f::Symbol) where {L<:Layer}
    if f ∈ (:descriptor, :graph, :v_V_associations)
        Base.getfield(layer, f) 
    elseif f ∈ (:name, :null_graph, :default_vertex_metadata, :default_edge_weight, :default_edge_metadata)
        Base.getfield(layer.descriptor, f)
    end
end