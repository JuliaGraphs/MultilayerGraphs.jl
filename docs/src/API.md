# API

This page lists all exported methods, organizing them by topic (whether it is a method that acts on vertices, edges, layers, etc) and by audience: due to how it was possible to integrate multilayer graphs within the `Graphs.jl` ecosystem, some methods are intended for developers who wish to use this library just as any other `Graphs.jl` package in their code, while others are meant to be employed by the end-user.


## End-user

### Nodes

```@docs
Node

```

### Vertices


```@docs
eltype
MultilayerVertex
MV
node
layer
metadata(mv::MultilayerVertex)
MissingVertex
```

### Edges


```@docs
MultilayerEdge
ME
weight(e::AbstractMultilayerEdge)
metadata(e::AbstractMultilayerEdge)
```

### Subgraphs


```@docs
nodes(subgraph::AbstractSubGraph)
has_vertex(layer::Layer, mv::MultilayerVertex)
has_vertex(interlayer::Interlayer, mv::MultilayerVertex)
nv(subgraph::AbstractSubGraph)
mv_vertices(subgraph::AbstractSubGraph)
mv_inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
mv_outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
mv_neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
has_edge(subgraph::AbstractSubGraph,me::MultilayerEdge)
has_edge( subgraph::AbstractSubGraph, s::MultilayerVertex, d::MultilayerVertex)
ne(subgraph::AbstractSubGraph)
edges(subgraph::S) where {T,U,S<:AbstractSubGraph{T,U}} 
add_edge!( subgraph::S, me::E) where {T,U<:Real,S<:AbstractSubGraph{T,U},E<:MultilayerEdge{ <: Union{U, Nothing}}}
rem_edge!(subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex)
rem_edge!(subgraph::AbstractSubGraph, me::MultilayerEdge)
get_metadata(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex)
get_metadata(subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex)
get_weight(subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex) 
is_directed(subgraph::AbstractSubGraph)
is_directed(::Type{S}) where {T,U,G,S <: AbstractSubGraph{T,U,G}}
adjacency_matrix(subgraph::AbstractSubGraph)
MultilayerGraphs.weights(subgraph::S) where {T,U,S<:AbstractSubGraph{T,U}}
name(subgraph::AbstractSubGraph)
Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}}
Layer(name::Symbol, vertices::Vector{<: MultilayerVertex}, edge_list::Vector{ <: MultilayerEdge}, null_graph::G, weighttype::Type{U};  default_vertex_metadata::Function = mv -> NamedTuple(), default_edge_weight::Function = (src, dst) -> one(U), default_edge_metadata::Function = (src, dst) -> NamedTuple()) where {T <: Integer, U <: Real,  G <: AbstractGraph{T}}

Layer(
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

has_node(layer::Layer, n::Node)
rem_vertex!(layer::Layer, mv::MultilayerVertex)
rem_vertex!(layer::Layer, n::Node)

Interlayer{T<:Integer,U<:Real,G<:AbstractGraph{T}}

Interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G,
    edge_list::Vector{ <: MultilayerEdge{<: Union{U, Nothing}}};
    default_edge_weight::Function = (x,y) -> nothing,
    default_edge_metadata::Function = (x,y) -> NamedTuple(),
    transfer_vertex_metadata::Bool = false,
    name::Symbol = Symbol("interlayer_$(layer_1.name)_$(layer_2.name)")
) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}


Interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64,
    null_graph::G;
    default_edge_weight::Function = (x,y) -> nothing,
    default_edge_metadata::Function = (x,y) -> NamedTuple(),
    name::Symbol = Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
    transfer_vertex_metadata::Bool = false
) where {T<:Integer, U <: Union{Nothing, <: Real},  G<:AbstractGraph{T}}


multiplex_interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G;
    default_edge_weight::Function = (x,y) -> nothing,
    default_edge_metadata::Function = (x,y) -> NamedTuple(),
    transfer_vertex_metadata::Bool = false,
    name::Symbol = Symbol("interlayer_$(layer_1.name)_$(layer_2.name)")
) where {T<:Integer, U <: Real, G<:AbstractGraph{T}} =  _multiplex_interlayer(collect(mv_vertices(layer_1)), collect(mv_vertices(layer_2)),  null_graph, U; default_edge_weight = default_edge_weight, default_edge_metadata = default_edge_metadata, transfer_vertex_metadata = transfer_vertex_metadata , name = name)


empty_interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G;
    default_edge_weight::Function = (x,y) -> nothing,
    default_edge_metadata::Function = (x,y) -> NamedTuple(),
    name::Symbol = Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
    transfer_vertex_metadata::Bool = false
) where {T<:Integer, U <: Real, G<:AbstractGraph{T}} 

is_multiplex_interlayer(interlayer::Interlayer)

get_symmetric_interlayer(
    interlayer::In;
    symmetric_interlayer_name::String = String(interlayer.name) * "_rev"
) where {T,U,G,In<:Interlayer{T,U,G}}


```

### Multilayer-specific methods


```@docs
nodes(mg::MultilayerGraph)
has_vertex(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
nv(mg::M) where {M <: AbstractMultilayerGraph }
mv_vertices(mg::AbstractMultilayerGraph)
mv_inneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
mv_outneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
mv_neighbors( mg::AbstractMultilayerGraph, mv::MultilayerVertex)
has_edge(mg::AbstractMultilayerGraph, edge::MultilayerEdge) 
has_edge( subgraph::AbstractMultilayerGraph, s::MultilayerVertex, d::MultilayerVertex)
ne(mg::AbstractMultilayerGraph)
edges(mg::AbstractMultilayerUGraph)
edges(mg::M) where {T,U,M<:AbstractMultilayerUGraph{T,U}}
edges(mg::M) where {T,U,M<:AbstractMultilayerDiGraph{T,U}}
add_edge!(mg::M, src::V, dst::V; weight::Union{Nothing, U} = one(U), metadata::Union{Tuple,NamedTuple} = NamedTuple() ) where {T,U, M <: AbstractMultilayerGraph{T,U}, V <: MultilayerVertex}
add_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerUGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}
add_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerDiGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}
rem_edge!(mg::AbstractMultilayerGraph, me::MultilayerEdge)
rem_edge!(mg::MultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)
rem_edge!(mg::MultilayerDiGraph, src::MultilayerVertex, dst::MultilayerVertex)
get_metadata(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
get_metadata(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)
get_weight(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex
is_directed(mg::AbstractMultilayerUGraph)
is_directed(m::M) where { M <: Type{ <: AbstractMultilayerUGraph}}
is_directed(mg::AbstractMultilayerDiGraph)
is_directed(m::M) where { M <: Type{ <: AbstractMultilayerDiGraph}}
has_node(mg::AbstractMultilayerGraph, n::Node)
rem_vertex!(mg::AbstractMultilayerUGraph, V::MultilayerVertex)
rem_vertex!(mg::AbstractMultilayerDiGraph, V::MultilayerVertex)
```

----------------------------------------------------
####################################################
----------------------------------------------------


## Developer

### Nodes

```@docs
AbstractNode
```



### Vertices

```@docs
AbstractVertex
AbstractMultilayerVertex
```


### Edges


```@docs
AbstractMultilayerEdge
metadata(he::MultilayerGraphs.HalfEdge)
weight(he::MultilayerGraphs.HalfEdge)
```

# Subgraphs

```@docs
has_vertex(subgraph::S, v::T ) where {T,S<:AbstractSubGraph{T}}
vertices(subgraph::AbstractSubGraph)
inneighbors(subgraph::S, v::T) where {T, S <: AbstractSubGraph{T}}
inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
mv_inneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}
outneighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
neighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
edgetype(::S) where {T,U,S<:AbstractSubGraph{T,U}}
has_edge(subgraph::S, s::MultilayerVertex, d::MultilayerVertex) where { T, S <: AbstractSubGraph{T}}
add_edge!(subgraph::S, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U<: Real, W<:Union{ U, Nothing},G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} 
rem_edge!(subgraph::S, src::T, dst::T) where {T, S<:AbstractSubGraph{T}}
AbstractLayer
Layer(descriptor::MultilayerGraphs.LayerDescriptor{T}, vertices::Vector{<: MultilayerVertex}, edge_list::Vector{<:MultilayerEdge}) where {T <: Integer}
rem_vertex!(layer::L, v::T) where {T, L <: Layer{T}} 
AbstractInterlayer
```

### Traits
```@docs
```

### Multilayer-specific methods

```@docs
has_vertex(mg::M, v::T) where {T, M <: AbstractMultilayerGraph{T}}
vertices(mg::AbstractMultilayerGraph)
inneighbors(mg::AbstractMultilayerUGraph, v::T) where {T <: Integer, M <: AbstractMultilayerUGraph{T}}
inneighbors(mg::M, v::T) where {T <: Integer, M <: AbstractMultilayerDiGraph{T}}
inneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
outneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}
outneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}
outneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}
neighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
edgetype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}
has_edge(mg::M, src::T, dst::T) where { T, M <: AbstractMultilayerUGraph{T}}
has_edge(mg::M, src::T, dst::T) where { T, M <: AbstractMultilayerDiGraph{T}}
rem_edge!(mg::M, src::T, dst::T) where {T, M <: AbstractMultilayerGraph{T}
```




```@docs
```
