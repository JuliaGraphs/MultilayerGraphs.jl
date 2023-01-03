# API

This page provides a list of exported methods organized by topic and audience. Methods that act on vertices, edges, and layers are grouped together. Some methods are intended for developers who want to use the `Graphs.jl` library as part of their code, while others are meant for end-users.

## End-User

### [Nodes](@id nodes_eu)

```@docs
Node
id
```

### [Vertices](@id vertices_eu)

```@docs
eltype
MultilayerVertex
MV
node
layer
metadata(mv::MultilayerVertex)
MissingVertex
```

### [Edges](@id edges_eu)

```@docs
MultilayerEdge
ME
weight(e::AbstractMultilayerEdge)
metadata(e::AbstractMultilayerEdge)
```

### [Subgraphs](@id subgraphs_eu)

```@docs
Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}}
Layer(
    name::Symbol, 
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}, 
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}}, 
    null_graph::G, 
    weighttype::Type{U}; 
    default_vertex_metadata::Function = mv -> NamedTuple(), 
    default_edge_weight::Function = (src, dst) -> one(U), 
    default_edge_metadata::Function = (src, dst) -> NamedTuple()
) where {T <: Integer, U <: Real, G <: AbstractGraph{T}}

Layer(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64,
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function = mv -> NamedTuple(),
    default_edge_weight::Function = (src, dst) -> nothing,
    default_edge_metadata::Function = (src, dst) -> NamedTuple(),
    allow_self_loops::Bool = false
) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

layer_simplegraph
layer_simpledigraph
layer_simpleweightedgraph
layer_simpleweighteddigraph
layer_metadigraph
layer_valgraph
layer_valoutdigraph
layer_valdigraph
layer_metagraph

has_node(layer::Layer, n::Node)
add_vertex!(layer::Layer, mv::MultilayerVertex)
add_vertex!(layer::L, n::Node, args...; kwargs...) where {T, L <: Layer{T}} 
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

interlayer_simplegraph
interlayer_simpleweightedgraph
interlayer_metagraph
interlayer_valgraph
interlayer_simpledigraph
interlayer_simpleweighteddigraph
interlayer_metadigraph
interlayer_valoutdigraph
interlayer_valdigraph


multiplex_interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G;
    default_edge_weight::Function = (x,y) -> nothing,
    default_edge_metadata::Function = (x,y) -> NamedTuple(),
    transfer_vertex_metadata::Bool = false,
    name::Symbol = Symbol("interlayer_$(layer_1.name)_$(layer_2.name)")
) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}


empty_interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G;
    default_edge_weight::Function = (x,y) -> nothing,
    default_edge_metadata::Function = (x,y) -> NamedTuple(),
    name::Symbol = Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
    transfer_vertex_metadata::Bool = false
) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

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
has_edge( layer::Layer, s::MultilayerVertex{nothing}, d::MultilayerVertex{nothing})
ne(subgraph::AbstractSubGraph)
edges(subgraph::S) where {T,U,S<:AbstractSubGraph{T,U}} 
add_edge!( subgraph::S, me::E) where {T,U<:Real,S<:AbstractSubGraph{T,U},E<:MultilayerEdge{ <: Union{U, Nothing}}}
add_edge!(layer::Layer, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...)
add_edge!(interlayer::Interlayer, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...)
rem_edge!(subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex)
rem_edge!(subgraph::AbstractSubGraph, me::MultilayerEdge)
get_metadata(subgraph::AbstractSubGraph, bare_mv::MultilayerVertex)
get_metadata(subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex)
get_weight(subgraph::AbstractSubGraph, src::MultilayerVertex, dst::MultilayerVertex) 
is_directed(subgraph::AbstractSubGraph)
is_directed(::Type{S}) where {T,U,G,S <: AbstractSubGraph{T,U,G}}
adjacency_matrix(subgraph::AbstractSubGraph)
MultilayerGraphs.weights(subgraph::S) where {T,U,S<:AbstractSubGraph{T,U}}

is_multiplex_interlayer(interlayer::Interlayer)

get_symmetric_interlayer(
    interlayer::In;
    symmetric_interlayer_name::String = String(interlayer.name) * "_rev"
) where {T,U,G,In<:Interlayer{T,U,G}}

name(subgraph::AbstractSubGraph)
graph(subgraph::AbstractSubGraph)
```

### [Multilayer-Specific Methods](@id msm_eu)

```@docs
MultilayerGraph{T,U}

MultilayerDiGraph{T,U}

MultilayerGraph(T::Type{<:Number}, U::Type{<:Number})

MultilayerDiGraph(T::Type{<:Number}, U::Type{<:Number})

MultilayerGraph(
    layers::Vector{<:Layer{T,U}},
    specified_interlayers::Vector{<:Interlayer{T,U}};
    default_interlayers_null_graph::H = SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U, H <: AbstractGraph{T}}

MultilayerDiGraph(
    layers::Vector{<:Layer{T,U}},
    specified_interlayers::Vector{<:Interlayer{T,U}};
    default_interlayers_null_graph::H = SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U, H <: AbstractGraph{T}}

MultilayerGraph(
    empty_layers::Vector{<:Layer{T,U}},
    empty_interlayers::Vector{<:Interlayer{T,U}},
    degree_distribution::UnivariateDistribution;
    allow_self_loops::Bool = false,
    default_interlayers_null_graph::H = SimpleGraph{T}(),
) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

MultilayerDiGraph(
    empty_layers::Vector{<:Layer{T,U}},
    empty_interlayers::Vector{<:Interlayer{T,U}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution;
    allow_self_loops::Bool = false,
    default_interlayers_null_graph::H = SimpleGraph{T}(),
) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

MultilayerGraph(
    empty_multilayergraph::MultilayerGraph{T,U}, 
    degree_sequence::Vector{<:Integer}; 
    allow_self_loops::Bool = false, 
    perform_checks::Bool = true
) where {T,U}

MultilayerDiGraph(
    empty_multilayerdigraph::MultilayerDiGraph{T,U}, 
    indegree_sequence::Vector{<:Integer},
    outdegree_sequence::Vector{<:Integer};
    allow_self_loops::Bool = false,
     perform_checks::Bool = false
) where {T,U}

nodes(mg::AbstractMultilayerGraph)
nn(mg::AbstractMultilayerGraph) 
add_node!(mg::AbstractMultilayerGraph, n::Node)
rem_node!(mg::AbstractMultilayerGraph, n::Node)
has_vertex(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
nv(mg::M) where {M <: AbstractMultilayerGraph }
mv_vertices(mg::AbstractMultilayerGraph)
mv_inneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
mv_outneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
mv_neighbors( mg::AbstractMultilayerGraph, mv::MultilayerVertex)
add_vertex!(mg::AbstractMultilayerGraph, mv::MultilayerVertex; add_node::Bool = true)
rem_vertex!(mg::AbstractMultilayerUGraph, V::MultilayerVertex)
rem_vertex!(mg::AbstractMultilayerDiGraph, V::MultilayerVertex)
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
get_weight(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)
set_weight!(mg::M, src::MultilayerVertex, dst::MultilayerVertex, weight::U) where { T,U, M <: AbstractMultilayerUGraph{T,U}}
set_weight!(mg::M, src::MultilayerVertex, dst::MultilayerVertex, weight::U) where { T,U, M <: AbstractMultilayerDiGraph{T,U}}
is_directed(mg::AbstractMultilayerUGraph)
is_directed(m::M) where { M <: Type{ <: AbstractMultilayerUGraph}}
is_directed(mg::AbstractMultilayerDiGraph)
is_directed(m::M) where { M <: Type{ <: AbstractMultilayerDiGraph}}
has_node(mg::AbstractMultilayerGraph, n::Node)
set_metadata!(mg::AbstractMultilayerGraph, mv::MultilayerVertex, metadata::Union{Tuple, NamedTuple}) 
set_metadata!(mg::AbstractMultilayerDiGraph, src::MultilayerVertex, dst::MultilayerVertex, metadata::Union{Tuple, NamedTuple})
set_metadata!(mg::AbstractMultilayerUGraph, src::MultilayerVertex, dst::MultilayerVertex, metadata::Union{Tuple, NamedTuple})
nl(mg::AbstractMultilayerGraph)
nIn(mg::AbstractMultilayerGraph)
has_layer(mg::AbstractMultilayerGraph, layer_name::Symbol)

add_layer!(mg::M, new_layer::L; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String ="multiplex"
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerUGraph{T,U},L<:Layer{T,U,G}, H <: AbstractGraph{T}}

add_layer!(
    mg::M, new_layer::L; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String ="multiplex"
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerDiGraph{T,U},L<:Layer{T,U,G}, H <: AbstractGraph{T}}

specify_interlayer!(
    mg::M,
    new_interlayer::In
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerDiGraph{T,U},In<:Interlayer{T,U,G}}

specify_interlayer!(
    mg::M,
    new_interlayer::In
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerUGraph{T,U},In<:Interlayer{T,U,G}}

get_interlayer(
    mg::AbstractMultilayerGraph, layer_1_name::Symbol, layer_2_name::Symbol
)


indegree( mg::AbstractMultilayerGraph, v::MultilayerVertex) 
indegree(mg::AbstractMultilayerGraph, vs::AbstractVector{<:MultilayerVertex}=vertices(mg))

outdegree(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
outdegree(mg::AbstractMultilayerGraph, vs::AbstractVector{<:MultilayerVertex}=vertices(mg))

degree(mg::AbstractMultilayerGraph, vs::AbstractVector{<:MultilayerVertex}=vertices(mg)) 

mean_degree(mg::AbstractMultilayerGraph)

degree_second_moment(mg::AbstractMultilayerGraph) 

degree_variance(mg::AbstractMultilayerGraph)

MultilayerGraphs.weighttype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}

multilayer_global_clustering_coefficient(
    mg::AbstractMultilayerGraph, norm_factor::Union{Float64,Symbol}=:max
)

overlay_clustering_coefficient(
    mg::AbstractMultilayerGraph,
    norm_factor::Union{Float64,Symbol}=:max
)

eigenvector_centrality(
    mg::M; weighted::Bool = true,  norm::String="1", tol::Float64=1e-6, maxiter::Int64=2000
) where {T,U,M<:AbstractMultilayerGraph{T,U}}

modularity(
    mg::M, c::Matrix{Int64}; null_model::Union{String,Array{U,4}}="degree"
) where {T,U,M<:AbstractMultilayerGraph{T,U}}


von_neumann_entropy(mg::M) where {T,U,M<:AbstractMultilayerUGraph{T,U}}
```

### [Representations](@id representations_eu)
```@docs
array(atr::AbstractTensorRepresentation)
WeightTensor{U}
weight_tensor(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}
MetadataTensor{U}
metadata_tensor(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}
array(amr::AbstractMatrixRepresentation)
SupraWeightMatrix{T,U}
supra_weight_matrix(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}
```

### [Traits](@id traits_eu)
```@docs
is_weighted(g::G) where { G <: AbstractGraph}
is_weighted(g::G) where {G<:Type{<:AbstractGraph}}

is_meta(g::G) where {G <: AbstractGraph}
is_meta(g::G) where {G<:Type{<:AbstractGraph}}
```

### [Utilities](@id utilities_eu)
```@docs
multilayer_kronecker_delta(dims::NTuple{4,Int64})
δk{T}
δk(N::Int64)
δ_1{T<: Number}
δ_2{T<:Number}
δ_3{T<:Number}
δ_Ω{T}
```


## Developer

### [Nodes](@id nodes_dev)

```@docs
AbstractNode
```



### [Vertices](@id vertices_dev)

```@docs
AbstractVertex
AbstractMultilayerVertex
```


### [Edges](@id edges_dev)


```@docs
AbstractMultilayerEdge
metadata(he::MultilayerGraphs.HalfEdge)
weight(he::MultilayerGraphs.HalfEdge)
```

### [Subgraphs](@id subgraphs_dev)

```@docs
has_vertex(subgraph::S, v::T ) where {T,S<:AbstractSubGraph{T}}
vertices(subgraph::AbstractSubGraph)
inneighbors(subgraph::S, v::T) where {T, S <: AbstractSubGraph{T}}
inneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
outneighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
outneighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
neighbors(subgraph::S, v::T) where {T,S<:AbstractSubGraph{T}}
neighbors(subgraph::AbstractSubGraph, mv::MultilayerVertex)
edgetype(::S) where {T,U,S<:AbstractSubGraph{T,U}}
has_edge( subgraph::S, s::T, d::T) where {T,S<:AbstractSubGraph{T}}
add_edge!(subgraph::S, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U<: Real, W<:Union{ U, Nothing},G<:AbstractGraph{T},S<:AbstractSubGraph{T,U,G}} 
rem_edge!(subgraph::S, src::T, dst::T) where {T, S<:AbstractSubGraph{T}}
AbstractLayer
Layer(descriptor::MultilayerGraphs.LayerDescriptor{T}, vertices::Vector{<: MultilayerVertex}, edge_list::Vector{<:MultilayerEdge}) where {T <: Integer}
rem_vertex!(layer::L, v::T) where {T, L <: Layer{T}} 
AbstractInterlayer
```


### [Multilayer-Specific Methods](@id msm_dev)

```@docs
AbstractMultilayerGraph{T <: Integer, U <: Real}
has_vertex(mg::M, v::T) where {T, M <: AbstractMultilayerGraph{T}}
vertices(mg::AbstractMultilayerGraph)
inneighbors(mg::M, v::T) where {T,M<:AbstractMultilayerUGraph{T,<:Real}}
inneighbors(mg::M, v::T) where {T, M<:AbstractMultilayerGraph{T,<:Real}}
inneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
outneighbors(mg::M, v::T) where {T, M<:AbstractMultilayerGraph{T}}
neighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
edgetype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}
has_edge(mg::M, src::T, dst::T) where { T, M <: AbstractMultilayerUGraph{T}}
has_edge(mg::M, src::T, dst::T) where { T, M <: AbstractMultilayerDiGraph{T}}
rem_edge!(mg::M, src::T, dst::T) where {T, M <: AbstractMultilayerGraph{T}
AbstractMultilayerUGraph{T,U}
AbstractMultilayerDiGraph{T,U}
```

### [Representations](@id representations_dev)
```@docs
AbstractTensorRepresentation{U}
AbstractMatrixRepresentation{T,U}
```

### [Traits](@id traits_dev)
```@docs
IsWeighted{X}
IsMeta{X}
```