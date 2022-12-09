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
nodes(::AbstractSubGraph)
has_vertex(layer::L, v::MultilayerVertex) where { T,U,G, L <: Layer{T,U,G}}
has_vertex(interlayer::In, v::MultilayerVertex) where { T,U,G, In <: Interlayer{T,U,G}}
nv(subgraph::S) where { S <: AbstractSubGraph}
mv_vertices(subgraph::S) where {S <: AbstractSubGraph{ <: Integer, <: AbstractSimpleGraph}}
```

### Multilayer-specific methods


```@docs
nodes(::MultilayerGraph)
has_vertex(mg::M, mv::MultilayerVertex) where {T,U, M <: AbstractMultilayerGraph{T,U}}
nv(mg::M) where {M <: AbstractMultilayerGraph }
 mv_vertices(mg::AbstractMultilayerGraph)
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
metadata(he::HalfEdge)
weight(he::HalfEdge)
```

# Subgraphs

```@docs
has_vertex(subgraph::S, v::T ) where {T,U,G,S<:AbstractSubGraph{T,U,G}}
vertices(subgraph::S) where {S <: AbstractSubGraph{ <: Integer, <: AbstractSimpleGraph}}
inneighbors(subgraph::S, v::T) where {T,U,G, S <: AbstractSubGraph{T,U,G}}
inneighbors(subgraph::S, mv::MultilayerVertex)  where {T,U,G, S <: AbstractSubGraph{T,U,G}}
```

### Traits
```@docs
```

### Multilayer-specific methods

```@docs
has_vertex(mg::M, v::T) where {T,U, M <: AbstractMultilayerGraph{T,U}}
vertices(mg::M) where {M<:AbstractMultilayerGraph}
inneighbors(mg::M, v::T) where {M <: AbstractMultilayerUGraph{T} } where { T <: Integer}
inneighbors(mg::M, v::T) where {M <: AbstractMultilayerDirDiGraph{T} } where { T <: Integer}
inneighbors( mg::M, mv::V ) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}
```




```@docs
```
