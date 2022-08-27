```@meta
CurrentModule = MultilayerGraphs
```

```@raw html
<div style="width:100%; height:150px;border-width:4px;border-style:solid;padding-top:25px;
        border-color:#000;border-radius:10px;text-align:center;background-color:#B3D8FF;
        color:#000">
    <h3 style="color: black;">Star us on GitHub!</h3>
    <a class="github-button" href="https://github.com/JuliaGraphs/MultilayerGraphs.jl" data-icon="octicon-star" data-size="large" data-show-count="true" aria-label="Star JuliaGraphs/MultilayerGraphs.jl on GitHub" style="margin:auto">Star</a>
    <script async defer src="https://buttons.github.io/buttons.js"></script>
</div>
```

# MultilayerGraphs.jl

**MultilayerGraphs.jl** is a Julia package for the construction, manipulation and analysis of multilayer graphs [extending Graphs.jl](https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/).

## Overview

**MultilayerGraphs.jl** implements the mathematical formulation of multilayer graphs proposed by [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). It mainly revolves around two custom types, [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref), encoding undirected and directed multilayer graphs respectively.

Roughly speaking, a multilayer graph is a collection of *layers*, i.e. graphs whose vertices are representations of the same set of nodes, and *interlayers*, i.e the [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) whose vertices are those of any two layers and whose edges are those between vertices of the same two layers. See below for the distinction between ***nodes*** and ***vertices***.

[`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) are fully-fledged [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) extensions. Both structs are designed so that their layers and interlayers can be of any type (as long as they are Graphs.jl extensions themselves) and they need not be all of the same type. It is anyway required that all layers and interlayers of [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) are respectively undirected and directed. Directedness is checked via the `IsDirected` trait defined in Graphs.jl adopting [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl). Since the layers' and interlayers' graph types don't need to be the same, multilayer graph types are considered weighted graphs by default, and thus are assigned the trait `IsWeighted`.

## Installation

Press `]` in the Julia REPL and then

```julia
pkg> add MultilayerGraphs
```

## Tutorial

Here we illustrate how to define, handle and analyse a [`MultilayerGraph`](@ref) (the directed version is completely analogous).

### Layers and Interlayers

Let's import some necessary packages

```julia
# Import necessary dependencies
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs
using MultilayerGraphs
```

We define some methods and constants that will prove useful later in the tutorial

```julia
# Set the number of nodes, minimum and maximum number of edges for random graphs
const n_nodes   = 5
const min_edges = n_nodes
const max_edges = 10
```

As said before, to define a multilayer graph we need to specify its layers and interlayers. We proceed by constructing a layer (see [`Layer`](@ref))

```julia
# Construct a layer
layer = Layer(:layer_1, SimpleGraph(n_nodes, rand(min_edges:max_edges)); U = Float64)
```

A `Layer` has a name (here `:layer_1`), an underlying graph (`SimpleGraph(n_nodes, rand(min_edges:max_edges))`) and a weight matrix `eltype` `U` (it defaults to the adjacency matrix's `eltype` if the graph is unweighted). To correctly specify a multilayer graph all layers and interlayers must have the same `U`, otherwise the multilayers's adjacency tensor would be poorly specified.

Notice that `U` does not need to coincide with the `eltype` of the adjacency matrix of the underlying graph: as far as we know, there is no way to set it explicitly for all Graphs.jl extensions, nor it is required for extensions to implement such feature, so our package converts to `U` the `eltype` of `Layer`s and `Interlayer`s weight (/adjacency) matrices every time they are invoked

```julia
adjacency_matrix(layer)
```
```nothing
5×5 SparseMatrixCSC{Float64, Int64} with 12 stored entries:
  ⋅    ⋅   1.0   ⋅    ⋅ 
  ⋅    ⋅   1.0  1.0   ⋅ 
 1.0  1.0   ⋅   1.0  1.0
  ⋅   1.0  1.0   ⋅   1.0
  ⋅    ⋅   1.0  1.0   ⋅ 
```

We may define more `Layer`s for future use

```julia
layers = [Layer(n_nodes, :layer_1, SimpleGraph{Int64}, rand(min_edges:max_edges); U = Float64), 
          Layer(n_nodes, :layer_2, SimpleWeightedGraph{Int64}, rand(min_edges:max_edges); U = Float64),
          Layer(n_nodes, :layer_3, MetaGraph{Int64, Float64}, rand(min_edges:max_edges); U = Float64),
          Layer(:layer_4, ValGraph( SimpleGraph{Int64}(5, rand(min_edges:max_edges)), edgeval_types=(Int64, ), edgeval_init=(s, d) -> (s + d, ), vertexval_types=(String, ), vertexval_init=undef); U = Float64)
]
```

There are other constructors for the `Layer` struct you may want to consult via `?Layer`.

We similarly define an interlayer (see [`Interlayer`](@ref))

```julia
interlayer = Interlayer(n_nodes, :interlayer_layer_1_layer_2 , :layer_1, :layer_2, SimpleGraph{Int64}, rand(min_edges:max_edges); U = Float64)
```

Here we used a constructor that returns a random `Interlayer`. Its arguments are the number of nodes `n_nodes`, the name of the `Interlayer` `interlayer_layer_1_layer_2`, the name of the `Layers` that it connects (`:layer_1` and `:layer_2`), the underlying graph type `SimpleGraph{Int64}`, and the number of edges `rand(min_edges:max_edges)` and the weight/adjacency matrix `eltype` again needs to be specified (although it may be left blank and the constructor will default to `eltype(adjacency_matrix(SimpleGraph{Int64}))`).

The adjacency matrix of an `Interlayer` is that of a bipartite graph

```julia
adjacency_matrix(interlayer)
```
```nothing
10×10 SparseMatrixCSC{Float64, Int64} with 14 stored entries:
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅ 
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0  1.0
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅   1.0   ⋅   1.0
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0
  ⋅    ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
 1.0   ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅   1.0   ⋅   1.0  1.0   ⋅    ⋅    ⋅    ⋅    ⋅ 
```

It is a 4 block matrix where the first `n_nodes` rows and columns refer to `:layer_1`'s vertices, while the last `n_nodes` rows and columns refer to `:layer_2`'s vertices.

We may define more `Interlayer`s for future use:

```julia
interlayers = [
                Interlayer(n_nodes, :interlayer_layer_1_layer_2, :layer_1, :layer_2, SimpleGraph{Int64}, rand(min_edges:max_edges); U = Float64),
                Interlayer(n_nodes, :interlayer_layer_1_layer_3, :layer_1, :layer_3, SimpleWeightedGraph{Int64,Float64}, rand(min_edges:max_edges); U = Float64),
                Interlayer(n_nodes, :interlayer_layer_1_layer_4, :layer_2, :layer_3, MetaGraph{Int64, Float64}, rand(min_edges:max_edges); U = Float64),
                
]
```

There are other constructors for the `Interlayer` struct you may want to consult via `?Interlayer`.

`Layer`s and `Interlayer`s and complete extensions of Graphs.jl, so all methods in Graphs.jl should just work. The explicitly extended methods are [`edges`](@ref), [`eltype`](@ref), [`edgetype`](@ref), [`has_edge`](@ref), [`has_vertex`](@ref), [`inneighbors`](@ref), [`ne`](@ref), [`nv`](@ref), [`outneighbors`](@ref), [`vertices`](@ref), [`is_directed`](@ref), [`add_edge!`](@ref), [`rem_edge!`](@ref).

### Instantiation and Handling of `MultilayerGraph`

We can define a `MultilayerGraph` by specifying its layers and interlayers:

```julia
multilayergraph = MultilayerGraph(layers, interlayers; default_interlayer = "multiplex")
```
```nothing
MultilayerGraph{Int64, Float64}([0.0 1.0 … 1.0 1.0; 1.0 0.0 … 0.0 1.0; … ; 1.0 0.0 … 0.0 0.0; 1.0 1.0 … 0.0 0.0;;; 0.0 0.0 … 0.0 0.0; 0.0 1.0 …],...)
```

It is not important that the `interlayers` array contains all the interlayers needed to specify the multilayer graph: the unspecified interlayers are automatically generated according to the `default_interlayer` argument. Right now only the `"multiplex"` value is supported, and will generate interlayers that have edges between pair of vertices of the two layers that represent the same node.

Notice that in the output the signature of `MultilayerGraph` has two parametric types, namely `MultilayerGraph{Int64, Float64}`. The first is referred to the node type, that just as in every `Graphs.jl` extension it is a subtype of `Integer`. The second parameter is instead the `eltype` of the equivalent of the adjacency matrix for multilayer graphs: the adjacency tensor (see below for more).

You may also specify random `MultilayerGraph`

```julia
random_multilayergraph = MultilayerGraph( 3,         # Number of layers
                                          n_nodes,   # Number of nodes
                                          min_edges, # Minimum number of edges in each layer/interlayer
                                          max_edges, # Maximum number of edges in each layer/interlayer
                                          [
                                             SimpleGraph{Int64},
                                             SimpleWeightedGraph{Int64, Float64},
                                          ]          # The set of graph types to random draw from when constructing layers and interlayers
                                        )
```

Alternatively, one may add layers and interlayers calling the [`add_layer!`](@ref) and [`specify_interlayer!`](@ref) functions respectively, perhaps starting with an empty multilayer graph created with the constructor:

```julia
empty_multilayergraph  = MultilayerGraph( n_nodes, # Number of nodes
                                          Int64,   # Node type
                                          Float64  # Adjacency tensor's eltype, see below
                                        )
```

Let's explore some properties of the `MultilayerGraph` struct.

#### Layers

It is an `OrderedDict` where the keys are the layers' indexes within the multilayer graph (the index is repeated twice in a tuple to be consistent with `multilayergraph.interlayers`' keys, see below) and the values are the actual `Layer`s.

```julia
multilayergraph.layers
```
```nothing
OrderedCollections.OrderedDict{Tuple{Int64, Int64}, Layer{Int64, U, G} where {U<:Real, G<:AbstractGraph{Int64}}} with 4 entries:
  (1, 1) => Layer{Int64, Float64, SimpleGraph{Int64}}(:layer_1, SimpleGraph{Int…
  (2, 2) => Layer{Int64, Float64, SimpleWeightedGraph{Int64, Float64}}(:layer_2…
  (3, 3) => Layer{Int64, Float64, MetaGraph{Int64, Float64}}(:layer_3, {5, 10} …
  (4, 4) => Layer{Int64, Float64, ValGraph{Int64, Tuple{String}, Tuple{Int64}, …
```

#### Interlayers

It is an `OrderedDict` where each key is the pair of indexes of the layers that the corresponding value, i.e. the interlayer, connects within the multilayer graph

```julia
multilayergraph.interlayers
```
```nothing
OrderedCollections.OrderedDict{Tuple{Int64, Int64}, Interlayer{Int64, U, G} where {U<:Real, G<:AbstractGraph{Int64}}} with 12 entries:
  (2, 1) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_2_…
  (1, 2) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_1_…
  (3, 1) => Interlayer{Int64, Float64, SimpleWeightedGraph{Int64, Float64}}(:in…
  (1, 3) => Interlayer{Int64, Float64, SimpleWeightedGraph{Int64, Float64}}(:in…
  (3, 2) => Interlayer{Int64, Float64, MetaGraph{Int64, Float64}}(:interlayer_l…
  (2, 3) => Interlayer{Int64, Float64, MetaGraph{Int64, Float64}}(:interlayer_l…
  (4, 1) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_4_…
  (1, 4) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_1_…
  (4, 2) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_4_…
  (2, 4) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_2_…
  (4, 3) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_4_…
  (3, 4) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_3_…
```

Note that the (1,2) interlayer (i.e. the interlayer between layer (1,1) and layer (2,2)) is very similar to interlayer (2,1), but not identical: its adjacency matrix rows and columns are reordered. One may get interlayer (2,1) from interlayer (1,2) (i.e. one may get the *symmetric* interlayer of (1,2)) as follows

```julia
symmetric_interlayer = get_symmetric_interlayer(multilayergraph.interlayers[(1,2)])
symmetric_interlayer == multilayergraph.interlayers[(2,1)]
```
```nothing
true
```

You may access individual layers and interlayers with the "dot" notation:

```julia
multilayergraph.layer_1
```
```nothing
Layer{Int64, Float64, SimpleGraph{Int64}}(:layer_1, SimpleGraph{Int64}(5, [[2, 3, 5], [1, 3], [1, 2], [5], [1, 4]]), MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVertex{Int64}}[])
```

#### Adjacency Tensor

The adjacency tensor is a 4-dimensional array

```julia
multilayergraph.adjacency_tensor
```
```nothing
5×5×4×4 Array{Float64, 4}:
...
```

To understand its indexing, consider the following example

```julia
multilayergraph.adjacency_tensor[1,5,2,3]
```
```nothing
0.0
```

This means that there is an edge of zero weight between the vertex representing node 1 in layer 5 and the vertex representing node 2 in layer 3. It is a good time to note the difference between *nodes* and *vertices*. In the context of multilayer graphs, the vertices of every layer and interlayer represent the same set of nodes. That is, vertex 1 in layer (1,1) represents the same node as vertex 1 in layer (2,2) and so on. To make this distinction clearer the package implements the `MultilayerVertex` type, that represents vertices within the multilayer graph. The implementation of `MultilayerVertex` is

```julia
struct MultilayerVertex{T <: Integer} <: AbstractMultilayerVertex{T}
    node::T        # The node  the vertex represents
    layer::Symbol  # The layer the vertex belongs to
end
```

To get the vertices of a `Layer` or an `Interlayer`, one may use the Graphs.jl APIs

```julia
vertices(multilayergraph.layers[(1,1)])
```
```nothing
5-element Vector{MultilayerVertex{Int64}}:
 MultilayerVertex{Int64}(1, :layer_1)
 MultilayerVertex{Int64}(2, :layer_1)
 MultilayerVertex{Int64}(3, :layer_1)
 MultilayerVertex{Int64}(4, :layer_1)
 MultilayerVertex{Int64}(5, :layer_1)
```

```julia
vertices(multilayergraph.interlayers[(1,2)])
```
```nothing
10-element Vector{Any}:
 MultilayerVertex{Int64}(1, :layer_1)
 MultilayerVertex{Int64}(2, :layer_1)
 MultilayerVertex{Int64}(3, :layer_1)
 MultilayerVertex{Int64}(4, :layer_1)
 MultilayerVertex{Int64}(5, :layer_1)
 MultilayerVertex{Int64}(1, :layer_2)
 MultilayerVertex{Int64}(2, :layer_2)
 MultilayerVertex{Int64}(3, :layer_2)
 MultilayerVertex{Int64}(4, :layer_2)
 MultilayerVertex{Int64}(5, :layer_2)
```

To get a specific layer of a multilayer graph from its name, one may also write

```julia
layer_1 = get_layer(multilayergraph, :layer_1)
```

Same for interlayers

```julia
interlayer_2_1 = get_interlayer(multilayergraph, :layer_2, :layer_1)
```

Both `MultilayerGraph` and `MultilayerDiGraph` fully extend `Graphs.jl`, so they have access to Graphs.jl API as one would expect, just keeping in mind that vertices are `MultilayerVertex`s and not subtypes of `Integer` (`MultilayerVertex` is actually a subtype of `AbstractVertex` that this package defines, see [Future Developments](#Future-Developments)), and that edges are `MultilayerEdge`s, which subtype `AbstractEdge`.

Some notable examples are

```julia
edges(multilayergraph)
```
```nothing
73-element Vector{MultilayerEdge}:
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_1), MultilayerVertex{Int64}(2, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_1), MultilayerVertex{Int64}(3, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_1), MultilayerVertex{Int64}(5, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(2, :layer_1), MultilayerVertex{Int64}(3, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(4, :layer_1), MultilayerVertex{Int64}(5, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(1, :layer_2), MultilayerVertex{Int64}(2, :layer_2), 0.32349997890438353)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(1, :layer_2), MultilayerVertex{Int64}(4, :layer_2), 0.3023269958672971)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(3, :layer_2), MultilayerVertex{Int64}(4, :layer_2), 0.2660051156748072)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(1, :layer_2), MultilayerVertex{Int64}(5, :layer_2), 0.23249716071401488)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(2, :layer_2), MultilayerVertex{Int64}(5, :layer_2), 0.05755829337786644)
 ⋮
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(2, :layer_4), MultilayerVertex{Int64}(2, :layer_2), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(3, :layer_4), MultilayerVertex{Int64}(3, :layer_2), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(4, :layer_4), MultilayerVertex{Int64}(4, :layer_2), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(5, :layer_4), MultilayerVertex{Int64}(5, :layer_2), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_4), MultilayerVertex{Int64}(1, :layer_3), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(2, :layer_4), MultilayerVertex{Int64}(2, :layer_3), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(3, :layer_4), MultilayerVertex{Int64}(3, :layer_3), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(4, :layer_4), MultilayerVertex{Int64}(4, :layer_3), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(5, :layer_4), MultilayerVertex{Int64}(5, :layer_3), 1)
```

The implementation of `MultilayerEdge` is

```julia
struct MultilayerEdge{ T <: MultilayerVertex, U <: Union{ <: Real, Nothing}} <: AbstractMultilayerEdge{T} # AbstractMultilayerEdge{T} subtypes AbstractEdge
    src::T    # The source vertex
    dst::T    # The destination vertex
    weight::U # The edge weight. Can be `nothing` to signify an unweighted edge, or a Real
end
```

#### Other Example APIs

Let's showcase some other key functionalities.

##### Get the node type

```julia
eltype(multilayergraph)
```
```nothing
Int64
```

##### Get the edge type

```julia
edgetype(multilayergraph)
```
```nothing
MultilayerEdge{MultilayerVertex{Int64}, Float64}
```

##### Check whether an edge exists

```julia
has_edge(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(4, :layer_2))
```
```nothing
false
```

##### Remove an edge

`rem_edge!` mimics the behaviour of the analogous function in Graphs.jl

```julia
rem_edge!(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2))
```
```nothing
The multilayer doesn't have any edge between MultilayerVertex{Int64}(1, :layer_1) and MultilayerVertex{Int64}(2, :layer_2)
false
```

The message tells us that the edge was already non existent. In fact, if we check the `adjacency_tensor` in the corresponding entry, we see that

```julia
multilayergraph.adjacency_tensor[1,2,1,2]
```
```nothing
0.0
```

##### Add an edge

We may add an edge using the `add_edge!` function. Since the interlayer we are adding the edge has an unweighted underlying graph (we will say that the interlayer is unweighted), we have to add an unweighted edge, so we don't specify the weight after the vertices. The adjacency tensor will be updated with a `one(U)` in the correct position. `add_edge!` mimics the behaviour of the analogous function in Graphs.jl

```julia
add_edge!(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2))
```
```nothing
true
```

To add a weighted edge we just need to write

```julia
add_edge!(multilayergraph,  MultilayerEdge(MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3), 3.14))
```
```nothing
true
```

##### Get the inneighbors

To get all the inneighbors of a vertex we just need to write

```julia
inneighbors(multilayergraph, MultilayerVertex(1, :layer_1))
```
```nothing
7-element Vector{MultilayerVertex{Int64}}:
 MultilayerVertex{Int64}(2, :layer_1)
 MultilayerVertex{Int64}(3, :layer_1)
 MultilayerVertex{Int64}(5, :layer_1)
 MultilayerVertex{Int64}(2, :layer_2)
 MultilayerVertex{Int64}(2, :layer_3)
 MultilayerVertex{Int64}(5, :layer_3)
 MultilayerVertex{Int64}(1, :layer_4)
```

[`outneighbors`](@ref) would be analogous.

##### Get the global clustering coefficient

```julia
multilayer_global_clustering_coefficient(multilayergraph)
```
```nothing
0.10125075759461827
```

Since our implementation of the global clustering coefficient follows [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022) rather than `Graphs.jl`'s implementation, we did not override `Graphs.jl`'s `global_clustering_coefficient`, which works on `MultilayerGraph` and `MultilayerDiGraph` but yields different results. For details, consult `?multilayer_global_clustering_coefficient` or read the comments in the source code.

### Multilayer-Specific Functions and Analysis

The following functions are specific to multilayer graphs or their implementations radically differ from their monoplex counterparts. For more information on every function, please refer to [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022) or consult the associated docstrings.

#### Overlay Monoplex Graph

Get the overlay monoplex graph: the monoplex graph whose nodes are the nodes of the multilayer graph and the edge between node $i$ and node $j$ has weight equal to the sum of all the weights of the edges between all vertex representations of $i$ and $j$ that belong to the same layer, for all the layers in the multilayer

```julia
get_overlay_monoplex_graph(multilayergraph)
```
```nothing
{5, 10} undirected simple Int64 graph with Float64 weights
```

#### Depth-Weighted Clustering Coefficient

Get the global clustering coefficient where triplets are weighted by how many layers they span

```julia
w =  [1/3, 1/3, 1/3]
multilayer_weighted_global_clustering_coefficient(multilayergraph,w)
```
```nothing
0.10125075759461745
```

The first component of `w` is the weight associated to triplets that are contained in one layer, the second component to triplets whose vertices are spread across exactly two layers, the third to triplets whose vertices are spread across exactly three layers. Weights must sum to `1.0`. When they are all equal (like in this example), the weighted global clustering coefficient coincides with the global clustering coefficient.

#### Eigenvector Centrality

Calculated via an iterative algorithm, its normalization is different from the Graphs.jl implementation. See `?eigenvector_centrality` for further details and context.

```julia
# The returned values are: the eigenvector centrality and the relative error at each iteration, that is, the summed absolute values of the componentwise differences between the centrality computed at the current iteration minus the centrality computed at the previous iteration.
eig_centrality, errs = eigenvector_centrality(multilayergraph; norm = "n", tol = 1e-3)
```
```nothing
([0.38448513173718 0.27376350094456875 0.29826894940433096 0.27928525958391875; 0.2095929988838825 0.212196610759091 0.440289042003606 0.2674763179954188; … ; 0.16489215196660842 0.17144436260745544 0.27826601575999654 0.23633134587445048; 0.27317146954683846 0.09008404085312828 0.3071561889077506 0.24335745518752538], [15.0, 0.2594289747318097, 0.11219926189825861, 0.045817552094948685, 0.02592915947732642, 0.011130317397910233, 0.006662267409348527, 0.003338105729180016, 0.0019678266671620814, 0.001239872171004075, 0.0007450702254126473])
```

#### Modularity

Compute the modularity of the multilayer graph. The signature mimics the Graphs.jl `modularity` implementation

```julia
modularity(multilayergraph,
          rand([1, 2, 3, 4], length(nodes(multilayergraph)),length(multilayergraph.layers)) # communities
          )
```
```nothing
-0.08991619711269938
```

#### Von Neumann Entropy

Compute the Von Neumann entropy as presented in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022)

```julia
von_neumann_entropy(multilayergraph)
```
```nothing
4.126274547913075
```

The Von Neumann entropy is currently available only for undirected multilayer graphs.

Other extended functions are: [`is_directed`](@ref), [`has_vertex`](@ref), [`ne`](@ref), [`nv`](@ref), [`outneighbors`](@ref), [`indegree`](@ref), [`outdegree`](@ref), [`degree`](@ref), [`mean_degree`](@ref), [`degree_second_moment`](@ref), [`degree_variance`](@ref), [`nn`](@ref), [`nodes`](@ref).

### Multiplex graphs

Special support has been given to multiplex graphs, i.e. multilayer graphs whose interlayer links are only allowed between vertices representing the same node. Such graphs come both in undirected `MultiplexGraph` and directed `MultiplexDiGraph` form. Their usage is completely analogous to `MultilayerGraph` and `MultilayerDiGraph`, with the few exceptions listed below.

Multiplex graphs need only the layers to be specified:

```julia
multiplexgraph = MultiplexGraph(layers)
```

They have been given a specialized random constructor:

```julia
multiplexgraph_random = MultiplexGraph( 4, n_nodes, min_edges, max_edges, [SimpleGraph{Int64}, SimpleWeightedGraph{Int64, Float64}, MetaGraph{Int64, Float64}])
```

One may not add edges between vertices  belonging to different layers:

```julia
add_edge!(multiplexgraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2), 3.14 )
```
```nothing
Adding an edge between vertices of different layers is not allowed within a multiplex graph.
```

## How to Contribute

The package is currently under development and further steps would benefit enormously from the precious feedback of the [JuliaGraph people](https://github.com/orgs/JuliaGraphs/people), graph theorists, network scientists and all the users who might have general questions or suggestions. 

Therefore feel free to open [discussions](https://github.com/JuliaGraphs/MultilayerGraphs.jl/discussions), [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues) or [PRs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/pulls). They are very welcome!   

## How to Cite

If you use this package in your work, please cite this repository using the metadata in [`CITATION.bib`](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/CITATION.bib).

## References

De Domenico et al. (2013) [Mathematical Formulation of Multilayer Networks](https://doi.org/10.1103/PhysRevX.3.041022). *Physical Review X*.
