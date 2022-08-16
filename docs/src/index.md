```@meta
CurrentModule = MultilayerGraphs
```

```@raw html
<div style="width:100%; height:150px;border-width:4px;border-style:solid;padding-top:25px;
        border-color:#000;border-radius:10px;text-align:center;background-color:#B3D8FF;
        color:#000">
    <h3 style="color: black;">Star us on GitHub!</h3>
    <a class="github-button" href="https://github.com/InPhyT/MultilayerGraphs.jl" data-icon="octicon-star" data-size="large" data-show-count="true" aria-label="Star InPhyT/MultilayerGraphs.jl on GitHub" style="margin:auto">Star</a>
    <script async defer src="https://buttons.github.io/buttons.js"></script>
</div>
```

# MultilayerGraphs.jl

**MultilayerGraphs.jl** is a Julia package for the construction, manipulation and analysis of multilayer graphs [extending Graphs.jl](https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/).

## Overview 

**MultilayerGraphs.jl** implements the mathematical formulation of multilayer graphs proposed by [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). It mainly revolves around two custom types, [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref), encoding undirected and directed multilayer graphs respectively. 

Roughly speaking, a multilayer graph is a collection of ***layers***, i.e. graphs whose vertices are representations of the same set of nodes, and ***interlayers***, i.e the [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) whose vertices are those of any two layers and whose edges are those between vertices of the same two layers. See below for the distinction between ***nodes*** and ***vertices***.

[`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) are fully-fledged [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) extensions. Both structs are designed so that their layers and interlayers can be of any type (as long as they are Graphs.jl extensions themselves) and they need not be all of the same type. It is anyway required that all layers and interlayers of [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) are respectively undirected and directed. Directedness is checked via the `IsDirected` trait defined in Graphs.jl adopting [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl). Since the layers' and interlayers' graph types don't need to be the same, multilayer graph types are considered weighted graphs by default, and thus are assigned the trait `IsWeighted`.

## Installation

Press `]` in the Julia REPL and then

```julia
pkg> add https://github.com/InPhyT/MultilayerGraphs.jl
```
[Registration](https://github.com/JuliaRegistries/General/pull/66311) is in progress.

## Tutorial 

Here we illustrate how to define, handle and analyse a [`MultilayerGraph`](@ref) (the directed version is completely analogous).

### Layers and Interlayers

Let's import some necessary packages

```julia
# Import necessary dependencies
using Graphs, SimpleWeightedGraphs, MultilayerGraphs
```

We define some methods and constants that will prove useful later in the tutorial

```julia
# Set the number of nodes, minimum and maximum number of edges for random graphs
const n_nodes   = 5
const min_edges = n_nodes
const max_edges = 10

# Define methods generating random graphs
get_SimpleGraph()   = SimpleGraph(n_nodes, rand(min_edges:max_edges))   # Undirected graph
get_SimpleDiGraph() = SimpleDiGraph(n_nodes, rand(min_edges:max_edges)) # Directed graph

# Define variables for random weighted graphs
const simpleweightedgraph_sources      = 1:n_nodes
const simpleweightedgraph_destinations = rand(1:n_nodes, n_nodes)
const simpleweightedgraph_weights      = rand(n_nodes)

# Define methods generating random weighted graphs
get_SimpleWeightedGraph()   = SimpleWeightedGraph(simpleweightedgraph_sources, rand(1:n_nodes, n_nodes), rand(n_nodes))    # Undirected graph
get_SimpleWeightedDiGraph() = SimpleWeightedDiGraph(simpleweightedgraph_sources, rand(1:n_nodes, n_nodes), rand(n_nodes))  # Directed graph
```

As said before, to define a multilayer graph we need to specify its layers and interlayers. We proceed by constructing a layer (see [`Layer`](@ref))

```julia
# Construct a layer 
layer = Layer(:layer_1, SimpleGraph(n_nodes, rand(min_edges:max_edges)); U = Float64)
```

A `Layer` has a name (here `:layer_1`), an underlying graph (`SimpleGraph(n_nodes, rand(min_edges:max_edges))`) and a weight matrix `eltype` `U` (it defaults to the adjacency matrix's `eltype` if the graph is unweighted). To correctly specify a multilayer graph all layers and interlayers must have the same `U`, otherwise the multilayers's adjacency tensor would be poorly specified. 

Notice that `U` does not need to coincide with the `eltype` of the adjacency matrix of the underlying graph: as far as we know, there is no way to set it explicitly for all Graphs.jl extensions, nor it is required for extensions to implement such feature, so our package converts to `U` the `eltype` of `Layer`s and `Interlayer`s weight (/adjacency) matrices every time they are invoked

```julia
adjacency_matrix(layers[1])
```
```nothing
5×5 SparseMatrixCSC{Float64, Int64} with 16 stored entries:
  ⋅   1.0   ⋅   1.0  1.0
 1.0   ⋅   1.0  1.0   ⋅
  ⋅   1.0   ⋅   1.0  1.0
 1.0  1.0  1.0   ⋅   1.0
 1.0   ⋅   1.0  1.0   ⋅
```

We may define more `Layer`s for future use

```julia
layers = [
            Layer(:layer_1, get_SimpleGraph(); U = Float64),
            Layer(:layer_2, get_SimpleWeightedGraph(); U = Float64),
            Layer(:layer_3, get_SimpleWeightedGraph(); U = Float64),
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
10×10 SparseMatrixCSC{Float64, Int64} with 18 stored entries:
  ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0
  ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅   1.0  1.0   ⋅ 
  ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅ 
  ⋅    ⋅    ⋅    ⋅    ⋅   1.0  1.0   ⋅    ⋅   1.0
  ⋅    ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅ 
  ⋅   1.0  1.0  1.0  1.0   ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅    ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
  ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
 1.0   ⋅    ⋅   1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
```

It is a 4 block matrix where the first `n_nodes` rows and columns refer to `:layer_1`'s vertices, while the last `n_nodes` rows and columns refer to `:layer_2`'s vertices.

We may define more `Interlayer`s for future use:

```julia
interlayers = [ 
                Interlayer(n_nodes, :interlayer_layer_1_layer_2, :layer_1, :layer_2, SimpleGraph{Int64}, rand(min_edges:max_edges); U = Float64), 
                Interlayer(n_nodes, :interlayer_layer_1_layer_3, :layer_1, :layer_3, SimpleWeightedGraph{Int64,Float64}, rand(min_edges:max_edges)) 
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
OrderedDict{Tuple{Int64, Int64}, Layer{Int64, U, G} where {U<:Real, G<:AbstractGraph{Int64}}} with 3 entries:
  (1, 1) => Layer{Int64, Float64, SimpleGraph{Int64}}(:layer_1, SimpleGraph{Int64}(5, [[2, 4], [1, 3, 4], [2], [1, 2, 5], [4]]), MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVertex{Int64}}[])
  (2, 2) => Layer{Int64, Float64, SimpleWeightedGraph{Int64, Float64}}(:layer_2, {5, 4} undirected simple Int64 graph with Float64 weights, MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVertex{Int64}}[])
  (3, 3) => Layer{Int64, Float64, SimpleWeightedGraph{Int64, Float64}}(:layer_3, {5, 4} undirected simple Int64 graph with Float64 weights, MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVertex{Int64}}[])
```

#### Interlayers

It is an `OrderedDict` where each key is the pair of indexes of the layers that the corresponding value, i.e. the interlayer, connects within the multilayer graph

```julia
multilayergraph.interlayers
```
```nothing
OrderedDict{Tuple{Int64, Int64}, Interlayer{Int64, U, G} where {U<:Real, G<:AbstractGraph{Int64}}} with 6 entries:
  (2, 1) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_2_layer_1, :layer_2, :layer_1, SimpleGraph{Int64}(9, [[7, 8, 9, 10], [9], [7], [7], [6, 9], [5], [1, 3, 4], [1], [1, 2, 5], [1]]), MultilayerVertex{Int64}[], Tuple{MultilayerVerte…  
  (1, 2) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_1_layer_2, :layer_1, :layer_2, SimpleGraph{Int64}(9, [[10], [6, 8, 9], [6], [6, 7, 10], [6], [2, 3, 4, 5], [4], [2], [2], [1, 4]]), MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64},…  
  (3, 1) => Interlayer{Int64, Float64, SimpleWeightedGraph{Int64, Float64}}(:interlayer_layer_3_layer_1, :layer_3, :layer_1, {10, 7} undirected simple Int64 graph with Float64 weights, MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVer…  
  (1, 3) => Interlayer{Int64, Float64, SimpleWeightedGraph{Int64, Float64}}(:interlayer_layer_1_layer_3, :layer_1, :layer_3, {10, 4} undirected simple Int64 graph with Float64 weights, MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVertex{Int64}…  
  (3, 2) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_3_layer_2, :layer_3, :layer_2, SimpleGraph{Int64}(5, [[6], [7], [8], [9], [10], [1], [2], [3], [4], [5]]), MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVert…  
  (2, 3) => Interlayer{Int64, Float64, SimpleGraph{Int64}}(:interlayer_layer_2_layer_3, :layer_2, :layer_3, SimpleGraph{Int64}(5, [[6], [7], [8], [9], [10], [1], [2], [3], [4], [5]]), MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVert…
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
Layer{Int64, Float64, SimpleGraph{Int64}}(:layer_1, SimpleGraph{Int64}(7, [[2, 3, 4, 5], [1, 4], [1, 5], [1, 2, 5], [1, 3, 4]]), MultilayerVertex{Int64}[], Tuple{MultilayerVertex{Int64}, MultilayerVertex{Int64}}[])
```

#### Adjacency Tensor

The adjacency tensor is a 4-dimensional array

```julia
multilayergraph.adjacency_tensor
```
```nothing
5×5×3×3 Array{Float64, 4}:
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
interlayer_2_1 = get_layer(multilayergraph, :interlayer_layer_2_layer_1)
```

Both `MultilayerGraph` and `MultilayerDiGraph` fully extend `Graphs.jl`, so they have access to Graphs.jl API as one would expect, just keeping in mind that vertices are `MultilayerVertex`s and not subtypes of `Integer` (`MultilayerVertex` is actually a subtype of `AbstractVertex` that this package defines, see [Future Developments](#Future-Developments)), and that edges are `MultilayerEdge`s, which subtype `AbstractEdge`.

Some notable examples are

```julia
edges(multilayergraph)
```
```nothing
34-element Vector{MultilayerEdge}:
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_1), MultilayerVertex{Int64}(2, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_1), MultilayerVertex{Int64}(3, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_1), MultilayerVertex{Int64}(4, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(1, :layer_1), MultilayerVertex{Int64}(5, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(2, :layer_1), MultilayerVertex{Int64}(3, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(2, :layer_1), MultilayerVertex{Int64}(5, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Int64}(MultilayerVertex{Int64}(3, :layer_1), MultilayerVertex{Int64}(5, :layer_1), 1)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(1, :layer_2), MultilayerVertex{Int64}(1, :layer_2), 0.7188425521261754)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(2, :layer_2), MultilayerVertex{Int64}(3, :layer_2), 0.9012061650463197) 
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(2, :layer_2), MultilayerVertex{Int64}(4, :layer_2), 0.6163304419976594) 
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(3, :layer_2), MultilayerVertex{Int64}(5, :layer_2), 1.0046265072746847) 
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(1, :layer_3), MultilayerVertex{Int64}(2, :layer_3), 0.2819477742859873) 
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(2, :layer_3), MultilayerVertex{Int64}(4, :layer_3), 0.40111133874926597)
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(3, :layer_3), MultilayerVertex{Int64}(4, :layer_3), 0.9498077050078636) 
 MultilayerEdge{MultilayerVertex{Int64}, Float64}(MultilayerVertex{Int64}(1, :layer_3), MultilayerVertex{Int64}(5, :layer_3), 0.9618455695308973) 
 ⋮
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
8-element Vector{MultilayerVertex{Int64}}:
 MultilayerVertex{Int64}(2, :layer_1)
 MultilayerVertex{Int64}(3, :layer_1)
 MultilayerVertex{Int64}(4, :layer_1)
 MultilayerVertex{Int64}(5, :layer_1)
 MultilayerVertex{Int64}(1, :layer_3)
 MultilayerVertex{Int64}(2, :layer_3)
 MultilayerVertex{Int64}(4, :layer_3)
 MultilayerVertex{Int64}(5, :layer_3)
```

[`outneighbors`](@ref) would be analogous.

##### Get the global clustering coefficient

```julia
multilayer_global_clustering_coefficient(multilayergraph)
```
```nothing
0.12667622867320932
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
{5, 8} undirected simple Int64 graph with Float64 weights
```

#### Depth-Weighted Clustering Coefficient

Get the global clustering coefficient where triplets are weighted by how many layers they span

```julia
w =  [1/3, 1/3, 1/3] 
multilayer_weighted_global_clustering_coefficient(multilayergraph,w)
```
```nothing 
0.12667622867320916
```

The first component of `w` is the weight associated to triplets that are contained in one layer, the second component to triplets whose vertices are spread across exactly two layers, the third to triplets whose vertices are spread across exactly three layers. Weights must sum to `1.0`. When they are all equal (like in this example), the weighted global clustering coefficient coincides with the global clustering coefficient. 

#### Eigenvector Centrality

Calculated via an iterative algorithm, its normalization is different from the Graphs.jl implementation. See `?eigenvector_centrality` for further details and context.

```julia
# The returned values are: the eigenvector centrality and the relative error at each iteration, that is, the summed absolute values of the componentwise differences between the centrality computed at the current iteration minus the centrality computed at the previous iteration.
eig_centrality, errs = eigenvector_centrality(multilayergraph; norm = "n", tol = 1e-3)
```
```nothing 
([0.260450377858897 0.02358226618172732 0.08408641909534659; 0.6517125919818912 0.646849109975291 0.25212976452498587; … ; 0.21160486454350666 0.24850764907380082 0.2809259415084613; 0.5089709123112656 0.27775286552140954 0.22663014660085468], [10.000000000000004, 0.7095737447337795, 0.329314405320397, 0.16256401870083886, 0.07860800793625616, 0.04125218689190476, 0.021191358681432643, 0.011056758224944913, 0.006831972142635288, 0.0035133270128280235, 0.0024675719418019056, 0.0013403531942171656, 0.0009739276503972424])
```

#### Modularity

Compute the modularity of the multilayer graph. The signature mimics the Graphs.jl `modularity` implementation

```julia
modularity(multilayergraph,
          rand([1, 2, 3, 4], length(nodes(multilayergraph)),length(multilayergraph.layers)) # communities
          )
```
```nothing
-0.039890139283044884
```

#### Von Neumann Entropy

Compute the Von Neumann entropy as presented in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022)

```julia
von_neumann_entropy(multilayergraph)
```
```nothing
3.3980014398404834
```

The Von Neumann entropy is currently available only for undirected multilayer graphs.

Other extended functions are: [`is_directed`](@ref), [`has_vertex`](@ref), [`ne`](@ref), [`nv`](@ref), [`outneighbors`](@ref), [`indegree`](@ref), [`outdegree`](@ref), [`degree`](@ref), [`mean_degree`](@ref), [`degree_second_moment`](@ref), [`degree_variance`](@ref), [`nn`](@ref), [`nodes`](@ref).

## How to Contribute

If you wish to change or add some functionality, please file an [issue](https://github.com/InPhyT/MultilayerGraphs.jl/issues). 

## How to Cite 

If you use this package in your work, please cite this repository using the metadata in [`CITATION.bib`](https://github.com/InPhyT/MultilayerGraphs.jl/blob/main/CITATION.bib).

## References 

De Domenico et al. (2013) [Mathematical Formulation of Multilayer Networks](https://doi.org/10.1103/PhysRevX.3.041022). *Physical Review X*