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

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs [extending Graphs.jl](https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/).

## Overview

**MultilayerGraphs.jl** implements the mathematical formulation of multilayer graphs proposed by [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022) together with insights from [Kivelä et al. (2014)](https://doi.org/10.1093/comnet/cnu016) and [Bianconi (2018)](https://global.oup.com/academic/product/multilayer-networks-9780192865540). It mainly revolves around two custom types, [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref), encoding undirected and directed multilayer graphs respectively.

Roughly speaking, a multilayer graph is a collection of ***layers***, i.e. graphs whose vertices are representations of the same set of nodes (not all nodes have to be represented in every layer), and ***interlayers***, i.e the [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) whose two sets of vertices are those of any two layers. A vertex of a multilayer graph will be represented via a [`MultilayerVertex`](@ref) struct, and nodes via a [`Node`](@ref) struct.

[`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) are fully-fledged [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) extensions. Both structs are designed so that their layers and interlayers can be of any type (as long as they are Graphs.jl extensions themselves) and they can be of different types. It is anyway required that all layers and interlayers of [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) are respectively undirected and directed.

Both [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) allow for vertex and edge metadata, provided that the layer (or interlayer) vertex or the layer (or interlayer) edge belongs to supports metadata.

## Installation

Press `]` in the Julia REPL and then

```nothing
pkg> add MultilayerGraphs
```

## Tutorial

Here we illustrate how to define, handle and analyse a [`MultilayerGraph`](@ref) (the directed version is completely analogous).

```julia
using Revise
using StatsBase, Distributions
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs
using MultilayerGraphs
```

Define some constants that will prove useful later in the tutorial:

```julia
# Set the minimum and maximum number of all_nodes and edges for random graphs
const vertextype   = Int64
const _weighttype  = Float64
const min_vertices = 5
const max_vertices = 7
const min_edges    = 1
const max_edges    = max_vertices*(max_vertices-1)
const n_nodes      = max_vertices
; #hide
```

Next we define nodes:

```julia
# The constructor for nodes (which are immutable) only requires a name (`id`) for the node
const all_nodes = [Node("node_$i") for i in 1:n_nodes]
```
```nothing
7-element Vector{Node}:
 Node("node_1")
 Node("node_2")
 Node("node_3")
 Node("node_4")
 Node("node_5")
 Node("node_6")
 Node("node_7")
```

And construct `MultilayerVertex`s from these nodes:

```julia
## Convert nodes to multilayer vertices without metadata
const multilayervertices = MV.(all_nodes)
## Convert nodes multilayer vertices with metadata
const multilayervertices_meta  = [MV(node, ("I'm node $(node.id)",)) for node in all_nodes]
```
```nothing
7-element Vector{MultilayerVertex{nothing}}:
 MV(Node("node_1"), :nothing, ("I'm node node_1",))
 MV(Node("node_2"), :nothing, ("I'm node node_2",))
 MV(Node("node_3"), :nothing, ("I'm node node_3",))
 MV(Node("node_4"), :nothing, ("I'm node node_4",))
 MV(Node("node_5"), :nothing, ("I'm node node_5",))
 MV(Node("node_6"), :nothing, ("I'm node node_6",))
 MV(Node("node_7"), :nothing, ("I'm node node_7",))
```

This conversion is done since it is logical to add vertices to a graph, not nodes, and also for consistency reasons with the ecosystem.

Printing a `MultilayerVertex` returns:

```julia
multilayervertices_meta[1]
```
```nothing 
MV(Node("node_1"), :nothing, ("I'm node node_1",))
```

Where `MV` is an alias for `MultilayerVertex`. The first field is the `Node` being represented, the second the (name of) the layer the vertex is represented in (here it is set to `nothing`, since these vertices are yet to be assigned), and the metadata associated to the vertex (no metadata are currently represented via an empty `NamedTuple`). `MultilayerVertex` metadata can be represented via a `Tuple` or a `NamedTuple` (see below for examples).

### Layers

As said before, to define a multilayer graph we need to specify its layers and interlayers. Layers and the individual graphs that make up multilayer graphs. We proceed by constructing a [`Layer`](@ref) using the constructor that randomly specifies the edges:

```julia
Layer(
    name::Symbol,                                                 # The name of the layer
    vertices::Vector{ <: MultilayerVertex},                       # The `MultilayerVertex`s of the Layer
    ne::Int64,                                                    # The number of edges of the Layer
    null_graph::G,                                                # The Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
    weighttype::Type{U};                                          # The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
    default_vertex_metadata::Function = mv -> NamedTuple(),       # Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
    default_edge_weight::Function = (src, dst) -> nothing,        #  Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
    default_edge_metadata::Function = (src, dst) -> NamedTuple(), # Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
    allow_self_loops::Bool = false                                # whether to allow self loops to be generated or not. Defaults to `false`.
)
```

A `Layer` is considered "weighted" if its underlying graph (`null_graph` argument) has been given the `IsWeighted` trait (traits throughout this package are implemented via from [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl), just like Graphs.jl does). Since one may at any moment add a new weighted `Layer` to a `MultilayerGraph` (see below for details), the latter is always considered a "weighted graph", so it is given the `IsWeighted` trait. Thus, all `Layer`s and `Interlayer`s (collectively named "subgraphs" hereafter) must specify their `weighttype` as the last argument of their constructor, so the user may debug their weight matrices immediately after construction. As better specified below, all subgraphs that are meant to be part of the same `MultilayerGraph` must have the same `weighttype`.   

Before instantiating `Layer`s, we define an utility function to ease randomization:

```julia
# Utility function that returns a random number of vertices and edges each time it is called:
function rand_nv_ne_layer(min_vertices, max_vertices)
    _nv = rand(min_vertices:max_vertices)
    _ne = rand(1:(_nv*(_nv-1)) ÷ 2 )
    return (_nv,_ne)
end

# Utility function that returns two vertices of a Layer that are not adjacent.
function _get_srcmv_dstmv_layer(layer::Layer)

    mvs = get_bare_mv.(collect(mv_vertices(layer)))

    src_mv = nothing    
    _collection = []

    while isempty(_collection)
        src_mv = rand(mvs)
        _collection = setdiff(Set(mvs), Set(vcat(get_bare_mv.(mv_outneighbors(layer, src_mv)), src_mv ) ) )  
    end

    dst_mv = get_bare_mv(rand(_collection))

    return mvs, src_mv, dst_mv
end
; #hide
```

We are now are ready to define some `Layer`s. Every type of graph from the Graphs.jl ecosystem may underlie a `Layer` (or an `Interlayer`). We will construct a few of them, each time with a different number of vertices and edges.

```julia
# An unweighted simple layer:
_nv, _ne  = rand_nv_ne_layer(min_vertices,max_vertices)
layer_sg = Layer(   :layer_sg,
                    sample(multilayervertices, _nv, replace = false),
                    _ne, 
                    SimpleGraph{vertextype}(),
                    _weighttype
)

# A weighted `Layer`
_nv, _ne  = rand_nv_ne_layer(min_vertices,max_vertices)
layer_swg = Layer(  :layer_swg, 
                    sample(multilayervertices, _nv, replace = false),
                    _ne, 
                    SimpleWeightedGraph{vertextype, _weighttype}(),
                    _weighttype; 
                    default_edge_weight = (src,dst) -> rand()
)
# A `Layer` with an underlying `MetaGraph`: 
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
layer_mg = Layer(   :layer_mg, 
                    sample(multilayervertices_meta, _nv, replace = false), 
                    _ne, 
                    MetaGraph{vertextype, _weighttype}(),
                    _weighttype; 
                    default_edge_metadata = (src,dst) -> (from_to = "from_$(src)_to_$(dst)",)
)
# `Layer` with an underlying `ValGraph` from `SimpleValueGraphs.jl`
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
layer_vg = Layer(   :layer_vg, 
                    sample(multilayervertices_meta, _nv, replace = false), 
                    _ne,
                    MultilayerGraphs.ValGraph{vertextype}(;edgeval_types=(Float64, String, ),
                                            edgeval_init=(s, d) -> (s+d, "hi"),
                                            vertexval_types=(String,),
                                            vertexval_init=v -> ("$v",),),
                    _weighttype;
                    default_edge_metadata = (src,dst) -> (rand(), "from_$(src)_to_$(dst)",),
                    default_vertex_metadata = mv -> ("This metadata had been generated via the default_vertex_metadata method",)
)

# Collect all layers in an ordered list. Order will be recorded when instantiating the multilayer graph.
layers = [layer_sg, layer_swg, layer_mg, layer_vg]
; #hide
```

The API that inspects and modifies `Layer`s will be shown below together with that of `Interlayer`s, since they are usually the same. There are of course other constructors that you may discover by typing `?Layer` in the console.

### Interlayers

Now we turn to defining `Interlayer`s. Interlayers are the graphs containing all the edges between vertices is two distinct layers. As before, we need an utility to ease randomization:


```julia
# Utilities for Interlayer
## Utility function that returns two vertices of an Interlayer that are not adjacent.
function _get_srcmv_dstmv_interlayer(interlayer::Interlayer)

    mvs = get_bare_mv.(collect(mv_vertices(interlayer)))

    src_mv = nothing    
    _collection = []

    while isempty(_collection)
        src_mv = rand(mvs)
        _collection = setdiff(Set(mvs), Set(vcat(get_bare_mv.(mv_outneighbors(interlayer, src_mv)), src_mv, get_bare_mv.(mv_vertices( eval(src_mv.layer) ))) ) )  
    end

    dst_mv = get_bare_mv(rand(_collection))

    return mvs, src_mv, dst_mv
end

## Utility function that returns a random number edges between its arguments `layer_1` and `layer_2`:
function rand_ne_interlayer(layer_1, layer_2)
    _nv = nv(layer_1) + nv(layer_2)
    _ne = rand(_nv:(_nv*(_nv-1)) ÷ 2 )
    return _ne
end
; #hide
```

An `Interlayer` is constructed by passing its name, the two `Layer`s it should connect, and the other parameters just like the `Layer`'s constructor. The random constructor reads:
```julia
Interlayer(
    layer_1::Layer{T,U},                                                 # One of the two layers connected by the Interlayer
    layer_2::Layer{T,U},                                                 # One of the two layers connected by the Interlayer  
    ne::Int64,                                                           # The number of edges of the Interlayer
    null_graph::G;                                                       # the Interlayer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
    default_edge_weight::Function = (x,y) -> nothing,                    # Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
    default_edge_metadata::Function = (x,y) -> NamedTuple(),             # Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
    name::Symbol = Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"), # The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");
    transfer_vertex_metadata::Bool = false                               # if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors.;
) 
```
We will build a few of random `Interlayer`s:

```julia
# Define the random undirected simple Interlayer
_ne = rand_ne_interlayer(layer_sg, layer_swg)
interlayer_sg_swg = Interlayer( layer_sg,                  # The first layer to be connected
                                layer_swg,                 # The second layer to be connected
                                _ne,                       # The number of edges to randomly generate
                                SimpleGraph{vertextype}(), # The underlying graph, passed as a null graph
                                name = :random_interlayer  # The name of the interlayer. We will be able to access it as a property of the multilayer graph via its name. This kwarg's default value is given by a combination of the two layers' names.
)
# Define a weighted `Interlayer`
_ne = rand_ne_interlayer(layer_swg, layer_mg)
interlayer_swg_mg = Interlayer( layer_swg, 
                                layer_mg,
                                _ne,
                                SimpleWeightedGraph{vertextype, _weighttype}(); 
                                default_edge_weight = (x,y) -> rand() # Arguments follow the same rules as in Layer
) 
# Define an `Interlayer` with an underlying `MetaGraph`
_ne = rand_ne_interlayer(layer_mg, layer_vg)
interlayer_mg_vg = Interlayer(  layer_mg,
                                layer_vg, 
                                _ne, 
                                MetaGraph{vertextype, _weighttype}(); 
                                default_edge_metadata = (x,y) -> (mymetadata = rand(),), 
                                transfer_vertex_metadata = true # This boolean kwarg controls whether vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors.
)
# Define an `Interlayer` with an underlying `ValGraph` from `SimpleValueGraphs.jl`, with diagonal couplings only:
interlayer_multiplex_sg_mg = multiplex_interlayer(  layer_sg, 
                                                    layer_mg, 
                                                    ValGraph{vertextype}(; edgeval_types=(from_to = String,), edgeval_init=(s, d) -> (from_to = "from_$(s)_to_$(d)")); 
                                                    default_edge_metadata = (x,y) -> (from_to = "from_$(src)_to_$(dst)",)
) 
# Finally, An `Interlayer` with no couplings (an "empty" interlayer):
interlayer_empty_sg_vg = empty_interlayer(  layer_sg, 
                                            layer_vg, 
                                            SimpleGraph{vertextype}()
)

# Collect all interlayers. Even though the list is ordered, order will not matter when instantiating the multilayer graph.
interlayers = [interlayer_sg_swg, interlayer_swg_mg, interlayer_mg_vg, interlayer_multiplex_sg_mg, interlayer_empty_sg_vg]
; #hide
```

Next, we explore the API associated to modify and analyze `Layer`s and `Interlayer`s.

### Subgraphs API

API for  `Layer`s and `Interlayer`s (collectively, "subgraphs") are very similar, so we will just show them for the `Layer` case, pointing out differences to the `Interlayer` scenario whenever they occur.

Subgraphs extend the Graphs.jl's interface, so one may expect every method from Graphs.jl to apply. Anyway, the output and signature is slightly different and thus worth pointing out below.

#### Nodes

One may retrieve the `Node`s that a `Layer` represents via:

```julia
layer_sg_nodes = nodes(layer_sg)
```
```nothing 
6-element Vector{Node}:
 Node("node_2")
 Node("node_3")
 Node("node_4")
 Node("node_6")
 Node("node_5")
 Node("node_7")
```

The same would be for `Interlayer`s. In this case, the union of the set of nodes represented by the two layers the interlayer connects is returned:

```julia
interlayer_sg_swg_nodes  = nodes(interlayer_sg_swg)
```
```nothing 
7-element Vector{Node}:
 Node("node_2")
 Node("node_3")
 Node("node_4")
 Node("node_6")
 Node("node_5")
 Node("node_7")
 Node("node_1")
``` 

One may check for the existence of a node within a layer (or interlayer) via:

```julia
has_node(layer_sg, layer_sg_nodes[1])
```
```nothing 
true
```

#### Vertices

One may retrieve the `MultilayerVertex`s of a layer by calling:

```julia
layer_sg_vertices = mv_vertices(layer_sg)
```
```nothing 
6-element Vector{MultilayerVertex{:layer_sg}}:
 MV(Node("node_2"), :layer_sg, NamedTuple())
 MV(Node("node_3"), :layer_sg, NamedTuple())
 MV(Node("node_4"), :layer_sg, NamedTuple())
 MV(Node("node_6"), :layer_sg, NamedTuple())
 MV(Node("node_5"), :layer_sg, NamedTuple())
 MV(Node("node_7"), :layer_sg, NamedTuple())
```

While vertices with metadata would look like:

```julia
mv_vertices(layer_mg)
```
```nothing 
6-element Vector{MultilayerVertex{:layer_mg}}:
 MV(Node("node_7"), :layer_mg, (var"1" = "I'm node node_7",))
 MV(Node("node_6"), :layer_mg, (var"1" = "I'm node node_6",))
 MV(Node("node_2"), :layer_mg, (var"1" = "I'm node node_2",))
 MV(Node("node_4"), :layer_mg, (var"1" = "I'm node node_4",))
 MV(Node("node_1"), :layer_mg, (var"1" = "I'm node node_1",))
 MV(Node("node_5"), :layer_mg, (var"1" = "I'm node node_5",))
```

The vertices of an interlayer are the union of the sets of vertices of the two layers it connects:

```julia
interlayer_sg_swg_vertices = mv_vertices(interlayer_sg_swg)
```
```nothing 
11-element Vector{MultilayerVertex}:
 MV(Node("node_2"), :layer_sg, NamedTuple())
 MV(Node("node_3"), :layer_sg, NamedTuple())
 MV(Node("node_4"), :layer_sg, NamedTuple())
 MV(Node("node_6"), :layer_sg, NamedTuple())
 MV(Node("node_5"), :layer_sg, NamedTuple())
 MV(Node("node_7"), :layer_sg, NamedTuple())
 MV(Node("node_1"), :layer_swg, NamedTuple())
 MV(Node("node_4"), :layer_swg, NamedTuple())
 MV(Node("node_5"), :layer_swg, NamedTuple())
 MV(Node("node_2"), :layer_swg, NamedTuple())
 MV(Node("node_3"), :layer_swg, NamedTuple())
```

The `vertices` command would return an internal representation of the `MultilayerVertex`s. This method, together with others, serves to make `MultilayerGraphs.jl` compatible with the Graphs.jl ecosystem, but it is not meant to be called by the end user. It is, anyway, thought to be used by developers who wish to interface their packages with `MultilayerGraphs.jl` just as with other packages of the `Graphs.jl` ecosystem: a developer-oriented guide will be compiled if there is the need. 

In the [API](@ref) page the intended usage of all methods (*end-user* or *developer*) is highlighted.

To add a vertex, simply use [`add_vertex!`](@ref). Let us define a vertex with metadata to add. Since nodes may not be represented more than once in layers, we have to define a new node to:

```julia
new_node     = Node("missing_node")
new_metadata = (meta = "my_metadata",)
new_vertex   = MV(new_node, new_metadata)
```
```nothing 
MV(Node("missing_node"), :nothing, (meta = "my_metadata",))
```

Of course, to be able to add a vertex with metadata to a layer, one must make sure that the underlying graph supports vertex-level metadata. Should one try to add a vertex with metadata different from an empty `NamedTuple` (i.e. no metadata) to a layer whose underlying graph does not support metadata, a warning is issued and the metadata are discarded.

Thus, if we consider a layer whose underlying graph is a `MetaGraph`, the following three syntaxes would be equivalent.

- The *standard* interface:
```julia
add_vertex!(layer_mg, new_vertex)
```
- The *uniform* interface. This signature has one keyword argument, `metadata`:
```julia
add_vertex!(layer_mg, new_node, metadata = new_metadata)
```
- The *transparent* interface. After you pass to `add_vertex` the `Layer` and the `Node` you wish to add, you  may pass the same `args` and `kwargs`  that you would pass to the `add_vertex!` dispatch that acts on the underlying graph (after the graph argument). This is a way to let the user directly exploit the API of the underlying graph package, which could be useful for two reasons:
    1. They may be more convenient;
    1. They should work even if we are not able to integrate the *standard* and the *uniform* interface with a particular `Graphs.jl`'s extension.
  
    Here is an example on how to use it:
```julia
add_vertex!(layer_mg, new_node, Dict(pairs(new_metadata)))
```
where `Dict(pairs(new_metadata))` is exactly what you would pass to the `add_vertex!` method that acts on `MetaGraphs`:
```julia
metagraph = MetaGraph()
add_vertex!(metagraph,  Dict(pairs(new_metadata)))
```

If an underlying graph has an `add_vertex!` interface whose signature overlaps with that of the uniform interface, the uniform interface will be prevail.

If, using the *transparent* interface, one does not specify any `metadata`, the `default_vertex_metadata` function passed to the `Layer`'s constructor is called to provide `metadata` to the vertex (type `?Layer` in the REPL for more information).

To remove the vertex, simply do:
```julia
rem_vertex!(layer_sg, new_vertex) # Returns true if succeeds
```s
To extract metadata:
```julia
get_metadata(layer_mg, MV(new_node))
```

By design, one may not add nor remove vertices to `Interlayer`s.

Please refer to the Vertex section of the API page ([end-user]() and [developer]()) to discover more methods related to `MultilayerVertex`s. 

### Edges

The edge type for multilayer graphs (and thus for this subgraphs) is `MultilayerEdge`, which has a type parameter corresponding to the chosen weight type:

```julia
edgetype(layer_sg)
```
```nothing 
MultilayerEdge{Float64}
```

The `MultilayerEdge`s of an unweighted simple layer are:

```julia
collect(edges(layer_sg))
```
```nothing
9-element Vector{MultilayerEdge{Float64}}:
 ME(MV(Node("node_2"), :layer_sg, NamedTuple()) --> MV(Node("node_4"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_2"), :layer_sg, NamedTuple()) --> MV(Node("node_6"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_2"), :layer_sg, NamedTuple()) --> MV(Node("node_5"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_3"), :layer_sg, NamedTuple()) --> MV(Node("node_4"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_3"), :layer_sg, NamedTuple()) --> MV(Node("node_6"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_4"), :layer_sg, NamedTuple()) --> MV(Node("node_6"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_4"), :layer_sg, NamedTuple()) --> MV(Node("node_5"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_6"), :layer_sg, NamedTuple()) --> MV(Node("node_5"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
 ME(MV(Node("node_6"), :layer_sg, NamedTuple()) --> MV(Node("node_7"), :layer_sg, NamedTuple()),	weight = 1.0,	metadata = NamedTuple())
```

Where `ME` is a shorthand for `MultilayerEdge`. Besides the two vertices connected, each `MultilayerEdge` carries the information about its `weight` and `metadata`. For unweighted subgraphs, the weight is just `one(weighttype)` and for non-meta subgraphs the metadata are an empty `NamedTuple`s. See `?MultilayerEdge` for additional information.

The `add_edge` function has the standard, uniform and transparent interfaces too. To understand how they work, let's define a weighted edge:

```julia
# Define a weighted edge for the layer_swg
## Define the weight
_weight = rand()
## Select two non-adjacent vertices in layer_swg
_, src_w, dst_w  = _get_srcmv_dstmv_layer(layer_swg)
## Construct a weighted MultilayerEdge
me_w = ME(src_w, dst_w, _weight) # ME is an alias for MultilayerEdge
```
```nothing
ME(MV(Node("node_4"), :layer_swg, NamedTuple()) --> MV(Node("node_1"), :layer_swg, NamedTuple()),	weight = 0.6369546116248217,	metadata = NamedTuple())
```

Of course, to be able to add a weighted edge to a subgraph, one must make sure that the underlying graph supports edge weights. Should one try to add a weight different from `one(weighttype)` or `nothing` to an edge of a subgraph whose underlying graph does not support edge weights, a warning is issued and the weight is discarded.

Thus, if we consider a layer whose underlying graph is a `SimpleWeightedGraph`, the following three syntaxes would be equivalent.

- The *standard* interface:
```julia
add_edge!(layer_swg, me_w)
```
- The *uniform* interface. This signature has two keyword arguments, `weight` and `metadata` that could be used exclusively (if, respectively, the underlying graph is weighted or supports edge-level metadata) or in combination (if the underlying graph supports both edge weights and edge-level metadata):
```julia
add_edge!(layer_swg, src_w, dst_w, weight = _weight)
```
- The *transparent* interface. After you pass to `add_edge!` the `Layer` and the two vertices you wish to connect, you  may pass the same `args` and `kwargs`  that you would pass to the `add_edge!` dispatch that acts on the underlying graph (after the graph and vertices arguments). This is done for the same reasons explained above.
  Here is an example on how to use it:
```julia
add_edge!(layer_swg, src_w, dst_w, _weight)
```
where `_weight` is exactly what you would pass to the `add_edge!` method that acts on `SimpleWeightedGraph` after:
```julia
simpleweightedgraph = SimpleWeightedGraph(5, 0)
add_edge!(simpleweightedgraph, 1, 2, _weight)
```

If an underlying graph has an `add_edge!` interface whose signature overlaps with that of the uniform interface, the uniform interface will prevail.

If, using the *transparent* interface, one does not specify any `weight` or (inclusively) `metadata` keyword argument, the `default_edge_weight` or (inclusively) the `default_edge_metadata` function passed to the `Layer`'s constructor will be called to provide `weight` or `metadata` to the edge (type `?Layer` in the REPL for more information).

To remove the edge, simply do:
```julia
rem_edge!(layer_swg, src_w, dst_w) # Returns true if succeeds
```

To extract weight:
```julia
get_weight(layer_swg, src_w, dst_w)
```

For an edge with metadata, it would be analogous. Let's define an edge with metadata:

```julia
# Define an edge with metadata for the layer_mg
## Define the metadata
_metadata  = (meta = "mymetadata",)
## Select two non-adjacent vertices in layer_mg
_, src_m, dst_m  = _get_srcmv_dstmv_layer(layer_mg)
## Construct a MultilayerEdge with metadata
me_m = ME(src_m, dst_m, _metadata)
```
```nothing 
ME(MV(Node("node_6"), :layer_mg, NamedTuple()) --> MV(Node("node_5"), :layer_mg, NamedTuple()),	weight = nothing,	metadata = (meta = "mymetadata",))
```

Then the following three signatures would be equivalent:

- *standard* interface:
```julia
add_edge!(layer_mg, me_m)
```
- *uniform* interface:
```julia
add_edge!(layer_mg, src_m, dst_m, metadata = _metadata)
```
- *transparent* interface
```julia
add_edge!(layer_mg, src_m, dst_m, Dict(pairs(_metadata)))
```

To extract metadata:
```julia
get_metadata(layer_mg, src_m, dst_m)
```

Please refer to the Vertex section of the API page ([end-user]() and [developer]()) to discover more methods related to `MultilayerEdges`s.

For the `layer_swg`, the following three signatures would be equivalent:

- *standard* interface:
```julia
add_edge!(layer_swg, me_w)
```
- *uniform* interface:
```julia
add_edge!(layer_swg, src_w, dst_w, weight = _weight)
```
- *transparent* interface
```julia
add_edge!(layer_swg, src_w, dst_w, _weight)
```

The uniform interface of `add_edge!` works so that the user may specify the keyword `weight` and/or the keyword `metadata`. If an underlying subgraph has a transparent interface whose signature overlaps with that of the uniform interface, the uniform interface will be prevail.

The edge may be removed via 

```julia
rem_edge!(layer_swg, src_w, dst_w)
```

Please refer to the `MultilayerEdge` section of the API page ([end-user]() and [developer]()) to discover more methods related to `MultilayerEdge`s.

### Multilayer Graphs

Given all the `Layer`s and the `Interlayer`s, let's instantiate a multilayer graph as follows:

```julia
multilayergraph = MultilayerGraph(  layers,                                                 # The (ordered) list of layers the multilayer graph will have
                                    interlayers;                                            # The list of interlayers specified by the user. Note that the user does not need to specify all interlayers, as the unspecified ones will be automatically constructed using the indications given by the `default_interlayers_null_graph` and `default_interlayers_structure` keywords.
                                    default_interlayers_null_graph = SimpleGraph{vertextype}(), # Sets the underlying graph for the interlayers that are to be automatically specified.  Defaults to `SimpleGraph{T}()`. See the `Layer` constructors for more information.
                                    default_interlayers_structure = "multiplex" # Sets the structure of the interlayers that are to be automatically specified. May be "multiplex" for diagonally coupled interlayers, or "empty" for empty interlayers (no edges).  "multiplex". See the `Interlayer` constructors for more information.
);
```

Keep in mind that `Multilayer(Di)Graph` only supports uniform and standard interface for both `add_vertex!` and `add_edge!`.

As already stated, a `MultilayerGraph` is an object made of `Layer`s and `Interlayer`s whose collections of vertices each represents a subset of the set of nodes, here being `nodes`.

*Adding* a `Node` to a `MultilayerGraph` will enable its `Layer`s  (and thus its `Interlayer`s) to represent it i.e. you will be able to add `MultilayerVertex`s that represent that `Node` to the multilayer graph.

Another constructor allows for a limited configuration model-like specification. It allows to generate a multilayer graph with a specific degree distribution. Since edges are created according to the provided distribution, it is necessary that the layers and interlayers specified are empty (i.e. they have no edges). Notice that layers and interlayers whose underlying graph is a `SimpleWeighted(Di)Graph` may not be used until [this PR](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl/pull/14) is merged.

It is used as:


```julia
# First, we need to empty the above layers and interlayers, and remove the ones having a `SimpleWeightedGraph`s. These lines are not necessary to comprehend the tutorial, they may be skipped. Just know that the variables `empty_layers` and `empty_interlayers` are two lists of. respectively, empty layers and interlayers that do not have `SimpleWeightedGraph`s as their underlying graphs

empty_layers =  deepcopy([layer for layer in layers if !(layer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph)])

empty_layers_names = name.(empty_layers)

empty_interlayers =  deepcopy([interlayer for interlayer in interlayers if all(in.(interlayer.layers_names, Ref(empty_layers_names))) && !(interlayer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph) ])

for layer in empty_layers
    for edge in edges(layer)
        rem_edge!(layer, edge)
    end
end

for interlayer in empty_interlayers
    for edge in edges(interlayer)
        rem_edge!(interlayer, edge)
    end
end

# Construct a multilayer graph that has a normal degree distribution. The support of the distribution must be positive, since negative degrees are not possible
configuration_multilayergraph = MultilayerGraph(empty_layers, empty_interlayers, truncated(Normal(10), 0.0, 20.0));
```

    ┌ Info: Looping through wirings to find one that works...
    └ @ MultilayerGraphs E:\other_drives\developer07\package_development\MultilayerGraphs\dev\MultilayerGraphs\src\utilities.jl:441
    

Note that this is not an implementation of a fully-fledged configuration model, which would require to be able to specify a degree distribution for every dimension of multiplexity. Please refer to [#future-developments]().

There is a similar constructor for `MultilayerDiGraph` which requires both the indegree distribution and the outdegree distribution. Anyway due to current performance limitations in the graph realization algorithms, it is suggested to provide two "similar" distributions (similar mean or location parameter, similar variance or shape parameter), in order not to incur in lengthy computational times.  

#### Nodes

You may add a node via `add_node`:

```julia
new_node = Node("new_node")
add_node!(multilayergraph, new_node) # Return true if succeeds
```
```nothing
true
```

Now one may add vertices that represent that node, e.g.:

```julia
new_vertex = MV(new_node, :layer_sg)
add_vertex!(multilayergraph, new_vertex)
rem_vertex!(multilayergraph, new_vertex) # hide
```
```nothing 
true
```

And remove the node via `rem_node!`:

```julia
rem_node!(multilayergraph, new_node) # Return true if succeeds
```
```nothing 
true
```

#### Modifying edge weight and metadata and vertex metadata

One may modify the weight of the edge of a multilayer graph via the `set_weight!` function. The call will succeed only if the edge that is acted upon exists and belongs to a weighted subgraph:


```julia
# This will succeed
random_weighted_edge = rand(collect(edges(multilayergraph.layer_swg)))
set_weight!(multilayergraph, src(random_weighted_edge), dst(random_weighted_edge), rand())
```
```nothing
true
```

```julia
# This will not succeed
random_unweighted_edge = rand(collect(edges(multilayergraph.layer_sg)))
set_weight!(multilayergraph, src(random_unweighted_edge), dst(random_unweighted_edge), rand())
```
```nothing 
false
```

Equivalent arguments can be made for [`set_metadata!`](@ref) (both vertex and edge dispatches).

#### Adding, Removing, Modifying and Accessing layers and interlayers
One may of course add layers on the fly:

```julia
# Instantiate a new Layer
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
new_layer = Layer(  :new_layer,
                    sample(multilayervertices, _nv, replace = false),
                    _ne, 
                    SimpleGraph{vertextype}(),
                    _weighttype
)

# Add the Layer
add_layer!(
            multilayergraph,                                # the `Multilayer(Di)Graph` which the new layer will be added to;
            new_layer;                                      # the new `Layer` to add to the `multilayergraph`
            default_interlayers_null_graph = SimpleGraph{vertextype}(), # upon addition of a new `Layer`, all the `Interlayer`s between the new and the existing `Layer`s are immediately created. This keyword argument specifies their `null_graph` See the `Layer` constructor for more information. Defaults to `SimpleGraph{T}()`
            default_interlayers_structure = "empty"         # The structure of the `Interlayer`s created by default. May either be "multiplex" to have diagonally-coupled only interlayers, or "empty" for empty interlayers. Defaults to "multiplex".
)

# Check that the new layer now exists within the multilayer graph
has_layer(multilayergraph, :new_layer)
```
```nothing 
true
```

The `add_layer!` function will automatically instantiate all the `Interlayer`s between the newly added `Layer` and the `Layer`s already present in the multilayer graph.

If you wish to manually specify an interlayer, just do:

```julia
# Instantiate a new Interlayer. Notice that its name will be given by default as 
_ne = rand_ne_interlayer(layer_sg, new_layer)
new_interlayer = Interlayer(    layer_sg,                
                                new_layer,               
                                _ne,                     
                                SimpleGraph{vertextype}(),
                                name = :new_interlayer
)

# Modify an existing interlayer with the latter i.e. specify the latter interlayer:
specify_interlayer!( multilayergraph,
                     new_interlayer)

# Now the interlayer between `layer_sg` and `new_layer` is `new_interlayer`
```
```nothing
true
```

Suppose that, after some modifications of `multilayergraph`, you would like to inspect a particular slice (or subgraph) of it (i.e. a `Layer` or an `Interlayer`). You may get both layers and interlayers as properties of the multilayer graph itself.

```julia
# Get a layer by name 
multilayergraph.new_layer
```
```nothing
Layer{Int64, Float64, SimpleGraph{Int64}}(LayerDescriptor{Int64, Float64, SimpleGraph{Int64}}(:new_layer, SimpleGraph{Int64}(0, Vector{Int64}[]), MultilayerGraphs.var"#54#60"(), MultilayerGraphs.var"#55#61"(), MultilayerGraphs.var"#56#62"()), SimpleGraph{Int64}(10, [[3, 5, 7], [4, 5, 6], [1, 6], [2, 7], [1, 2, 6], [2, 3, 5, 7], [1, 4, 6]]), Bijection{Int64,MultilayerVertex{:new_layer}} (with 7 pairs))
```

```julia
# Get an Interlayer by name
multilayergraph.new_interlayer
```
```nothing 
Interlayer{Int64, Float64, SimpleGraph{Int64}}(InterlayerDescriptor{Int64, Float64, SimpleGraph{Int64}}(:new_interlayer, :layer_sg, :new_layer, SimpleGraph{Int64}(0, Vector{Int64}[]), MultilayerGraphs.var"#92#96"(), MultilayerGraphs.var"#93#97"(), false), SimpleGraph{Int64}(19, [[7, 12], [7, 8, 10, 11, 13], [7, 10], [7, 13], [7, 9, 11, 13], [7, 8, 12, 13], [1, 2, 3, 4, 5, 6], [2, 6], [5], [2, 3], [2, 5], [1, 6], [2, 4, 5, 6]]), Bijection{Int64,MultilayerVertex} (with 13 pairs))
```

`Interlayer`s may also be accessed by remembering the names of the `Layer`s they connect:

```julia
# Get an Interlayer from the names of the two layers that it connects
get_interlayer(multilayergraph, :new_layer, :layer_sg )
```
```nothing 
Interlayer{Int64, Float64, SimpleGraph{Int64}}(InterlayerDescriptor{Int64, Float64, SimpleGraph{Int64}}(:new_interlayer_rev, :new_layer, :layer_sg, SimpleGraph{Int64}(0, Vector{Int64}[]), MultilayerGraphs.var"#92#96"(), MultilayerGraphs.var"#93#97"(), false), SimpleGraph{Int64}(19, [[8, 9, 10, 11, 12, 13], [9, 13], [12], [9, 10], [9, 12], [8, 13], [9, 11, 12, 13], [1, 6], [1, 2, 4, 5, 7], [1, 4], [1, 7], [1, 3, 5, 7], [1, 2, 6, 7]]), Bijection{Int64,MultilayerVertex} (with 13 pairs))
```

**NB:** Although the interlayer from an arbitrary `layer_1` to `layer_2` is the same mathematical object as the interlayer from `layer_2` to `layer_1`, their representations as `Interlayer`s differ in the internals, and most notably in the order of the vertices. The `Interlayer` from `layer_1` to `layer_2` orders its internal vertices label so that the `MultilayerVertex`s of `layer_1` (in the order they were in `layer_1` when the `Interlayer` was instantiated) come before the `MultilayerVertex`s of `layer_2` (in the order they were in `layer_2` when the `Interlayer` was instantiated).

When calling `get_interlayer(multilayergraph, :layer_1, :layer_2)` it is returned the `Interlayer` from `layer_1` to `layer_2`. If the Interlayer from `layer_2` to `layer_1` was manually specified or automatically generated during  during the instantiation of the multilayer graph with name, say, `"some_interlayer"`, then the returned `Interlayer` will be named `"some_interlayer_rev"`.

To remove a layer:

```julia
# Remove the layer. This will also remove all the interlayers associated to it.
rem_layer!( multilayergraph,
            :new_layer;
            remove_nodes = false # Whether to also remove all nodes represented in the to-be-removed layer from the multilayer graph
)
```
```nothing
true
```

Visit the **Layers and Interlayers** subsection of the [end-user]() and [developer]() APIs to discover more useful methods.

#### Weight/Adjacency Tensor, Metadata Tensor and Supra Weight/Adjacency Matrix

One may extract the weight tensor of a `multilayergraph` via:

```julia
wgt = weight_tensor(multilayergraph)
```
```nothing
WeightTensor{Float64}([0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0;;; 0.0 1.0 … 0.0 0.0; 1.0 1.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 1.0 … 1.0 0.0;;; 1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 1.0 0.0; 0.0 0.0 … 0.0 0.0;;; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0;;;; 0.0 1.0 … 0.0 0.0; 1.0 1.0 … 0.0 1.0; … ; 0.0 0.0 … 0.0 1.0; 0.0 0.0 … 0.0 0.0;;; 0.0 0.9828581516714545 … 0.0 0.7736606234481569; 0.9828581516714545 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.7736606234481569 0.0 … 0.0 0.0;;; 0.19989407135094706 0.0 … 0.0 0.16650315003660054; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.6643136923736794 … 0.0 0.0; 0.40879510776523964 0.7458197468092816 … 0.0 0.0;;; 1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 1.0;;;; 1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 1.0 0.0; 0.0 0.0 … 0.0 0.0;;; 0.19989407135094706 0.0 … 0.0 0.40879510776523964; 0.0 0.0 … 0.6643136923736794 0.7458197468092816; … ; 0.0 0.0 … 0.0 0.0; 0.16650315003660054 0.0 … 0.0 0.0;;; 0.0 0.0 … 0.0 1.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 1.0 0.0 … 0.0 0.0;;; 1.0 0.0 … 1.0 1.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 1.0 0.0; 0.0 0.0 … 1.0 1.0;;;; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0;;; 1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 1.0;;; 1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 1.0 0.0 … 1.0 1.0; 1.0 0.0 … 0.0 1.0;;; 0.0 0.0 … 1.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 1.0 0.0 … 0.0 1.0; 0.0 0.0 … 1.0 0.0], [:layer_sg, :layer_swg, :layer_mg, :layer_vg], Bijection{Int64,Node} (with 7 pairs))
```

Note that `wgt` is an object of type [`WeightTensor`](@ref). You may access its array representation using:

```julia
array(wgt)
; # hide
```

Also, you may index it using `MultilayerVertex`s:

```julia
# Get two random vertices from the MultilayerGraph
mv1, mv2 = rand(mv_vertices(multilayergraph), 2)

# Get the strength of the edge between them (0 for no edge):
wgt[mv1, mv2]
```
```nothing
0.0
```

Similarly, there is a [`MetadataTensor`](@ref), that may be created via `metadata_tensor(multilayergraph)`

The package also exports a [`SupraWeightMatrix`](@ref) which is a supra (weighted) adjacency matrix with the same indexing functionality as above. You may instantiate it via `supra_weight_matrix(multilayergraph)`.

#### Multilayer-specific analytical tools

Read a complete list of analytical methods exclusive to multilayer graphs in the dedicated [API section](@ref API) (here "exclusive" means that wither those methods do not exists for standard graphs, or that they had to be reimplemented and so may present some caveats). Refer to their docstrings for more information.

### Future Developments

- [Implement faster graph realization algorithms](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/32);
- [Implement more general configuration models / graph generators](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/33);
- [Implement graph of layers](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/34);
- [Implement projected monoplex and overlay graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/35);
- [Implement more default multilayer graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/36) (e.g. multiplex graphs).

## How to Contribute

The package is currently under development and further steps would benefit enormously from the precious feedback of the [JuliaGraph people](https://github.com/orgs/JuliaGraphs/people), graph theorists, network scientists and all the users who might have general questions or suggestions. 

Therefore feel free to open [discussions](https://github.com/JuliaGraphs/MultilayerGraphs.jl/discussions), [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues) or [PRs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/pulls). They are very welcome!   

## How to Cite

If you use this package in your work, please cite this repository using the metadata in [`CITATION.bib`](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/CITATION.bib).

## Announcements 

- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-a-package-to-construct-handle-and-analyse-multilayer-graphs/85988)
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-a-package-to-construct-handle-and-analyse-multilayer-graphs-3k22)
- [Twitter](https://twitter.com/In_Phy_T/status/1560594513189638146)

## References

1. De Domenico et al. (2013) [Mathematical Formulation of Multilayer Networks](https://doi.org/10.1103/PhysRevX.3.041022). *Physical Review X*; 
2. Kivelä et al. (2014) [Multilayer networks](https://doi.org/10.1093/comnet/cnu016). *Journal of Complex Networks*; 
3. Bianconi (2018) [Multilayer Networks: Structure and Function](https://global.oup.com/academic/product/multilayer-networks-9780192865540). *Oxford University Press*.