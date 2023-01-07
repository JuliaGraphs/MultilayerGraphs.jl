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

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs. 

## Overview

A multilayer graph is a graph consisting of multiple standard subgraphs called *layers* which can be interconnected through [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) called *interlayers* composed of the vertex sets of two different layers and the edges between them. The vertices in each layer represent a single set of nodes, although not all nodes have to be represented in every layer. 

Formally, a multilayer graph can be defined as a triple $G=(V,E,L)$, where:

- $V$ is the set of vertices;
- $E$ is the set of edges, pairs of nodes $(u, v)$ representing a connection, relationship or interaction between the nodes $u$ and $v$;
- $L$ is a set of layers, which are subsets of $V$ and $E$ encoding the nodes and edges within each layer.

Each layer $\ell$ in $L$ is a tuple $(V_\ell, E_\ell)$, where $V_\ell$ is a subset of $V$ that represents the vertices within that layer, and $E_\ell$ is a subset of $E$ that represents the edges within that layer.

Multiple theoretical frameworks have been proposed to formally subsume all instances of multilayer graphs [^1]. 

Multilayer graphs have been adopted to model the structure and dynamics of a wide spectrum of high-dimensional, non-linear, multi-scale, time-dependent complex systems including physical, chemical, biological, neuronal, socio-technical, epidemiological, ecological and economic networks [^2]. 

MultilayerGraphs.jl is an integral part of the [JuliaGraphs](https://github.com/JuliaGraphs) ecosystem extending [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) so all the methods and metrics exported by Graphs.jl work for multilayer graphs, but due to the special nature of multilayer graphs the package features a peculiar implementation that maps a standard integer-labelled vertex representation to a more user-friendly framework exporting all the objects an experienced practitioner would expect such as nodes ([`Node`](@ref)), vertices ([`MultilayerVertex`](@ref)), layers ([`Layer`](@ref)), interlayers ([`Interlayer`](@ref)), etc.

`MultilayerGraph` and `MultilayerDiGraph` are fully-fledged Graphs.jl extensions. Both structs are designed to allow for layers and interlayers of any type (as long as they are Graphs.jl extensions themselves) and to permit layers and interlayers of different types. However, it is required that all layers and interlayers in `MultilayerGraph` are undirected, and all layers and interlayers in `MultilayerDiGraph` are directed.

`MultilayerGraph` and `MultilayerDiGraph` support the specification of vertex and edge metadata, provided that the underlying layer or interlayer also supports metadata.

The documentation is organized as follows: you will find a comprehensive [Tutorial](@ref) below, complemented by an [API](@ref) page. The API page is organized in two sections: the [End-User](@ref) section lists all the methods intended for the user who does not need to write code that is also compatible with other libraries in the Graphs.jl's ecosystem, while the [Developer](@ref) section contains methods that allow MultilayerGraphs.jl to be used as any package that extend Graphs.jl . Bot section are further stratified by topic.
The tutorial below will be focused on the end-used experience, as developer methods often have very similar signature and will be better addressed in a future developer-oriented guide, should the community manifest the need of it.

## Installation

To install MultilayerGraphs.jl it is sufficient to activate the `pkg` mode by pressing `]` in the Julia REPL and then run the following command:

```nothing
pkg> add MultilayerGraphs
```

## Tutorial

Here we illustrate how to define, handle and analyse a `MultilayerGraph` (the directed version is completely analogous).

```julia
using Revise
using StatsBase, Distributions
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs
using MultilayerGraphs
```

Define some constants that will prove useful later in the tutorial:

```julia
# Set the minimum and maximum number of nodes_list and edges for random graphs
const vertextype   = Int64
const _weighttype  = Float64
const min_vertices = 5
const max_vertices = 7
const n_nodes      = max_vertices
```

Next we define the list of immutable objects that are represented (through vertices, see below) in the various layers and interlayers of a multilayer graph. These objects are called [`Node`](@ref)s. The constructor for a `Node` reads:

```julia
Node(
    id::String
)
```

Where `id` is a `String` that is the name of what the `Node` stands for (could be cities in a transportation network, users in a social network, etc.). Let's construct a list of `Node`s to use in the remainder of the tutorial:

```julia
# The constructor for nodes (which are immutable) only requires a name (`id`) for the node
const nodes_list = [Node("node_$i") for i in 1:n_nodes]
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

You may access (but not modify) the `id` of a `Node` via the [`id`](@ref) function. `Node`s are represented throughout layers and interlayers via a struct named [`MultilayerVertex`](@ref). It has several convenience constructors, the most complete of them reads:

```julia
MultilayerVertex( 
                node::Node,                            # The Node that the vertex will represent     
                layer::Union{Nothing,Symbol},          # The layer which the `Node` will be represented in. Should be set to `nothing` when constructing layers.
                metadata::Union{<:NamedTuple,<:Tuple} # The metadata associated to this vertex
)
```

Let's contruct a list of `MultilayerVertex`s to use in the remainder of the tutorial:

```julia
## Convert nodes to multilayer vertices without metadata
const multilayervertices = MV.(nodes_list)
## Convert nodes multilayer vertices with metadata
const multilayervertices_meta  = [MV(node, ("I'm node $(node.id)",)) for node in nodes_list] # `MV` is an alias for `MultilayerVertex`
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

This conversion from `Node`s to `MultilayerVertex`s is performed since it is logical to add vertices to a graph, not nodes, and also for consistency reasons with the ecosystem. Regarding layers, adding/removing nodes or vertices allows for selecting the most comfortable interface. A similar mechanism is implemented for edges (see below).

Printing a `MultilayerVertex` returns:

```julia
multilayervertices_meta[1]
```
```nothing
MV(Node("node_1"), :nothing, ("I'm node node_1",))
```

Where `MV` is an alias for `MultilayerVertex`. The first field is the `Node` being represented (accessible via the [`node`](@ref) function), the second the (name of) the layer the vertex is represented in (accessible via the [`layer`](@ref) function, here it is set to `nothing`, since these vertices are yet to be assigned), and the metadata associated to the vertex (accessible via the [`metadata`](@ref) function, no metadata are currently represented via an empty `NamedTuple`). `MultilayerVertex` metadata can be represented via a `Tuple` or a `NamedTuple` (see below for examples). For a complete list of methods applicable to `MultilayerVertices`, please refer to the [Vertices](@ref vertices_eu) of the API.

### Layers

As said before, to define a multilayer graph we need to specify its layers and interlayers. Layers and the individual graphs that make up multilayer graphs. We proceed by constructing a [`Layer`](@ref) using the constructor that randomly specifies the edges:

```julia
Layer(
    name::Symbol,                                                     # The name of the layer
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},  # The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
    ne::Int64,                                                        # The number of edges of the Layer
    null_graph::AbstractGraph{T},                                     # The Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
    weighttype::Type{U};                                              # The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
    default_vertex_metadata::Function = mv -> NamedTuple(),           # Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. Defaults to `mv -> NamedTuple()`;
    default_edge_weight::Function = (src, dst) -> nothing,            #  Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
    default_edge_metadata::Function = (src, dst) -> NamedTuple(),     # Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
    allow_self_loops::Bool = false                                    # whether to allow self loops to be generated or not. Defaults to `false`.
) where {T<:Integer, U<: Real}
```

A `Layer` is considered "weighted" if its underlying graph (`null_graph` argument) has been given the `IsWeighted` trait (traits throughout this package are implemented via [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl), just like Graphs.jl does). Since one may at any moment add a new weighted `Layer` to a `MultilayerGraph` (see below for details), the latter is always considered a "weighted graph", so it is given the `IsWeighted` trait. Thus, all `Layer`s and `Interlayer`s (collectively named "subgraphs" hereafter) must specify their `weighttype` as the last argument of their constructor, so the user may debug their weight matrices ([`weights(subgraph::AbstractSubGraph)`](@ref)) immediately after construction. As better specified below, all subgraphs that are meant to be part of the same `MultilayerGraph` must have the same `weighttype`. Moreover, also the vertex type `T` (i.e. the internal representation of vertices) should be the same.

Before instantiating `Layer`s, we define an utility function to ease randomisation:

```julia
# Utility function that returns a random number of vertices and edges each time it is called:
function rand_nv_ne_layer(min_vertices, max_vertices)
    _nv = rand(min_vertices:max_vertices)
    _ne = rand(1:(_nv*(_nv-1)) ÷ 2 )
    return (_nv,_ne)
end

# Utility function that returns two vertices of a Layer that are not adjacent.
function _get_srcmv_dstmv_layer(layer::Layer)
    mvs = MultilayerGraphs.get_bare_mv.(collect(mv_vertices(layer)))

    src_mv_idx = findfirst(mv -> !isempty(setdiff(
        Set(mvs),
        Set(
            vcat(MultilayerGraphs.get_bare_mv.(mv_outneighbors(layer, mv)), mv)
        ),
    )), mvs)

    src_mv = mvs[src_mv_idx]

    _collection = setdiff(
        Set(mvs),
        Set(
            vcat(MultilayerGraphs.get_bare_mv.(mv_outneighbors(layer, src_mv)), src_mv)
        ),
    )
    
    dst_mv = MultilayerGraphs.get_bare_mv(rand(_collection))

    return mvs, src_mv, dst_mv
end
```

We are now are ready to define some `Layer`s. Every type of graph from the Graphs.jl ecosystem may underlie a `Layer` (or an `Interlayer`). We will construct a few of them, each time with a different number of vertices and edges.

```julia
# An unweighted simple layer:
_nv, _ne  = rand_nv_ne_layer(min_vertices,max_vertices)
layer_sg = Layer(   :layer_sg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    SimpleGraph{vertextype}(),
                    _weighttype
)

# A weighted `Layer`
_nv, _ne  = rand_nv_ne_layer(min_vertices,max_vertices)
layer_swg = Layer(  :layer_swg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    SimpleWeightedGraph{vertextype, _weighttype}(),
                    _weighttype;
                    default_edge_weight = (src,dst) -> rand()
)
# A `Layer` with an underlying `MetaGraph`:
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
layer_mg = Layer(   :layer_mg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    MetaGraph{vertextype, _weighttype}(),
                    _weighttype;
                    default_edge_metadata = (src,dst) -> (from_to = "from_$(src)_to_$(dst)",)
)
# `Layer` with an underlying `ValGraph` from `SimpleValueGraphs.jl`
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
layer_vg = Layer(   :layer_vg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    MultilayerGraphs.ValGraph(SimpleGraph{vertextype}();
                                                edgeval_types=(Float64, String, ),
                                                edgeval_init=(s, d) -> (s+d, "hi"),
                                                vertexval_types=(String,),
                                                vertexval_init=v -> ("$v",),),
                                                _weighttype;
                                                default_edge_metadata = (src,dst) -> (rand(), "from_$(src)_to_$(dst)",),
                                                default_vertex_metadata = mv -> ("This metadata had been generated via the default_vertex_metadata method",)
)

# Collect all layers in an ordered list. Order will be recorded when instantiating the multilayer graph.
layers = [layer_sg, layer_swg, layer_mg, layer_vg]
```

The API that inspects and modifies `Layer`s will be shown below together with that of `Interlayer`s, since they are usually the same.  There are of course other constructors that you may discover by reading the [API](@ref subgraphs_eu). They include:

1. Constructors that exempt the user from having to explicitly specify the `null_graph`, at the cost of some flexibility;
2. Constructors that allow for a configuration model-like specifications.

### Interlayers

Now we turn to defining [`Interlayer`](@ref)s. Interlayers are the graphs containing all the edges between vertices is two distinct layers. As before, we need an utility to ease randomization:


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
    nv_1 = nv(layer_1)
    nv_2 = nv(layer_2)
    _ne = rand(1: (nv_1 * nv_2 - 1) )
    return _ne
end
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
    transfer_vertex_metadata::Bool = false                               # If true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors.;
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
                                interlayer_name = :random_interlayer  # The name of the interlayer. We will be able to access it as a property of the multilayer graph via its name. This kwarg's default value is given by a combination of the two layers' names.
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
                                                    ValGraph(SimpleGraph{vertextype}(); edgeval_types=(from_to = String,), edgeval_init=(s, d) -> (from_to = "from_$(s)_to_$(d)"));
                                                    default_edge_metadata = (x,y) -> (from_to = "from_$(src)_to_$(dst)",)
)
# Finally, An `Interlayer` with no couplings (an "empty" interlayer):
interlayer_empty_sg_vg = empty_interlayer(  layer_sg,
                                            layer_vg,
                                            SimpleGraph{vertextype}()
)

# Collect all interlayers. Even though the list is ordered, order will not matter when instantiating the multilayer graph.
interlayers = [interlayer_sg_swg, interlayer_swg_mg, interlayer_mg_vg, interlayer_multiplex_sg_mg, interlayer_empty_sg_vg]
```

There are of course other constructors that you may discover by reading the [API](@ref subgraphs_eu). They include constructors that exempt the user from having to explicitly specify the `null_graph`, at the cost of some flexibility;

Next, we explore the API associated to modify and analyze `Layer`s and `Interlayer`s.

### Subgraphs API

API for  `Layer`s and `Interlayer`s (collectively, "subgraphs") are very similar, so we will just show them for the `Layer` case, pointing out differences to the `Interlayer` scenario whenever they occur.

Subgraphs extend the Graphs.jl's interface, so one may expect every method from Graphs.jl to apply. Anyway, the output and signature is slightly different and thus worth pointing out below.

#### [Nodes](@ref nodes_tut_subg)

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

#### [Vertices](@ref vertices_tut_subg)

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

The [`vertices(subgraph::AbstractSubGraph)`](@ref) command would return an internal representation of the `MultilayerVertex`s (that of type `vertextype`). This method, together with others, serves to make `MultilayerGraphs.jl` compatible with the Graphs.jl ecosystem, but it is not meant to be called by the end user. It is, anyway, thought to be used by developers who wish to interface their packages with `MultilayerGraphs.jl` just as with other packages of the `Graphs.jl` ecosystem: as said above, a developer-oriented guide will be compiled if there is the need, although docstrings are already completed.

To add a vertex, simply use `add_vertex!`. Let us define a vertex with metadata to add. Since nodes may not be represented more than once in layers, we have to define a new node too:

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
    2. They should work even if we are not able to integrate the *standard* and the *uniform* interface with a particular `Graphs.jl`'s extension.

    Here is an example on how to use it:
```julia
add_vertex!(layer_mg, new_node, Dict(pairs(new_metadata)))
```
where `Dict(pairs(new_metadata))` is exactly what you would pass to the `add_vertex!` method that acts on `MetaGraphs`:
```julia
metagraph = MetaGraph()
add_vertex!(metagraph,  Dict(pairs(new_metadata))) # true
```

If an underlying graph has an `add_vertex!` interface whose signature overlaps with that of the uniform interface, the uniform interface will be prevail.

If, using the *transparent* interface, one does not specify any `metadata`, the `default_vertex_metadata` function passed to the `Layer`'s constructor is called to provide `metadata` to the vertex (type `?Layer` in the REPL for more information).

To remove the vertex, simply do:
```julia
rem_vertex!(layer_sg, new_vertex) # Returns true if succeeds
```

To extract metadata:
```julia
get_metadata(layer_mg, MV(new_node))
```

By design, one may not add nor remove vertices to `Interlayer`s.


### [Edges](@ref edges_tut_subg)

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

Where `ME` is a shorthand for `MultilayerEdge`. Besides the two vertices connected, each `MultilayerEdge` carries the information about its `weight` and `metadata`. For unweighted subgraphs, the weight is just `one(weighttype)` and for non-meta subgraphs the metadata are an empty `NamedTuple`s. See `?MultilayerEdge` for additional information, or refer to the [Edges](@ref edges_eu) to discover more methods related to `MultilayerEdges`s.

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

If, using the *transparent* interface, one does not specify any `weight` or (inclusively) `metadata` keyword argument, the `default_edge_weight` or (inclusively) the `default_edge_metadata` function passed to the subgraph's constructor will be called to provide `weight` or `metadata` to the edge (type `?Layer` in the REPL for more information).

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

A complete list of methods relating to subgraphs can be found [here](@ref subgraphs_eu).

### Multilayer Graphs

Given all the `Layer`s and the `Interlayer`s, let's instantiate a multilayer graph as follows:

```julia
multilayergraph = MultilayerGraph(  layers,                                                 # The (ordered) list of layers the multilayer graph will have
                                    interlayers;                                            # The list of interlayers specified by the user. Note that the user does not need to specify all interlayers, as the unspecified ones will be automatically constructed using the indications given by the `default_interlayers_null_graph` and `default_interlayers_structure` keywords.
                                    default_interlayers_null_graph = SimpleGraph{vertextype}(), # Sets the underlying graph for the interlayers that are to be automatically specified.  Defaults to `SimpleGraph{T}()`, where `T` is the `T` of all the `layers` and `interlayers`. See the `Layer` constructors for more information.
                                    default_interlayers_structure = "multiplex" # Sets the structure of the interlayers that are to be automatically specified. May be "multiplex" for diagonally coupled interlayers, or "empty" for empty interlayers (no edges).  "multiplex". See the `Interlayer` constructors for more information.
)
```
```nothing
`MultilayerGraph` with vertex type `Int64` and weight type `Float64`.

### LAYERS
┌───────────┬────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│   NAME    │                                                                UNDERLYING GRAPH                                                                │
├───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ layer_sg  │                                                               SimpleGraph{Int64}                                                               │
├───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ layer_swg │                                                      SimpleWeightedGraph{Int64, Float64}                                                       │
├───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ layer_mg  │                                                           MetaGraph{Int64, Float64}                                                            │
├───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ layer_vg  │ ValGraph{Int64, Tuple{String}, Tuple{Float64, String}, Tuple{}, Tuple{Vector{String}}, Tuple{Vector{Vector{Float64}}, Vector{Vector{String}}}} │
└───────────┴────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘


### INTERLAYERS
┌───────────────────────────────┬───────────┬───────────┬────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┬──────────────────────────┐
│             NAME              │  LAYER 1  │  LAYER 2  │                                                              UNDERLYING GRAPH                                                              │ TRANSFER VERTEX METADATA │
├───────────────────────────────┼───────────┼───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼──────────────────────────┤
│ interlayer_layer_sg_layer_swg │ layer_sg  │ layer_swg │                                                             SimpleGraph{Int64}                                                             │          false           │
├───────────────────────────────┼───────────┼───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼──────────────────────────┤
│ interlayer_layer_sg_layer_mg  │ layer_sg  │ layer_mg  │ ValGraph{Int64, Tuple{}, NamedTuple{(:from_to,), Tuple{String}}, Tuple{}, Tuple{}, NamedTuple{(:from_to,), Tuple{Vector{Vector{String}}}}} │          false           │
├───────────────────────────────┼───────────┼───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼──────────────────────────┤
│ interlayer_layer_swg_layer_mg │ layer_swg │ layer_mg  │                                                    SimpleWeightedGraph{Int64, Float64}                                                     │          false           │
├───────────────────────────────┼───────────┼───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼──────────────────────────┤
│ interlayer_layer_sg_layer_vg  │ layer_sg  │ layer_vg  │                                                             SimpleGraph{Int64}                                                             │          false           │
├───────────────────────────────┼───────────┼───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼──────────────────────────┤
│ interlayer_layer_vg_layer_swg │ layer_vg  │ layer_swg │                                                             SimpleGraph{Int64}                                                             │          false           │
├───────────────────────────────┼───────────┼───────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┼──────────────────────────┤
│ interlayer_layer_mg_layer_vg  │ layer_mg  │ layer_vg  │                                                         MetaGraph{Int64, Float64}                                                          │           true           │
└───────────────────────────────┴───────────┴───────────┴────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┴──────────────────────────┘
```

**Keep in mind that `Multilayer(Di)Graph` only supports uniform and standard interface for both `add_vertex!` and `add_edge!`.**

As already stated, a `MultilayerGraph` is an object made of `Layer`s and `Interlayer`s whose collections of vertices each represents a subset of the set of nodes, here being the variable `nodes`.

*Adding* a `Node` to a `MultilayerGraph` will enable its `Layer`s  (and thus its `Interlayer`s) to represent it i.e. you will be able to add `MultilayerVertex`s that represent that `Node` to the multilayer graph.

Another constructor allows for a limited configuration model-like specification. It allows to generate a multilayer graph with a specific degree distribution. Since edges are created according to the provided distribution, it is necessary that the layers and interlayers specified are empty (i.e. they have no edges). Notice that layers and interlayers whose underlying graph is a `SimpleWeighted(Di)Graph` may not be used until [this PR](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl/pull/14) is merged.

It is used as:


```julia
# The configuration model-like constructor will be responsible for creating the edges, so we need to provide it with empty layers and interlayers.
# To create empty layers and interlayers, we will empty the above subgraphs, and, for compatibility reasons, we'll remove the ones having a `SimpleWeightedGraph`s. These lines are not necessary to comprehend the tutorial, they may be skipped. Just know that the variables `empty_layers` and `empty_interlayers` are two lists of, respectively, empty layers and interlayers that do not have `SimpleWeightedGraph`s as their underlying graphs

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

Note that this is not an implementation of a fully-fledged configuration model, which would require to be able to specify a degree distribution for every dimension of multiplexity. Moreover, it still lacks the ability to specify a minimum discrepancy (w.r.t. a yet-to-be-chosen metric) between the empirical distributions of the sampled sequence and the provided theoretical distribution. Please refer to [Future Developments](@ref).

There is a similar constructor for `MultilayerDiGraph` which requires both the indegree distribution and the outdegree distribution. Anyway due to current performance limitations in the graph realization algorithms, it is suggested to provide two "similar" distributions (similar mean or location parameter, similar variance or shape parameter), in order not to incur in lengthy computational times.  

#### [Nodes](@ref nodes_tut_multig)

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

One may modify the weight of the edge of a multilayer graph via the [`set_weight!`](@ref) function. The call will succeed only if the edge that is acted upon exists and belongs to a weighted subgraph:


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
                    sample(nodes_list, _nv, replace = false),
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

The [`add_layer!`](@ref) function will automatically instantiate all the `Interlayer`s between the newly added `Layer` and the `Layer`s already present in the multilayer graph, according to its kwargs `default_interlayers_null_graph` and `default_interlayers_structure`.

If you wish to manually specify an interlayer, just do:

```julia
# Instantiate a new Interlayer. Notice that its name will be given by default as
_ne = rand_ne_interlayer(layer_sg, new_layer)
new_interlayer = Interlayer(    layer_sg,                
                                new_layer,               
                                _ne,                     
                                SimpleGraph{vertextype}(),
                                interlayer_name = :new_interlayer
)

# Modify an existing interlayer with the latter i.e. specify the latter interlayer:
specify_interlayer!( multilayergraph,
                     new_interlayer)

# Now the interlayer between `layer_sg` and `new_layer` is `new_interlayer`
```
```nothing
true
```

Suppose that, after some modifications of `multilayergraph`, you would like to inspect a particular slice (or subgraph) of it (i.e. a `Layer` or an `Interlayer`). You may use both layers and interlayers names as properties of the multilayer graph itself.

```julia
# Get a layer by name
multilayergraph.new_layer
```
```nothing
Layer   new_layer
underlying_graph: SimpleGraph{Int64}
vertex type: Int64
weight type: Float64
nv = 7
ne = 11
```

```julia
# Get an Interlayer by name
multilayergraph.new_interlayer
```
```nothing
Interlayer      new_interlayer
layer_1: layer_sg
layer_2: new_layer
underlying graph: SimpleGraph{Int64}
vertex type : Int64
weight type : Float64
nv : 14
ne : 46
```

`Interlayer`s may also be accessed by remembering the names of the `Layer`s they connect:

```julia
# Get an Interlayer from the names of the two layers that it connects
get_interlayer(multilayergraph, :new_layer, :layer_sg )
```
```nothing
Interlayer      new_interlayer_rev
layer_1: new_layer
layer_2: layer_sg
underlying graph: SimpleGraph{Int64}
vertex type : Int64
weight type : Float64
nv : 14
ne : 46
```

**NB:** Although the interlayer from an arbitrary `layer_1` to `layer_2` is the same mathematical object as the interlayer from `layer_2` to `layer_1`, their representations as `Interlayer`s differ in the internals, and most notably in the order of the vertices. The `Interlayer` from `layer_1` to `layer_2` orders its vertices so that the `MultilayerVertex`s of `layer_1` (in the order they were in `layer_1` when the `Interlayer` was instantiated) come before the `MultilayerVertex`s of `layer_2` (in the order they were in `layer_2` when the `Interlayer` was instantiated).

When calling `get_interlayer(multilayergraph, :layer_1, :layer_2)` it is returned the `Interlayer` from `layer_1` to `layer_2`. If the Interlayer from `layer_2` to `layer_1` was manually specified or automatically generated during the instantiation of the multilayer graph with name, say, `"some_interlayer"`, then the returned `Interlayer` will be named `"some_interlayer_rev"`.

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

Visit the [multilayer graph](@ref msm_eu) subsection of the edn-user APIs to discover more useful methods.

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

Read a complete list of analytical methods exclusive to multilayer graphs in the dedicated [API section](@ref msm_eu) (here "exclusive" means that wither those methods do not exists for standard graphs, or that they had to be reimplemented and so may present some caveats). Refer to their docstrings for more information.

#### Compatibility with Agents.jl

`Multilayer(Di)Graph`s may be used as an argument to `GraphSpace` in [Agents.jl](https://github.com/JuliaDynamics/Agents.jl). A complete compatibility example may be found in [this test](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/test/agents_jl_integration.jl).

### Future Developments

- [Implement graph of layers](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/34);
- [Implement projected monoplex and overlay graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/35);
- [Implement more default multilayer graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/36) (e.g. multiplex graphs);
- [Implement configuration models / graph generators for interlayers](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/46);
- [Implement a fully-fledged multilayer configuration model / graph generator](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/48);
- [Relax the requirement of same `T` and `U` for all `Layer`s and `Interlayer`s that are meant to constitute a `Multilayer(Di)Graph`](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/53);
- [Implement multilayer graph data visualisation functionalities](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/54);
- [Infer `weighttype` from `default_edge_weight`](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/58);
- [Improve error explanations](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/59); 
- [Improve integration with Agents.jl](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/61);
- [Allow configuration models to specify a minimum discrepancy between the sampled (di)graphical sequence(s) and the provided distribution](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/62);
- [Add to `add_layer!` a kwarg that allows the user to specify some new interlayers, skipping the instantiation of the default ones.](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/63).

## How to Contribute

The ongoing development of this package would greatly benefit from the valuable feedback of the esteemed members of the [JuliaGraph](https://github.com/orgs/JuliaGraphs/people) community, as well as from graph theorists, network scientists, and any users who may have general questions or suggestions.

We therefore encourage you to participate in [discussions](https://github.com/JuliaGraphs/MultilayerGraphs.jl/discussions), raise [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues), or submit [pull requests](https://github.com/JuliaGraphs/MultilayerGraphs.jl/pulls). Your contributions are most welcome!

## How to Cite

If you utilize this package in your project, please consider citing this repository using the citation information provided in [`CITATION.bib`](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/CITATION.bib). 

This will help to give appropriate credit to the [contributors](https://github.com/JuliaGraphs/MultilayerGraphs.jl/graphs/contributors) and support the continued development of the package.

## Announcements

The package and its features were announced on the following platforms:

- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-a-package-to-construct-handle-and-analyse-multilayer-graphs/85988)
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-a-package-to-construct-handle-and-analyse-multilayer-graphs-3k22)
- [Twitter](https://twitter.com/In_Phy_T/status/1560594513189638146)

## Related Packages 

### R 

Here is a list of software packages for the creation, manipulation, analysis and visualisation of multilayer graphs implemented in the [R language](https://www.r-project.org): 

- [`muxViz`](https://github.com/manlius/muxViz) implements functions to perform multilayer correlation analysis, multilayer centrality analysis, multilayer community structure detection, multilayer structural reducibility, multilayer motifs analysis and utilities to statically and dynamically visualise multilayer graphs;
- [`multinet`](https://github.com/cran/multinet) implements functions to import, export, create and manipulate multilayer graphs, several state-of-the-art multiplex graph analysis algorithms for centrality measures, layer comparison, community detection and visualization;
- [`mully`](https://github.com/frankkramer-lab/mully) implements functions to import, export, create, manipulate and merge multilayer graphs and utilities to visualise multilayer graphs in 2D and 3D;
- [`multinets`](https://github.com/neylsoncrepalde/multinets) implements functions to import, export, create, manipulate multilayer graphs and utilities to visualise multilayer graphs.

### Python

Here is a list of software packages for the creation, manipulation, analysis and visualisation of multilayer graphs implemented in the [Python language](https://www.python.org): 

- [`MultiNetX`](https://github.com/nkoub/multinetx) implements methods to create undirected networks with weighted or unweighted links, to analyse the spectral properties of adjacency or Laplacian matrices and to visualise multilayer graphs and dynamical processes by coloring the nodes and links accordingly;
- [`PyMNet`](https://github.com/bolozna/Multilayer-networks-library) implements data structures for multilayer graphs and multiplex graphs, methods to import, export, create, manipulate multilayer graphs and for the rule-based generation and lazy-evaluation of coupling edges and utilities to visualise multilayer graphs.

### Julia 

At the best of our knowledge there are currently no software packages dedicated to the creation, manipulation and analysis of multilayer graphs implemented in the [Julia language](https://julialang.org) apart from MultilayerGraphs.jl itself.

## References

1. De Domenico et al. (2013) [Mathematical Formulation of Multilayer Networks](https://doi.org/10.1103/PhysRevX.3.041022). *Physical Review X*; 
2. Kivelä et al. (2014) [Multilayer networks](https://doi.org/10.1093/comnet/cnu016). *Journal of Complex Networks*; 
3. Boccaletti et al. (2014) [The structure and dynamics of multilayer networks](https://doi.org/10.1016/j.physrep.2014.07.001). *Physics Reports*; 
4. Lee et al. (2015) [Towards real-world complexity: an introduction to multiplex networks](https://doi.org/10.1140/epjb/e2015-50742-1). *The European Physical Journal B*; 
5. Bianconi (2018) [Multilayer Networks: Structure and Function](https://global.oup.com/academic/product/multilayer-networks-9780192865540). *Oxford University Press*;
6. Cozzo et al. (2018) [Multiplex Networks: Basic Formalism and Structural Properties](https://doi.org/10.1007/978-3-319-92255-3). *SpringerBriefs in Complexity*; 
7. Aleta and Moreno (2019) [Multilayer Networks in a Nutshell](https://doi.org/10.1146/annurev-conmatphys-031218-013259). *Annual Review of Condensed Matter Physics*; 
8. Artime et al. (2022) [Multilayer Network Science: From Cells to Societies](https://doi.org/10.1017/9781009085809). *Cambridge University Press*; 
9. De Domenico (2022) [Multilayer Networks: Analysis and Visualization](https://doi.org/10.1007/978-3-030-75718-2). *Springer Cham*. 

[^1]: [De Domenico  et al. (2013)](https://doi.org/10.1103/physrevx.3.041022); [Kivelä et al. (2014)](https://doi.org/10.1093/comnet/cnu016); [Boccaletti et al. (2014)](https://doi.org/10.1016/j.physrep.2014.07.001); [Lee et al. (2015)](https://doi.org/10.1140/epjb/e2015-50742-1); [Aleta and Moreno (2019)](https://doi.org/10.1146/annurev-conmatphys-031218-013259); [Bianconi (2018)](https://doi.org/10.1093/oso/9780198753919.001.0001); [Cozzo et al. (2018)](https://doi.org/10.1007/978-3-319-92255-3); [Artime et al. (2022)](https://doi.org/10.1017/9781009085809); [De Domenico (2022)](https://doi.org/10.1007/978-3-030-75718-2).

[^2]: [Cozzo et al. (2013)](https://doi.org/10.1103/physreve.88.050801); [Granell et al. (2013)](https://doi.org/10.1103/physrevlett.111.128701); [Massaro and Bagnoli (2014)](https://doi.org/10.1103/physreve.90.052817); [Estrada and Gomez-Gardenes (2014)](https://doi.org/10.1103/physreve.89.042819); [Azimi-Tafreshi (2016)](https://doi.org/10.1103/physreve.93.042303); [Baggio et al. (2016)](https://doi.org/10.1073/pnas.1604401113); [DeDomenico et al. (2016)](https://doi.org/10.1038/nphys3865); [Amato et al. (2017)](https://doi.org/10.1038/s41598-017-06933-2); [DeDomenico (2017)](https://doi.org/10.1093/gigascience/gix004); [Pilosof et al. (2017)](https://doi.org/10.1038/s41559-017-0101); [de Arruda et al. (2017)](https://doi.org/10.1103/physrevx.7.011014); [Gosak et al. (2018)](https://doi.org/10.1016/j.plrev.2017.11.003); [Soriano-Panos et al. (2018)](https://doi.org/10.1103/physrevx.8.031039); [Timteo et al. (2018)](https://doi.org/10.1038/s41467-017-02658-y); [Buldú et al. (2018)](https://doi.org/10.1162/netn_a_00033); [Lim et al. (2019)](https://doi.org/10.1038/s41598-019-39243-w); [Mangioni et al. (2020)](https://doi.org/10.1109/tnse.2018.2871726); [Aleta et al. (2020)](https://doi.org/10.1038/s41562-020-0931-9); [Aleta et al. (2022)](https://doi.org/10.1073/pnas.2112182119)).
