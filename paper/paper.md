---
title: 'MultilayerGraphs.jl: Multilayer Network Science in Julia'
tags:
  - Julia
  - Graphs
  - Networks
  - Complexity
  - Complex-Networks
  - Complex-Systems
  - Network-Science
  - Network-Analysis
  - Multilayer-Graphs
  - Multilayer-Networks
authors:
  - name: Claudio Moroni
    orcid: 0000-0003-1274-6937
    equal-contrib: true
    affiliation: "1, 2" 
  - name: Pietro Monticone
    orcid: 0000-0002-2731-9623
    equal-contrib: true
    affiliation: "1, 2" 
affiliations:
 - name: University of Turin, Italy
   index: 1
 - name: Interdisciplinary Physics Team, Italy
   index: 2
date: 30 December 2022
bibliography: paper.bib
---

# Summary

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs.

A multilayer graph consists of multiple subgraphs called *layers* which can be interconnected through [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) called *interlayers* composed of the vertex sets of two different layers and the edges between them. The vertices in each layer represent a single set of nodes, although not all nodes have to be represented in every layer. 

Formally, a multilayer graph can be defined as a triple $G=(V,E,L)$, where:

- $V$ is the set of vertices;
- $E$ is the set of edges, pairs of nodes $(u, v)$ representing a connection, relationship or interaction between the nodes $u$ and $v$;
- $L$ is a set of layers, which are subsets of $V$ and $E$ encoding the nodes and edges within each layer.

Each layer $\ell$ in $L$ is a tuple $(V_\ell, E_\ell)$, where $V_\ell$ is a subset of $V$ that represents the vertices within that layer, and $E_\ell$ is a subset of $E$ that represents the edges within that layer.

MultilayerGraphs.jl is an integral part of the [JuliaGraphs](https://github.com/JuliaGraphs) ecosystem extending Graphs.jl [@Graphs2021] so all the methods and metrics exported by Graphs.jl work for multilayer graphs, but due to the special nature of multilayer graphs the package features a peculiar implementation that maps a standard integer-labelled vertex representation to a more user-friendly framework exporting all the objects an experienced practitioner would expect such as nodes (`Node`), vertices (`MultilayerVertex`), layers (`Layer`), interlayers (`Interlayer`), etc.

MultilayerGraphs.jl features multilayer-specific methods and metrics including the global clustering coefficient, the overlay clustering coefficient, the multilayer eigenvector centrality, the multilayer modularity and the Von Neumann entropy.

Finally, MultilayerGraphs.jl has been integrated within the [JuliaDynamics](https://github.com/JuliaDynamics) ecosystem so that any `Multilayer(Di)Graph` can be utilised as an argument to the `GraphSpace` constructor in Agents.jl [@Datseris2022]. 

# Statement of Need

Several theoretical frameworks have been proposed to formally subsume all instances of multilayer graphs [@DeDomenico2013; @Kivela2014; @Boccaletti2014; @Lee2015; @Aleta2019; @Bianconi2018; @Cozzo2018; @Artime2022; @DeDomenico2022]. 

Multilayer graphs have been adopted to model the structure and dynamics of a wide spectrum of high-dimensional, non-linear, multi-scale, time-dependent complex systems including physical, chemical, biological, neuronal, socio-technical, epidemiological, ecological and economic networks [@Cozzo2013; @Granell2013; @Massaro2014; @Estrada2014; @AzimiTafreshi2016; @Baggio2016; @DeDomenico2016; @Amato2017; @DeDomenico2017; @Pilosof2017; @deArruda2017; @Gosak2018; @SorianoPaos2018; @Timteo2018; @Buldu2018; @Lim2019; @Mangioni2020; @Aleta2020; @Aleta2022]. 

We have chosen the [Julia language](https://julialang.org) for this software package because it is a modern, open-source, high-level, high-performance dynamic language for technical computing [@Bezanson2017]. At the best of our knowledge there are currently no software packages dedicated to the creation, manipulation and analysis of multilayer graphs implemented in the Julia language apart from MultilayerGraphs.jl itself [@Moroni_Monticone_MultilayerGraphs_2022]. 

# Main Features 

The two main data structures are `MultilayerGraph` and `MultilayerDiGraph`: collections of layers connected through interlayers. 

The **vertices** of a multilayer graph are representations of one set of distinct objects called `Node`s. Each layer may represent all the node set or just a subset of it. The vertices of `Multilayer(Di)Graph` are implemented via the `MultilayerVertex` custom type. Each `MultilayerVertex` encodes information about the node it represents, the layer it belongs to and its metadata. 

Both the **intra-layer** and **inter-layer edges** are embedded in the `MultilayerEdge` struct, whose arguments are the two connected multilayer vertices, the edge weight and its metadata. It's important to highlight that `Multilayer(Di)Graph`s are weighted and able to store metadata by default (i.e. they have been assigned the `IsWeighted` and `IsMeta` traits from [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl)). 

The **layers** are implemented via the `Layer` struct composed of an underlying graph and a mapping from its integer-labelled vertices to the collection of `MultilayerVertex`s the layer represents. **Interlayers** are similarly implemented via the `Interlayer` mutable struct, and they are generally constructed by providing the two connected layers, the (multilayer) edge list between them and a graph. This usage of underlying graphs allows for an easier debugging procedure during construction and a more intuitive analysis afterwards allowing the package to leverage all the features of the JuliaGraphs ecosystem so that it can be effectively considered as a real proving ground of its internal consistency.

The `Multilayer(Di)Graph` structs are weighted and endowed with the functionality to store both vertex-level and edge-level metadata by default so that at any moment the user may add or remove a `Layer` or specify an `Interlayer` and since different layers and interlayers could be better represented by graphs that are weighted or unweighted and with or without metadata, it was crucial for us to provide the most general and adaptable structure. A `Multilayer(Di)Graph` is instantiated by providing the ordered list of layers and the list of interlayers to the constructor. The latter are automatically specified, so there is no need to instantiate all of them. 

Alternatively, it is possible to construct a `Multilayer(Di)Graph` making use of a graph generator-like signature allowing the user to set the degree distribution or the degree sequence and employs graph realisation methods such as the Havel-Hakimi algorithm for undirected graphs [@Hakimi1962] and the Kleitman-Wang algorithm for directed ones [@Kleitman1973]. 

`Multilayer(Di)Graph`s structure may be represented via dedicated `WeightTensor`, `MetadataTensor` and `SupraWeightMatrix` structs, all of which support indexing with `MultilayerVertex`s. Once a `Multilayer(Di)Graph` has been instantiated, its layers and interlayers can be accessed as their properties.

# Installation and Usage 

To install MultilayerGraphs.jl it is sufficient to activate the `pkg` mode by pressing `]` in the Julia REPL and then run the following command:

```nothing
pkg> add MultilayerGraphs
```

In the following code chunks we synthetically illustrate some of the main features outlined in the previous section.

Let's begin by importing the necessary dependencies and setting the relevant constants: 

```julia
using Distributions, Graphs, SimpleValueGraphs
using MultilayerGraphs

# Set the number of nodes: objects represented by multilayer vertices
const n_nodes = 100 
# Create a list of nodes
const node_list = [Node("node_$i") for i in 1:n_nodes]
```

## Layers and Interlayers 

We will instantiate layers and interlayers with randomly-selected edges and vertices adopting a variety of techniques.

Here we define a layer with an underlying simple directed graph using a graph generator-like (or "configuration model"-like) constructor which allows us to specify both the **indegree** and the **outdegree sequences**. Before instantiating each layer we sample the number of its vertices and, optionally, of its edges.

```julia
n_vertices = rand(1:100)                          # Number of vertices 
layer_simple_directed = layer_simpledigraph(      # Layer constructor 
    :layer_simpledigraph,                         # Layer name
    sample(node_list, n_vertices; replace=false), # Nodes represented in the layer
    Truncated(Normal(5, 5), 0, 20), # Indegree sequence distribution 
    Truncated(Normal(5, 5), 0, 20)  # Outdegree sequence distribution
)
```

Then we define a layer with an underlying simple weighted directed graph. This is another kind of constructor that allows the user to specify the number of edges to be randomly distributed among the vertices. 

```julia
n_vertices = rand(1:n_nodes)                                   # Number of vertices 
n_edges = rand(n_vertices:(n_vertices * (n_vertices - 1) - 1)) # Number of edges 
layer_simple_directed_weighted = layer_simpleweighteddigraph(  # Layer constructor 
    :layer_simpleweighteddigraph,                              # Layer name
    sample(node_list, n_vertices; replace=false), # Nodes represented in the layer
    n_edges;                                      # Number of randomly distributed edges
    default_edge_weight=(src, dst) -> rand()      # Function assigning weights to edges 
)
```

Similar constructors, more flexible at the cost of ease of use, enable a finer tuning and should be necessary only in rare circumstances, e.g. if the equivalent simplified constructor `layer_simple_directed_value` is not able to infer the correct return types of `default_vertex_metadata` or `default_edge_metadata`, or to use and underlying graph structure that isn't currently supported.

```julia
n_vertices = rand(1:n_nodes)                                   # Number of vertices 
n_edges = rand(n_vertices:(n_vertices * (n_vertices - 1) - 1)) # Number of edges 
default_vertex_metadata = v -> ("vertex_$(v)_metadata")        # Vertex metadata 
default_edge_metadata = (s, d) -> (rand(),)                    # Edge metadata 
layer_simple_directed_value = Layer(                           # Layer constructor
    :layer_simplevaldigraph,                                   # Layer name
    sample(node_list, n_vertices; replace=false),              # Nodes represented in the layer
    n_edges,                                                   # Number of randomly distributed edges
    ValDiGraph(                                                
        SimpleDiGraph{Int64}(); 
        vertexval_types=(String,),
        vertexval_init=default_vertex_metadata,
        edgeval_types=(Float64,),
        edgeval_init=default_edge_metadata,
    ),
    Float64;
    default_vertex_metadata=default_vertex_metadata,           # Vertex metadata 
    default_edge_metadata=default_edge_metadata                # Edge metadata 
)

layers = [layer_simple_directed, layer_simple_directed_weighted, layer_simple_directed_value]
```

There are many more constructors the user is encouraged to explore in the package [documentation](https://juliagraphs.org/MultilayerGraphs.jl).
We may now move to Interlayers. Note that, in order to define a `Multilayer(Di)Graph`, 
interlayers do not need to be explicitly constructed by the user, 
since they are automatically specified by the `Multilayer(Di)Graph` constructor. 
Anyway, more complex interlayers need to be manually instantiated. 
The interface is very similar to the layers.

```julia
#=
There are many more constructors that the user is encouraged to explore in the package documentation.
We may now move to Interlayers. Note that, in order to define a `Multilayer(Di)Graph`, 
interlayers do not need to be explicitly constructed by the user, 
since they are automatically specified by the `Multilayer(Di)Graph` constructor. 
Anyway, more complex interlayers need to be manually instantiated. 
The interface is very similar to the layers.
=#

# Interlayer with an underlying simple directed graph and `n_edges` edges
n_vertices_1 = nv(layer_simple_directed)               # Number of vertices of layer 1
n_vertices_2 = nv(layer_simple_directed_weighted)      # Number of vertices of layer 2
n_edges = rand(1:(n_vertices_1 * n_vertices_2 - 1))    # Number of interlayer edges 
interlayer_simple_directed = interlayer_simpledigraph( # Interlayer constructor 
    layer_simple_directed,                             # Layer 1 
    layer_simple_directed_weighted,                    # Layer 2 
    n_edges                                            # Number of edges 
)

## The interlayer exports a more flexible constructor too. 
n_vertices_1 = nv(layer_simple_directed_weighted)         # Number of vertices of layer 1
n_vertices_2 = nv(layer_simple_directed_value)            # Number of vertices of layer 2
n_edges = rand(1:(n_vertices_1 * n_vertices_2 - 1))       # Number of interlayer edges 
interlayer_simple_directed_meta = interlayer_metadigraph( # Interlayer constructor
    layer_simple_directed_weighted,                       # Layer 1 
    layer_simple_directed_value,                          # Layer 2
    n_edges;                                              # Number of edges
    default_edge_metadata=(src, dst) ->                   # Edge metadata 
        (edge_metadata="metadata_of_edge_from_$(src)_to_$(dst)"),
    transfer_vertex_metadata=true                         # Boolean deciding layer vertex metadata inheritance
)

interlayers = [interlayer_simple_directed, interlayer_simple_directed_meta]

# Layers and Interlayers are not immutable, and mostly behave like normal graphs. The reader is invited to consult the API for more information.
# A MultilayerDiGraph (i.e. a directed multilayer graph, following the naming convention of the JuliaGraph ecosystem) may now be specified 

multilayerdigraph = MultilayerDiGraph(
    layers,                                    # The (ordered) collection of layers
    interlayers;                               # The manually specified interlayers
                                               # The interlayers that are left unspecified 
                                               # will be automatically inserted according 
                                               # to the keyword argument below
    default_interlayers_structure="multiplex"  # The automatically specified interlayers will have only diagonal couplings
)

# Layers and interlayer may be accessed as properties using their names
multilayerdigraph.layer_simple_directed_value

# We proceed to show some basic functionality.
## Add a vertex
### The user may add vertices that do or do not represent node_list already present in the multilayergraph. In the latter case, we have a new node:
new_node_1 = Node("new_node_1")
### Before adding any vertex representing such node to the multilayer graph, the user should first add the Node:
add_node!(multilayerdigraph, new_node_1)
### Define a vertex that represents that node
new_vertex_1 = MV(           # Constructor (alias for "MultilayerVertex")
    new_node_1,              # Node represented by the vertex
    :layer_simplevaldigraph, # Layer containing the vertex 
    ("new_metadata")         # Vertex metadata 
)
### Add the vertex
add_vertex!(
    multilayerdigraph, # MultilayerDiGraph the vertex will be added to
    new_vertex_1       # MultilayerVertex to add
)
# NB: The vertex-level metadata are considered iff the graph underlying the layer/interlayer which the vertex belongs to supports them, otherwise they are discarded.

### `add_vertex!` implements multiple interfaces.

## Add an edge
### Let's represent another node in another layer
new_node_2 = Node("new_node_2")
new_vertex_2 = MV(new_node_2, :layer_simpledigraph)
add_vertex!(
    multilayerdigraph,
    new_vertex_2;
    add_node=true # Add the associated node before adding the vertex
)
### Construct a new edge
new_edge = MultilayerEdge( # Constructor 
    new_vertex_1,          # Source vertex
    new_vertex_2,          # Destination vertex 
    ("some_edge_metadata") # Edge metadata 
)
### Add the edge
add_edge!(
    multilayerdigraph, # MultilayerDiGraph the edge will be added to
    new_edge           # MultilayerVertex to add
)
# NB: The edge-level metadata and/or weight are considered iff the graph underlying the layer/interlayer which the edge belongs to supports them, otherwise they are discarded.

#=
Using the provided `add_layer!`, `rem_layer!` and `specify_interlayer!`, 
Layers and Interlayers may be added, removed or specified on the fly. 
Since MultilayerGraphs.jl extends Graphs.jl, all metrics from the JuliaGraphs ecosystem 
should be available by default. Anyway, some multilayer-specific metrics have been implemented, 
and others required to be re-implemented. 

We showcase a few of them here:
=#

## Compute the global clustering coefficient as in @DeDomenico2014
multilayer_global_clustering_coefficient(multilayerdigraph) # A weighted version `multilayer_weighted_global_clustering_coefficient`
## Compute the overlay clustering coefficient as in @DeDomenico2013
overlay_clustering_coefficient(multilayerdigraph)
## Compute the eigenvector centrality (the implementation is sp that it coincides with Graphs.jl's `eigenvector_centrality` on monoplex graphs)
eigenvector_centrality(multilayerdigraph)
## Compute the multilayer modularity as in @DeDomenico2013
modularity(
    multilayerdigraph,
    rand([1, 2, 3, 4], length(nodes(multilayerdigraph)), length(multilayerdigraph.layers))
)

## Currently, Von Neumann entropy is available only for undirected multilayer graphs.
# NB: this brief script is far from complete: many more features and functionalities are detailed in the documentation.

```

For a more comprehensive exploration of the package features and functionalities we strongly recommend consulting the [tutorial](https://juliagraphs.org/MultilayerGraphs.jl/stable/#Tutorial) included in the package documentation.  

# Related Packages 

## R 

Here is a list of software packages for the creation, manipulation, analysis and visualisation of multilayer graphs implemented in the [R language](https://www.r-project.org): 

- [`muxViz`](https://github.com/manlius/muxViz) implements functions to perform multilayer correlation analysis, multilayer centrality analysis, multilayer community structure detection, multilayer structural reducibility, multilayer motifs analysis and utilities to statically and dynamically visualise multilayer graphs [@DeDomenico2014];
- [`multinet`](https://github.com/cran/multinet) implements functions to import, export, create and manipulate multilayer graphs, several state-of-the-art multiplex graph analysis algorithms for centrality measures, layer comparison, community detection and visualization [@Magnani2021];
- [`mully`](https://github.com/frankkramer-lab/mully) implements functions to import, export, create, manipulate and merge multilayer graphs and utilities to visualise multilayer graphs in 2D and 3D [@Hammoud2018];
- [`multinets`](https://github.com/neylsoncrepalde/multinets) implements functions to import, export, create, manipulate multilayer graphs and utilities to visualise multilayer graphs [@Lazega2008].

## Python

Here is a list of software packages for the creation, manipulation, analysis and visualisation of multilayer graphs implemented in the [Python language](https://www.python.org): 

- [`MultiNetX`](https://github.com/nkoub/multinetx) implements methods to create undirected networks with weighted or unweighted links, to analyse the spectral properties of adjacency or Laplacian matrices and to visualise multilayer graphs and dynamical processes by coloring the nodes and links accordingly;
- [`PyMNet`](https://github.com/bolozna/Multilayer-networks-library) implements data structures for multilayer graphs and multiplex graphs, methods to import, export, create, manipulate multilayer graphs and for the rule-based generation and lazy-evaluation of coupling edges and utilities to visualise multilayer graphs [@Kivela2014].

# Acknowledgements

This open-source research software project received no financial support.

# References