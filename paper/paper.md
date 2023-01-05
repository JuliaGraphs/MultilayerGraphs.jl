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

A multilayer graph is a graph consisting of multiple standard subgraphs called *layers* which can be interconnected through [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) called *interlayers* composed of the vertex sets of two different layers and the edges between them. The vertices in each layer represent a single set of nodes, although not all nodes have to be represented in every layer. 

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

In the following code chunks we synthetically illustrate how to define, handle and analyse a `MultilayerGraph` in order to showcase some of the main features outlined in the previous section.

```julia
using Revise
using Graphs, SimpleValueGraphs
using Distributions
using MultilayerGraphs

# The objects that `MultilayerVertex`s represent are `Node`s
const n_nodes = 100
const nodes = [Node("node_$i") for i in 1:n_nodes]

#= # Every layer and Interlayer that is involved in a `Multilayer(Di)Graph` must have the same vertex type (i.e. the type of the internal representation of vertices) and (edge) weight type. We proceed to define them hereafter:
const vertex_type = Int64
const weight_type = Float64 =#

# We will be instantiating layers and interlayers with randomly-chosen edges and vertices, using different techniques with the intent to showcase them.
## Here, we define a layer with an underlying simple directed graph, using a configuration model-like constructor that allows for specifying both the indegree and the outdegree sequences. Before instantiating each layer we will sample its number of vertices and, optionally, of edges.
_nv = rand(1:100)
_layer_simpledigraph = layer_simpledigraph(
                                                    :layer_simpledigraph,           # The name of the layer_metagraph
                                                    sample(nodes, _nv, replace=false),                          # The nodes that the layer will represent
                                                    Truncated(Normal(5,5), 0, 20),  # The distribution from which the indegree sequence will be sampled from
                                                    Truncated(Normal(5,5), 0, 20)   # The distribution from which the outdegree sequence will be sampled from
)

## Next, we define a layer with an underliyng simple weighted directed graph. This time we show another kind of constructor that allows the user to specify the number of edges to be randomly distributed among the vertices. The keyword argument `default_edge_weight` will assign a weight to such edges before they are added to the layer
_nv = rand(1:n_nodes)
_ne = rand(_nv:(_nv* (_nv-1) - 1))
_layer_simpleweighteddigraph = layer_simpleweighteddigraph(
                                                                    :layer_simpleweighteddigraph,           
                                                                    sample(nodes, _nv, replace=false),                     
                                                                    _ne;                                      # The number of edges that will be randomly distributed among the vertices
                                                                    default_edge_weight = (src,dst) -> rand() # The function that assigns a weight to such edges before they are added to the layer
)


## Similar constructors, more flexible at the cost of ease of use, allows for finer tuning:
### NB: This constructor should be necessary only in rare circumstances, where e.g. the equivalent simplified constructor `layer_simplevaldigraph` is not able to infer the correct return types of the `default_vertex/edge_metadata`s, or to use underlying graph that aren't currently supported.
_nv = rand(1:n_nodes)
_ne = rand(_nv:(_nv* (_nv-1) - 1))
default_vertex_metadata = v -> ("vertex_$(v)_metadata",)
default_edge_metadata  =  (s,d) -> (rand(),)
_layer_valdigraph = Layer(
                                        :layer_simplevaldigraph, 
                                        sample(nodes, _nv, replace=false),
                                        _ne,
                                        ValDiGraph(
                                                    SimpleDiGraph{Int64}();
                                                    vertexval_types = (String,),
                                                    vertexval_init = default_vertex_metadata,
                                                    edgeval_types = (Float64,),
                                                    edgeval_init =   default_edge_metadata 
                                        ),
                                        Float64;
                                        default_vertex_metadata = default_vertex_metadata,
                                        default_edge_metadata = default_edge_metadata
)

layers = [_layer_simpledigraph, _layer_simpleweighteddigraph, _layer_valdigraph]

# There are many more constructors that the user is encouraged to explore in the package documentation.

# We may now move to Interlayers. Note that, in order to define a `Multilayer(Di)Graph`, interlayers do not need to be explicitly constructed by the user, since they are automatically specified by the `Multilayer(Di)Graph` constructor. Anyway, more complex interlayers need to be manually instantiated. The interface is very similar to the layers.

## Interlayer with an underlying simple directed graph and `_ne` edges
nv_1 = nv(_layer_simpledigraph) 
nv_2 = nv(_layer_simpleweighteddigraph)
_ne = rand(1:(nv_1*nv_2 -1))             # The interlayer is a bipartite graph between the vertices of the two layers
_interlayer_simpledigraph = interlayer_simpledigraph(
                                                        _layer_simpledigraph,          # One of the two layers connected by this interlayer
                                                        _layer_simpleweighteddigraph,  # One of the two layers connected by this interlayer
                                                        _ne                            # # The number of edges that will be randomly distributed among the vertices

)

## The interlayer exports a more flexible constructor too. 
nv_1 = nv(_layer_simpledigraph) 
nv_2 = nv(_layer_simpleweighteddigraph)
_ne = rand(1:(nv_1*nv_2 -1))
_interlayer_metadigraph = interlayer_metadigraph(
                                                    _layer_simpleweighteddigraph,
                                                    _layer_valdigraph,
                                                    _ne;
                                                    default_edge_metadata = (src,dst) -> (edge_metadata = "metadata_of_edge_from_$(src)_to_$(dst)",),
                                                    transfer_vertex_metadata = true # Whether to have the vertices of the interlayer endowed with the metadata of the vertices of the layers it connects. By default it is set to false.

)

interlayers = [_interlayer_simpledigraph, _interlayer_metadigraph]

# Layers and Interlayers are not immutable, and mostly behave like normal graphs. The reader is invited to consult the API for more information.
# A MultilayerDiGraph (i.e. a directed multilayer graph, following the naming convention of the JuliaGraph ecosystem) may now be specified 

multilayerdigraph = MultilayerDiGraph(
                                        layers,     # The (ordered) collection of layers that constitute the multilayer graph
                                        interlayers, # The manually-specified interlayers. The interlayers that are left unspecified (in this example, the interlayer between `_layer_simpledigraph` and `_layer_valdigraph`), will be automatically inserted according to the keyword argument below
                                        default_interlayers_structure = "multiplex" # The automatically-specified interlayers will have only diagonal couplings

)
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