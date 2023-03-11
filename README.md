# MultilayerGraphs.jl 

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/LICENSE)
[![Docs: Stable](https://img.shields.io/badge/Docs-Stable-blue.svg)](https://juliagraphs.org/MultilayerGraphs.jl/stable)
[![Docs: Dev](https://img.shields.io/badge/Docs-Dev-lightblue.svg)](https://juliagraphs.org/MultilayerGraphs.jl/dev)
[![CI](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CI.yml)
[![Compat Helper](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CompatHelper.yml)
[![Format Check](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/FormatCheck.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/FormatCheck.yml)
[![Coverage: Codecov](https://codecov.io/gh/JuliaGraphs/MultilayerGraphs.jl/branch/main/graph/badge.svg?token=Z758JuxDJX)](https://codecov.io/gh/JuliaGraphs/MultilayerGraphs.jl)
[![Coverage: Coveralls](https://coveralls.io/repos/github/JuliaGraphs/MultilayerGraphs.jl/badge.svg?branch=main)](https://coveralls.io/github/JuliaGraphs/MultilayerGraphs.jl?branch=main)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/MultilayerGraphs&label=Downloads)](https://pkgs.genieframework.com?packages=MultilayerGraphs)
[![DOI: Zenodo](https://zenodo.org/badge/490352002.svg)](https://zenodo.org/badge/latestdoi/490352002)
[![JOSS](https://joss.theoj.org/papers/7d1b94ff7a1a2ebd5e86075e86fc62fb/status.svg)](https://joss.theoj.org/papers/7d1b94ff7a1a2ebd5e86075e86fc62fb)

<img align="right" width="220" height="220" src="https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/docs/src/assets/logo.png?raw=true">

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs. 

## 🌐 Overview

A multilayer graph is a graph consisting of multiple standard subgraphs called *layers* which can be interconnected through [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) called *interlayers* composed of the vertex sets of two different layers and the edges between them. The vertices in each layer represent a single set of nodes, although not all nodes have to be represented in every layer. 

Formally, a multilayer graph can be defined as a triple $G=(V,E,L)$, where:

- $V$ is the set of vertices;
- $E$ is the set of edges, pairs of nodes $(u, v)$ representing a connection, relationship or interaction between the nodes $u$ and $v$;
- $L$ is a set of layers, which are subsets of $V$ and $E$ encoding the nodes and edges within each layer.

Each layer $\ell$ in $L$ is a tuple $(V_\ell, E_\ell)$, where $V_\ell$ is a subset of $V$ that represents the vertices within that layer, and $E_\ell$ is a subset of $E$ that represents the edges within that layer.

Multiple theoretical frameworks have been proposed to formally subsume all instances of multilayer graphs ([De Domenico  et al. (2013)](https://doi.org/10.1103/physrevx.3.041022); [Kivelä et al. (2014)](https://doi.org/10.1093/comnet/cnu016); [Boccaletti et al. (2014)](https://doi.org/10.1016/j.physrep.2014.07.001); [Lee et al. (2015)](https://doi.org/10.1140/epjb/e2015-50742-1); [Aleta and Moreno (2019)](https://doi.org/10.1146/annurev-conmatphys-031218-013259); [Bianconi (2018)](https://doi.org/10.1093/oso/9780198753919.001.0001); [Cozzo et al. (2018)](https://doi.org/10.1007/978-3-319-92255-3); [Artime et al. (2022)](https://doi.org/10.1017/9781009085809); [De Domenico (2022)](https://doi.org/10.1007/978-3-030-75718-2)). 

Multilayer graphs have been adopted to model the structure and dynamics of a wide spectrum of high-dimensional, non-linear, multi-scale, time-dependent complex systems including physical, chemical, biological, neuronal, socio-technical, epidemiological, ecological and economic networks ([Cozzo et al. (2013)](https://doi.org/10.1103/physreve.88.050801); [Granell et al. (2013)](https://doi.org/10.1103/physrevlett.111.128701); [Massaro and Bagnoli (2014)](https://doi.org/10.1103/physreve.90.052817); [Estrada and Gomez-Gardenes (2014)](https://doi.org/10.1103/physreve.89.042819); [Azimi-Tafreshi (2016)](https://doi.org/10.1103/physreve.93.042303); [Baggio et al. (2016)](https://doi.org/10.1073/pnas.1604401113); [DeDomenico et al. (2016)](https://doi.org/10.1038/nphys3865); [Amato et al. (2017)](https://doi.org/10.1038/s41598-017-06933-2); [DeDomenico (2017)](https://doi.org/10.1093/gigascience/gix004); [Pilosof et al. (2017)](https://doi.org/10.1038/s41559-017-0101); [de Arruda et al. (2017)](https://doi.org/10.1103/physrevx.7.011014); [Gosak et al. (2018)](https://doi.org/10.1016/j.plrev.2017.11.003); [Soriano-Panos et al. (2018)](https://doi.org/10.1103/physrevx.8.031039); [Timteo et al. (2018)](https://doi.org/10.1038/s41467-017-02658-y); [Buldú et al. (2018)](https://doi.org/10.1162/netn_a_00033); [Lim et al. (2019)](https://doi.org/10.1038/s41598-019-39243-w); [Mangioni et al. (2020)](https://doi.org/10.1109/tnse.2018.2871726); [Aleta et al. (2020)](https://doi.org/10.1038/s41562-020-0931-9); [Aleta et al. (2022)](https://doi.org/10.1073/pnas.2112182119)). 

MultilayerGraphs.jl is an integral part of the [JuliaGraphs](https://github.com/JuliaGraphs) ecosystem extending [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) so all the methods and metrics exported by Graphs.jl work for multilayer graphs, but due to the special nature of multilayer graphs the package features a peculiar implementation that maps a standard integer-labelled vertex representation to a more user-friendly framework exporting all the objects an experienced practitioner would expect such as nodes (`Node`), vertices (`MultilayerVertex`), layers (`Layer`), interlayers (`Interlayer`), etc.

MultilayerGraphs.jl features multilayer-specific methods and metrics including the global clustering coefficient, the overlay clustering coefficient, the multilayer eigenvector centrality, the multilayer modularity and the Von Neumann entropy.

Finally, MultilayerGraphs.jl has been integrated within the [JuliaDynamics](https://github.com/JuliaDynamics) ecosystem so that any `Multilayer(Di)Graph` can be utilised as an argument to the `GraphSpace` constructor in [Agents.jl](https://github.com/JuliaDynamics/Agents.jl). 

## 🔰 Installation

To install the latest stable release of MultilayerGraphs.jl, make sure you have [installed](https://julialang.org/downloads/) Julia v1.8 or later and run the following command:

``` julia
using Pkg
Pkg.add("MultilayerGraphs")
```

The development version can be installed as follows:

``` julia
using Pkg
Pkg.add(url="https://github.com/JuliaGraphs/MultilayerGraphs.jl")
```

## 🔍 Usage

Let's begin by importing the necessary dependencies and setting the relevant constants.

```julia
# Import necessary dependencies
using Distributions, Graphs, SimpleValueGraphs
using MultilayerGraphs
# Set the number of nodes
const n_nodes = 100 
# Create a list of nodes
const node_list = [Node("node_$i") for i in 1:n_nodes]
```

### Layers and Interlayers 

We will instantiate layers and interlayers with randomly-selected edges and vertices adopting a variety of techniques. Layers and Interlayers are not immutable, and mostly behave like normal graphs. The user is invited to consult the [API](https://juliagraphs.org/MultilayerGraphs.jl/stable/API/) for further details.

Here we define a layer with an underlying simple directed graph using a graph generator-like (or "configuration model"-like) constructor which allows us to specify both the **indegree** and the **outdegree sequences**. Before instantiating each layer we sample the number of its vertices and, optionally, of its edges.

```julia
# Create a simple directed layer
n_vertices = rand(1:100)                          # Number of vertices 
layer_simple_directed = layer_simpledigraph(      # Layer constructor 
    :layer_simple_directed,                       # Layer name
    sample(node_list, n_vertices; replace=false), # Nodes represented in the layer
    Truncated(Normal(5, 5), 0, 20), # Indegree sequence distribution 
    Truncated(Normal(5, 5), 0, 20)  # Outdegree sequence distribution
)
```

Then we define a layer with an underlying simple weighted directed graph. This is another kind of constructor that allows the user to specify the number of edges to be randomly distributed among vertices. 

```julia
# Create a simple directed weighted layer
n_vertices = rand(1:n_nodes)                                   # Number of vertices 
n_edges = rand(n_vertices:(n_vertices * (n_vertices - 1) - 1)) # Number of edges 
layer_simple_directed_weighted = layer_simpleweighteddigraph(  # Layer constructor 
    :layer_simple_directed_weighted,                           # Layer name
    sample(node_list, n_vertices; replace=false), # Nodes represented in the layer
    n_edges;                                 # Number of randomly distributed edges
    default_edge_weight=(src, dst) -> rand() # Function assigning weights to edges 
)
```

Similar constructors, more flexible at the cost of ease of use, enable a finer tuning. The constructor we use below should be necessary only in rare circumstances, e.g. if the equivalent simplified constructor `layer_simplevaldigraph` is not able to infer the correct return types of `default_vertex_metadata` or `default_edge_metadata`, or to use and underlying graph structure that isn't currently supported.

```julia
# Create a simple directed value layer
n_vertices = rand(1:n_nodes)                                   # Number of vertices 
n_edges = rand(n_vertices:(n_vertices * (n_vertices - 1) - 1)) # Number of edges 
default_vertex_metadata = v -> ("vertex_$(v)_metadata",)       # Vertex metadata 
default_edge_metadata = (s, d) -> (rand(),)                    # Edge metadata 
layer_simple_directed_value = Layer(                           # Layer constructor
    :layer_simple_directed_value,                              # Layer name
    sample(node_list, n_vertices; replace=false), # Nodes represented in the layer
    n_edges,                                      # Number of randomly distributed edges
    ValDiGraph(                                                
        SimpleDiGraph{Int64}(); 
        vertexval_types=(String,),
        vertexval_init=default_vertex_metadata,
        edgeval_types=(Float64,),
        edgeval_init=default_edge_metadata,
    ),
    Float64;
    default_vertex_metadata=default_vertex_metadata, # Vertex metadata 
    default_edge_metadata=default_edge_metadata      # Edge metadata 
)

# Create a list of layers 
layers = [layer_simple_directed, layer_simple_directed_weighted, layer_simple_directed_value]
```

There are many more constructors the user is encouraged to explore in the package [documentation](https://juliagraphs.org/MultilayerGraphs.jl).

The interface of interlayers is very similar to that of layers. It is very important to notice that, in order to define a `Multilayer(Di)Graph`, interlayers don't need to be explicitly constructed by the user since they are automatically identified by the `Multilayer(Di)Graph` constructor, but for more complex interlayers the manual instantiation is required.

Here we define an interlayer with an underlying simple directed graph.

```julia
# Create a simple directed interlayer
n_vertices_1 = nv(layer_simple_directed)               # Number of vertices of layer 1
n_vertices_2 = nv(layer_simple_directed_weighted)      # Number of vertices of layer 2
n_edges = rand(1:(n_vertices_1 * n_vertices_2 - 1))    # Number of interlayer edges 
interlayer_simple_directed = interlayer_simpledigraph( # Interlayer constructor 
    layer_simple_directed,                             # Layer 1 
    layer_simple_directed_weighted,                    # Layer 2 
    n_edges                                            # Number of edges 
)
```

The interlayer exports a more flexible constructor too.

```julia
# Create a simple directed meta interlayer 
n_vertices_1 = nv(layer_simple_directed_weighted)   # Number of vertices of layer 1
n_vertices_2 = nv(layer_simple_directed_value)      # Number of vertices of layer 2
n_edges = rand(1:(n_vertices_1 * n_vertices_2 - 1)) # Number of interlayer edges 
interlayer_simple_directed_meta = interlayer_metadigraph( # Interlayer constructor
    layer_simple_directed_weighted,                       # Layer 1 
    layer_simple_directed_value,                          # Layer 2
    n_edges;                                              # Number of edges
    default_edge_metadata=(src, dst) ->                   # Edge metadata 
        (edge_metadata="metadata_of_edge_from_$(src)_to_$(dst)",),
    transfer_vertex_metadata=true # Boolean deciding layer vertex metadata inheritance
)

# Create a list of interlayers 
interlayers = [interlayer_simple_directed, interlayer_simple_directed_meta]
```

### Multilayer Graphs

Let's construct a directed multilayer graph (`MultilayerDiGraph`).

```julia
# Create a simple directed multilayer graph
multilayerdigraph = MultilayerDiGraph( # Constructor 
    layers,                     # The (ordered) collection of layers
    interlayers;                # The manually specified interlayers
                                # The interlayers that are left unspecified 
                                # will be automatically inserted according 
                                # to the keyword argument below
    default_interlayers_structure="multiplex" 
    # The automatically specified interlayers will have only diagonal couplings
)

# Layers and interlayer can be accessed as properties using their names
multilayerdigraph.layer_simple_directed_value
```

Then we proceed by showing how to add nodes, vertices and edges to a directed multilayer graph. The user may add vertices that do or do not represent nodes which are already present in the multilayer graph. In the latter case, we have to create a node first and then add the vertex representing such node to the multilayer graph. The vertex-level metadata are effectively considered only if the graph underlying the relevant layer or interlayer supports them, otherwise they are discarded. The same holds for edge-level metadata and/or weight. 

```julia
# Create a node 
new_node_1 = Node("new_node_1")
# Add the node to the multilayer graph 
add_node!(multilayerdigraph, new_node_1)
# Create a vertex representing the node 
new_vertex_1 = MV(                # Constructor (alias for "MultilayerVertex")
    new_node_1,                   # Node represented by the vertex
    :layer_simple_directed_value, # Layer containing the vertex 
    ("new_metadata",)             # Vertex metadata 
)
# Add the vertex 
add_vertex!(
    multilayerdigraph, # MultilayerDiGraph the vertex will be added to
    new_vertex_1       # MultilayerVertex to add
)

# Create another node in another layer 
new_node_2 = Node("new_node_2")
# Create another vertex representing the new node
new_vertex_2 = MV(new_node_2, :layer_simple_directed_value)
# Add the new vertex
add_vertex!(
    multilayerdigraph,
    new_vertex_2;
    add_node=true # Add the associated node before adding the vertex
)
# Create an edge 
new_edge = MultilayerEdge(  # Constructor 
    new_vertex_1,           # Source vertex
    new_vertex_2,           # Destination vertex 
    ("some_edge_metadata",) # Edge metadata 
)
# Add the edge 
add_edge!(
    multilayerdigraph, # MultilayerDiGraph the edge will be added to
    new_edge           # MultilayerVertex to add
)
```

Finally we illustrate how to compute a few multilayer metrics such as the global clustering coefficient, the overlay clustering coefficient, the multilayer eigenvector centrality, and the multilayer modularity as defined in [De Domenico  et al. (2013)](https://doi.org/10.1103/physrevx.3.041022). 

```julia
# Compute the global clustering coefficient
multilayer_global_clustering_coefficient(multilayerdigraph) 
# Compute the overlay clustering coefficient
overlay_clustering_coefficient(multilayerdigraph)
# Compute the multilayer eigenvector centrality 
eigenvector_centrality(multilayerdigraph)
# Compute the multilayer modularity 
modularity(
    multilayerdigraph,
    rand([1, 2, 3, 4], length(nodes(multilayerdigraph)), length(multilayerdigraph.layers))
)
```

## 🎯 Future Developments 

All the information regarding the future developments of MultilayerGraphs.jl can be found in the [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues).

## 🛠 How to Contribute 

The ongoing development of this package would greatly benefit from the valuable feedback of the esteemed members of the [JuliaGraph](https://github.com/orgs/JuliaGraphs/people) community, as well as from graph theorists, network scientists, and any users who may have general questions or suggestions. 

We therefore encourage you to participate in [discussions](https://github.com/JuliaGraphs/MultilayerGraphs.jl/discussions), raise [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues), or submit [pull requests](https://github.com/JuliaGraphs/MultilayerGraphs.jl/pulls). Your contributions are welcome!

## 🎓 How to Cite

If you utilize this package in your project, please consider citing this repository using the citation information provided in [`CITATION.bib`](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/CITATION.bib). 

This will help to give appropriate credit to the [contributors](https://github.com/JuliaGraphs/MultilayerGraphs.jl/graphs/contributors) and support the continued development of the package.

## 📢 Announcements 

### v0.1

MultilayerGraphs.jl (v0.1) and its features were announced on the following platforms:

- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-a-package-to-construct-handle-and-analyse-multilayer-graphs/85988);
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-a-package-to-construct-handle-and-analyse-multilayer-graphs-3k22);
- [Twitter](https://twitter.com/In_Phy_T/status/1560594513189638146).

### v1.1

MultilayerGraphs.jl (v1.1) and its features were announced on the following platforms:

- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-v1-1-multilayer-network-science-in-julia/92680);
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-v11-multilayer-network-science-in-julia-2oa3);
- [Twitter](https://twitter.com/In_Phy_T/status/1612460371939581955).

## 📦 Related Packages 

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

## 📚 References

1. De Domenico et al. (2013) [Mathematical Formulation of Multilayer Networks](https://doi.org/10.1103/PhysRevX.3.041022). *Physical Review X*; 
2. Kivelä et al. (2014) [Multilayer networks](https://doi.org/10.1093/comnet/cnu016). *Journal of Complex Networks*; 
3. Boccaletti et al. (2014) [The structure and dynamics of multilayer networks](https://doi.org/10.1016/j.physrep.2014.07.001). *Physics Reports*; 
4. Lee et al. (2015) [Towards real-world complexity: an introduction to multiplex networks](https://doi.org/10.1140/epjb/e2015-50742-1). *The European Physical Journal B*; 
5. Bianconi (2018) [Multilayer Networks: Structure and Function](https://global.oup.com/academic/product/multilayer-networks-9780192865540). *Oxford University Press*;
6. Cozzo et al. (2018) [Multiplex Networks: Basic Formalism and Structural Properties](https://doi.org/10.1007/978-3-319-92255-3). *SpringerBriefs in Complexity*; 
7. Aleta and Moreno (2019) [Multilayer Networks in a Nutshell](https://doi.org/10.1146/annurev-conmatphys-031218-013259). *Annual Review of Condensed Matter Physics*; 
8. Artime et al. (2022) [Multilayer Network Science: From Cells to Societies](https://doi.org/10.1017/9781009085809). *Cambridge University Press*; 
9. De Domenico (2022) [Multilayer Networks: Analysis and Visualization](https://doi.org/10.1007/978-3-030-75718-2). *Springer Cham*. 