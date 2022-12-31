---
title: 'MultilayerGraphs.jl: A Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs'
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

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs extending Graphs.jl [@Graphs2021] and fully integrating with the [JuliaGraphs](https://github.com/JuliaGraphs) ecosystem.

A multilayer graph is a graph consisting of multiple standard subgraphs called *layers* which can be interconnected through bipartite graphs called *interlayers* composed of the vertex sets of two different layers and the edges between them. The vertices in each layer represent a single set of nodes, although not all nodes have to be represented in every layer. 

Formally, a multilayer graph can be defined as a triple $G=(V,E,L)$, where:

- $V$ is the set of vertices;
- $E$ is the set of edges, pairs of nodes $(u, v)$ representing a connection, relationship or interaction between the nodes $u$ and $v$;
- $L$ is a set of layers, which are subsets of $V$ and $E$ encoding the nodes and edges within each layer.

Each layer $\ell$ in $L$ is a tuple $(V_\ell, E_\ell)$, where $V_\ell$ is a subset of $V$ that represents the vertices within that layer, and $E_\ell$ is a subset of $E$ that represents the edges within that layer.

[A FEW WORDS ABOUT THE MAIN FEATURES, POSSIBLY EXTRACTED FROM TUTORIAL / README]

# Statement of Need

Multiple theoretical frameworks have been proposed to formally incorporate all instances of multilayer graphs [@DeDomenico2013; @Kivela2014; @Boccaletti2014; @Aleta2019; @Bianconi2018; @Cozzo2018; @Artime2022; @DeDomenico2022]. 

Multilayer graphs have been adopted to model the structure and dynamics of a wide spectrum of high-dimensional and heterogeneous complex systems, including physical, chemical, biological, neuronal, socio-technical, ecological and economic networks [@Baggio2016; @LazegaSnijders2016; @Dickison2016; @Timteo2018; @Buldu2018; @Lim2019]. 

At the best of our knowledge there are currently no software packages dedicated to the creation, manipulation and analysis of multilayer graphs implemented in the [Julia language](https://julialang.org) apart from MultilayerGraphs.jl itself [@Moroni_Monticone_MultilayerGraphs_2022].

# Main Features 

- Main structs 
- Different formalisms 
- Main methods and metrics 
- Extension of Graphs.jl [@Graphs2021], fully integrated within the [JuliaGraphs](https://github.com/JuliaGraphs) ecosystem
- Integration with Agents.jl [@Datseris2022], fully integrated within the [JuliaDynamics](https://github.com/JuliaDynamics) ecosystem

Although being part of the `Graphs.jl`'s ecosystem, due to the special nature of multilayer graphs this package features a peculiar implementation that maps a standard integer-labelled vertex representation to a more user friendly framework that exports all the objects a practitioner would expect (`Node`s, `MultilayerVertex`s, `Layer`s, `Interlayer`s, etc). The details are briefly described hereafter. 
The package revolves around two data structures, `MultilayerGraph` and `MultilayerDiGraph`. As said above, they are collection of layers whose couplings form the edge sets of the so-called interlayers. The vertices of a multilayer graph are representations of one set of distinct objects named `Node`s. Each layer may represent all or just part of such set. The vertices of `Multilayer(Di)Graph` are implemented via the `MultilayerVertex` custom type. Each `MultilayerVertex` carries information about the `node` it represents, the `layer` it belongs to and its `metadata`. Edges, both intra- and inter-layer, are embodied in the `MultilayerEdge` struct, whose fields are the two `MultilayerVertex`s involved, the edge `weight` and its `metadata`. Note that `Multilayer(Di)Graph`s are weighted and able to carry metadata by default (i.e. they are given the `IsWeighted` and `IsMeta` traits from [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl)). Layers are implemented via the `Layer` struct, which is constituted by an underlying graph from the `Graphs.jl` ecosystem and a mapping from its integer-labelled vertices to the collection of `MultilayerVertex`s the layer represents. Interlayers are similarly implemented via the `Interlayer` mutable struct, and they are generally constructed by providing the two `Layers`s involved, the (multilayer) edge list between them and an underlying graph. This usage of underlying graphs allows for easy debugging during construction and more intuitive analysis afterwards. It also allows the package to leverage all the features of the ecosystem, and acts as a proving ground of its consistency and coherence.  Now we may understand why `Multilayer(Di)Graph` are weighted and able to carry both vertex and edge-level metadata by default: since they are designed so that at any moment the user may add or remove a `Layer` or specify an `Interlayer`, and since it could be that different layers and interlayers are better substantiated by graphs that are weighted or unweighted and with or without metadata, it was necessary to provide a structure capable to adapt to the most general scenario. As specified in the [Future Developments](https://github.com/JuliaGraphs/MultilayerGraphs.jl#future-developments) section of the package README, future enhancements may provide more stringent multilayer graphs data structures, by restricting to specific traits, types and/or special cases defined in the literature.
A `Multilayer(Di)Graph` is instantiated by providing to the constructor the ordered list of layers and the list of interlayers. The latter are automatically specified, so there is no need to instantiate all of them.
Another way of constructing a `Multilayer(Di)Graph` uses a configuration model-like signature: it allows to select the degree distribution or the degree sequence (indegree and outdegree distributions or sequences may be provided separately for the `MultilayerDiGraph`) and uses the Havel-Hakimi algorithm from {Hakimi (1962)}[@havelhakimi] (or {Kleitman and Wang (1973)}[Kleitman1973] for the directed `MultilayerDiGraph`). Please note that, although inspired from BIANCONI??, this is not a complete implementation of a multilayer configuration model: it lacks the capability to specify a different distributions for different groups of layers and/or interlayers (aspects).
Once specified, the full API of `Graphs.jl` works on `Multilayer(Di)Graph`s as they were ordinary extensions of the ecosystem. Moreover, multilayer-specific methods or implementations thereof have been developed, mainly drawing from {De Domenico et al, 2013}[@DeDomenico]. They include:
- Global Clustering Coefficient
- Overlay Clustering Coefficient
- (Multilayer) Eigenvector Centrality
- (Multilayer) Modularity
- Von Neumann Entropy

`Multilayer(Di)Graph`s structure may be represented via dedicated `WeightTensor`, `MetadataTensor` and `SupraWeightMatrix` structs, all of which support indexing with `MultilayerVertex`s.

Once  a `Multilayer(Di)Graph` has been instantiated, its layers and interlayers may be accessed as they where its properties. In order to simplify the code and improve performance, `Layer`s and `Interlayers`s are not fully stored within `Multilayer(Di)Graph`s, only enough information to reconstruct them when accessed as properties is saved, in the form of `LayerDescriptor` and `InterlayerDescriptor`s.

# Installation and Usage 

To install MultilayerGraphs.jl it's sufficient to activate the `pkg` mode by pressing `]` in the Julia REPL and then run the following command:

```nothing
pkg> add MultilayerGraphs
```

[HERE WE SHOULD INSERT A FEW LINES OF CODE SHOWACASING THE MAIN FEATURES WRITTEN ABOVE]

In the package documentation you can find a comprehensive [tutorial](https://juliagraphs.org/MultilayerGraphs.jl/stable/#Tutorial) that illustrates all its main features and functionalities.

# Related Packages 

## R 

Here is a list of software packages for the creation, manipulation, analysis and visualisation of multilayer graphs implemented in the [R language](https://www.r-project.org): 

- [`muxViz`](https://github.com/manlius/muxViz) implements functions to perform multilayer correlation analysis, multilayer centrality analysis, multilayer community structure detection, multilayer structural reducibility, multilayer motifs analysis and utilities to statically and dynamically visualise multilayer graphs [@DeDomenico2014];
- [`multinet`](https://github.com/cran/multinet) implements functions to import, export, create and manipulate multilayer graphs, several state-of-the-art multiplex graph analysis algorithms for centrality measures, layer comparison, community detection and visualization [@Magnani2021];
- [`mully`](https://github.com/frankkramer-lab/mully) implements functions to import, export, create, manipulate and merge multilayer graphs and utilities to visualise multilayer graphs in 2D and 3D [@Hammoud2018];
- [`multinets`](https://github.com/neylsoncrepalde/multinets) implements functions to import/export, create, manipulate multilayer graphs and utilities to visualise multilayer graphs [@Lazega2008].

## Python

Here is a list of software packages for the creation, manipulation, analysis and visualisation of multilayer graphs implemented in the [Python language](https://www.python.org): 

- [`MultiNetX`](https://github.com/nkoub/multinetx) implements methods to create undirected networks with weighted or unweighted links, to analyse the spectral properties of adjacency or Laplacian matrices and to visualise multilayer graphs and dynamical processes by coloring the nodes and links accordingly;
- [`PyMNet`](https://github.com/bolozna/Multilayer-networks-library) implements data structures for multilayer graphs and multiplex graphs, methods to import/export, create, manipulate multilayer graphs and for the rule-based generation and lazy-evaluation of coupling edges and utilities to visualise multilayer graphs [@Kivela2014].

## Julia 

At the best of our knowledge there are currently no software packages dedicated to the creation, manipulation and analysis of multilayer graphs implemented in the [Julia language](https://julialang.org) apart from MultilayerGraphs.jl itself [@Moroni_Monticone_MultilayerGraphs_2022].

# Acknowledgements

This open-source research software project received no financial support.

# References