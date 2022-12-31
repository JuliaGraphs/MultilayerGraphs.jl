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

A multilayer graph may also be formally defined as a collection of individual graphs $G = (G_1, G_2, ..., G_n)$, where each graph represents a layer and the edges between the layers represent inter-layer connections. 

Multiple theoretical frameworks have been proposed to formally incorporate all instances of multilayer graphs [@DeDomenico2013; @Kivela2014; @Boccaletti2014; @Aleta2019; @Bianconi2018; @Cozzo2018; @Artime2022; @DeDomenico2022]. 

Multilayer graphs have been adopted to model the structure and dynamics of a wide range of high-dimensional and heterogeneous complex systems, including physical, chemical, biological, socio-technical and economic networks [@LazegaSnijders2016; @Dickison2016]. 

# Statement of Need

- Highlight the importance of multilayer graphs in the modern computational modelling of high-dimensional, non-linear and highly heterogeneous phenomena both in the natural and in the social sciences.
- At the best of our knowledge there are currently no software packages dedicated to the creation, manipulation and analysis of multilayer graphs implemented in the [Julia language](https://julialang.org) apart from MultilayerGraphs.jl itself [@Moroni_Monticone_MultilayerGraphs_2022].

# Main Features 

- Main structs 
- Different formalisms 
- Main methods and metrics 
- Extension of Graphs.jl [@Graphs2021], fully integrated within the [JuliaGraphs](https://github.com/JuliaGraphs) ecosystem
- Integration with Agents.jl [@Datseris2022], fully integrated within the [JuliaDynamics](https://github.com/JuliaDynamics) ecosystem

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