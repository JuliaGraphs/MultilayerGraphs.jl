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



## To-Write
- One or two sentences on the mathematical formulation of graphs (with LaTeX) and scientific applications citing the relevant scientific literature
- One or two sentences on the mathematical formulation of multilayer graphs (with LaTeX) and scientific applications citing the relevant scientific literature.
- Highlight the importance of multilayer graphs in the modern computational modelling of high-dimensional, non-linear and highly heterogeneous phenomena both in the natural and in the social sciences.

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