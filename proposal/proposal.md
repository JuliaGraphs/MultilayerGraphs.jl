# JuliaCon 2023 Proposal 

**NB**:
- Session type is pretty flexible: it can be changed during the review process either by authors' or reviewers' choice;
- Make sure to read proposals and watch talks from previous JuliaCons to see what's expected; 
- It's possible to submit the same talk in multiple languages: you just need to write it in the "Note" section.
- It's not mandatory to present just a solution to a problem: pointing to a problem to solve would be highly appreciated by the community. It might be very useful to highlight the limitations of the current development of your package (e.g. features you aren't able implement, etc.). 

## Title 

MultilayerGraphs.jl: Multilayer Network Science in Julia

## Session Type 

Lightning talk (10 minutes). 

## Track 

JuliaCon 

## Abstract 

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of multilayer graphs, which have been adopted to model of a wide range of complex systems from bio-chemical to socio-technical networks.

We will synthetically introduce multilayer graph theory and applications, illustrate some of the main features of the current version of the package and talk about its future developments.

## Description 

[**MultilayerGraphs.jl**](https://github.com/JuliaGraphs/MultilayerGraphs.jl) is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs. 

A multilayer graph consists of multiple subgraphs called *layers* which can be interconnected through [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) called *interlayers*.

In order to formally represent multilayer networks, several theoretical paradigms have been proposed (e.g. see [Bianconi (2018)](https://doi.org/10.1093/oso/9780198753919.001.0001) and [De Domenico (2022)](https://doi.org/10.1007/978-3-030-75718-2)) and adopted to model the structure and dynamics of a wide spectrum of high-dimensional, multi-scale, time-dependent complex systems including molecular, neuronal, social, ecological and economic networks (e.g. see [Amato et al. (2017)](https://doi.org/10.1038/s41598-017-06933-2), [DeDomenico (2017)](https://doi.org/10.1093/gigascience/gix004), [Timteo et al. (2018)](https://doi.org/10.1038/s41467-017-02658-y), [Aleta et al. (2020)](https://doi.org/10.1038/s41562-020-0931-9), [Aleta et al. (2022)](https://doi.org/10.1073/pnas.2112182119)).

The package features an implementation that maps a standard integer-labelled vertex representation to a more user-friendly framework exporting all the objects a practitioner would expect such as nodes, vertices, layers, interlayers, etc.

MultilayerGraphs.jl has been integrated with the [JuliaGraphs](https://github.com/JuliaGraphs) and the [JuliaDynamics](https://github.com/JuliaDynamics) ecosystems through: 

- the extension of [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) with several methods and metrics including the multilayer eigenvector centrality, the multilayer modularity and the Von Newman entropy; 
- the compatibility with [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) allowing for agent-based modelling on general multilayer networks. 

We will synthetically introduce multilayer graph theory and applications, illustrate some of the main features of the current version of the package and talk about its future developments.

In our talk we will briefly introduce the theory and applications of multilayer graphs and showcase some of the main features of the current version of the package through a quick tutorial including: 

- how to install the package;
- how to define layers and interlayers with a variety of constructors and underlying graphs;
- how to construct a directed multilayer graph with those layers and interlayers;
- how to add nodes, vertices and edges to the multilayer graph;
- how to compute some standard multilayer metrics.

For a more comprehensive exploration of the package functionalities and further details on the future developments the user is invited to consult the package [README](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/README.md), [documentation](https://juliagraphs.org/MultilayerGraphs) and [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues).

## Notes 

### Multiple Languages 

We're available to submit the same proposal / talk in two languages: 
- **English**: "MultilayerGraphs.jl: Multilayer Network Science in Julia"; 
- **Italian**: "MultilayerGraphs.jl: Scienza delle Reti Multistrato in Julia". 

### Statement of Need 

At the best of our knowledge there are currently no software packages dedicated to the creation, manipulation and analysis of multilayer graphs implemented in the Julia language apart from MultilayerGraphs.jl itself.

### Future Developments 

See the relevant [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues). 

### Publication

We have recently submitted a [paper](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/JOSS/paper/paper.pdf) to [JOSS](https://joss.theoj.org) presenting MultilayerGraphs.jl. [Here](https://github.com/openjournals/joss-reviews/issues/5055) you can follow the open review process. 

MultilayerGraphs.jl has been published on the following open access websites: 
- [GitHub](https://github.com/JuliaGraphs/MultilayerGraphs.jl);
- [Zenodo](https://doi.org/10.5281/zenodo.7009172);
- [Research Software Directory](https://research-software-directory.org/software/multilayergraphs).

### Announcements 

#### v0.1

MultilayerGraphs.jl (v0.1) and its features were announced on the following platforms:

- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-a-package-to-construct-handle-and-analyse-multilayer-graphs/85988);
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-a-package-to-construct-handle-and-analyse-multilayer-graphs-3k22);
- [Twitter](https://twitter.com/In_Phy_T/status/1560594513189638146).

#### v1.1

MultilayerGraphs.jl (v1.1) and its features were announced on the following platforms:

- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-v1-1-multilayer-network-science-in-julia/92680);
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-v11-multilayer-network-science-in-julia-2oa3);
- [Twitter](https://twitter.com/In_Phy_T/status/1612460371939581955).

## Session Image 

![logo](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/docs/src/assets/logo.png?raw=true)

