# MultilayerGraphs.jl 

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/LICENSE)
[![Stable](https://img.shields.io/badge/Docs-Stable-blue.svg)](https://juliagraphs.org/MultilayerGraphs.jl/stable)
[![Dev](https://img.shields.io/badge/Docs-Dev-lightblue.svg)](https://juliagraphs.org/MultilayerGraphs.jl/dev)
[![CI](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CI.yml)
[![Compat Helper](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CompatHelper.yml)
[![Format Check](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/FormatCheck.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/FormatCheck.yml)
[![Coverage: Codecov](https://codecov.io/gh/JuliaGraphs/MultilayerGraphs.jl/branch/main/graph/badge.svg?token=Z758JuxDJX)](https://codecov.io/gh/JuliaGraphs/MultilayerGraphs.jl)
[![Coverage: Coveralls](https://coveralls.io/repos/github/JuliaGraphs/MultilayerGraphs.jl/badge.svg?branch=main)](https://coveralls.io/github/JuliaGraphs/MultilayerGraphs.jl?branch=main)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![DOI: Zenodo](https://zenodo.org/badge/490352002.svg)](https://zenodo.org/badge/latestdoi/490352002)

<img align="right" width="220" height="220" src="https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/docs/src/assets/logo.png?raw=true">

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs [extending Graphs.jl](https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/).

## Overview

**MultilayerGraphs.jl** provides an implementation of the mathematical formulation of multilayer graphs as proposed by [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022) and incorporates insights from [Kivelä et al. (2014)](https://doi.org/10.1093/comnet/cnu016) and [Bianconi (2018)](https://global.oup.com/academic/product/multilayer-networks-9780192865540). The package focuses on two custom types, [`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref), which represent undirected and directed multilayer graphs, respectively.

A multilayer graph is composed of ***layers***, i.e. graphs whose vertices represent the same set of nodes (not all nodes need to be represented in every layer), and ***interlayers***, i.e. the [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) that connect vertices in two different layers. Vertices in a multilayer graph are represented using the [`MultilayerVertex`](@ref) struct, while nodes are represented using the [`Node`](@ref) struct.

[`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) are fully-fledged [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) extensions. Both structs are designed to allow for layers and interlayers of any type (as long as they are Graphs.jl extensions themselves) and to permit layers and interlayers of different types. However, it is required that all layers and interlayers in [`MultilayerGraph`](@ref) are undirected, and all layers and interlayers in [`MultilayerDiGraph`](@ref) are directed.

[`MultilayerGraph`](@ref) and [`MultilayerDiGraph`](@ref) support the specification of vertex and edge metadata, provided that the underlying layer or interlayer also supports metadata.

## Installation

Press `]` in the Julia REPL and then

```nothing
pkg> add MultilayerGraphs
```

## Usage

In the package documentation you can find a [tutorial](https://juliagraphs.org/MultilayerGraphs.jl/stable/#Tutorial) that illustrates all its main features and functionalities.

## Future Developments 

- [ ] [Implement more general configuration models / graph generators](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/33);
- [ ] [Implement graph of layers](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/34);
- [ ] [Implement projected monoplex and overlay graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/35);
- [ ] [Implement more default multilayer graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/36) (e.g. multiplex graphs).

## How to Contribute 

The ongoing development of this package would greatly benefit from the valuable feedback of the esteemed members of the [JuliaGraph](https://github.com/orgs/JuliaGraphs/people) community, as well as from graph theorists, network scientists, and any users who may have general questions or suggestions. 

We therefore encourage you to participate in [discussions](https://github.com/JuliaGraphs/MultilayerGraphs.jl/discussions), raise [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues), or submit [pull requests](https://github.com/JuliaGraphs/MultilayerGraphs.jl/pulls). Your contributions are most welcome!

## How to Cite

If you utilize this package in your project, please consider citing this repository using the citation information provided in [`CITATION.bib`](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/CITATION.bib). This will help to give appropriate credit to the [contributors](https://github.com/JuliaGraphs/MultilayerGraphs.jl/graphs/contributors) and support the continued development of the package.

## Announcements 

The package and its features were announced on the following platforms:

- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-a-package-to-construct-handle-and-analyse-multilayer-graphs/85988)
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-a-package-to-construct-handle-and-analyse-multilayer-graphs-3k22)
- [Twitter](https://twitter.com/In_Phy_T/status/1560594513189638146)

## References

1. De Domenico et al. (2013) [Mathematical Formulation of Multilayer Networks](https://doi.org/10.1103/PhysRevX.3.041022). *Physical Review X*; 
2. Kivelä et al. (2014) [Multilayer networks](https://doi.org/10.1093/comnet/cnu016). *Journal of Complex Networks*; 
3. Bianconi (2018) [Multilayer Networks: Structure and Function](https://global.oup.com/academic/product/multilayer-networks-9780192865540). *Oxford University Press*.
 
