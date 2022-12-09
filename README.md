# MultilayerGraphs.jl 

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/LICENSE)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagraphs.org/MultilayerGraphs.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://juliagraphs.org/MultilayerGraphs.jl/dev)
[![CI](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CI.yml)
[![Compat Helper](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/CompatHelper.yml)
[![Format Check](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/FormatCheck.yml/badge.svg)](https://github.com/JuliaGraphs/MultilayerGraphs.jl/actions/workflows/FormatCheck.yml)
[![Coverage: Codecov](https://codecov.io/gh/JuliaGraphs/MultilayerGraphs.jl/branch/main/graph/badge.svg?token=Z758JuxDJX)](https://codecov.io/gh/JuliaGraphs/MultilayerGraphs.jl)
[![Coverage: Coveralls](https://coveralls.io/repos/github/JuliaGraphs/MultilayerGraphs.jl/badge.svg?branch=main)](https://coveralls.io/github/JuliaGraphs/MultilayerGraphs.jl?branch=main)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![DOI: Zenodo](https://zenodo.org/badge/490352002.svg)](https://zenodo.org/badge/latestdoi/490352002)

<img align="right" width="215" height="215" src="https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/docs/src/assets/logo.png?raw=true">

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

## Usage

In the package documentation you can find a [tutorial](https://juliagraphs.org/MultilayerGraphs.jl/stable/#Tutorial) that illustrates all its main features.

## Future Developments 

- [ ] [Implement faster graph realization algorithms](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/32);
- [ ] [Implement more general configuration models / graph generators](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/33);
- [ ] [Implement graph of layers](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/34);
- [ ] [Implement projected monoplex and overlay graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/35);
- [ ] [Implement more default multilayer graphs](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues/36) (e.g. multiplex graphs).

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
 
