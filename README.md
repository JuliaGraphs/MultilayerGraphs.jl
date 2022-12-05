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

**MultilayerGraphs.jl** is a Julia package for the construction, manipulation and analysis of multilayer graphs [extending Graphs.jl](https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/).

## Overview

**MultilayerGraphs.jl** implements the mathematical formulation of multilayer graphs proposed by [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). It mainly revolves around two custom types, [`MultilayerGraph`](https://juliagraphs.org/MultilayerGraphs.jl/dev/internals/#MultilayerGraphs.MultilayerGraph) and [`MultilayerDiGraph`](https://juliagraphs.org/MultilayerGraphs.jl/dev/internals/#MultilayerGraphs.MultilayerDiGraph), encoding undirected and directed multilayer graphs respectively.

Roughly speaking, a multilayer graph is a collection of ***layers***, i.e. graphs whose vertices are representations of the same set of nodes, and ***interlayers***, i.e the [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) whose vertices are those of any two layers and whose edges are those between vertices of the same two layers.

[`MultilayerGraph`](https://juliagraphs.org/MultilayerGraphs.jl/dev/internals/#MultilayerGraphs.MultilayerGraph) and [`MultilayerDiGraph`](https://juliagraphs.org/MultilayerGraphs.jl/dev/internals/#MultilayerGraphs.MultilayerDiGraph) are fully-fledged [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) extensions. Both structs are designed so that their layers and interlayers can be of any type (as long as they are Graphs.jl extensions themselves) and they need not be all of the same type. It is anyway required that all layers and interlayers of [`MultilayerGraph`](https://juliagraphs.org/MultilayerGraphs.jl/dev/internals/#MultilayerGraphs.MultilayerGraph) and [`MultilayerDiGraph`](https://juliagraphs.org/MultilayerGraphs.jl/dev/internals/#MultilayerGraphs.MultilayerDiGraph) are respectively undirected and directed. Directedness is checked via the `IsDirected` trait defined in Graphs.jl adopting [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl). Since the layers' and interlayers' graph types don't need to be the same, multilayer graph types are considered weighted graphs by default, and thus are assigned the trait `IsWeighted`.

## Installation

Press `]` in the Julia REPL and then

```julia
pkg> add MultilayerGraphs
```

## Tutorial

In the package documentation we have prepared a [tutorial](https://juliagraphs.org/MultilayerGraphs.jl/stable/#Tutorial) to illustrate how to define, handle and analyse a [`MultilayerGraph`](https://juliagraphs.org/MultilayerGraphs.jl/dev/internals/#MultilayerGraphs.MultilayerGraph) (the directed version is completely analogous).

## Future Developments

Here we highlight the major future developments we have currently identified:

- [ ] Better integration with [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl). We've recently opened an [issue](https://github.com/JuliaGraphs/Graphs.jl/issues/165) on Graphs.jl illustrating our proposals;
- [ ] Better integration with [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl) and [SimpleValueGraphs.jl](https://github.com/simonschoelly/SimpleValueGraphs.jl). Although it is possible to specify a `MetaGraph` and `SimpleValueGraph` as layer and/or interlayer, they are not yet fully supported (i.e. API may be a little unfit for them). An example using MetaGraphs, SimpleValueGraphs can be found at our announcement post [here]();
- [ ] Optimise the adjacency tensor;
- [ ] More intuitive constructor for `Interlayer`;
- [ ] Implement specialised and simplified API for `MultiplexGraph`;
- [ ] Implement visualisation functionalities;
- [ ] Implement other features and methods for the analysis of multilayer graphs following the scientific literature:
  - Kivelä et al. (2014) [Multilayer networks](https://doi.org/10.1093/comnet/cnu016). *Journal of Complex Networks*
  - Cozzo et al. (2015) [Structure of triadic relations in multiplex networks](https://doi.org/10.1088/1367-2630/17/7/073029). *New Journal of Physics*
  - De Domenico et al. (2015) [MuxViz: a tool for multilayer analysis and visualization of networks](https://doi.org/10.1093/comnet/cnu038). *Journal of Complex Networks*
  - De Domenico et al. (2015) [Ranking in interconnected multilayer networks reveals versatile nodes](https://doi.org/10.1038/ncomms7868). *Nature Communications*
  - De Domenico (2022) [Multilayer Networks: Analysis and Visualization](https://doi.org/10.1007/978-3-030-75718-2). *Springer Cham*
  - Artime et al. (2022) [Multilayer Network Science: From Cells to Societies](https://doi.org/10.1017/9781009085809). *Cambridge Elements*

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

De Domenico et al. (2013) [Mathematical Formulation of Multilayer Networks](https://doi.org/10.1103/PhysRevX.3.041022). *Physical Review X*.
