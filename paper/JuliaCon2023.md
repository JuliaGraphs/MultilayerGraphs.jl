# JuliaCon 2023

## Title 

MultilayerGraphs.jl: Multilayer Network Science in Julia

## Session Type 

Lightning talk (10 minutes) 

## Abstract 

**MultilayerGraphs.jl** is a Julia package for the creation, manipulation and analysis of multilayer graphs, which have been adopted to model of a wide range of complex systems from bio-chemical to socio-technical networks.

It is integrated within the **JuliaGraphs** ecosystem extending **Graphs.jl** with several multilayer-specific methods and metrics including the multilayer eigenvector centrality, the multilayer modularity and the Von Newman entropy. 

## Description 

[**MultilayerGraphs.jl**](https://github.com/JuliaGraphs/MultilayerGraphs.jl) is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs. 

A multilayer graph is a graph consisting of multiple subgraphs called **layers** which can be interconnected through [bipartite graphs](https://en.wikipedia.org/wiki/Bipartite_graph) called **interlayers**. 

Multiple theoretical frameworks have been proposed to formally subsume all instances of multilayer graphs ([De Domenico  et al. (2013)](https://doi.org/10.1103/physrevx.3.041022); [Kivelä et al. (2014)](https://doi.org/10.1093/comnet/cnu016); [Bianconi (2018)](https://doi.org/10.1093/oso/9780198753919.001.0001); [De Domenico (2022)](https://doi.org/10.1007/978-3-030-75718-2)). 

Multilayer graphs have been adopted to model the structure and dynamics of a wide spectrum of high-dimensional, non-linear, multi-scale, time-dependent complex systems including physical, chemical, biological, neuronal, socio-technical, epidemiological, ecological and economic networks ([Cozzo et al. (2013)](https://doi.org/10.1103/physreve.88.050801); [Estrada and Gomez-Gardenes (2014)](https://doi.org/10.1103/physreve.89.042819); [Baggio et al. (2016)](https://doi.org/10.1073/pnas.1604401113); [DeDomenico et al. (2016)](https://doi.org/10.1038/nphys3865); [Amato et al. (2017)](https://doi.org/10.1038/s41598-017-06933-2); [DeDomenico (2017)](https://doi.org/10.1093/gigascience/gix004); [Timteo et al. (2018)](https://doi.org/10.1038/s41467-017-02658-y); [Aleta et al. (2020)](https://doi.org/10.1038/s41562-020-0931-9); [Aleta et al. (2022)](https://doi.org/10.1073/pnas.2112182119)). 

MultilayerGraphs.jl is an integral part of the [JuliaGraphs](https://github.com/JuliaGraphs) ecosystem extending [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) so all the methods and metrics exported by Graphs.jl work for multilayer graphs, but due to the special nature of multilayer graphs the package features an implementation that maps a standard integer-labelled vertex representation to a more user-friendly framework exporting all the objects a practitioner would expect such as nodes (`Node`), vertices (`MultilayerVertex`), layers (`Layer`), interlayer (`Interlayer`), etc.

The two main data structures are `MultilayerGraph` and `MultilayerDiGraph`: collections of layers connected through interlayers. 

The **vertices** of a multilayer graph are representations of one set of distinct objects called `Node`s. Each layer may represent all the node set or just a subset of it. The vertices of `Multilayer(Di)Graph` are implemented via the `MultilayerVertex` custom type. Each `MultilayerVertex` encodes information about the node it represents, the layer it belongs to and its metadata. 

Both the **intra-layer** and **inter-layer edges** are embedded in the `MultilayerEdge` struct, whose arguments are the two connected multilayer vertices, the edge weight and its metadata. It's important to highlight that `Multilayer(Di)Graph`s are weighted and able to store metadata by default.

The **layers** are implemented via the `Layer` struct composed of an underlying graph and a mapping from its integer-labelled vertices to the collection of `MultilayerVertex`s the layer represents. **Interlayers** are similarly implemented via the `Interlayer` mutable struct, and they are generally constructed by providing the two connected layers, the (multilayer) edge list between them and a graph. This usage of underlying graphs allows for an easier debugging procedure during construction and a more intuitive analysis afterwards allowing the package to leverage all the features of the JuliaGraphs ecosystem so that it can be effectively considered as a real proving ground of its internal consistency.

The `Multilayer(Di)Graph` structs are weighted and endowed with the functionality to store both vertex-level and edge-level metadata by default so that at any moment the user may add or remove a `Layer` or specify an `Interlayer` and since different layers and interlayers could be better represented by graphs that are weighted or unweighted and with or without metadata, it was crucial for us to provide the most general and adaptable structure. A `Multilayer(Di)Graph` is instantiated by providing the ordered list of layers and the list of interlayers to the constructor. The latter are automatically specified, so there is no need to instantiate all of them. 

Finally, MultilayerGraphs.jl has been integrated within the [JuliaDynamics](https://github.com/JuliaDynamics) ecosystem so that any `Multilayer(Di)Graph` can be an argument to the `GraphSpace` constructor in [Agents.jl](https://github.com/JuliaDynamics/Agents.jl). 

For a more comprehensive exploration of the package features and functionalities we strongly recommend consulting the [tutorial](https://juliagraphs.org/MultilayerGraphs.jl/stable/#Tutorial) included in the package documentation. 

## Notes 

### Statement of Need 
At the best of our knowledge there are currently no software packages dedicated to the creation, manipulation and analysis of multilayer graphs implemented in the Julia language apart from MultilayerGraphs.jl itself.

### Future Developments 
See the relevant [issues](https://github.com/JuliaGraphs/MultilayerGraphs.jl/issues). 

### Publication
We [are are going to submit / have submitted] a [paper] to [JOSS](https://joss.theoj.org) presenting MultilayerGraphs.jl. 

MultilayerGraphs.jl has been published on the following open access websites: 
- [GitHub](https://github.com/JuliaGraphs/MultilayerGraphs.jl);
- [Zenodo](https://doi.org/10.5281/zenodo.7009172);
- [Research Software Directory](https://research-software-directory.org/software/multilayergraphs).

MultilayerGraphs.jl and its features were announced on the following platforms:
- [Discourse](https://discourse.julialang.org/t/ann-multilayergraphs-jl-a-package-to-construct-handle-and-analyse-multilayer-graphs/85988);
- [Forem](https://forem.julialang.org/inphyt/ann-multilayergraphsjl-a-package-to-construct-handle-and-analyse-multilayer-graphs-3k22);
- [Twitter](https://twitter.com/In_Phy_T/status/1560594513189638146).

## Session Image

Social preview or logo. 

## Additional Speaker 

Claudio's email address (the same used for GitHub?). 