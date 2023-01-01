# Summary

A multilayer graph is, loosely speaking, a collection of "layers" (represented by regular graphs) 
that also allows for links between the vertices of different layers. The bipartite graphs constituted 
by the two sets of vertices of two different layers and the edges between them are called "interlayers". 
The vertices in each layer represent a single set of nodes, although not all nodes have to be represented in every layer. 
There are multiple special cases of multilayer graphs, and multiple frameworks have been proposed to explain them all (see @Kivela2014). 
Common application of multilayer graphs are social network and epidemiological modeling.

From muxViz: "Multilayer networks are a class of models widely adopted to represent empirical complex system, including biomolecular networks (eg, interactomes, metabolomes), neuronal networks (eg, connectomes), information and communication networks, social/socio-technical/socio-ecological networks, economic and financial networks, urban and transportation networks."

From mully: "Network theory has been used for many years in the modelling and analysis of complex systems, as epidemiology, biology and biomedicine. As the data evolves and become more heterogeneous and complex, monoplex networks become an oversimplification of the corresponding systems. This imposes a need to go beyond traditional networks into a richer framework capable of hosting objects and relations of different scales and attributes called Multilayered Network."

# Statement of Need

The Julia graph ecosystem, which gravitates around the {Graph.jl}[@Graphs2021] package, was lacking an implementation of general multilayer graphs, particularly one that was simultaneously integrated with the main agent-based modeling library, {Agents.jl}[@Agents.jl]. 
Great care has been devoted to seamlessly integrate the package with the existing ecosystem, filling gaps in the latter where necessary:
- [Implementation of `isdigraphical` and fix of `isgraphical`](https://github.com/JuliaGraphs/Graphs.jl/pull/186);
- [Implementation of Havel-Hakimi and Kleitman-Wang algorithm for simple graph realization](https://github.com/JuliaGraphs/Graphs.jl/pull/202));
- [Better integration of Agents.jl with the graph ecosystem](https://github.com/JuliaDynamics/Agents.jl/pull/693);
- [Feedback on the state of the graph ecosystem](https://github.com/JuliaGraphs/Graphs.jl/issues/165).

This resulted in the creation of two API sets: one meant for the end-user, and the other for the developer.

# Main Features 

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