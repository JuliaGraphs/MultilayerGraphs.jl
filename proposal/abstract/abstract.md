# Abstract

## Title 

Multilayer Network Science in Julia with MultilayerGraphs.jl

## Speakers 

- Claudio Moroni (University of Turin);
- Pietro Monticone (University of Turin).

## Text 

[**MultilayerGraphs.jl**](https://github.com/JuliaGraphs/MultilayerGraphs.jl) is a Julia package for the creation, manipulation and analysis of the structure, dynamics and functions of multilayer graphs.

A multilayer graph consists of multiple subgraphs called **layers** which can be interconnected through bipartite graphs called **interlayers** composed of the sets of vertices of two different layers and the edges between them.

In order to formally represent multilayer networks, multiple theoretical paradigms 
have been proposed and adopted to model a wide spectrum of high-dimensional, 
multi-scale, time-dependent complex systems including molecular,
neuronal, social, ecological and economic networks.

The package features an implementation that maps a standard integer-labelled vertex representation to a more user-friendly framework exporting all the objects a practitioner would expect such as nodes, vertices, layers, interlayers, etc.

MultilayerGraphs.jl is integrated within the [**JuliaGraphs**](https://github.com/JuliaGraphs) ecosystem extending [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) with several multilayer-specific methods and metrics and within the [**JuliaDynamics**](https://github.com/JuliaDynamics) ecosystem allowing for agent-based modelling on general multilayer networks via [Agents.jl](https://github.com/JuliaDynamics/Agents.jl).

MultilayerGraphs.jl has been integrated within the [JuliaGraph](https://github.com/JuliaGraphs) and the [JuliaDynamics](https://github.com/JuliaDynamics) ecosystems through: 

- the extension of [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) with several methods and metrics including the multilayer eigenvector centrality, the multilayer modularity and the Von Newman entropy; 
- the compatibility with [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) allowing for agent-based modelling on general multilayer networks. 

For a comprehensive exploration of the package features and functionalities the readers is invited to consult the [README](https://github.com/JuliaGraphs/MultilayerGraphs.jl/blob/main/README.md) and [documentation](https://juliagraphs.org/MultilayerGraphs).