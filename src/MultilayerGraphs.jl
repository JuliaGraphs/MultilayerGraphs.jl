module MultilayerGraphs

export getindex,
    δ_Ω,
    tensoreig,
    AbstractMultilayerGraph,
    AbstractMultilayerUGraph,
    AbstractMultilayerDiGraph,
    IsWeighted,
    MultilayerGraph,
    MultilayerDiGraph,
    MultiplexGraph,
    AbstractVertex,
    AbstractMultilayerVertex,
    MultilayerVertex,
    AbstractMultilayerEdge,
    MultilayerEdge,
    MultilayerWeightedEdge,
    GraphOfGraphs,
    DiGraphOfGraphs,
    AbstractBipartiteGraph,
    SimpleBipartiteGraph,
    add_edge!,
    Layer,
    Interlayer,
    get_symmetric_interlayer,
    add_layer!,
    get_layer,
    specify_interlayer!,
    get_interlayer,
    get_subgraph,
    indegree,
    outdegree,
    degree,
    neighbors,
    edges,
    is_directed,
    eltype,
    edgetype,
    has_edge,
    has_vertex,
    inneighbors,
    ne,
    nv,
    nn,
    outneighbors,
    vertices,
    nodes,
    mean_degree,
    degree_second_moment,
    degree_variance,
    multilayer_global_clustering_coefficient,
    multilayer_weighted_global_clustering_coefficient,
    overlay_clustering_coefficient,
    eigenvector_centrality,
    eigenvector_centrality_2,
    modularity,
    von_neumann_entropy,
    get_projected_monoplex_graph,
    get_overlay_monoplex_graph,
    get_graph_of_layers,
    get_oom

using Base, InteractiveUtils, IterTools, SimpleTraits
using LinearAlgebra, Statistics, OMEinsum, TensorOperations
using DataStructures, SparseArrays
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs

include("utilities.jl")
include("tensorsfactorizations.jl")
include("traits.jl")
include("multilayervertex.jl")
include("multilayeredge.jl")
include("layer.jl")
include("interlayer.jl")
include("graphs_extensions.jl")
include("graph_of_graphs.jl")
include("abstractmultilayergraph.jl")
include("abstractmultilayerugraph.jl")
include("abstractmultilayerdigraph.jl")
include("abstractmultiplexugraph.jl")
# include("abstractmultilayerdigraph.jl")
include("multilayergraph.jl")
include("multiplexgraph.jl")
include("multilayerdigraph.jl")


end
