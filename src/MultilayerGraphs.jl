# TODO: FIXME:

# SimpleWeightedGraphs cannot be used for the configuration/random multilayergraph until this PR: https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl/pull/14 is merged .

# Default value for  metadata should be `nothing` and not `NamedTuple()`? This would be more consistent with how edge weights behave.

# In a future release we may implement Multiplex(Di)Graph

# Implement https://en.wikipedia.org/wiki/Kleitman%E2%80%93Wang_algorithms and https://en.wikipedia.org/wiki/Havel%E2%80%93Hakimi_algorithm for wiring (un)directed configuration models.

# Double check modularity. Also compare with https://juliagraphs.org/Graphs.jl/v1.5/community/#Graphs.modularity-Tuple{AbstractGraph,%20AbstractVector{%3C:Integer}}. We should (probably) correct its implementation or at least compare it with a simpler one made in terms of unpadded supra adjacency matrix (that we have yet to implement).

# GraphOfGraphs and DiGraphOfGraphs could be improved/redesigned. Also, they don't yet extend Graphs.jl

# The δ_Ω implementation could be moved to a separate file.

# The usage of mutable `MissingVertex`s (although limited to SupraWeightMatrix) points to the fact that Bijections (at least as they are used right now) are not the best way to represent integer label-MultilayerVertex associations. We my implement our own container object or use the exixtsing ons differently.

# We need a quick Multilayer(Di)Graph constructor of the form Multilayer(Di)Graph(nn, nl; nv = rand(0:nn*nl), ne = rand(0:nv*(nv-1)) kwargs...) where kwargs may be used to further specify it.

module MultilayerGraphs

export 
    getindex,
    δ_Ω,
    tensoreig,
    AbstractMultilayerGraph,
    AbstractMultilayerUGraph,
    AbstractMultilayerDiGraph,
    IsWeighted,
    MultilayerGraph,
    MultilayerDiGraph,
    # MultiplexGraph,
    # MultiplexDiGraph,
    AbstractNode,
    Node,
    AbstractVertex,
    AbstractMultilayerVertex,
    MultilayerVertex,
    MV,
    ME,
    AbstractMultilayerEdge,
    MultilayerEdge,
    MultilayerWeightedEdge,
    AbstractIntraLayerEdge,
    IntraLayerEdge,
    # GraphOfGraphs,
    # DiGraphOfGraphs,
    add_edge!,
    LayerDescriptor,
    Layer,
    add_node!,
    rem_node!,
    set_prop!,
    get_prop,
    InterlayerDescriptor,
    Interlayer,
    get_symmetric_interlayer,
    multiplex_interlayer,
    is_multiplex_interlayer,
    has_node,
    add_vertex!,
    rem_vertex!,
    add_layer!,
    # get_layer,
    specify_interlayer!,
    get_interlayer,
    get_subgraph,
    indegree,
    outdegree,
    degree,
    neighbors,
    mv_neighbors,
    edges,
    is_directed,
    eltype,
    edgetype,
    has_edge,
    has_vertex,
    inneighbors,
    node_inneighbors,
    inneighbors_mv,
    ne,
    nv,
    nn,
    outneighbors,
    node_outneighbors,
    outneighbors_mv,
    nodes,
    mv_vertices,
    mean_degree,
    degree_second_moment,
    degree_variance,
    multilayer_global_clustering_coefficient,
    multilayer_weighted_global_clustering_coefficient,
    overlay_clustering_coefficient,
    eigenvector_centrality,
    modularity,
    von_neumann_entropy,
    get_projected_monoplex_graph,
    get_overlay_monoplex_graph,
    # get_graph_of_layers,
    get_oom,
    adjacency_matrix, 
    weights,
    get_bare_mv,
    get_rich_mv,
    mv_inneighbors,
    mv_outneighbors,
    getproperty,
    isequal,
    weighttype,
    MissingVertex,
    nv_withmissing,
    nl,
    nIn,
    get_supra_weight_matrix_from_weight_tensor,
    get_weight_tensor_from_supra_weight_matrix,
    rem_layer!,
    has_layer,
    weight_tensor,
    weight,
    supra_weight_matrix,
    empty_interlayer,
    isdigraphical,
    get_prop,
    name,
    set_weight!,
    set_metadata!,
    get_metadata,
    get_weight,
    MetadataTensor,
    metadata_tensor,
    WeightTensor,
    node,
    array

using Base, InteractiveUtils, IterTools, SimpleTraits, Bijections
using Distributions: Uniform
using LinearAlgebra, Statistics, OMEinsum, TensorOperations, Distributions
using DataStructures, SparseArrays
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs
using Agents

include("traits.jl")
include("tensorsfactorizations.jl")
include("node.jl")
include("vertices/abstractvertex.jl")
include("vertices/multilayervertex.jl")
include("vertices/missingvertex.jl")
include("multilayeredge.jl")
include("halfedge.jl")
include("subgraphs/layerdescriptor.jl")
include("subgraphs/interlayerdescriptor.jl")
include("subgraphs/abstractsubgraph.jl")
include("subgraphs/layer.jl")
include("subgraphs/interlayer.jl")
include("graphs_extensions/graphs_extensions.jl")
include("abstracttensorrepresentation.jl")
include("abstractmatrixrepresentation.jl")
include("weighttensor.jl")
include("metadatatensor.jl")
include("supraweightmatrix.jl")
include("abstractmultilayergraph.jl")
include("abstractmultilayerugraph.jl")
include("abstractmultilayerdigraph.jl")
include("multilayergraph.jl")
include("multilayerdigraph.jl")
include("utilities.jl")

end