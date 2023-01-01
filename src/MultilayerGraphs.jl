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
    # Node.jl 
    AbstractNode,
    Node,
    id,
    # abstractvertex.jl
    AbstractVertex,
    # multilayervertex.jl 
    AbstractMultilayerVertex,
    MultilayerVertex,
    MV,
    node,
    layer,
    metadata,
    # missingvertex.jl
    MissingVertex,
    # multilayeredge.jl
    AbstractMultilayerEdge,
    MultilayerEdge,
    ME,
    weight,
    metadata,
    # halfedge.jl
    # layerdescriptor.jl
    # interlayerdescriptor.jl
    # abstractsubgraph.jl
    AbstractSubGraph,
    nodes,
    eltype,
    has_vertex,
    nv,
    vertices,
    mv_vertices,
    inneighbors,
    mv_inneighbors,
    outneighbors,
    mv_outneighbors,
    neighbors,
    mv_neighbors,
    get_v,
    edgetype,
    has_edge,
    ne,
    edges,
    add_edge!,
    rem_edge!,
    get_metadata,
    get_weight,
    is_directed,
    adjacency_matrix,
    weights,
    name,
    graph,
    # layer.jl
    AbstractLayer,
    Layer,
    layer_simplegraph,
    layer_simpledigraph,
    layer_simpleweightedgraph,
    layer_simpleweighteddigraph,
    layer_metadigraph,
    layer_valgraph,
    layer_valoutdigraph,
    layer_valdigraph,
    layer_metagraph,
    has_node,
    add_vertex!,
    rem_vertex!,
    # interlayer.jl
    AbstractInterlayer,
    Interlayer,
    interlayer_simplegrah,
    multiplex_interlayer,
    empty_interlayer,
    is_multiplex_interlayer,
    get_symmetric_interlayer,
    # abstracttensorrepresentation.jl
    AbstractTensorRepresentation,
    getindex,
    array,
    # abstractmatrixrrepresentation.jl
    AbstractMatrixRepresentation,
    # weighttensor.jl
    WeightTensor,
    # metadatatensor.jl
    MetadataTensor,
    # supraweightmatrix.jl
    SupraWeightMatrix,
    # traits.jl
    IsWeighted,
    is_weighted,
    IsMeta,
    is_meta,
    # abstractmultilayergraph.jl
    AbstractMultilayerGraph,
    nn,
    add_node!,
    rem_node!,
    set_metadata!,
    nl,
    nIn,
    has_layer,
    rem_layer!,
    get_interlayer,
    indegree,
    outdegree,
    degree,
    weighttype,
    weight_tensor,
    supra_weight_matrix,
    metadata_tensor,
    mean_degree,
    degree_second_moment,
    degree_variance,
    multilayer_global_clustering_coefficient,
    multilayer_weighted_global_clustering_coefficient,
    overlay_clustering_coefficient,
    eigenvector_centrality,
    modularity,
    # abstractmultilayerugraph.jl
    AbstractMultilayerUGraph,
    set_weight!,
    add_layer!,
    specify_interlayer!,
    von_neumann_entropy,
    # abstractmultilayerdigraph.jl
    AbstractMultilayerDiGraph,
    # multilayergraph.jl
    MultilayerGraph,
    # multilayerdigraph.jl
    MultilayerDiGraph,
    # utilities
    multilayer_kronecker_delta,
    δk,
    size,
    δ_1,
    δ_2,
    δ_3,
    δ_Ω,
    havel_hakimi_graph_generator,
    kleitman_wang_graph_generator
# tensorfacoriazations.jl

using Base, InteractiveUtils, IterTools, SimpleTraits, Bijections, PrettyTables
using Distributions: Uniform
using LinearAlgebra, Statistics, OMEinsum, TensorOperations, Distributions
using DataStructures, SparseArrays
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs

include("node.jl")
include("vertices/abstractvertex.jl")
include("vertices/multilayervertex.jl")
include("vertices/missingvertex.jl")
include("multilayeredge.jl")
include("halfedge.jl")
include("subgraphs/abstractdescriptor.jl")
include("subgraphs/layerdescriptor.jl")
include("subgraphs/interlayerdescriptor.jl")
include("subgraphs/abstractsubgraph.jl")
include("subgraphs/layer.jl")
include("subgraphs/interlayer.jl")
include("graphs_extensions/graphs_extensions.jl")
include("representations/abstracttensorrepresentation.jl")
include("representations/abstractmatrixrepresentation.jl")
include("representations/weighttensor.jl")
include("representations/metadatatensor.jl")
include("representations/supraweightmatrix.jl")
include("traits.jl")
include("abstractmultilayergraph.jl")
include("abstractmultilayerugraph.jl")
include("abstractmultilayerdigraph.jl")
include("multilayergraph.jl")
include("multilayerdigraph.jl")
include("utilities.jl")
include("tensorsfactorizations.jl")

end
