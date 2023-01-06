using Revise
using Graphs, SimpleValueGraphs
using Distributions
using MultilayerGraphs

# The objects that `MultilayerVertex`s represent are `Node`s
const n__nodes = 100
const _nodes = [Node("node_$i") for i in 1:n__nodes]

#= # Every layer and Interlayer that is involved in a `Multilayer(Di)Graph` must have the same vertex type (i.e. the type of the internal representation of vertices) and (edge) weight type. We proceed to define them hereafter:
const vertex_type = Int64
const weight_type = Float64 =#

# We will be instantiating layers and interlayers with randomly-chosen edges and vertices, using different techniques with the intent to showcase them.
## Here, we define a layer with an underlying simple directed graph, using a configuration model-like constructor that allows for specifying both the indegree and the outdegree sequences. Before instantiating each layer we will sample its  nuber of vertices and, optionally, of edges.
_nv = rand(1:100)
_layer_simpledigraph = layer_simpledigraph(
    :layer_simpledigraph,           # The name of the layer_metagraph
    sample(_nodes, _nv; replace=false),                          # The _nodes that the layer will represent
    Truncated(Normal(5, 5), 0, 20),  # The distribution from which the indegree sequence will be sampled from
    Truncated(Normal(5, 5), 0, 20),   # The distribution from which the outdegree sequence will be sampled from
)

## Next, we define a layer with an underliyng simple weighted directed graph. This time we show another kind of constructor that allows the user to specify the number of edges to be randomly distributed among the vertices. The keyword argument `default_edge_weight` will assign a weight to such edges before they are added to the layer
_nv = rand(1:n__nodes)
_ne = rand(_nv:(_nv * (_nv - 1) - 1))
_layer_simpleweighteddigraph = layer_simpleweighteddigraph(
    :layer_simpleweighteddigraph,
    sample(_nodes, _nv; replace=false),
    _ne;                                      # The number of edges that will be randomly distributed among the vertices
    default_edge_weight=(src, dst) -> rand(), # The function that assigns a weight to such edges before they are added to the layer
)

## Similar constructors, more flexible at the cost of ease of use, allows for finer tuning:
### NB: This constructor should be necessary only in rare circumstances, where e.g. the equivalent simplified constructor `layer_simplevaldigraph` is not able to infer the correct return types of the `default_vertex/edge_metadata`s, or to use underlying graph that aren't currently supported.
_nv = rand(1:n__nodes)
_ne = rand(_nv:(_nv * (_nv - 1) - 1))
default_vertex_metadata = v -> ("vertex_$(v)_metadata",)
default_edge_metadata = (s, d) -> (rand(),)
_layer_valdigraph = Layer(
    :layer_simplevaldigraph,
    sample(_nodes, _nv; replace=false),
    _ne,
    ValDiGraph(
        SimpleDiGraph{Int64}();
        vertexval_types=(String,),
        vertexval_init=default_vertex_metadata,
        edgeval_types=(Float64,),
        edgeval_init=default_edge_metadata,
    ),
    Float64;
    default_vertex_metadata=default_vertex_metadata,
    default_edge_metadata=default_edge_metadata,
)

layers = [_layer_simpledigraph, _layer_simpleweighteddigraph, _layer_valdigraph]

# There are many more constructors that the user is encouraged to explore in the package documentation.

# We may now move to Interlayers. Note that, in order to define a `Multilayer(Di)Graph`, interlayers do not need to be explicitly constructed by the user, since they are automatically specified by the `Multilayer(Di)Graph` constructor. Anyway, more complex interlyaers need to be manually instantiated. The interface is very similar to the layers.

## Interlayer with an underlying simple directed graph and `_ne` edges
nv_1 = nv(_layer_simpledigraph)
nv_2 = nv(_layer_simpleweighteddigraph)
_ne = rand(1:(nv_1 * nv_2 - 1))             # The interlayer is a bipartite graph between the vertices of the two layers
_interlayer_simpledigraph = interlayer_simpledigraph(
    _layer_simpledigraph,          # One of the two layers connected by this interlayer
    _layer_simpleweighteddigraph,  # One of the two layers connected by this interlayer
    _ne,                            # # The number of edges that will be randomly distributed among the vertices
)

## The interlayer exports a more flexible constructor too. 
nv_1 = nv(_layer_simpleweighteddigraph)
nv_2 = nv(_layer_valdigraph)
_ne = rand(1:(nv_1 * nv_2 - 1))
_interlayer_metadigraph = interlayer_metadigraph(
    _layer_simpleweighteddigraph,
    _layer_valdigraph,
    _ne;
    default_edge_metadata=(src, dst) ->
        (edge_metadata="metadata_of_edge_from_$(src)_to_$(dst)",),
    transfer_vertex_metadata=true, # Whether to have the vertices of the interlayer endowed with the metadata of the vertices of the layers it connects. By default it is set to false.
)

interlayers = [_interlayer_simpledigraph, _interlayer_metadigraph]

# Layers and Interlayers are not immutable, and mostly behave like normal graphs. The reader is invited to consult the API for more information.
# A MultilayerDiGraph (i.e. a directed multilayer graph, following the naming convention of the JuliaGraph ecosystem) may now be specified 

multilayerdigraph = MultilayerDiGraph(
    layers,     # The (ordered) collection of layers that consistute the multilayer graph
    interlayers; # The manually-specified interlayers. The interlayers that are left unspecified (in this example, the interlayer between `_layer_simpledigraph` and `_layer_valdigraph`), will be automatically inserted according to the keyword argument below
    default_interlayers_structure="multiplex", # The automatically-specified interlayers will have only diagonal couplings
)

# Layers and interlayer may be accessed as properties using their names
multilayerdigraph.layer_simplevaldigraph

# We proceed to show some basic funcionality.
## Add a vertex
### The user may add vertices that do or do not represent _nodes already present in the multilayergraph. In the latter case, we have a new node:
new_node_1 = Node("new_node_1")
### Before adding any vertex representing such node to the multilayer graph, the user should first add the Node:
add_node!(multilayerdigraph, new_node_1)
### Define a vertex that represents that node
new_vertex_1 = MV(                         # MV is an alias for MultilayerVertex
    new_node_1,              # The Node that the vetex represents 
    :layer_simplevaldigraph, # The layer which the vertex belongs to
    ("new_metadata",),        # Metadata to associate to the vertex.
)
### Add the vertex
add_vertex!(
    multilayerdigraph,  # The MultilayerDiGraph which the vertex will be added to
    new_vertex_1,       # The MultilayerVertex to add. Note that vertex-level metadata are considered iff the graph underlying the layer/interlayer which the vertex belongs to supports them, otherwise they are discarded.
)
### `add_vertex!` implements multiple interfaces.

## Add an edge
### Let's represent another node in another layer
new_node_2 = Node("new_node_2")
new_vertex_2 = MV(new_node_2, :layer_simpledigraph)
add_vertex!(
    multilayerdigraph,
    new_vertex_2;
    add_node=true,                                # Also perform `add_node!(new_node_2)` before adding the vertex
)
### Construct a new edge
new_edge = MultilayerEdge(
    new_vertex_1,           # The source of the edge
    new_vertex_2,           # The destination of the edge
    ("some_edge_metadata",), # Metadata to associate to the edge. 
)
### Add the edge
add_edge!(
    multilayerdigraph, # The MultilayerDiGraph which the edge will be added to
    new_edge,           # The MultilayerVertex to add. Note that edge-level metadata and/or weight are considered iff the graph underlying the layer/interlayer which the edge belongs to supports them, otherwise they are discarded.
)

# Using the provided `add_layer!`, `rem_layer!` and `specify_interlayer!`, Layers and Interlayers may be added, removed or specified on the fly. 

# Since MultilayerGraphs.jl extends Graphs.jl, all metrics from the JuliaGraphs ecosystem should be avaliable by default. Anyway, some  Multilayer graph-specific metrics have been implemented, and others required to be re-implemented.  A few of them are shown  hereafter

## Compute the global clustering coefficient as in @DeDomenico2014
multilayer_global_clustering_coefficient(multilayerdigraph) # A weighted version `multilayer_weighted_global_clustering_coefficient`, from the same authors ad been implemented

## Compute the overlay clustering coefficient as in @DeDomenico2014
overlay_clustering_coefficient(multilayerdigraph)

## Compute the eigenvector centality (the implementaation is sp that it coincides with Graphs.jl's `eigenvector_centality` on monoplex graphs)
eigenvector_centrality(multilayerdigraph)

## Compute the multilayer modularity as in @DeDomenico2014
modularity(
    multilayerdigraph,
    rand([1, 2, 3, 4], length(nodes(multilayerdigraph)), length(multilayerdigraph.layers)),
)

## Currently, Von Neumann entropy is available only for undirected multilayer graphs.

# NB: this brief script is far from complete: many more features and funcionalities are detailed in the documentation.
