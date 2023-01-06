using Revise
using Graphs, SimpleValueGraphs
using Distributions
using MultilayerGraphs

# The objects represented by `MultilayerVertex`s are `Node`s
const n_nodes = 100 
const node_list = [Node("node_$i") for i in 1:n_nodes]

#= 
Each layer and Interlayer that is involved in a `Multilayer(Di)Graph` must have the same 
vertex type (i.e. the type of the internal representation of vertices) and (edge) weight type. 
We proceed by defining them:
const vertex_type = Int64
const weight_type = Float64 
=#

#= 
We will instantiate layers and interlayers with randomly-selected edges and vertices adopting a portfolio of techniques.
Here, we define a layer with an underlying simple directed graph
using a graph generator-like (or "configuration model"-like) constructor which allows 
for specifying both the indegree and the outdegree sequences. 
    
Before instantiating each layer we sample its number of vertices and, optionally, of edges.
=#
n_vertices = rand(1:100)                          # Number of vertices 
layer_simple_directed = layer_simpledigraph(      # Layer constructor 
    :layer_simpledigraph,                         # Layer name
    sample(node_list, n_vertices; replace=false), # Nodes represented in the layer
    Truncated(Normal(5, 5), 0, 20),               # Indegree sequence sample distribution 
    Truncated(Normal(5, 5), 0, 20)                # Outdegree sequence sample distribution
)

#= Next, we define a layer with an underlying simple weighted directed graph. 
This time we show another kind of constructor that allows the user to specify 
the number of edges to be randomly distributed among the vertices. 
The keyword argument `default_edge_weight` will assign a weight to such edges before they are added to the layer
=#
n_vertices = rand(1:n_nodes)                                   # Number of vertices 
n_edges = rand(n_vertices:(n_vertices * (n_vertices - 1) - 1)) # Number of edges 
layer_simple_directed_weighted = layer_simpleweighteddigraph(  # Layer constructor 
    :layer_simpleweighteddigraph,                              # Layer name
    sample(node_list, n_vertices; replace=false),              # Nodes represented in the layer
    n_edges;                                                   # Number of randomly distributed edges
    default_edge_weight=(src, dst) -> rand()                   # Function assigning weights to edges 
)

#=
## Similar constructors, more flexible at the cost of ease of use, allows for finer tuning:
## NB: This constructor should be necessary only in rare circumstances, where e.g. 
the equivalent simplified constructor `layer_simplevaldigraph` is not able to infer
the correct return types of the `default_vertex/edge_metadata`s, or to use underlying 
graph that aren't currently supported.
=#
n_vertices = rand(1:n_nodes)                                   # Number of vertices 
n_edges = rand(n_vertices:(n_vertices * (n_vertices - 1) - 1)) # Number of edges 
default_vertex_metadata = v -> ("vertex_$(v)_metadata")        # Vertex metadata 
default_edge_metadata = (s, d) -> (rand(),)                    # Edge metadata 
layer_simple_directed_value = Layer(                           # Layer constructor
    :layer_simplevaldigraph,                                   # Layer name
    sample(node_list, n_vertices; replace=false),              # Nodes represented in the layer
    n_edges,                                                   # Number of randomly distributed edges
    ValDiGraph(                                                
        SimpleDiGraph{Int64}(); 
        vertexval_types=(String,),
        vertexval_init=default_vertex_metadata,
        edgeval_types=(Float64,),
        edgeval_init=default_edge_metadata,
    ),
    Float64;
    default_vertex_metadata=default_vertex_metadata,           # Vertex metadata 
    default_edge_metadata=default_edge_metadata                # Edge metadata 
)

layers = [layer_simple_directed, layer_simple_directed_weighted, layer_simple_directed_value]

#=
There are many more constructors that the user is encouraged to explore in the package documentation.
We may now move to Interlayers. Note that, in order to define a `Multilayer(Di)Graph`, 
interlayers do not need to be explicitly constructed by the user, 
since they are automatically specified by the `Multilayer(Di)Graph` constructor. 
Anyway, more complex interlayers need to be manually instantiated. 
The interface is very similar to the layers.
=#

# Interlayer with an underlying simple directed graph and `n_edges` edges
n_vertices_1 = nv(layer_simple_directed)               # Number of vertices of layer 1
n_vertices_2 = nv(layer_simple_directed_weighted)      # Number of vertices of layer 2
n_edges = rand(1:(n_vertices_1 * n_vertices_2 - 1))    # Number of interlayer edges 
interlayer_simple_directed = interlayer_simpledigraph( # Interlayer constructor 
    layer_simple_directed,                             # Layer 1 
    layer_simple_directed_weighted,                    # Layer 2 
    n_edges                                            # Number of edges 
)

## The interlayer exports a more flexible constructor too. 
n_vertices_1 = nv(layer_simple_directed_weighted)         # Number of vertices of layer 1
n_vertices_2 = nv(layer_simple_directed_value)            # Number of vertices of layer 2
n_edges = rand(1:(n_vertices_1 * n_vertices_2 - 1))       # Number of interlayer edges 
interlayer_simple_directed_meta = interlayer_metadigraph( # Interlayer constructor
    layer_simple_directed_weighted,                       # Layer 1 
    layer_simple_directed_value,                          # Layer 2
    n_edges;                                              # Number of edges
    default_edge_metadata=(src, dst) ->                   # Edge metadata 
        (edge_metadata="metadata_of_edge_from_$(src)_to_$(dst)"),
    transfer_vertex_metadata=true                         # Boolean deciding layer vertex metadata inheritance
)

interlayers = [interlayer_simple_directed, interlayer_simple_directed_meta]

# Layers and Interlayers are not immutable, and mostly behave like normal graphs. The reader is invited to consult the API for more information.
# A MultilayerDiGraph (i.e. a directed multilayer graph, following the naming convention of the JuliaGraph ecosystem) may now be specified 

multilayerdigraph = MultilayerDiGraph(
    layers,                                    # The (ordered) collection of layers
    interlayers;                               # The manually specified interlayers
                                               # The interlayers that are left unspecified 
                                               # will be automatically inserted according 
                                               # to the keyword argument below
    default_interlayers_structure="multiplex"  # The automatically specified interlayers will have only diagonal couplings
)

# Layers and interlayer may be accessed as properties using their names
multilayerdigraph.layer_simple_directed_value

# We proceed to show some basic functionality.
## Add a vertex
### The user may add vertices that do or do not represent node_list already present in the multilayergraph. In the latter case, we have a new node:
new_node_1 = Node("new_node_1")
### Before adding any vertex representing such node to the multilayer graph, the user should first add the Node:
add_node!(multilayerdigraph, new_node_1)
### Define a vertex that represents that node
new_vertex_1 = MV(           # Constructor (alias for "MultilayerVertex")
    new_node_1,              # Node represented by the vertex
    :layer_simplevaldigraph, # Layer containing the vertex 
    ("new_metadata")         # Vertex metadata 
)
### Add the vertex
add_vertex!(
    multilayerdigraph, # MultilayerDiGraph the vertex will be added to
    new_vertex_1       # MultilayerVertex to add
)
# NB: The vertex-level metadata are considered iff the graph underlying the layer/interlayer which the vertex belongs to supports them, otherwise they are discarded.

### `add_vertex!` implements multiple interfaces.

## Add an edge
### Let's represent another node in another layer
new_node_2 = Node("new_node_2")
new_vertex_2 = MV(new_node_2, :layer_simpledigraph)
add_vertex!(
    multilayerdigraph,
    new_vertex_2;
    add_node=true # Add the associated node before adding the vertex
)
### Construct a new edge
new_edge = MultilayerEdge( # Constructor 
    new_vertex_1,          # Source vertex
    new_vertex_2,          # Destination vertex 
    ("some_edge_metadata") # Edge metadata 
)
### Add the edge
add_edge!(
    multilayerdigraph, # MultilayerDiGraph the edge will be added to
    new_edge           # MultilayerVertex to add
)
# NB: The edge-level metadata and/or weight are considered iff the graph underlying the layer/interlayer which the edge belongs to supports them, otherwise they are discarded.

#=
Using the provided `add_layer!`, `rem_layer!` and `specify_interlayer!`, 
Layers and Interlayers may be added, removed or specified on the fly. 
Since MultilayerGraphs.jl extends Graphs.jl, all metrics from the JuliaGraphs ecosystem 
should be available by default. Anyway, some multilayer-specific metrics have been implemented, 
and others required to be re-implemented. 

We showcase a few of them here:
=#

## Compute the global clustering coefficient as in @DeDomenico2014
multilayer_global_clustering_coefficient(multilayerdigraph) # A weighted version `multilayer_weighted_global_clustering_coefficient`
## Compute the overlay clustering coefficient as in @DeDomenico2013
overlay_clustering_coefficient(multilayerdigraph)
## Compute the eigenvector centrality (the implementation is sp that it coincides with Graphs.jl's `eigenvector_centrality` on monoplex graphs)
eigenvector_centrality(multilayerdigraph)
## Compute the multilayer modularity as in @DeDomenico2013
modularity(
    multilayerdigraph,
    rand([1, 2, 3, 4], length(nodes(multilayerdigraph)), length(multilayerdigraph.layers))
)

## Currently, Von Neumann entropy is available only for undirected multilayer graphs.
# NB: this brief script is far from complete: many more features and functionalities are detailed in the documentation.
