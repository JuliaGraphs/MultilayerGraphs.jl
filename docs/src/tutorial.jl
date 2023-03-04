# Install necessary tutorial dependencies 
using Pkg
Pkg.add(["Revise", "Distributions", "Graphs", "SimpleValueGraphs", 
         "LoggingExtras", "StatsBase", "SimpleWeightedGraphs", 
         "MetaGraphs", "Agents", "MultilayerGraphs"])

# Import necessary tutorial dependencies
using Revise
using StatsBase, Distributions
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs
using MultilayerGraphs

# Set the minimum and maximum number of nodes_list and edges for random graphs
const vertextype   = Int64
const _weighttype  = Float64
const min_vertices = 5
const max_vertices = 7
const n_nodes      = max_vertices

# The constructor for nodes (which are immutable) only requires a name (`id`) for the node
const nodes_list = [Node("node_$i") for i in 1:n_nodes]

## Convert nodes to multilayer vertices without metadata
const multilayervertices = MV.(nodes_list)
## Convert nodes multilayer vertices with metadata
const multilayervertices_meta  = [MV(node, ("I'm node $(node.id)",)) for node in nodes_list] # `MV` is an alias for `MultilayerVertex`

multilayervertices_meta[1]

# Utility function that returns a random number of vertices and edges each time it is called:
function rand_nv_ne_layer(min_vertices, max_vertices)
    _nv = rand(min_vertices:max_vertices)
    _ne = rand(1:(_nv*(_nv-1)) ÷ 2 )
    return (_nv,_ne)
end

# Utility function that returns two vertices of a Layer that are not adjacent.
function _get_srcmv_dstmv_layer(layer::Layer)
    mvs = MultilayerGraphs.get_bare_mv.(collect(mv_vertices(layer)))

    src_mv_idx = findfirst(mv -> !isempty(setdiff(
        Set(mvs),
        Set(
            vcat(MultilayerGraphs.get_bare_mv.(mv_outneighbors(layer, mv)), mv)
        ),
    )), mvs)

    src_mv = mvs[src_mv_idx]

    _collection = setdiff(
        Set(mvs),
        Set(
            vcat(MultilayerGraphs.get_bare_mv.(mv_outneighbors(layer, src_mv)), src_mv)
        ),
    )
    
    dst_mv = MultilayerGraphs.get_bare_mv(rand(_collection))

    return mvs, src_mv, dst_mv
end

# An unweighted simple layer:
_nv, _ne  = rand_nv_ne_layer(min_vertices,max_vertices)
layer_sg = Layer(   :layer_sg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    SimpleGraph{vertextype}(),
                    _weighttype
)

# A weighted `Layer`
_nv, _ne  = rand_nv_ne_layer(min_vertices,max_vertices)
layer_swg = Layer(  :layer_swg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    SimpleWeightedGraph{vertextype, _weighttype}(),
                    _weighttype;
                    default_edge_weight = (src,dst) -> rand()
)
# A `Layer` with an underlying `MetaGraph`:
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
layer_mg = Layer(   :layer_mg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    MetaGraph{vertextype, _weighttype}(),
                    _weighttype;
                    default_edge_metadata = (src,dst) -> (from_to = "from_$(src)_to_$(dst)",)
)
# `Layer` with an underlying `ValGraph` from `SimpleValueGraphs.jl`
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
layer_vg = Layer(   :layer_vg,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    MultilayerGraphs.ValGraph(SimpleGraph{vertextype}();
                                                edgeval_types=(Float64, String, ),
                                                edgeval_init=(s, d) -> (s+d, "hi"),
                                                vertexval_types=(String,),
                                                vertexval_init=v -> ("$v",),),
                                                _weighttype;
                                                default_edge_metadata = (src,dst) -> (rand(), "from_$(src)_to_$(dst)",),
                                                default_vertex_metadata = mv -> ("This metadata had been generated via the default_vertex_metadata method",)
)

# Collect all layers in an ordered list. Order will be recorded when instantiating the multilayer graph.
layers = [layer_sg, layer_swg, layer_mg, layer_vg]

# Utilities for Interlayer
## Utility function that returns two vertices of an Interlayer that are not adjacent.
function _get_srcmv_dstmv_interlayer(interlayer::Interlayer)

    mvs = get_bare_mv.(collect(mv_vertices(interlayer)))

    src_mv = nothing    
    _collection = []

    while isempty(_collection)
        src_mv = rand(mvs)
        _collection = setdiff(Set(mvs), Set(vcat(get_bare_mv.(mv_outneighbors(interlayer, src_mv)), src_mv, get_bare_mv.(mv_vertices( eval(src_mv.layer) ))) ) )  
    end

    dst_mv = get_bare_mv(rand(_collection))

    return mvs, src_mv, dst_mv
end

## Utility function that returns a random number edges between its arguments `layer_1` and `layer_2`:
function rand_ne_interlayer(layer_1, layer_2)
    nv_1 = nv(layer_1)
    nv_2 = nv(layer_2)
    _ne = rand(1: (nv_1 * nv_2 - 1) )
    return _ne
end

# Define the random undirected simple Interlayer
_ne = rand_ne_interlayer(layer_sg, layer_swg)
interlayer_sg_swg = Interlayer( layer_sg,                  # The first layer to be connected
                                layer_swg,                 # The second layer to be connected
                                _ne,                       # The number of edges to randomly generate
                                SimpleGraph{vertextype}(), # The underlying graph, passed as a null graph
                                interlayer_name = :random_interlayer  # The name of the interlayer. We will be able to access it as a property of the multilayer graph via its name. This kwarg's default value is given by a combination of the two layers' names.
)
# Define a weighted `Interlayer`
_ne = rand_ne_interlayer(layer_swg, layer_mg)
interlayer_swg_mg = Interlayer( layer_swg,
                                layer_mg,
                                _ne,
                                SimpleWeightedGraph{vertextype, _weighttype}();
                                default_edge_weight = (x,y) -> rand() # Arguments follow the same rules as in Layer
)
# Define an `Interlayer` with an underlying `MetaGraph`
_ne = rand_ne_interlayer(layer_mg, layer_vg)
interlayer_mg_vg = Interlayer(  layer_mg,
                                layer_vg,
                                _ne,
                                MetaGraph{vertextype, _weighttype}();
                                default_edge_metadata = (x,y) -> (mymetadata = rand(),),
                                transfer_vertex_metadata = true # This boolean kwarg controls whether vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors.
)
# Define an `Interlayer` with an underlying `ValGraph` from `SimpleValueGraphs.jl`, with diagonal couplings only:
interlayer_multiplex_sg_mg = multiplex_interlayer(  layer_sg,
                                                    layer_mg,
                                                    ValGraph(SimpleGraph{vertextype}(); edgeval_types=(from_to = String,), edgeval_init=(s, d) -> (from_to = "from_$(s)_to_$(d)"));
                                                    default_edge_metadata = (x,y) -> (from_to = "from_$(src)_to_$(dst)",)
)
# Finally, An `Interlayer` with no couplings (an "empty" interlayer):
interlayer_empty_sg_vg = empty_interlayer(  layer_sg,
                                            layer_vg,
                                            SimpleGraph{vertextype}()
)

# Collect all interlayers. Even though the list is ordered, order will not matter when instantiating the multilayer graph.
interlayers = [interlayer_sg_swg, interlayer_swg_mg, interlayer_mg_vg, interlayer_multiplex_sg_mg, interlayer_empty_sg_vg]

# Nodes
layer_sg_nodes = nodes(layer_sg)
interlayer_sg_swg_nodes  = nodes(interlayer_sg_swg)
has_node(layer_sg, layer_sg_nodes[1])

# Vertices
layer_sg_vertices = mv_vertices(layer_sg)
mv_vertices(layer_mg)
interlayer_sg_swg_vertices = mv_vertices(interlayer_sg_swg)

new_node     = Node("missing_node")
new_metadata = (meta = "my_metadata",)
new_vertex   = MV(new_node, new_metadata)

add_vertex!(layer_mg, new_vertex)
add_vertex!(layer_mg, new_node, metadata = new_metadata)
add_vertex!(layer_mg, new_node, Dict(pairs(new_metadata)))

metagraph = MetaGraph()
add_vertex!(metagraph,  Dict(pairs(new_metadata))) # true
rem_vertex!(layer_sg, new_vertex) # Returns true if succeeds

get_metadata(layer_mg, MV(new_node))

# Edges 
edgetype(layer_sg)
collect(edges(layer_sg))
# Define a weighted edge for the layer_swg
## Define the weight
_weight = rand()
## Select two non-adjacent vertices in layer_swg
_, src_w, dst_w  = _get_srcmv_dstmv_layer(layer_swg)
## Construct a weighted MultilayerEdge
me_w = ME(src_w, dst_w, _weight) # ME is an alias for MultilayerEdge

add_edge!(layer_swg, me_w)
add_edge!(layer_swg, src_w, dst_w, weight = _weight)
add_edge!(layer_swg, src_w, dst_w, _weight)

simpleweightedgraph = SimpleWeightedGraph(SimpleGraph(5, 0))
add_edge!(simpleweightedgraph, 1, 2, _weight)

rem_edge!(layer_swg, src_w, dst_w) # Returns true if succeeds

get_weight(layer_swg, src_w, dst_w)

# Define an edge with metadata for the layer_mg
## Define the metadata
_metadata  = (meta = "mymetadata",)
## Select two non-adjacent vertices in layer_mg
_, src_m, dst_m  = _get_srcmv_dstmv_layer(layer_mg)
## Construct a MultilayerEdge with metadata
me_m = ME(src_m, dst_m, _metadata)

add_edge!(layer_mg, me_m)
add_edge!(layer_mg, src_m, dst_m, metadata = _metadata)
add_edge!(layer_mg, src_m, dst_m, Dict(pairs(_metadata)))
get_metadata(layer_mg, src_m, dst_m)
add_edge!(layer_swg, me_w)
add_edge!(layer_swg, src_w, dst_w, weight = _weight)
add_edge!(layer_swg, src_w, dst_w, _weight)

rem_edge!(layer_swg, src_w, dst_w)

# Multilayer Graphs
multilayergraph = MultilayerGraph(  layers,                                                 # The (ordered) list of layers the multilayer graph will have
                                    interlayers;                                            # The list of interlayers specified by the user. Note that the user does not need to specify all interlayers, as the unspecified ones will be automatically constructed using the indications given by the `default_interlayers_null_graph` and `default_interlayers_structure` keywords.
                                    default_interlayers_null_graph = SimpleGraph{vertextype}(), # Sets the underlying graph for the interlayers that are to be automatically specified.  Defaults to `SimpleGraph{T}()`, where `T` is the `T` of all the `layers` and `interlayers`. See the `Layer` constructors for more information.
                                    default_interlayers_structure = "multiplex" # Sets the structure of the interlayers that are to be automatically specified. May be "multiplex" for diagonally coupled interlayers, or "empty" for empty interlayers (no edges).  "multiplex". See the `Interlayer` constructors for more information.
)

# The configuration model-like constructor will be responsible for creating the edges, so we need to provide it with empty layers and interlayers.
# To create empty layers and interlayers, we will empty the above subgraphs, and, for compatibility reasons, we'll remove the ones having a `SimpleWeightedGraph`s. These lines are not necessary to comprehend the tutorial, they may be skipped. Just know that the variables `empty_layers` and `empty_interlayers` are two lists of, respectively, empty layers and interlayers that do not have `SimpleWeightedGraph`s as their underlying graphs

empty_layers =  deepcopy([layer for layer in layers if !(layer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph)])

empty_layers_names = name.(empty_layers)

empty_interlayers =  deepcopy([interlayer for interlayer in interlayers if all(in.(interlayer.layers_names, Ref(empty_layers_names))) && !(interlayer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph) ])

for layer in empty_layers
    for edge in edges(layer)
        rem_edge!(layer, edge)
    end
end

for interlayer in empty_interlayers
    for edge in edges(interlayer)
        rem_edge!(interlayer, edge)
    end
end

# Construct a multilayer graph that has a normal degree distribution. The support of the distribution must be positive, since negative degrees are not possible
configuration_multilayergraph = MultilayerGraph(empty_layers, empty_interlayers, truncated(Normal(10), 0.0, 20.0));

new_node = Node("new_node")
add_node!(multilayergraph, new_node) # Return true if succeeds
new_vertex = MV(new_node, :layer_sg)
add_vertex!(multilayergraph, new_vertex)
rem_node!(multilayergraph, new_node) # Return true if succeeds

# This will succeed
random_weighted_edge = rand(collect(edges(multilayergraph.layer_swg)))
set_weight!(multilayergraph, src(random_weighted_edge), dst(random_weighted_edge), rand())

# This will not succeed
random_unweighted_edge = rand(collect(edges(multilayergraph.layer_sg)))
set_weight!(multilayergraph, src(random_unweighted_edge), dst(random_unweighted_edge), rand())

# Instantiate a new Layer
_nv, _ne = rand_nv_ne_layer(min_vertices,max_vertices)
new_layer = Layer(  :new_layer,
                    sample(nodes_list, _nv, replace = false),
                    _ne,
                    SimpleGraph{vertextype}(),
                    _weighttype
)

# Add the Layer
add_layer!(
            multilayergraph,                                # the `Multilayer(Di)Graph` which the new layer will be added to;
            new_layer;                                      # the new `Layer` to add to the `multilayergraph`
            default_interlayers_null_graph = SimpleGraph{vertextype}(), # upon addition of a new `Layer`, all the `Interlayer`s between the new and the existing `Layer`s are immediately created. This keyword argument specifies their `null_graph` See the `Layer` constructor for more information. Defaults to `SimpleGraph{T}()`
            default_interlayers_structure = "empty"         # The structure of the `Interlayer`s created by default. May either be "multiplex" to have diagonally-coupled only interlayers, or "empty" for empty interlayers. Defaults to "multiplex".
)

# Check that the new layer now exists within the multilayer graph
has_layer(multilayergraph, :new_layer)

# Instantiate a new Interlayer. Notice that its name will be given by default as
_ne = rand_ne_interlayer(layer_sg, new_layer)
new_interlayer = Interlayer(    layer_sg,                
                                new_layer,               
                                _ne,                     
                                SimpleGraph{vertextype}(),
                                interlayer_name = :new_interlayer
)

# Modify an existing interlayer with the latter i.e. specify the latter interlayer:
specify_interlayer!( multilayergraph,
                     new_interlayer)

# Now the interlayer between `layer_sg` and `new_layer` is `new_interlayer`

# Get a layer by name
multilayergraph.new_layer

# Get an Interlayer by name
multilayergraph.new_interlayer

# Get an Interlayer from the names of the two layers that it connects
get_interlayer(multilayergraph, :new_layer, :layer_sg )

# Remove the layer. This will also remove all the interlayers associated to it.
rem_layer!( multilayergraph,
            :new_layer;
            remove_nodes = false # Whether to also remove all nodes represented in the to-be-removed layer from the multilayer graph
)

wgt = weight_tensor(multilayergraph)
array(wgt)

# Get two random vertices from the MultilayerGraph
mv1, mv2 = rand(mv_vertices(multilayergraph), 2)

# Get the strength of the edge between them (0 for no edge):
wgt[mv1, mv2]