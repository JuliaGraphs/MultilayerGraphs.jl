using Revise, Test, Logging, LoggingExtras, Beep
using DataStructures
using StatsBase, Distributions #, Bijections, SparseArrays
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs
using Agents
using MultilayerGraphs

degree_sequence = degree(SimpleGraph(10, 20))

havel_hakimi_graph_generator_unsafe(degree_sequence)

graph = SimpleGraph(length(degree_sequence))
v_ds = OrderedDict(i => deg for (i, deg) in enumerate(degree_sequence))

while (any(values(v_ds) .!= 0))
    v_ds = OrderedDict(sort(collect(v_ds); by=last, rev=true))
    S, s = popfirst!(v_ds)
    @debug "" S, s
    all(collect(values(v_ds))[1:s] .> 0) ||
        throw(ErrorException("The provided `degree_sequence` is not graphical"))
    # println("here")
    for v in collect(keys(v_ds))[1:s]
        add_edge!(graph, S, v)
        @debug "", v, length(v_ds), collect(keys(v_ds))
        v_ds[v] -= 1
        println("here")
    end
end

degree(graph) .== degree_sequence

# cd("dev/MultilayerGraphs")
# cd("dev/MultilayerGraphs/test")

# Create the logger ushc that it prints debug messages only coming from our package. Taken from https://julialogging.github.io/how-to/filter-messages/#How-to-filter-messages
logger = EarlyFilteredLogger(ConsoleLogger(stderr, Logging.Debug)) do args
    r = args._module === MultilayerGraphs || args._module === Main
    return r
end
global_logger(logger)

const vertextype = Int64
const _weighttype = Float64
const min_vertices = 5
const max_vertices = 7
const min_edges = 1
const max_edges = max_vertices * (max_vertices - 1)

const multilayer_nodes = [Node("node_$i") for i in 1:max_vertices]
const mvs_layers = MV.(multilayer_nodes)
# Multilayer vertices with metadata (to be used with MetaGraphs, SimpleValueGraphs, etc)
const mvs_metadata = [MV(node, ("I'm node $(node.id)",)) for node in multilayer_nodes]

function _get_srcmv_dstmv_layer(layer::Layer)
    mvs = MultilayerGraphs.get_bare_mv.(collect(mv_vertices(layer)))

    src_mv = nothing
    _collection = []

    while isempty(_collection)
        src_mv = rand(mvs)
        _collection = setdiff(
            Set(mvs),
            Set(
                vcat(MultilayerGraphs.get_bare_mv.(mv_outneighbors(layer, src_mv)), src_mv)
            ),
        )
    end

    dst_mv = MultilayerGraphs.get_bare_mv(rand(_collection))

    return mvs, src_mv, dst_mv
end

_nv = rand(min_vertices:max_vertices)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
layer_sg = Layer(
    :layer_sg,
    sample(mvs_layers, _nv; replace=false),
    _ne,
    SimpleGraph{vertextype}(),
    _weighttype,
)

_nv = rand(min_vertices:max_vertices)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
layer_swg = Layer(
    :layer_swg,
    sample(mvs_layers, _nv; replace=false),
    _ne,
    SimpleWeightedGraph{vertextype,_weighttype}(),
    _weighttype;
    default_edge_weight=(src, dst) -> rand(),
)

_nv = rand(min_vertices:max_vertices)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))

layer_mg = Layer(
    :layer_mg,
    sample(mvs_metadata, _nv; replace=false),
    _ne,
    MetaGraph{vertextype,_weighttype}(),
    _weighttype;
    default_edge_metadata=(src, dst) -> (from_to="from_$(src)_to_$(dst)",),
)

_nv = rand(min_vertices:max_vertices)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
layer_vg = Layer(
    :layer_vg,
    sample(mvs_metadata, _nv; replace=false),
    _ne,
    MultilayerGraphs.ValGraph(
        SimpleGraph{vertextype}();
        edgeval_types=(Float64, String),
        edgeval_init=(s, d) -> (s + d, "hi"),
        vertexval_types=(String,),
        vertexval_init=v -> ("$v",),
    ),
    _weighttype;
    default_edge_metadata=(src, dst) -> (rand(), "from_$(src)_to_$(dst)"),
    default_vertex_metadata=mv ->
        ("This metadata have been generated via the default_vertex_metadata method",),
)

all_layers_u = [layer_sg, layer_swg, layer_mg, layer_vg];

_nv = nv(layer_sg) + nv(layer_swg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_sg_swg = Interlayer(layer_sg, layer_swg, _ne, SimpleGraph{vertextype}())

_nv = nv(layer_swg) + nv(layer_mg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_swg_mg = Interlayer(
    layer_swg,
    layer_mg,
    _ne,
    SimpleWeightedGraph{vertextype,_weighttype}();
    default_edge_weight=(x, y) -> rand(),
)

_nv = nv(layer_mg) + nv(layer_vg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_mg_vg = Interlayer(
    layer_mg,
    layer_vg,
    _ne,
    MetaGraph{vertextype,_weighttype}();
    default_edge_metadata=(x, y) -> (mymetadata=rand(),),
    transfer_vertex_metadata=true,
)  #SimpleWeightedGraph{Int64, Float64}()

interlayer_multiplex_sg_mg = multiplex_interlayer(
    layer_sg,
    layer_mg,
    ValGraph(
        SimpleGraph{vertextype}();
        edgeval_types=(from_to=String,),
        edgeval_init=(s, d) -> (from_to = "from_$(s)_to_$(d)"),
    );
    default_edge_metadata=(x, y) -> (from_to="from_$(src)_to_$(dst)",),
)

interlayer_empty_sg_vg = empty_interlayer(layer_sg, layer_vg, SimpleGraph{vertextype}())

all_interlayers_u = [
    interlayer_sg_swg,
    interlayer_swg_mg,
    interlayer_mg_vg,
    interlayer_multiplex_sg_mg,
    interlayer_empty_sg_vg,
]

multilayergraph = MultilayerGraph(all_layers_u, all_interlayers_u) #, all_interlayers_u

# Test random multilayergraph (configuration model)

# We need to wait for this PR: https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl/pull/14 top be merged before being able to use SimpleWeightedGraphs in the configuration model. Also we cannot use AbstractValgraphs until we implement default weight and metadata functions for layers and interlayers to be kept into the descriptor.

layers_to_be_emptied = deepcopy([
    layer for layer in all_layers_u if
    !(layer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph)
])

layers_names_to_be_emptied = name.(layers_to_be_emptied)

interlayers_to_be_emptied = deepcopy([
    interlayer for interlayer in all_interlayers_u if
    all(in.(interlayer.layers_names, Ref(layers_names_to_be_emptied))) &&
    !(interlayer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph)
])

for layer in layers_to_be_emptied
    for edge in edges(layer)
        rem_edge!(layer, edge)
    end
end

for interlayer in interlayers_to_be_emptied
    for edge in edges(interlayer)
        rem_edge!(interlayer, edge)
    end
end

@test all(ne.(layers_to_be_emptied) .== 0)

@test all(ne.(interlayers_to_be_emptied) .== 0)

# Instantiate configuration-model multilayergraph

# configuration_multilayergraph = MultilayerGraph(layers_to_be_emptied, interlayers_to_be_emptied, truncated(Normal(10), 0.0, 20.0) );

# histogram(degree(configuration_multilayergraph), bins = 20)

# Test layers and interlayers
@test nl(multilayergraph) == length(all_layers_u)
_nl = nl(multilayergraph)
@test nIn(multilayergraph) == _nl * (_nl - 1) / 2

# Test get_interlayer
for (layer_1, layer_2) in
    Iterators.product(multilayergraph.layers_names, multilayergraph.layers_names)
    if layer_1 != layer_2
        interlayer = get_interlayer(multilayergraph, layer_1, layer_2)
        @test nv(interlayer) >= 0
    end
end

# Test add_layer! and rem_layer!
const n_vertices_missing = rand(min_vertices:max_vertices)
const n_edges_missing_u = rand(
    n_vertices_missing:((n_vertices_missing * (n_vertices_missing - 1)) ÷ 2)
)
missing_layer_u = Layer(
    :missing_layer_u,
    sample(mvs_metadata, n_vertices_missing; replace=false),
    n_edges_missing_u,
    MultilayerGraphs.ValGraph(
        SimpleGraph{vertextype}();
        edgeval_types=(Float64, String),
        edgeval_init=(s, d) -> (s + d, "missing vertex $(s+d)"),
        vertexval_types=(String,),
        vertexval_init=v -> ("$(v^2)",),
    ),
    _weighttype;
    default_edge_metadata=(src, dst) -> (rand(), "missing edge from_$(src)_to_$(dst)"),
) # SimpleGraph{vertextype}()

@test !has_layer(multilayergraph, missing_layer_u.name)
@test add_layer!(multilayergraph, missing_layer_u)
@test has_layer(multilayergraph, missing_layer_u.name)
@test rem_layer!(multilayergraph, missing_layer_u.name)
@test !has_layer(multilayergraph, missing_layer_u.name)

# Test nodes
@inferred(nodes(multilayergraph))

@inferred(nn(multilayergraph))

for node in nodes(multilayergraph)
    @test has_node(multilayergraph, node)
end

## Test MultilayerGraphs.add_node! and MultilayerGraphs.rem_node!
new_node = Node("new_node")
nv_prev = nv(multilayergraph)
ne_prev = ne(multilayergraph)
@test !has_node(multilayergraph, new_node)
@test MultilayerGraphs.add_node!(multilayergraph, new_node)
@test has_node(multilayergraph, new_node)
@test MultilayerGraphs.rem_node!(multilayergraph, new_node)
@test !has_node(multilayergraph, new_node)
### Test that nothing changed
@test nv_prev == nv(multilayergraph)
@test ne_prev == ne(multilayergraph)

# Test vertices
@test eltype(multilayergraph) == Int64
nv(multilayergraph)
@test length(multilayergraph.fadjlist) == length(vertices(multilayergraph)) # nv_withmissing(multilayergraph) ==

## Test that all multilayer vertices are present
for mv in vcat(mv_vertices.(all_layers_u)...)
    @test has_vertex(multilayergraph, mv)
end

for mv in mv_vertices(multilayergraph)
    # if !(mv isa MissingVertex)
    mv_inneighbors(multilayergraph, mv)
    mv_outneighbors(multilayergraph, mv)
    # end
end

## Test set_edge_weight!
_, rand_mv_1_weight, rand_mv_2_weight = _get_srcmv_dstmv_layer(layer_swg)
_weight = 3.14
@test !has_edge(multilayergraph, rand_mv_1_weight, rand_mv_2_weight)
@test add_edge!(multilayergraph, rand_mv_1_weight, rand_mv_2_weight, weight=_weight)
@test has_edge(multilayergraph, rand_mv_1_weight, rand_mv_2_weight)
wgt = weight_tensor(multilayergraph)
@test wgt[rand_mv_1_weight, rand_mv_2_weight] ==
    wgt[rand_mv_2_weight, rand_mv_1_weight] ==
    get_weight(multilayergraph, rand_mv_1_weight, rand_mv_2_weight) ==
    get_weight(multilayergraph, rand_mv_2_weight, rand_mv_1_weight) ==
    _weight
@test set_weight!(multilayergraph, rand_mv_1_weight, rand_mv_2_weight, _weight + 1)
wgt = weight_tensor(multilayergraph)
@test wgt[rand_mv_1_weight, rand_mv_2_weight] ==
    wgt[rand_mv_2_weight, rand_mv_1_weight] ==
    get_weight(multilayergraph, rand_mv_1_weight, rand_mv_2_weight) ==
    get_weight(multilayergraph, rand_mv_2_weight, rand_mv_1_weight) ==
    _weight + 1

## Test set_metadata!
_, rand_mv_1_meta, rand_mv_2_meta = _get_srcmv_dstmv_layer(layer_mg)
### On vertices
@test set_metadata!(multilayergraph, rand_mv_1_meta, (meta="new_metadata",))
@test get_metadata(multilayergraph, rand_mv_1_meta).meta == "new_metadata"
### On edges
@test !has_edge(multilayergraph, rand_mv_1_meta, rand_mv_2_meta)
@test add_edge!(multilayergraph, rand_mv_1_meta, rand_mv_2_meta, metadata=(meta="hello",))
@test has_edge(multilayergraph, rand_mv_1_meta, rand_mv_2_meta)
mt = metadata_tensor(multilayergraph)
@test mt[rand_mv_1_meta, rand_mv_2_meta].meta ==
    mt[rand_mv_2_meta, rand_mv_1_meta].meta ==
    get_metadata(multilayergraph, rand_mv_1_meta, rand_mv_2_meta).meta ==
    get_metadata(multilayergraph, rand_mv_2_meta, rand_mv_1_meta).meta ==
    "hello"
_metadata = (meta="byebye",)
@test set_metadata!(multilayergraph, rand_mv_1_meta, rand_mv_2_meta, _metadata)
mt = metadata_tensor(multilayergraph)
@test mt[rand_mv_1_meta, rand_mv_2_meta].meta ==
    mt[rand_mv_2_meta, rand_mv_1_meta].meta ==
    get_metadata(multilayergraph, rand_mv_1_meta, rand_mv_2_meta).meta ==
    get_metadata(multilayergraph, rand_mv_2_meta, rand_mv_1_meta).meta ==
    "byebye"

# Test edges
ne(multilayergraph)

## Test that all edges are present
for edge in vcat(collect.(edges.(all_layers_u))..., collect.(edges.(all_interlayers_u))...)
    @test has_edge(multilayergraph, edge)
end

@test MultilayerGraphs.weighttype(multilayergraph) == Float64
@test edgetype(multilayergraph) == MultilayerEdge{Float64}

# Test set_edge_weight!

## Test add_edge! and rem_edge!

# Test Graphs.jl's extra overrides
@test all(indegree(multilayergraph) .== degree(multilayergraph)) #.+ outdegree(multilayergraph)

@inferred(mean_degree(multilayergraph))
@inferred(degree_second_moment(multilayergraph))
@inferred(degree_variance(multilayergraph))

# Test multilayer-specific methods
@test all(
    get_supra_weight_matrix_from_weight_tensor(weight_tensor(multilayergraph).array) .==
    supra_weight_matrix(multilayergraph).array,
)
@test all(
    get_weight_tensor_from_supra_weight_matrix(
        multilayergraph, supra_weight_matrix(multilayergraph).array
    ) .== weight_tensor(multilayergraph).array,
)
@test_broken multilayer_global_clustering_coefficient(multilayergraph) .==
    global_clustering_coefficient(multilayergraph)

overlaygraph = get_overlay_monoplex_graph(multilayergraph)
@test_broken global_clustering_coefficient(overlaygraph) .==
    overlay_clustering_coefficient(multilayergraph)

@test multilayer_weighted_global_clustering_coefficient(
    multilayergraph, [1 / 3, 1 / 3, 1 / 3]
) .≈ multilayer_global_clustering_coefficient(multilayergraph)

eig_centr_u, errs_u = eigenvector_centrality(multilayergraph; norm="n", tol=1e-3)

modularity(
    multilayergraph,
    rand([1, 2, 3, 4], length(nodes(multilayergraph)), length(multilayergraph.layers)),
)

get_graph_of_layers(multilayergraph)

wgt = weight_tensor(multilayergraph)
sam = supra_weight_matrix(multilayergraph)
for edge in collect(edges(multilayergraph.layer_swg))
    @debug "" src(edge) dst(edge) edge
    @test wgt[src(edge), dst(edge)] ==
        wgt[dst(edge), src(edge)] ==
        MultilayerGraphs.weight(edge)
    @test sam[src(edge), dst(edge)] ==
        sam[dst(edge), src(edge)] ==
        MultilayerGraphs.weight(edge)
end

# Test that, given a 1-dimensional multilayergraph, we obtain the same metrics as we would by using Graphs.jl's utilities on the one and only layer
## unweighted and weighted case
for layer in all_layers_u
    if !(layer.graph isa SimpleValueGraphs.AbstractValGraph)
        monolayergraph = MultilayerGraph([layer])

        @test length(edges(monolayergraph)) == length(edges(layer.graph))

        @test eltype(monolayergraph) == eltype(layer.graph)

        @test ne(monolayergraph) == ne(layer.graph)

        @test length(nodes(monolayergraph)) == nv(layer.graph)

        @test nv(monolayergraph) .== nv(layer.graph)

        @test all(
            inneighbors.(Ref(monolayergraph), vertices(monolayergraph)) .==
            inneighbors.(Ref(layer.graph), vertices(layer.graph)),
        )

        @test all(indegree(monolayergraph) .== indegree(layer.graph))

        @test all(outdegree(monolayergraph) .== outdegree(layer.graph))

        @test all(degree(monolayergraph) .== degree(layer.graph))

        @test_broken vec(eigenvector_centrality(monolayergraph; norm="n", tol=1e-3)[1]) ==
            eigenvector_centrality(layer.graph)

        tests = Bool[]
        for i in 1:5
            clustering = rand(
                [1, 2, 3], length(nodes(monolayergraph)), length(monolayergraph.layers)
            )
            push!(
                tests,
                modularity(monolayergraph, clustering) ==
                modularity(layer.graph, vec(clustering)),
            )
        end
        @test_broken all(tests)

        for edge in edges(layer)
            @test has_edge(monolayergraph, edge)
        end
    end
end
