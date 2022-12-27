using Test, Logging, LoggingExtras
using StatsBase, Distributions
using Graphs, SimpleWeightedGraphs, MetaGraphs, SimpleValueGraphs
using Agents
using MultilayerGraphs

# Create the logger such that it prints debug messages only coming from our package. Taken from https://julialogging.github.io/how-to/filter-messages/#How-to-filter-messages
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

_nv = rand(min_vertices:max_vertices)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
layer_sg = Layer(
    :layer_sg,
    sample(mvs_layers, _nv; replace=false),
    _ne,
    SimpleGraph{vertextype}(),
    _weighttype,
)
_layer_simplegraph = layer_simplegraph(:layer_simplegraph, sample(mvs_layers, _nv, replace = false), Truncated(Normal(), 0, 10))


layer_sdg = Layer(
    :layer_sdg,
    sample(mvs_layers, _nv; replace=false),
    _ne,
    SimpleDiGraph{vertextype}(),
    _weighttype,
)
_layer_simpledigraph = layer_simpledigraph(:layer_simpledigraph, sample(mvs_layers, _nv, replace = false), Truncated(Normal(), 0, 10), Truncated(Normal(), 0, 10))

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
_layer_simpleweightedgraph = layer_simpleweightedgraph(:layer_simpleweightedgraph, sample(mvs_layers, _nv, replace = false), Truncated(Normal(), 0, 10), default_edge_weight = (src,dst) -> 3.0)

layer_swdg = Layer(
    :layer_swdg,
    sample(mvs_layers, _nv; replace=false),
    _ne,
    SimpleWeightedDiGraph{vertextype,_weighttype}(),
    _weighttype;
    default_edge_weight=(src, dst) -> rand(),
)
_layer_simpleweighteddigraph = layer_simpleweighteddigraph(:layer_simpleweighteddigraph, sample(mvs_layers, _nv, replace = false), Truncated(Normal(), 0, 10), Truncated(Normal(), 0, 10), default_edge_weight = (src,dst) -> 3.0)

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
_layer_metagraph = layer_metagraph(:layer_metagraph, sample(mvs_metadata, _nv, replace = false), Truncated(Normal(), 0, 10), default_vertex_metadata = mv -> (metamv = "metadata of $mv",), default_edge_metadata = (src,dst) -> (metaedge = "metadata of edge from $src to $dst",))


layer_mdg = Layer(
    :layer_mdg,
    sample(mvs_metadata, _nv; replace=false),
    _ne,
    MetaDiGraph{vertextype,_weighttype}(),
    _weighttype;
    default_edge_metadata=(src, dst) -> (from_to="from_$(src)_to_$(dst)",),
)
_layer_metadigraph = layer_metadigraph(:layer_metadigraph, sample(mvs_metadata, _nv, replace = false), Truncated(Normal(), 0, 10), Truncated(Normal(), 0, 10), default_vertex_metadata = mv -> (metamv = "metadata of $mv",), default_edge_metadata = (src,dst) -> (metaedge = "metadata of edge from $src to $dst",))

_nv = rand(min_vertices:max_vertices)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
layer_vg = Layer(
    :layer_vg,
    sample(mvs_metadata, _nv; replace=false),
    _ne,
    ValGraph(
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
_layer_valgraph = layer_valgraph(:layer_valgraph, sample(mvs_metadata, _nv, replace = false), Truncated(Normal(), 0, 10), default_vertex_metadata = mv -> ("metadata of $mv",), default_edge_metadata = (src,dst) -> (metaedge = "metadata of edge from $src to $dst",))

layer_vodg = Layer(
    :layer_vodg,
    sample(mvs_metadata, _nv; replace=false),
    _ne,
    ValOutDiGraph(
        SimpleDiGraph{vertextype}();
        edgeval_types=(a=Float64, b=String),
        edgeval_init=(s, d) -> (a=s + d, b="hi"),
        vertexval_types=(String,),
        vertexval_init=v -> ("$v",),
    ),
    _weighttype;
    default_edge_metadata=(src, dst) -> (a=rand(), b="from_$(src)_to_$(dst)"),
)
_layer_valoutdigraph = layer_valoutdigraph(:layer_valoutgraph, sample(mvs_metadata, _nv, replace = false), Truncated(Normal(), 0, 10),Truncated(Normal(), 0, 10), default_vertex_metadata = mv -> ("metadata of $mv",), default_edge_metadata = (src,dst) -> (metaedge = "metadata of edge from $src to $dst",))

layer_vdg = Layer(
    :layer_vdg,
    sample(mvs_metadata, _nv; replace=false),
    _ne,
    ValDiGraph(
        SimpleDiGraph{vertextype}();
        edgeval_types=(Float64, String),
        edgeval_init=(s, d) -> (s + d, "hi"),
        vertexval_types=(String,),
        vertexval_init=v -> ("$v",),
    ),
    _weighttype;
    default_edge_metadata=(src, dst) -> (rand(), "from_$(src)_to_$(dst)"),
)
_layer_valdigraph = layer_valdigraph(:layer_valdigraph, sample(mvs_metadata, _nv, replace = false), Truncated(Normal(), 0, 10), Truncated(Normal(), 0, 10),default_vertex_metadata = mv -> ("metadata of $mv",), default_edge_metadata = (src,dst) -> (metaedge = "metadata of edge from $src to $dst",))



collect(edges(_layer_simpledigraph))




collect(edges(_simple_weighted_layer))


collect(edges(_simple_weighted_dilayer))



collect(edges(_layer_metagraph))
mv_vertices(_layer_metagraph)


collect(edges(_layer_metagraph))
mv_vertices(_layer_metagraph)



collect(edges(_layer_valgraph))
mv_vertices(_layer_valgraph)



collect(edges(_layer_valoutdigraph))
mv_vertices(_layer_valoutdigraph)




collect(edges(_layer_valdigraph))
mv_vertices(_layer_valdigraph)




all_layers = [
    layer_sg,
    layer_sdg,
    layer_swg,
    layer_swdg,
    layer_mg,
    layer_mdg,
    layer_vg,
    layer_vodg,
    layer_vdg,
];

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
    default_edge_metadata=(x, y) -> (from_to="from_$(x)_to_$(y)",),
)

interlayer_empty_sg_vg = empty_interlayer(layer_sg, layer_vg, SimpleGraph{vertextype}())

_nv = nv(layer_sdg) + nv(layer_swdg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_sdg_swdg = Interlayer(layer_sdg, layer_swdg, _ne, SimpleDiGraph{vertextype}())

_nv = nv(layer_swdg) + nv(layer_mdg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_swdg_mdg = Interlayer(
    layer_swdg,
    layer_mdg,
    _ne,
    SimpleWeightedDiGraph{vertextype,_weighttype}();
    default_edge_weight=(x, y) -> rand(),
)

_nv = nv(layer_mdg) + nv(layer_vodg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_mdg_vodg = Interlayer(
    layer_mdg,
    layer_vodg,
    _ne,
    MetaDiGraph{vertextype,_weighttype}();
    default_edge_metadata=(x, y) -> (mymetadata=rand(),),
    transfer_vertex_metadata=true,
);

_nv = nv(layer_vodg) + nv(layer_vdg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_vodg_vdg = Interlayer(
    layer_vodg,
    layer_vdg,
    _ne,
    ValOutDiGraph(
        SimpleDiGraph{vertextype}();
        edgeval_types=(from_to=String,),
        edgeval_init=(s, d) -> (from_to="from_$(s)_to_$(d)",),
    );
    default_edge_metadata=(x, y) -> (from_to="from_$(x)_to_$(y)",),
);

_nv = nv(layer_sdg) + nv(layer_mdg)
_ne = rand(_nv:((_nv * (_nv - 1)) ÷ 2))
interlayer_sdg_mdg = Interlayer(
    layer_sdg,
    layer_mdg,
    _ne,
    ValDiGraph(
        SimpleDiGraph{vertextype}();
        edgeval_types=(from_to=String,),
        edgeval_init=(s, d) -> (from_to="from_$(s)_to_$(d)",),
    );
    default_edge_metadata=(x, y) -> (from_to="from_$(x)_to_$(y)",),
);

interlayer_multiplex_sdg_vodg = multiplex_interlayer(
    layer_sdg, layer_vodg, SimpleDiGraph{vertextype}()
);

interlayer_empty_sdg_vdg = empty_interlayer(
    layer_sdg,
    layer_vdg,
    SimpleWeightedDiGraph{vertextype,_weighttype}();
    default_edge_weight=(s, d) -> rand(),
)

all_interlayers = [
    interlayer_sg_swg,
    interlayer_swg_mg,
    interlayer_mg_vg,
    interlayer_multiplex_sg_mg,
    interlayer_empty_sg_vg,
    interlayer_sdg_swdg,
    interlayer_swdg_mdg,
    interlayer_mdg_vodg,
    interlayer_vodg_vdg,
    interlayer_sdg_mdg,
    interlayer_multiplex_sdg_vodg,
    interlayer_empty_sdg_vdg,
]

@testset verbose = true "MultilayerGraphs" begin
    @testset "layer" begin
        include("layer.jl")
    end
    @testset "interlayer" begin
        include("interlayer.jl")
    end
    @testset "abstractmultilayerugraph" begin
        include("abstractmultilayerugraph.jl")
    end
    @testset "abstractmultilayerdigraph" begin
        include("abstractmultilayerdigraph.jl")
    end
    @testset "agents_jl_integration" begin
        include("agents_jl_integration.jl")
    end
    @testset "utilities" begin
        include("utilities.jl")
    end
end
