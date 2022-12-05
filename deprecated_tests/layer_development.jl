
const min_vertices = 5
const n_vertices = 15
const min_edges  = 1
const max_edges  = n_vertices*(n_vertices-1)

const multilayer_nodes = [Node("node_$i") for i in 1:n_vertices]

const n_vertices_simple = rand(min_vertices:n_vertices)
const n_edges_simple =  rand(n_vertices_simple:(n_vertices_simple*(n_vertices_simple-1)) ÷ 2 )
layer_sg = Layer(:layer_sg, sample(multilayer_nodes, n_vertices_simple, replace = false), SimpleGraph(n_vertices_simple, n_edges_simple) )
layer_sdg = Layer(:layer_sdg, sample(multilayer_nodes, n_vertices_simple, replace = false), SimpleDiGraph(n_vertices_simple,n_edges_simple))

const n_vertices_weighted = rand(min_vertices:n_vertices)
const n_edges_weighted =  rand(n_vertices_weighted:(n_vertices_weighted*(n_vertices_weighted-1))÷ 2)
layer_swg = Layer(:layer_swg, sample(multilayer_nodes, n_vertices_weighted, replace = false), SimpleWeightedGraph(n_vertices_weighted,n_edges_weighted))
layer_swdg = Layer(:layer_swdg, sample(multilayer_nodes, n_vertices_weighted, replace = false), SimpleWeightedDiGraph(n_vertices_weighted,n_edges_weighted))

const n_vertices_meta = rand(min_vertices:n_vertices)
const n_edges_meta =  rand(n_vertices_meta:(n_vertices_meta*(n_vertices_meta-1))÷ 2)
layer_mg = Layer(:layer_mg, sample(multilayer_nodes, n_vertices_meta, replace = false), MetaGraph(n_vertices_meta,n_edges_meta))
layer_mdg = Layer(:layer_mdg, sample(multilayer_nodes, n_vertices_meta, replace = false), MetaDiGraph(n_vertices_meta, n_edges_meta))

const n_vertices_value = rand(min_vertices:n_vertices)
const n_edges_value =  rand(n_vertices_value:(n_vertices_value*(n_vertices_value-1))÷ 2)
layer_vg = Layer(:layer_vg, sample(multilayer_nodes, n_vertices_value, replace = false), MultilayerGraphs.ValGraph(n_vertices_value,n_edges_value;   edgeval_types=(Int64,String),
edgeval_init=(s, d) -> (s + d, "hi"),
vertexval_types=(String,),
vertexval_init=v -> ("$v",),))

# get_vertexval.(Ref(layer_vg.graph), vertices(layer_vg.graph), Ref(1))

layer_vodg = Layer(:layer_vodg, sample(multilayer_nodes, n_vertices_value, replace = false), ValOutDiGraph(n_vertices_value,n_edges_value;    edgeval_types=(a = Int64, b = String), edgeval_init=(s, d) -> (s + d, "hi"), vertexval_types=(String,), vertexval_init=v -> ("$v",),))
layer_vdg = Layer(:layer_vdg, sample(multilayer_nodes, n_vertices_value, replace = false), ValDiGraph(n_vertices_value,n_edges_value))



all_layers = [layer_sg, layer_sdg, layer_swg, layer_swdg, layer_mg, layer_mdg, layer_vg, layer_vodg, layer_vdg, ]

#= # test random constructor
Layer(:layer_random, n_vertices, rand(min_edges:max_edges), graph_type = ValOutDiGraph{Int64} , edgeval_types=(Int64,), edgeval_init=(s, d) -> (s + d,),vertexval_types=(String,),vertexval_init=undef) =#

edges.(all_layers)
eltype.(all_layers)
edgetype.(all_layers)
has_edge.(all_layers,  [layer.v_V_associations[1] for layer in all_layers],  [layer.v_V_associations[2] for layer in all_layers])
has_vertex.(all_layers,  Ref(1))

collect(edges(all_layers[8]))


mv_inneighbors.(all_layers,  [layer.v_V_associations[1] for layer in all_layers])
ne.(all_layers)
nv.(all_layers)
mv_outneighbors.(all_layers, [layer.v_V_associations[1] for layer in all_layers])
vertices.(all_layers)
is_directed.(all_layers)

rand_mv_1 =  rand(multilayer_vertices(layer_sg))
rand_mv_2 =  rand(multilayer_vertices(layer_sg))
layer = layer_sg
# test uniform add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test !has_edge(layer, rand_mv_1,rand_mv_2)
@test add_edge!(layer, rand_mv_1,rand_mv_2; weight = nothing, metadata = ())
@test has_edge(layer, rand_mv_1,rand_mv_2)
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
# test hybrid add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
add_edge!(layer,  rand_mv_1,rand_mv_2)
@test has_edge(layer, rand_mv_1,rand_mv_2)
# Test uniform add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1 )
@test has_vertex(layer, rand_mv_1)
# Test hybrid add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)


rand_mv_1 =  rand(multilayer_vertices(layer_sdg))
rand_mv_2 =  rand(multilayer_vertices(layer_sdg))
layer = layer_sdg
# test uniform add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test !has_edge(layer, rand_mv_1,rand_mv_2)
@test add_edge!(layer, rand_mv_1,rand_mv_2)
@test has_edge(layer, rand_mv_1,rand_mv_2)
# test hybrid add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
add_edge!(layer,  rand_mv_1,rand_mv_2)
@test has_edge(layer, rand_mv_1,rand_mv_2)
# Test uniform add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)
# Test hybrid add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)





rand_mv_1 =  rand(multilayer_vertices(layer_swg))
rand_mv_2 =  rand(multilayer_vertices(layer_swg))
layer = layer_swg
# test uniform add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test !has_edge(layer, rand_mv_1,rand_mv_2)
@test add_edge!(layer, rand_mv_1,rand_mv_2 ; weight = 3.14, metadata = ())
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test layer.graph.weights[layer.v_V_associations(rand_mv_1), layer.v_V_associations(rand_mv_2)] == layer.graph.weights[layer.v_V_associations(rand_mv_2), layer.v_V_associations(rand_mv_1)] == 3.14
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test add_edge!(layer,  rand_mv_1,rand_mv_2, 3.14)
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test layer.graph.weights[layer.v_V_associations(rand_mv_1), layer.v_V_associations(rand_mv_2)] == layer.graph.weights[layer.v_V_associations(rand_mv_2), layer.v_V_associations(rand_mv_1)] == 3.14
# Test uniform add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)
# Test hybrid add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)

layer = layer_swdg
rand_mv_1 =  rand(multilayer_vertices(layer))
rand_mv_2 =  rand(multilayer_vertices(layer))
# test uniform add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test !has_edge(layer, rand_mv_1,rand_mv_2)
@test add_edge!(layer, rand_mv_1,rand_mv_2 , weight = 3.14, metadata = ())
@test has_edge(layer, rand_mv_1,rand_mv_2)
# Why do I have to switch the vertices? 
@test layer.graph.weights[layer.v_V_associations(rand_mv_2), layer.v_V_associations(rand_mv_1)] == 3.14
# test hybrid add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test add_edge!(layer,  rand_mv_1,rand_mv_2, 3.14)
@test has_edge(layer, rand_mv_1,rand_mv_2)
# Why do I have to switch the vertices? 
@test layer.graph.weights[layer.v_V_associations(rand_mv_2), layer.v_V_associations(rand_mv_1)] == 3.14
# Test uniform add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)
# Test hybrid add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)

layer = layer_mg
rand_mv_1 =  rand(multilayer_vertices(layer))
rand_mv_2 =  rand(multilayer_vertices(layer))
# test uniform add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test !has_edge(layer, rand_mv_1,rand_mv_2)
@test add_edge!(layer, rand_mv_1,rand_mv_2 , weight = nothing,  metadata = (weight = 4, property_1 = "hello"))
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test get_prop(layer.graph, layer.v_V_associations(rand_mv_1), layer.v_V_associations(rand_mv_2), :weight  ) == 4
@test get_prop(layer.graph, layer.v_V_associations(rand_mv_1), layer.v_V_associations(rand_mv_2), :property_1  ) == "hello"
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test add_edge!(layer,  rand_mv_1,rand_mv_2)
@test has_edge(layer, rand_mv_1,rand_mv_2)
set_prop!(layer, rand_mv_1, rand_mv_2, :weight, 5)
set_prop!(layer, rand_mv_1, rand_mv_2, :property_1, "world")
@test  get_edge_metadata(layer, rand_mv_1, rand_mv_2).weight == 5 # get_prop(layer.graph, layer.v_V_associations(rand_mv_1), layer.v_V_associations(rand_mv_2), :weight  )
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2).property_1 == "world"
# Test uniform add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, MultilayerVertex(rand_mv_1.node, rand_mv_1.layer, (age = 28,)))
@test has_vertex(layer, rand_mv_1)
@test get_vertex_metadata(layer, rand_mv_1).age == 28
# Test hybrid add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1.node, Dict(:age => 28))
@test has_vertex(layer, rand_mv_1)
@test get_vertex_metadata(layer, rand_mv_1).age == 28

layer = layer_vg
rand_mv_1 =  rand(multilayer_vertices(layer))
rand_mv_2 =  rand(multilayer_vertices(layer))
# test uniform add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test !has_edge(layer, rand_mv_1,rand_mv_2)
@test add_edge!(layer, rand_mv_1,rand_mv_2 , weight = nothing, metadata = (4,"hello"))
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2)[1] == 4
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2)[2] == "hello"
@test add_edge!(layer, rand_mv_1,rand_mv_2 , weight = nothing,  metadata = (5,"world"))
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2)[1] == 5
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2)[2] == "world"
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test add_edge!(layer, rand_mv_1,rand_mv_2 , (4,"hello"))
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2)[1] == 4
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2)[2] ==  "hello"
# Test uniform add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test        has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata = (age = 28,))
#= @test has_vertex(layer, rand_mv_1)
@test get_prop(layer, rand_mv_1, :age) == 28 =#
# Test hybrid add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test_broken !has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata = (age = 28,))
#= @test has_vertex(layer, rand_mv_1)
@test get_prop(layer, rand_mv_1, :age) == 28 =#

layer = layer_vodg
rand_mv_1 =  rand(multilayer_vertices(layer))
rand_mv_2 =  rand(multilayer_vertices(layer))
# test uniform add_edge!
rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test !has_edge(layer, rand_mv_1,rand_mv_2)
@test add_edge!(layer, rand_mv_1,rand_mv_2 , weight = nothing, metadata = (a = 4, b = "hello"))
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2).a == 4
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2).b == "hello"
@test add_edge!(layer, rand_mv_1,rand_mv_2 , weight = nothing, metadata = (a = 5, b = "world"))
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2).a == 5
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2).b == "world"
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1 , rand_mv_2 )
@test add_edge!(layer, rand_mv_1,rand_mv_2 , (a = 4, b = "hello"))
@test has_edge(layer, rand_mv_1,rand_mv_2)
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2).a == 4
@test get_edge_metadata(layer, rand_mv_1, rand_mv_2).b == "hello"
# Test uniform add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test_broken !has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata = (age = 28,))
#= @test has_vertex(layer, rand_mv_1)
@test get_prop(layer, rand_mv_1, :age) == 28 =#
# Test hybrid add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test_broken !has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata = (age = 28,))
#= @test has_vertex(layer, rand_mv_1)
@test get_prop(layer, rand_mv_1, :age) == 28 =#

#= a=  sparse([20, 11, 17, 10, 19, 20, 19, 15, 6, 3, 9, 4, 6, 8, 2, 6], [2, 3, 4, 6, 6, 6, 8, 9, 10, 11, 15, 17, 19, 19, 20, 20], [0.2242017352436838, 0.07751430316809071, 0.0056804992063987925, 0.8508391801649005, 0.5428523433512523, 0.8200535617288671, 0.2666672566906638, 0.8700451363963257, 0.8508391801649005, 0.07751430316809071, 0.8700451363963257, 0.0056804992063987925, 0.5428523433512523, 0.2666672566906638, 0.2242017352436838, 0.8200535617288671], 20, 20)

b = sparse([20, 11, 17, 10, 19, 20, 19, 15, 6, 3, 9, 4, 6, 8, 2, 6], [2, 3, 4, 6, 6, 6, 8, 9, 10, 11, 15, 17, 19, 19, 20, 20], [0.2242017352436838, 0.07751430316809071, 0.0056804992063987925, 0.8508391801649005, 0.5428523433512523, 0.8200535617288671, 0.2666672566906638, 0.6191206440012805, 0.8508391801649005, 0.07751430316809071, 0.6191206440012805, 0.0056804992063987925, 0.5428523433512523, 0.2666672566906638, 0.2242017352436838, 0.8200535617288671], 20, 20)

_diffs = findall( x -> !x, a .== b)

a[_diffs]
b[_diffs]


edg_1= MultilayerEdge{MultilayerVertex, Float64}[MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_11"), :layer_mg, NamedTuple()), MultilayerVertex(Node{String}("node_10"), :layer_vg, NamedTuple()), 0.2666672566906638, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_1"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_5"), :layer_mg, NamedTuple()), 0.07556695303698568, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_9"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_3"), :layer_mg, NamedTuple()), 0.6191206440012805, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_1"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_5"), :layer_mg, NamedTuple()), 0.07751430316809071, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_4"), :layer_mg, NamedTuple()), MultilayerVertex(Node{String}("node_11"), :layer_vg, NamedTuple()), 0.2242017352436838, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_12"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_13"), :layer_mg, NamedTuple()), 0.8508391801649005, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_3"), :layer_mg, NamedTuple()), MultilayerVertex(Node{String}("node_9"), :layer_vg, NamedTuple()), 0.8700451363963257, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_10"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_13"), :layer_mg, NamedTuple()), 0.5428523433512523, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_15"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_14"), :layer_mg, NamedTuple()), 0.0056804992063987925, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_11"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_13"), :layer_mg, NamedTuple()), 0.8200535617288671, NamedTuple())]

edg_2 =  MultilayerEdge[MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_12"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_13"), :layer_mg, NamedTuple()), 0.8508391801649005, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_3"), :layer_mg, NamedTuple()), MultilayerVertex(Node{String}("node_9"), :layer_vg, NamedTuple()), 0.8700451363963257, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_11"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_13"), :layer_mg, NamedTuple()), 0.8200535617288671, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_10"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_13"), :layer_mg, NamedTuple()), 0.5428523433512523, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_15"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_14"), :layer_mg, NamedTuple()), 0.0056804992063987925, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_4"), :layer_mg, NamedTuple()), MultilayerVertex(Node{String}("node_11"), :layer_vg, NamedTuple()), 0.2242017352436838, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_1"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_5"), :layer_mg, NamedTuple()), 0.07751430316809071, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_11"), :layer_mg, NamedTuple()), MultilayerVertex(Node{String}("node_10"), :layer_vg, NamedTuple()), 0.2666672566906638, NamedTuple()), MultilayerEdge{MultilayerVertex, Float64}(MultilayerVertex(Node{String}("node_9"), :layer_vg, NamedTuple()), MultilayerVertex(Node{String}("node_3"), :layer_mg, NamedTuple()), 0.6191206440012805, NamedTuple())]

#= _compute_interlayer_graph; duplicate found: MultilayerVertex(Node{String}("node_1"), :layer_vg, NamedTuple()) and MultilayerVertex(Node{String}("node_5"), :layer_mg, NamedTuple())
_compute_interlayer_graph; duplicate found: MultilayerVertex(Node{String}("node_3"), :layer_mg, NamedTuple()) and MultilayerVertex(Node{String}("node_9"), :layer_vg, NamedTuple()) =#

union(Set(edg_1), Set(edg_2) )

wgt_1 = MultilayerGraphs.weight.(edg_1)

wgt_2 = MultilayerGraphs.weight.(edg_2)

wgt_1 .== wgt_2

setdiff(Set(wgt_1), Set(wgt_2) )


swg = SimpleWeightedGraph(10,20)
swg_dc = deepcopy(swg)

swg === swg_dc

which(==, typeof.((swg,swg_dc))) =#
# Test interlayer.
# Cosntructors already test funxtions: has_edge, has_vertex, 
nv_interlayer_sg_swg = nv(layer_sg) + nv(layer_swg)
interlayer_sg_swg = Interlayer(layer_sg, layer_swg, rand(nv_interlayer_sg_swg:(nv_interlayer_sg_swg*(nv_interlayer_sg_swg-1))÷ 2), SimpleGraph{Int64} )

nv_interlayer_sdg_swdg = nv(layer_sdg) + nv(layer_swdg)
interlayer_sdg_swdg = Interlayer(layer_sdg, layer_swdg, rand(nv_interlayer_sdg_swdg:(nv_interlayer_sdg_swdg*(nv_interlayer_sdg_swdg-1))÷ 2), SimpleDiGraph{Int64})


nv_interlayer_mg_vg = nv(layer_mg) + nv(layer_vg)
interlayer_mg_vg = Interlayer(layer_mg, layer_vg, rand(nv_interlayer_mg_vg:(nv_interlayer_mg_vg*(nv_interlayer_mg_vg-1))÷ 2), SimpleWeightedGraph{Int64, Float64}; edge_weight_function = (x,y) -> rand() ) # rand(nv_interlayer_mg_vg:(nv_interlayer_mg_vg*(nv_interlayer_mg_vg-1))÷ 2)


nv_interlayer_mdg_vodg = nv(layer_mdg) + nv(layer_vodg)
interlayer_mdg_vodg = Interlayer(layer_mdg, layer_vodg, rand(nv_interlayer_mdg_vodg:(nv_interlayer_mdg_vodg*(nv_interlayer_mdg_vodg-1))÷ 2), ValDiGraph{Int64}; edge_metadata_function = (src,dst) -> (weight = rand(),),  graph_kwargs = (edgeval_types=(weight = Float64,),
edgeval_init=(s, d) -> (s + d,) )); #


# nv_interlayer_vodg_vdg = nv(layer_mdg) + nv(layer_vodg)
interlayer_multiplex_mg_vg = multiplex_interlayer(layer_mg, layer_vg, ValOutDiGraph{Int64}; edge_metadata_function = node -> (weight = rand(),),  graph_kwargs = (edgeval_types=(weight = Float64,), edgeval_init=(s, d) -> (weight = s + d,))) # , vertexval_types=(String,), vertexval_init = v -> ("$v",) 

all_interlayers = [interlayer_sg_swg, interlayer_sdg_swdg, interlayer_mg_vg, interlayer_mdg_vodg, interlayer_multiplex_mg_vg]


@test is_multiplex_interlayer(interlayer_multiplex_mg_vg)


collect.(edges.(all_interlayers))
eltype.(all_interlayers)
edgetype.(all_interlayers)

for interlayer in all_interlayers
    @test all(has_edge.(Ref(interlayer), edges(interlayer)))
    @test all(has_vertex.(Ref(interlayer), multilayer_vertices(interlayer)))
    @test all(has_vertex.(Ref(interlayer), multilayer_vertices(interlayer)))
    mv = rand(collect(multilayer_vertices(interlayer)))
    mv_inneighbors(interlayer, mv)
    mv_outneighbors(interlayer, mv)
    adjacency_matrix(interlayer)
    MultilayerGraphs.weights(interlayer)
end

ne.(all_interlayers)
nv.(all_interlayers)


interlayer = interlayer_sg_swg
mvs = get_bare_mv.(collect(multilayer_vertices(interlayer)))
src_mv = rand(mvs)
dst_mv = get_bare_mv(rand(setdiff(Set(mvs), Set(vcat(get_bare_mv.(mv_outneighbors(interlayer, src_mv)), src_mv, get_bare_mv.(multilayer_vertices(Base.getproperty(interlayer, src_mv.layer)))) ) )  ))
missing_edge = ME(src_mv, dst_mv)
# test uniform add_edge!
@test !has_edge(interlayer,missing_edge)
@test !has_edge(interlayer,MultilayerGraphs.reverse(missing_edge))
@test add_edge!(interlayer,missing_edge)
@test has_edge(interlayer,missing_edge)
@test has_edge(interlayer,MultilayerGraphs.reverse(missing_edge))
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv , dst_mv )
@test add_edge!(interlayer,  src_mv,dst_mv)
@test has_edge(interlayer, src_mv,dst_mv)

missing_vertex = MV(Node("missing_node"), :layer_sg)
@test !has_vertex(interlayer, missing_vertex)
@test add_vertex!(interlayer, missing_vertex)
@test has_vertex(layer_sg, missing_vertex)
@test has_vertex(interlayer, missing_vertex)
@test rem_vertex!(interlayer,  missing_vertex)
@test !has_vertex(interlayer, missing_vertex)


interlayer = interlayer_sdg_swdg
mvs = get_bare_mv.(collect(multilayer_vertices(interlayer)))
src_mv = rand(mvs)
dst_mv = get_bare_mv(rand(setdiff(Set(mvs), Set(vcat(get_bare_mv.(mv_outneighbors(interlayer, src_mv)), src_mv, get_bare_mv.(multilayer_vertices(Base.getproperty(interlayer, src_mv.layer)))) ) )  ))
missing_edge = ME(src_mv, dst_mv)
# test uniform add_edge!
@test !has_edge(interlayer,missing_edge)
@test add_edge!(interlayer,missing_edge)
@test has_edge(interlayer,missing_edge)
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv , dst_mv )
@test add_edge!(interlayer,  src_mv,dst_mv)
@test has_edge(interlayer, src_mv,dst_mv)


interlayer = interlayer_mg_vg
mvs = get_bare_mv.(collect(multilayer_vertices(interlayer)))
src_mv = get_bare_mv(rand(mvs))
dst_mv = get_bare_mv(rand(setdiff(Set(mvs), Set(vcat(get_bare_mv.(mv_outneighbors(interlayer, src_mv)), src_mv, get_bare_mv.(multilayer_vertices(Base.getproperty(interlayer, src_mv.layer)))) ) )  ))
missing_edge = ME(src_mv, dst_mv, rand())
# test uniform add_edge!
@test !has_edge(interlayer,missing_edge)
@test add_edge!(interlayer,missing_edge)
@test has_edge(interlayer,missing_edge)
@test MultilayerGraphs.weights(interlayer)[interlayer.v_V_associations(src_mv),interlayer.v_V_associations(dst_mv)] == MultilayerGraphs.weights(interlayer)[interlayer.v_V_associations(dst_mv),interlayer.v_V_associations(src_mv)] != 0.0
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv , dst_mv )
@test add_edge!(interlayer,  src_mv,dst_mv, rand())
@test has_edge(interlayer, src_mv,dst_mv)
@test MultilayerGraphs.weights(interlayer)[interlayer.v_V_associations(src_mv),interlayer.v_V_associations(dst_mv)] == MultilayerGraphs.weights(interlayer)[interlayer.v_V_associations(dst_mv),interlayer.v_V_associations(src_mv)] != 0.0



interlayer = interlayer_mdg_vodg
mvs = get_bare_mv.(collect(multilayer_vertices(interlayer)))
src_mv = get_bare_mv(rand(mvs))
dst_mv = get_bare_mv(rand(setdiff(Set(mvs), Set(vcat(get_bare_mv.(mv_outneighbors(interlayer, src_mv)), src_mv, get_bare_mv.(multilayer_vertices(Base.getproperty(interlayer, src_mv.layer)))) ) )  ))
missing_edge = ME(src_mv, dst_mv, (weight = rand(),))
# test uniform add_edge!
@test !has_edge(interlayer,missing_edge)
@test add_edge!(interlayer,missing_edge)
@test has_edge(interlayer,missing_edge)
@test !isempty(get_edge_metadata(interlayer,src_mv, dst_mv )) && !isempty(get_edge_metadata(interlayer,src_mv, dst_mv ))
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv , dst_mv )
@test_throws ErrorException isempty(get_edge_metadata(interlayer,src_mv, dst_mv )) && isempty(get_edge_metadata(interlayer,src_mv, dst_mv ))
@test add_edge!(interlayer,  src_mv,dst_mv, (weight = rand(),))
@test has_edge(interlayer, src_mv,dst_mv)
@test !isempty(get_edge_metadata(interlayer,src_mv, dst_mv )) && !isempty(get_edge_metadata(interlayer,src_mv, dst_mv ))


get_symmetric_interlayer.(all_interlayers)