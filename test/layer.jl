edges.(all_layers)
eltype.(all_layers)
edgetype.(all_layers)
has_edge.(
    all_layers,
    [layer.v_V_associations[1] for layer in all_layers],
    [layer.v_V_associations[2] for layer in all_layers],
)
has_vertex.(all_layers, Ref(1))

collect(edges(all_layers[8]))

mv_inneighbors.(all_layers, [layer.v_V_associations[1] for layer in all_layers])
ne.(all_layers)
nv.(all_layers)
mv_outneighbors.(all_layers, [layer.v_V_associations[1] for layer in all_layers])
vertices.(all_layers)
is_directed.(all_layers)

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
@debug ""
#= rand_mv_1 =  rand(mv_vertices(layer_sg))
rand_mv_2 =  rand(mv_vertices(layer_sg)) =#

layer = layer_sg
_, rand_mv_1, rand_mv_2 = _get_srcmv_dstmv_layer(layer)
# test uniform add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test !has_edge(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2; weight=nothing, metadata=())
@test has_edge(layer, rand_mv_1, rand_mv_2)
rem_edge!(layer, rand_mv_1, rand_mv_2)
# test hybrid add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
add_edge!(layer, rand_mv_1, rand_mv_2)
@test has_edge(layer, rand_mv_1, rand_mv_2)
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

#= rand_mv_1 =  rand(mv_vertices(layer_sdg))
rand_mv_2 =  rand(mv_vertices(layer_sdg)) =#
layer = layer_sdg
_, rand_mv_1, rand_mv_2 = _get_srcmv_dstmv_layer(layer)
# test uniform add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test !has_edge(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2)
@test has_edge(layer, rand_mv_1, rand_mv_2)
# test hybrid add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
add_edge!(layer, rand_mv_1, rand_mv_2)
@test has_edge(layer, rand_mv_1, rand_mv_2)
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

@debug ""

#= rand_mv_1 =  rand(mv_vertices(layer_swg))
rand_mv_2 =  rand(mv_vertices(layer_swg)) =#
layer = layer_swg
_, rand_mv_1, rand_mv_2 = _get_srcmv_dstmv_layer(layer)
# test uniform add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test !has_edge(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2; weight=3.14, metadata=())
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test Graphs.weights(layer)[get_v(layer, rand_mv_1), get_v(layer, rand_mv_2) ] == Graphs.weights(layer)[get_v(layer, rand_mv_2), get_v(layer, rand_mv_1) ] == 3.14 
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, 3.14)
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test Graphs.weights(layer)[get_v(layer, rand_mv_1), get_v(layer, rand_mv_2) ] == Graphs.weights(layer)[get_v(layer, rand_mv_2), get_v(layer, rand_mv_1) ] == 3.14
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

@debug ""

layer = layer_swdg
#= rand_mv_1 =  rand(mv_vertices(layer))
rand_mv_2 =  rand(mv_vertices(layer)) =#
_, rand_mv_1, rand_mv_2 = _get_srcmv_dstmv_layer(layer)
# test uniform add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test !has_edge(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, weight=3.14, metadata=())
@test has_edge(layer, rand_mv_1, rand_mv_2)
# Why do I have to switch the vertices? 
@test Graphs.weights(layer)[get_v(layer, rand_mv_1), get_v(layer, rand_mv_2) ] == 3.14
# test hybrid add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, 3.14)
@test has_edge(layer, rand_mv_1, rand_mv_2)
# Why do I have to switch the vertices? 
@test Graphs.weights(layer)[get_v(layer, rand_mv_1), get_v(layer, rand_mv_2) ] == 3.14 
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
#= rand_mv_1 =  rand(mv_vertices(layer))
rand_mv_2 =  rand(mv_vertices(layer)) =#
_, rand_mv_1, rand_mv_2 = _get_srcmv_dstmv_layer(layer)
# test uniform add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test !has_edge(layer, rand_mv_1, rand_mv_2)
@test add_edge!(
    layer, rand_mv_1, rand_mv_2, weight=nothing, metadata=(weight=4, property_1="hello")
)
@test has_edge(layer, rand_mv_1, rand_mv_2)
#= @test get_prop(layer.graph, get_v(layer,rand_mv_1), get_v(layer, rand_mv_2), :weight  ) == 4
@test get_prop(layer.graph, get_v(layer,rand_mv_1), get_v(layer, rand_mv_2), :property_1  ) == "hello" =#
@test get_prop(layer, rand_mv_1, rand_mv_2, :weight) == 4
@test get_prop(layer, rand_mv_1, rand_mv_2, :property_1) == "hello"
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2)
@test has_edge(layer, rand_mv_1, rand_mv_2)
set_prop!(layer, rand_mv_1, rand_mv_2, :weight, 5)
set_prop!(layer, rand_mv_1, rand_mv_2, :property_1, "world")
@test get_metadata(layer, rand_mv_1, rand_mv_2).weight == 5
@test get_metadata(layer, rand_mv_1, rand_mv_2).property_1 == "world"
# Test uniform add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, MultilayerVertex(rand_mv_1.node, rand_mv_1.layer, (age=28,)))
@test has_vertex(layer, rand_mv_1)
@test get_metadata(layer, rand_mv_1).age == 28
# Test hybrid add_vertex!
@test rem_vertex!(layer, rand_mv_1)
@test !has_vertex(layer, rand_mv_1)
@test add_vertex!(layer, rand_mv_1.node, Dict(:age => 28))
@test has_vertex(layer, rand_mv_1)
@test get_metadata(layer, rand_mv_1).age == 28

layer = layer_vg
#= rand_mv_1 =  rand(mv_vertices(layer))
rand_mv_2 =  rand(mv_vertices(layer)) =#
_, rand_mv_1, rand_mv_2 = _get_srcmv_dstmv_layer(layer)
# test uniform add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test !has_edge(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, weight=nothing, metadata=(4, "hello"))
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test get_metadata(layer, rand_mv_1, rand_mv_2)[1] == 4
@test get_metadata(layer, rand_mv_1, rand_mv_2)[2] == "hello"
@test rem_edge!(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, weight=nothing, metadata=(5, "world"))
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test get_metadata(layer, rand_mv_1, rand_mv_2)[1] == 5
@test get_metadata(layer, rand_mv_1, rand_mv_2)[2] == "world"
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, (4, "hello"))
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test get_metadata(layer, rand_mv_1, rand_mv_2)[1] == 4
@test get_metadata(layer, rand_mv_1, rand_mv_2)[2] == "hello"
# Test uniform add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata=(age=28,))

# Test hybrid add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test_broken !has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata=(age=28,))

layer = layer_vodg
#= rand_mv_1 =  rand(mv_vertices(layer))
rand_mv_2 =  rand(mv_vertices(layer)) =#
_, rand_mv_1, rand_mv_2 = _get_srcmv_dstmv_layer(layer)
# test uniform add_edge!
rem_edge!(layer, rand_mv_1, rand_mv_2)
@test !has_edge(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, weight=nothing, metadata=(a=4, b="hello"))
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test get_metadata(layer, rand_mv_1, rand_mv_2).a == 4
@test get_metadata(layer, rand_mv_1, rand_mv_2).b == "hello"
@test rem_edge!(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, weight=nothing, metadata=(a=5, b="world"))
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test get_metadata(layer, rand_mv_1, rand_mv_2).a == 5
@test get_metadata(layer, rand_mv_1, rand_mv_2).b == "world"
# test hybrid add_edge!
@test rem_edge!(layer, rand_mv_1, rand_mv_2)
@test add_edge!(layer, rand_mv_1, rand_mv_2, (a=4, b="hello"))
@test has_edge(layer, rand_mv_1, rand_mv_2)
@test get_metadata(layer, rand_mv_1, rand_mv_2).a == 4
@test get_metadata(layer, rand_mv_1, rand_mv_2).b == "hello"
# Test uniform add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test_broken !has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata=(age=28,))

# Test hybrid add_vertex!
@test_throws MethodError rem_vertex!(layer, rand_mv_1)
@test_broken !has_vertex(layer, rand_mv_1)
@test_broken add_vertex!(layer, rand_mv_1; metadata=(age=28,))
