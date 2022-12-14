all_layers_d      = [layer for layer in all_layers if is_directed(layer)]
all_interlayers_d = [interlayer for interlayer in all_interlayers if is_directed(interlayer)]
multilayerdigraph = MultilayerDiGraph(all_layers_d, all_interlayers_d)
layers_to_be_emptied =  deepcopy([layer for layer in all_layers_d if !(layer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph)])
layers_names_to_be_emptied = name.(layers_to_be_emptied)
interlayers_to_be_emptied =  deepcopy([interlayer for interlayer in all_interlayers_d if all(in.(interlayer.layers_names, Ref(layers_names_to_be_emptied))) && !(interlayer.graph isa SimpleWeightedGraphs.AbstractSimpleWeightedGraph) ])

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

# Instantiate configuration-model multilayerdigraph
# configuration_multilayerDigraph = MultilayerDiGraph(layers_to_be_emptied, interlayers_to_be_emptied, truncated(Normal(10), 0.0, 20.0), truncated(Normal(11), 0.0, 22.0));

# Test get_interlayer
for (layer_1, layer_2) in Iterators.product(multilayerdigraph.layers_names, multilayerdigraph.layers_names)
    if layer_1 != layer_2
        interlayer = get_interlayer(multilayerdigraph,  layer_1, layer_2)
        @test nv(interlayer) >= 0
    end
end

# Test add_layer! and rem_layer!
const n_vertices_missing = rand(min_vertices:max_vertices)
const n_edges_missing_d =  rand(n_vertices_missing:(n_vertices_missing*(n_vertices_missing-1))÷ 2)
missing_layer_d = Layer(:missing_layer_d, 
sample(mvs_metadata, n_vertices_missing, replace = false), 
n_edges_missing_d,
MultilayerGraphs.ValDiGraph(SimpleDiGraph{vertextype}();edgeval_types=(Float64, String, ),
                           edgeval_init=(s, d) -> (s+d, "missing vertex $(s+d)"),
                           vertexval_types=(String,),
                           vertexval_init=v -> ("$(v^2)",),),
_weighttype;
default_edge_metadata = (src,dst) -> (rand(), "missing edge from_$(src)_to_$(dst)",)
) # SimpleGraph{vertextype}()

@test !has_layer(multilayerdigraph, missing_layer_d.name)
@test add_layer!(multilayerdigraph, missing_layer_d)
@test has_layer(multilayerdigraph, missing_layer_d.name)
@test rem_layer!(multilayerdigraph, missing_layer_d.name)
@test !has_layer(multilayerdigraph, missing_layer_d.name)

# Test nodes
@inferred(nodes(multilayerdigraph))

@inferred(nn(multilayerdigraph))

for node in nodes(multilayerdigraph)
    @test has_node(multilayerdigraph, node)
end

# Test MultilayerGraphs.add_node! and MultilayerGraphs.rem_node!
new_node = Node("new_node")
nv_prev = nv(multilayerdigraph)
ne_prev = ne(multilayerdigraph)
@test !has_node(multilayerdigraph, new_node)
@test MultilayerGraphs.add_node!(multilayerdigraph, new_node)
@test has_node(multilayerdigraph, new_node)
@test MultilayerGraphs.rem_node!(multilayerdigraph, new_node)
@test !has_node(multilayerdigraph, new_node)
# Test that nothing changed
@test nv_prev == nv(multilayerdigraph)
@test ne_prev == ne(multilayerdigraph)

# Test vertices
@test eltype(multilayerdigraph) == Int64
nv(multilayerdigraph)
@test length(multilayerdigraph.fadjlist) == length(vertices(multilayerdigraph)) # nv_withmissing(multilayerdigraph)

# Test that all multilayer vertices are present
for mv in vcat(mv_vertices.(all_layers_d)...)
    @test has_vertex(multilayerdigraph, mv)
end

for mv in mv_vertices(multilayerdigraph)
        mv_inneighbors(multilayerdigraph, mv)
        mv_outneighbors(multilayerdigraph, mv)
end

# Test add_vertex! and rem_vertex!

# Test edges
ne(multilayerdigraph)

## Test that all edges are present
for edge in vcat(collect.(edges.(all_layers_d))..., collect.(edges.(all_interlayers_d))... )
    @test has_edge(multilayerdigraph, edge)
end

@test MultilayerGraphs.weighttype(multilayerdigraph) == Float64
@test edgetype(multilayerdigraph) == MultilayerEdge{Float64}

## Test set_weight!
_, rand_mv_1_weight, rand_mv_2_weight = _get_srcmv_dstmv_layer(layer_swdg)
_weight = 3.14
@debug "" rand_mv_1_weight rand_mv_2_weight mv_vertices(multilayerdigraph)
@test !has_edge(multilayerdigraph, rand_mv_1_weight, rand_mv_2_weight)
@test add_edge!(multilayerdigraph, rand_mv_1_weight, rand_mv_2_weight, weight = _weight)
@test has_edge(multilayerdigraph, rand_mv_1_weight, rand_mv_2_weight)
wgt = weight_tensor(multilayerdigraph)
@test wgt[rand_mv_1_weight, rand_mv_2_weight] ==  get_weight(multilayerdigraph, rand_mv_1_weight, rand_mv_2_weight)  == _weight
@test set_weight!(multilayerdigraph , rand_mv_1_weight, rand_mv_2_weight, _weight + 1)
wgt = weight_tensor(multilayerdigraph)
@test wgt[rand_mv_1_weight, rand_mv_2_weight] == get_weight(multilayerdigraph, rand_mv_1_weight, rand_mv_2_weight) == _weight + 1

## Test set_metadata!
_, rand_mv_1_meta, rand_mv_2_meta = _get_srcmv_dstmv_layer(layer_mdg)
### On vertices
@test set_metadata!(multilayerdigraph, rand_mv_1_meta, (meta = "new_metadata",))
@test get_metadata(multilayerdigraph, rand_mv_1_meta).meta == "new_metadata"
## On edges
@test !has_edge(multilayerdigraph, rand_mv_1_meta, rand_mv_2_meta)
@test add_edge!(multilayerdigraph, rand_mv_1_meta, rand_mv_2_meta, metadata =  (meta = "hello",))
@test has_edge(multilayerdigraph, rand_mv_1_meta, rand_mv_2_meta)
mt = metadata_tensor(multilayerdigraph)
@test mt[rand_mv_1_meta, rand_mv_2_meta].meta == get_metadata(multilayerdigraph, rand_mv_1_meta, rand_mv_2_meta).meta == "hello"
_metadata = (meta = "bye",)
@test set_metadata!(multilayerdigraph , rand_mv_1_meta, rand_mv_2_meta, _metadata)
mt = metadata_tensor(multilayerdigraph)
@test mt[rand_mv_1_meta, rand_mv_2_meta].meta == get_metadata(multilayerdigraph, rand_mv_1_meta, rand_mv_2_meta).meta == "bye"

# Test Graphs.jl extra overrides
@test all(indegree(multilayerdigraph) .+ outdegree(multilayerdigraph) .== degree(multilayerdigraph))
@inferred(mean_degree(multilayerdigraph))
@inferred(degree_second_moment(multilayerdigraph))
@inferred(degree_variance(multilayerdigraph))

# Test multilayer-specific methods
@test all(MultilayerGraphs.get_supra_weight_matrix_from_weight_tensor(weight_tensor(multilayerdigraph).array) .== supra_weight_matrix(multilayerdigraph).array)
@test all(MultilayerGraphs.get_weight_tensor_from_supra_weight_matrix(multilayerdigraph, supra_weight_matrix(multilayerdigraph).array) .==  weight_tensor(multilayerdigraph).array)
@test_broken multilayer_global_clustering_coefficient(multilayerdigraph) .== global_clustering_coefficient(multilayerdigraph)

overlaygraph = MultilayerGraphs.get_overlay_monoplex_graph(multilayerdigraph)
@test_broken global_clustering_coefficient(overlaygraph) .== overlay_clustering_coefficient(multilayerdigraph)

@test multilayer_weighted_global_clustering_coefficient(multilayerdigraph, [1 / 3, 1 / 3, 1 / 3]) .≈ multilayer_global_clustering_coefficient(multilayerdigraph)

eig_centr_u, errs_u = eigenvector_centrality(multilayerdigraph; norm="n", tol=1e-3)

modularity(multilayerdigraph, rand([1, 2, 3, 4], length(nodes(multilayerdigraph)), length(multilayerdigraph.layers)), )

wgt = weight_tensor(multilayerdigraph)
sam = supra_weight_matrix(multilayerdigraph)
for edge in collect(edges(multilayerdigraph.layer_swdg))
    @debug "" src(edge) dst(edge) edge
    @test wgt[src(edge),dst(edge) ] == MultilayerGraphs.weight(edge)
    @test sam[src(edge),dst(edge)] == MultilayerGraphs.weight(edge)
end

# Test that, given a 1-dimensional multilayerdigraph, we obtain the same metrics as we would by using Graphs.jl utilities on the one and only layer
## unweighted and weighted case
for layer in all_layers_d
    if !(layer.graph isa SimpleValueGraphs.AbstractValGraph)
        monolayergraph = MultilayerDiGraph([layer])
        @test length(edges(monolayergraph)) == length(edges(layer.graph))
        @test eltype(monolayergraph) == eltype(layer.graph)
        @test ne(monolayergraph) == ne(layer.graph)
        @test length(nodes(monolayergraph)) == nv(layer.graph)
        @test nv(monolayergraph) .== nv(layer.graph)
        @test all(inneighbors.(Ref(monolayergraph), vertices(monolayergraph)) .== inneighbors.(Ref(layer.graph), vertices(layer.graph)))
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
                modularity(monolayergraph, clustering) == modularity(layer.graph, vec(clustering)),
            )
        end
        @test_broken all(tests)
        for edge in edges(layer)
            @test has_edge(monolayergraph, edge)
        end
    end
end