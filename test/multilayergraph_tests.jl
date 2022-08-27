# Specify layers
layers_u = [Layer(n_nodes, :layer_1, SimpleGraph{Int64}, rand(min_edges:max_edges); U = Float64), 
          Layer(n_nodes, :layer_2, SimpleWeightedGraph{Int64}, rand(min_edges:max_edges); U = Float64),
          Layer(n_nodes, :layer_3, MetaGraph{Int64, Float64}, rand(min_edges:max_edges); U = Float64),
          Layer(:layer_4, ValGraph( SimpleGraph{Int64}(5, rand(min_edges:max_edges)), edgeval_types=(Int64, ), edgeval_init=(s, d) -> (s + d, ), vertexval_types=(String, ), vertexval_init=undef); U = Float64)
]

const num_layers_u = length(layers_u)

# Specify interlayers
interlayers_u = [
    Interlayer(
        n_nodes,
        :interlayer_1_2,
        :layer_1,
        :layer_2,
        SimpleGraph{Int64},
        rand(min_edges:max_edges);
        U=Float64,
    ),
    Interlayer(
        n_nodes,
        :interlayer_1_3,
        :layer_1,
        :layer_3,
        SimpleWeightedGraph{Int64,Float64},
        rand(min_edges:max_edges),
    ),
]

# Test instantiation. This also tests add_layer! and specify_interlayer!
multilayergraph = MultilayerGraph(layers_u, interlayers_u)

# Test random multilayer
random_multilayergraph = MultilayerGraph(
    num_layers_u,
    n_nodes,
    min_edges,
    max_edges,
    [SimpleGraph{Int64}, SimpleWeightedGraph{Int64,Float64}],
)

# Test getproperty and getters
multilayergraph.interlayers
multilayergraph.interlayer_1_2
multilayergraph.layer_1
multilayergraph.layers
multilayergraph.graphs
get_layer(multilayergraph, :layer_1)
get_interlayer(multilayergraph, :layer_1, :layer_2)
get_subgraph(multilayergraph, :layer_1, :layer_2)

# Test Graphs.jl extensions
is_directed(multilayergraph)
@test !is_directed(multilayergraph)

edges(multilayergraph)
@inferred edges(multilayergraph)

eltype(multilayergraph)
@inferred eltype(multilayergraph)

edgetype(multilayergraph)
@inferred eltype(multilayergraph)

has_edge(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(4, :layer_2))
@inferred has_edge(
    multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(4, :layer_2)
)

rem_edge!(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2))
@test !has_edge(
    multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2)
)

rem_edge!(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3))
@test !has_edge(
    multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3)
)

rem_edge!(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1))
@test !has_edge(
    multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1)
)

rem_edge!(multilayergraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2))
@test !has_edge(
    multilayergraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2)
)

add_edge!(multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2))
@test has_edge(
    multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2)
)
@test multilayergraph.adjacency_tensor[1, 2, 1, 2] == 1.0

add_edge!(
    multilayergraph,
    MultilayerEdge(MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3), 3.14),
)
@test has_edge(
    multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3)
)
@test multilayergraph.adjacency_tensor[1, 2, 1, 3] == 3.14

add_edge!(
    multilayergraph,
    MultilayerEdge(MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1)),
)
@test has_edge(
    multilayergraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1)
)
@test multilayergraph.adjacency_tensor[1, 2, 1, 1] == 1.0

add_edge!(
    multilayergraph,
    MultilayerEdge(MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2), 3.14),
)
@test has_edge(
    multilayergraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2)
)
@test multilayergraph.adjacency_tensor[1, 2, 2, 2] == 3.14

has_vertex(multilayergraph, MultilayerVertex(1, :layer_1))
@inferred has_vertex(multilayergraph, MultilayerVertex(1, :layer_1))

inneighbors(multilayergraph, MultilayerVertex(1, :layer_1))
@inferred inneighbors(multilayergraph, MultilayerVertex(1, :layer_1))

ne(multilayergraph)
@inferred ne(multilayergraph)

nv(multilayergraph)
@inferred nv(multilayergraph)

outneighbors(multilayergraph, MultilayerVertex(1, :layer_1))
@inferred outneighbors(multilayergraph, MultilayerVertex(1, :layer_1))

vertices(multilayergraph)
@inferred vertices(multilayergraph)

indegree(multilayergraph)
@inferred indegree(multilayergraph)

outdegree(multilayergraph)
@inferred outdegree(multilayergraph)

# test degree centrality 
degree(multilayergraph)
@inferred degree(multilayergraph)

mean_degree(multilayergraph)
@inferred mean_degree(multilayergraph)

degree_second_moment(multilayergraph)
@inferred degree_second_moment(multilayergraph)

degree_variance(multilayergraph)
@inferred degree_variance(multilayergraph)

#
randoms_u = [
    MultilayerGraph(
        num_layers_u,
        n_nodes,
        min_edges,
        max_edges,
        [SimpleGraph{Int64}, SimpleWeightedGraph{Int64,Float64}, MetaGraph{Int64, Float64}],
    ) for i in 1:4
]
@test_broken multilayer_global_clustering_coefficient.(randoms_u) .==
    global_clustering_coefficient.(randoms_u)
overlays_randoms_u = get_overlay_monoplex_graph.(randoms_u)
@test_broken global_clustering_coefficient.(overlays_randoms_u) .==
    overlay_clustering_coefficient.(randoms_u)

multilayer_weighted_global_clustering_coefficient.(randoms_u, Ref([1/3,1/3,1/3])) .≈
multilayer_global_clustering_coefficient.(randoms_u)

eig_centr_u, errs_u = eigenvector_centrality(randoms_u[1]; norm="n", tol=1e-3)

modularity.(
    randoms_u,
    Ref(rand([1, 2, 3, 4], length(nodes(randoms_u[1])), length(randoms_u[1].layers))),
)

von_neumann_entropy.(randoms_u)

get_graph_of_layers.(randoms_u)

# Test that, given a 1-dimensional multilayergraph, we obtain the same metrics as we would by using Graphs.jl's utilities on the one and only layer

## unweighted case
layer_graph = SimpleGraph{Int64}(n_nodes, rand(min_edges:max_edges))
monolayergraph = MultilayerGraph([Layer(:layer_1, layer_graph)])

@test length(edges(monolayergraph)) == length(edges(layer_graph))

@test eltype(monolayergraph) == eltype(layer_graph)

@test ne(monolayergraph) == ne(layer_graph)

@test length(nodes(monolayergraph)) == nv(layer_graph)

@test length(vertices(monolayergraph)) .== length(vertices(monolayergraph))

@test all(indegree(monolayergraph) .== indegree(layer_graph))

@test all(outdegree(monolayergraph) .== outdegree(layer_graph))

@test all(degree(monolayergraph) .== degree(layer_graph))

@test_broken vec(eigenvector_centrality(monolayergraph; norm="n", tol=1e-3)[1]) ==
    eigenvector_centrality(layer_graph)

eigenvector_centrality(monolayergraph; norm="n", tol=1e-3)[1] ./
eigenvector_centrality(layer_graph)

tests = Bool[]
for i in 1:5
    clustering = rand(
        [1, 2, 3], length(nodes(monolayergraph)), length(monolayergraph.layers)
    )
    push!(
        tests,
        modularity(monolayergraph, clustering) == modularity(layer_graph, vec(clustering)),
    )
end
@test_broken all(tests)

layer_graph_vertices = vertices(layer_graph)
for (vert_1, vert_2) in Iterators.product(layer_graph_vertices, layer_graph_vertices)
    @test has_edge(
        monolayergraph,
        MultilayerVertex(vert_1, :layer_1),
        MultilayerVertex(vert_2, :layer_1),
    ) == has_edge(layer_graph, vert_1, vert_2)
end

for vertex in vertices(layer_graph)
    @test has_vertex(monolayergraph, MultilayerVertex(vertex, :layer_1)) ==
        has_vertex(layer_graph, vertex)

    @test all(
        [
            mv_neighbor.node for
            mv_neighbor in inneighbors(monolayergraph, MultilayerVertex(vertex, :layer_1))
        ] .== inneighbors(layer_graph, vertex),
    )

    @test all(
        [
            mv_neighbor.node for
            mv_neighbor in outneighbors(monolayergraph, MultilayerVertex(vertex, :layer_1))
        ] .== outneighbors(layer_graph, vertex),
    )
end

## weighted case
layer_w_graph = SimpleWeightedGraph{Int64}(n_nodes, rand(min_edges:max_edges))
monolayerweightedgraph = MultilayerGraph([Layer(:layer_1, layer_w_graph)])

@test length(edges(monolayerweightedgraph)) == length(edges(layer_w_graph))

@test eltype(monolayerweightedgraph) == eltype(layer_w_graph)

@test ne(monolayerweightedgraph) == ne(layer_w_graph)

@test length(nodes(monolayerweightedgraph)) == nv(layer_w_graph)

@test length(vertices(monolayerweightedgraph)) .== length(vertices(monolayerweightedgraph))

@test all(indegree(monolayerweightedgraph) .== indegree(layer_w_graph))

@test all(outdegree(monolayerweightedgraph) .== outdegree(layer_w_graph))

@test all(degree(monolayerweightedgraph) .== degree(layer_w_graph))

@test_broken vec(eigenvector_centrality(monolayerweightedgraph; norm="n", tol=1e-3)[1]) ==
    eigenvector_centrality(layer_w_graph)

tests = Bool[]
for i in 1:5
    clustering = rand(
        [1, 2, 3],
        length(nodes(monolayerweightedgraph)),
        length(monolayerweightedgraph.layers),
    )
    push!(
        tests,
        modularity(monolayerweightedgraph, clustering) ==
        modularity(layer_graph, vec(clustering)),
    )
end
@test_broken all(tests)

layer_w_graph_vertices = vertices(layer_w_graph)
for (vert_1, vert_2) in Iterators.product(layer_w_graph_vertices, layer_w_graph_vertices)
    @test has_edge(
        monolayerweightedgraph,
        MultilayerVertex(vert_1, :layer_1),
        MultilayerVertex(vert_2, :layer_1),
    ) == has_edge(layer_w_graph, vert_1, vert_2)
end

for vertex in vertices(layer_w_graph)
    @test has_vertex(monolayerweightedgraph, MultilayerVertex(vertex, :layer_1)) ==
        has_vertex(layer_w_graph, vertex)

    @test all(
        [
            mv_neighbor.node for mv_neighbor in
            inneighbors(monolayerweightedgraph, MultilayerVertex(vertex, :layer_1))
        ] .== inneighbors(layer_w_graph, vertex),
    )

    @test all(
        [
            mv_neighbor.node for mv_neighbor in
            outneighbors(monolayerweightedgraph, MultilayerVertex(vertex, :layer_1))
        ] .== outneighbors(layer_w_graph, vertex),
    )
end




# test multiplex graph

multiplexgraph = MultiplexGraph(layers_u)

multiplexgraph_random = MultiplexGraph( num_layers_u, n_nodes, min_edges, max_edges, [SimpleGraph{Int64}, SimpleWeightedGraph{Int64, Float64}, MetaGraph{Int64, Float64}])

add_edge!(multiplexgraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1) )
@test multiplexgraph.adjacency_tensor[1,2,1,1] == 1.0

@test_throws ErrorException add_edge!(multiplexgraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_1), 3.14 )


add_edge!(multiplexgraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2), 3.14 )
@test multiplexgraph.adjacency_tensor[1,2,2,2] == 3.14


@test rem_edge!(multiplexgraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2))
@test multiplexgraph.adjacency_tensor[1,2,2,2] == 0.0


get_graph_of_layers(multiplexgraph)
