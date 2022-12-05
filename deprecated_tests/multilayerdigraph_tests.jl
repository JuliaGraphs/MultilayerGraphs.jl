#= layers_d = [
    Layer(:layer_1, get_SimpleDiGraph(); U=Float64),
    Layer(:layer_2, get_SimpleWeightedDiGraph(); U=Float64),
    Layer(:layer_3, get_SimpleWeightedDiGraph(); U=Float64),
] =#

layers_d = [
    Layer(n_nodes, :layer_1, SimpleDiGraph{Int64}, rand(min_edges:max_edges); U=Float64),
    Layer(
        n_nodes,
        :layer_2,
        SimpleWeightedDiGraph{Int64},
        rand(min_edges:max_edges);
        U=Float64,
    ),
    Layer(
        n_nodes, :layer_3, MetaDiGraph{Int64,Float64}, rand(min_edges:max_edges); U=Float64
    ),
    Layer(
        :layer_4,
        ValOutDiGraph(
            SimpleDiGraph{Int64}(5, rand(min_edges:max_edges));
            edgeval_types=(Int64,),
            edgeval_init=(s, d) -> (s + d,),
            vertexval_types=(String,),
            vertexval_init=undef,
        );
        U=Float64,
    ),
    Layer(
        :layer_5,
        ValDiGraph(
            SimpleDiGraph{Int64}(5, rand(min_edges:max_edges));
            edgeval_types=(Int64,),
            edgeval_init=(s, d) -> (s + d,),
            vertexval_types=(String,),
            vertexval_init=undef,
        );
        U=Float64,
    ),
]

const num_layers_d = length(layers_d)

# Specify interlayers
interlayers_d = [
    Interlayer(
        n_nodes,
        :myinterlayer_1_2,
        :layer_1,
        :layer_2,
        SimpleDiGraph{Int64},
        rand(min_edges:max_edges);
        U=Float64,
    ),
    Interlayer(
        n_nodes,
        :myinterlayer_1_3,
        :layer_1,
        :layer_3,
        SimpleWeightedDiGraph{Int64,Float64},
        rand(min_edges:max_edges),
    ),
]

# Test instantiation. This also tests add_layer! and specify_interlayer!
multilayerdigraph = MultilayerDiGraph(layers_d, interlayers_d)

# Test random multilayer
random_multilayerdigraph = MultilayerDiGraph(
    num_layers_d,
    n_nodes,
    min_edges,
    max_edges,
    [
        SimpleDiGraph{Int64},
        SimpleWeightedDiGraph{Int64,Float64},
        MetaDiGraph{Int64,Float64},
    ],
)

# Test getproperty and getters
random_multilayerdigraph.interlayers
random_multilayerdigraph.interlayer_1_2
random_multilayerdigraph.layer_1
random_multilayerdigraph.layers
random_multilayerdigraph.graphs
get_layer(multilayerdigraph, :layer_1)
get_interlayer(multilayerdigraph, :layer_1, :layer_2)
get_subgraph(multilayerdigraph, :layer_1, :layer_2)

# Test Graphs.jl extensions
is_directed(multilayerdigraph)
@test is_directed(multilayerdigraph)

edges(multilayerdigraph)
@inferred edges(multilayerdigraph)

eltype(multilayerdigraph)
@inferred eltype(multilayerdigraph)

edgetype(multilayerdigraph)
@inferred eltype(multilayerdigraph)

has_edge(multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(4, :layer_2))
@inferred has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(4, :layer_2)
)

has_vertex(multilayerdigraph, MultilayerVertex(1, :layer_1))
@inferred has_vertex(multilayerdigraph, MultilayerVertex(1, :layer_1))

rem_edge!(multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2))
@test !has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2)
)

rem_edge!(multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3))
@test !has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3)
)

rem_edge!(multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1))
@test !has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1)
)

rem_edge!(multilayerdigraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2))
@test !has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2)
)

add_edge!(multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2))
@test has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_2)
)
@test multilayerdigraph.array[1, 2, 1, 2] == 1.0

add_edge!(
    multilayerdigraph,
    MultilayerEdge(MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3), 3.14),
)
@test has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_3)
)
@test multilayerdigraph.array[1, 2, 1, 3] == 3.14

add_edge!(
    multilayerdigraph,
    MultilayerEdge(MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1)),
)
@test has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1)
)
@test multilayerdigraph.array[1, 2, 1, 1] == 1.0

add_edge!(
    multilayerdigraph,
    MultilayerEdge(MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2), 3.14),
)
@test has_edge(
    multilayerdigraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2)
)
@test multilayerdigraph.array[1, 2, 2, 2] == 3.14

inneighbors(multilayerdigraph, MultilayerVertex(1, :layer_1))
@inferred inneighbors(multilayerdigraph, MultilayerVertex(1, :layer_1))

ne(multilayerdigraph)
@inferred ne(multilayerdigraph)

nv(multilayerdigraph)
@inferred nv(multilayerdigraph)

outneighbors(multilayerdigraph, MultilayerVertex(1, :layer_1))
@inferred outneighbors(multilayerdigraph, MultilayerVertex(1, :layer_1))

vertices(multilayerdigraph)
@inferred vertices(multilayerdigraph)

indegree(multilayerdigraph)
@inferred indegree(multilayerdigraph)

outdegree(multilayerdigraph)
@inferred outdegree(multilayerdigraph)

# test degree centrality 
degree(multilayerdigraph)
@inferred degree(multilayerdigraph)

mean_degree(multilayerdigraph)
@inferred mean_degree(multilayerdigraph)

degree_second_moment(multilayerdigraph)
@inferred degree_second_moment(multilayerdigraph)

degree_variance(multilayerdigraph)
@inferred degree_variance(multilayerdigraph)

#
randoms_d = [
    MultilayerDiGraph(
        4,
        n_nodes,
        min_edges,
        max_edges,
        [SimpleDiGraph{Int64}, SimpleWeightedDiGraph{Int64,Float64}],
    ) for i in 1:4
]
@test_broken multilayer_global_clustering_coefficient.(randoms_du) .==
    global_clustering_coefficient.(randoms_du)
overlays_randoms_du = get_overlay_monoplex_graph.(randoms_d)
@test_broken global_clustering_coefficient.(overlays_randoms_du) .==
    overlay_clustering_coefficient.(randoms_du)

multilayer_weighted_global_clustering_coefficient.(randoms_d, Ref([1 / 3, 1 / 3, 1 / 3])) .≈
multilayer_global_clustering_coefficient.(randoms_d)

eig_centr_d, errs_d = eigenvector_centrality(randoms_d[1]; norm="n", tol=1e-3)

modularity.(
    randoms_d,
    Ref(rand([1, 2, 3, 4], length(nodes(randoms_d[1])), length(randoms_d[1].layers))),
)

get_graph_of_layers.(randoms_d)

# Test that, given a 1-dimensional multilayerdigraph, we obtain the same metrics as we would by using Graphs.jl's utilities on the one and only layer

## unweighted case
layer_graph = SimpleDiGraph{Int64}(n_nodes, rand(min_edges:max_edges))
monolayerdigraph = MultilayerDiGraph([Layer(:layer_1, layer_graph)])

@test length(edges(monolayerdigraph)) == length(edges(layer_graph))

@test eltype(monolayerdigraph) == eltype(layer_graph)

@test ne(monolayerdigraph) == ne(layer_graph)

@test length(nodes(monolayerdigraph)) == nv(layer_graph)

@test length(vertices(monolayerdigraph)) .== length(vertices(monolayerdigraph))

@test all(indegree(monolayerdigraph) .== indegree(layer_graph))

@test all(outdegree(monolayerdigraph) .== outdegree(layer_graph))

@test all(degree(monolayerdigraph) .== degree(layer_graph))

@test_broken vec(eigenvector_centrality(monolayerdigraph; norm="n", tol=1e-3)[1]) ==
    eigenvector_centrality(layer_graph)

tests = Bool[]
for i in 1:5
    clustering = rand(
        [1, 2, 3], length(nodes(monolayerdigraph)), length(monolayerdigraph.layers)
    )
    push!(
        tests,
        modularity(monolayerdigraph, clustering) ==
        modularity(layer_graph, vec(clustering)),
    )
end
@test_broken all(tests)

layer_graph_vertices = vertices(layer_graph)
for (vert_1, vert_2) in Iterators.product(layer_graph_vertices, layer_graph_vertices)
    @test has_edge(
        monolayerdigraph,
        MultilayerVertex(vert_1, :layer_1),
        MultilayerVertex(vert_2, :layer_1),
    ) == has_edge(layer_graph, vert_1, vert_2)
end

for vertex in vertices(layer_graph)
    @test has_vertex(monolayerdigraph, MultilayerVertex(vertex, :layer_1)) ==
        has_vertex(layer_graph, vertex)

    @test all(
        [
            mv_neighbor.layer_vertex for
            mv_neighbor in inneighbors(monolayerdigraph, MultilayerVertex(vertex, :layer_1))
        ] .== inneighbors(layer_graph, vertex),
    )

    @test all(
        [
            mv_neighbor.layer_vertex for mv_neighbor in
            outneighbors(monolayerdigraph, MultilayerVertex(vertex, :layer_1))
        ] .== outneighbors(layer_graph, vertex),
    )
end

## weighted case
layer_w_graph = SimpleWeightedDiGraph{Int64}(n_nodes, rand(min_edges:max_edges))
monolayerweighteddigraph = MultilayerDiGraph([Layer(:layer_1, layer_w_graph)])

@test length(edges(monolayerweighteddigraph)) == length(edges(layer_w_graph))

@test eltype(monolayerweighteddigraph) == eltype(layer_w_graph)

@test ne(monolayerweighteddigraph) == ne(layer_w_graph)

@test length(nodes(monolayerweighteddigraph)) == nv(layer_w_graph)

@test length(vertices(monolayerweighteddigraph)) .==
    length(vertices(monolayerweighteddigraph))

@test all(indegree(monolayerweighteddigraph) .== indegree(layer_w_graph))

@test all(outdegree(monolayerweighteddigraph) .== outdegree(layer_w_graph))

@test all(degree(monolayerweighteddigraph) .== degree(layer_w_graph))

@test_broken vec(eigenvector_centrality(monolayerweighteddigraph; norm="n", tol=1e-3)[1]) ==
    eigenvector_centrality(layer_w_graph)

tests = Bool[]
for i in 1:5
    clustering = rand(
        [1, 2, 3],
        length(nodes(monolayerweighteddigraph)),
        length(monolayerweighteddigraph.layers),
    )
    push!(
        tests,
        modularity(monolayerweighteddigraph, clustering) ==
        modularity(layer_graph, vec(clustering)),
    )
end
@test_broken all(tests)

layer_w_graph_vertices = vertices(layer_w_graph)
for (vert_1, vert_2) in Iterators.product(layer_w_graph_vertices, layer_w_graph_vertices)
    @test has_edge(
        monolayerweighteddigraph,
        MultilayerVertex(vert_1, :layer_1),
        MultilayerVertex(vert_2, :layer_1),
    ) == has_edge(layer_w_graph, vert_1, vert_2)
end

for vertex in vertices(layer_w_graph)
    @test has_vertex(monolayerweighteddigraph, MultilayerVertex(vertex, :layer_1)) ==
        has_vertex(layer_w_graph, vertex)

    @test all(
        [
            mv_neighbor.layer_vertex for mv_neighbor in
            inneighbors(monolayerweighteddigraph, MultilayerVertex(vertex, :layer_1))
        ] .== inneighbors(layer_w_graph, vertex),
    )

    @test all(
        [
            mv_neighbor.layer_vertex for mv_neighbor in
            outneighbors(monolayerweighteddigraph, MultilayerVertex(vertex, :layer_1))
        ] .== outneighbors(layer_w_graph, vertex),
    )
end

# test multiplex graph

multiplexdigraph = MultiplexDiGraph(layers_d)

multiplexgraph_random = MultiplexDiGraph(
    num_layers_d,
    n_nodes,
    min_edges,
    max_edges,
    [
        SimpleDiGraph{Int64},
        SimpleWeightedDiGraph{Int64,Float64},
        MetaDiGraph{Int64,Float64},
    ],
)

add_edge!(multiplexdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1))
@test multiplexdigraph.array[1, 2, 1, 1] == 1.0

@test_throws ErrorException add_edge!(
    multiplexdigraph, MultilayerVertex(1, :layer_1), MultilayerVertex(2, :layer_1), 3.14
)

add_edge!(
    multiplexdigraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2), 3.14
)
@test multiplexdigraph.array[1, 2, 2, 2] == 3.14

@test rem_edge!(
    multiplexdigraph, MultilayerVertex(1, :layer_2), MultilayerVertex(2, :layer_2)
)
@test multiplexdigraph.array[1, 2, 2, 2] == 0.0

get_graph_of_layers(multiplexdigraph)
