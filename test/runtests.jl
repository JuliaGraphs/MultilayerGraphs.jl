using Test
using MultilayerGraphs
using Graphs, SimpleWeightedGraphs, MetaGraphs

# Test that methods work on a general multilayer graph
## Define variables and utilities
const n_nodes = 5
const min_edges = n_nodes
const max_edges = 10

get_SimpleGraph() = SimpleGraph(n_nodes, rand(min_edges:max_edges))
get_SimpleDiGraph() = SimpleDiGraph(n_nodes, rand(min_edges:max_edges))

const simpleweightedgraph_sources = 1:n_nodes
const simpleweightedgraph_destinations = rand(1:n_nodes, n_nodes)
const simpleweightedgraph_weights = rand(n_nodes)
simpleweightedgraph = SimpleWeightedGraph(simpleweightedgraph_sources,
                                          simpleweightedgraph_destinations,
                                          simpleweightedgraph_weights)
function get_SimpleWeightedGraph()
    SimpleWeightedGraph(simpleweightedgraph_sources, rand(1:n_nodes, n_nodes),
                        rand(n_nodes))
end
simpleweighteddigraph = SimpleWeightedDiGraph(simpleweightedgraph_sources,
                                              simpleweightedgraph_destinations,
                                              simpleweightedgraph_weights)
function get_SimpleWeightedDiGraph()
    SimpleWeightedDiGraph(simpleweightedgraph_sources, rand(1:n_nodes, n_nodes),
                          rand(n_nodes))
end

metadigraph = MetaDiGraph(simpleweighteddigraph)
metagraph = MetaGraph(simpleweightedgraph)

@testset "MultilayerGraphs" begin
    @testset "utilities" begin include("utilities_tests.jl") end

    @testset "multilayergraph" begin include("multilayergraph_tests.jl") end

    @testset "multilayerdigraph" begin include("multilayerdigraph_tests.jl") end
end
