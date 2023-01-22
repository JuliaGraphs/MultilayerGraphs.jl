
naecg = NodeAlignedEdgeColoredGraph(all_layers_u)
new_node = Node("new_node")
MultilayerGraphs.add_node!(naecg, new_node)
@test all([has_vertex(naecg, MV(new_node, name(layer))) for layer in naecg.layers])
MultilayerGraphs.rem_node!(naecg, new_node)
@test !has_node(naecg, new_node)
@test all([!has_vertex(naecg, MV(new_node, name(layer))) for layer in naecg.layers])


naecdg = NodeAlignedEdgeColoredDiGraph(all_layers_d)
new_node = Node("new_node")
MultilayerGraphs.add_node!(naecdg, new_node)
@test all([has_vertex(naecdg, MV(new_node, name(layer))) for layer in naecdg.layers])
MultilayerGraphs.rem_node!(naecdg, new_node)
@test !has_node(naecdg, new_node)
@test all([!has_vertex(naecdg, MV(new_node, name(layer))) for layer in naecdg.layers])