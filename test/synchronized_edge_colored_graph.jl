
secg = SynchronizedEdgeColoredGraph(all_layers_u)
new_node = Node("new_node")
MultilayerGraphs.add_node!(secg, new_node)
@test all([has_vertex(secg, MV(new_node, name(layer))) for layer in secg.layers])
MultilayerGraphs.rem_node!(secg, new_node)
@test !has_node(secg, new_node)
@test all([!has_vertex(secg, MV(new_node, name(layer))) for layer in secg.layers])


secdg = SynchronizedEdgeColoredDiGraph(all_layers_d)
new_node = Node("new_node")
MultilayerGraphs.add_node!(secdg, new_node)
@test all([has_vertex(secdg, MV(new_node, name(layer))) for layer in secdg.layers])
MultilayerGraphs.rem_node!(secdg, new_node)
@test !has_node(secdg, new_node)
@test all([!has_vertex(secdg, MV(new_node, name(layer))) for layer in secdg.layers])