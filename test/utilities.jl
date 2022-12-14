# Test `cartIndexTovecIndex`
const A = reshape(rand(1:9, 1000), 10, 10, 10)
carts_idxs = CartesianIndices(A)
vec_idxs = MultilayerGraphs.cartIndexTovecIndex.(carts_idxs, Ref(size(A)))
@test all(getindex.(Ref(A), carts_idxs) .== getindex.(Ref(A), vec_idxs)) 
directed_graphs = [SimpleDiGraph(10,50) for i in 1:5]
@test all(MultilayerGraphs.isdigraphical.(indegree.(directed_graphs), outdegree.(directed_graphs)))