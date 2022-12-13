# Test `cartIndexTovecIndex`
const A = reshape(rand(1:9, 1000), 10, 10, 10)
carts_idxs = CartesianIndices(A)
vec_idxs = MultilayerGraphs.cartIndexTovecIndex.(carts_idxs, Ref(size(A)))
@test all(getindex.(Ref(A), carts_idxs) .== getindex.(Ref(A), vec_idxs) ) 
# Test directed simple graphicality
sdg = SimpleDiGraph(10, 90)
@test @inferred(isdigraphical(indegree(sdg), outdegree(sdg)))
@test !@inferred(isdigraphical([1, 1, 1], [1, 1, 0]))
@test @inferred(isdigraphical(Integer[], Integer[]))