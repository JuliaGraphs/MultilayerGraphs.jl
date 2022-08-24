"""
AbstractMultiplexGraph{T}
"""
abstract type AbstractMultiplexUGraph{T,U} <: AbstractMultilayerUGraph{T,U} end


# Graphs.jl's internals extra overrides

"""
    add_edge!(mg::M, src::V, dst::V, weight::U) where { T, U <: Real, M <: AbstractMultiplexUGraph{T,U}, V <: MultilayerVertex{T}}

Add edge from `src` to `dst` with weight `weight` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, src::V, dst::V, weight::U
) where {T,U<:Real,M<:AbstractMultiplexUGraph{T,U},V<:MultilayerVertex{T}}

    if !(has_vertex(mg, src) && has_vertex(mg, dst))
        println("Either vertex $src or $dst does not exist within the multiplex graph")
        return false
    elseif src.layer != dst.layer
        throw(ErrorException("Adding an edge between vertices of different layers is not allowed within a multiplex graph."))
    end

    src_layer_cartIdxs, _src_layer = get_layer(mg, src.layer)
    dst_layer_cartIdxs, _dst_layer = get_layer(mg, dst.layer)
    added = addede_symmetric = false
    layerContainingEdge_cartIdxs, layerContainingEdge = get_layer(mg, _src_layer.name)
    istrait(IsWeighted{typeof(layerContainingEdge.graph)}) || throw(
        ErrorException(
            "You are trying to add a weighted edge to an unweighted layer. If you are convinced that the graph type of the layer to which the two vertices $(src) and $(dst) belong is weighted, then please submit an issue so that we may give the `IsWeighted` trait to the graph type.",
        ),
    )
    added = added_symmetric = add_edge!(layerContainingEdge, src, dst, weight)
    if added
        mg.adjacency_tensor[
            src.node, dst.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]
        ] = weight
        mg.adjacency_tensor[
            dst.node, src.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]
        ] = weight
        return true
    else
        return false
    end

    return added
end

"""
    add_edge!(mg::M, src::V, dst::V, weight::U) where { T, U <: Real, M <: AbstractMultiplexUGraph{T,U}, V <: MultilayerVertex{T}}

Add edge from `src` to `dst` with weight `one(U)` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, src::V, dst::V
) where {T,U<:Real,M<:AbstractMultiplexUGraph{T,U},V<:MultilayerVertex{T}}

    if !(has_vertex(mg, src) && has_vertex(mg, dst))
        println("Either vertex $src or $dst does not exist within the multiplex graph")
        return false
    elseif src.layer != dst.layer
        throw(ErrorException("Adding an edge between vertices of different layers is not allowed within a multiplex graph."))
    end


    src_layer_cartIdxs, _src_layer = get_layer(mg, src.layer)
    dst_layer_cartIdxs, _dst_layer = get_layer(mg, dst.layer)
    added = addede_symmetric = false
    layerContainingEdge_cartIdxs, layerContainingEdge = get_layer(mg, _src_layer.name)
    added = added_symmetric = add_edge!(layerContainingEdge, src, dst)
    if added
        mg.adjacency_tensor[src.node, dst.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]] = one(
            U
        )
        mg.adjacency_tensor[dst.node, src.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]] = one(
            U
        )
        return true
    else
        return false
    end

    return added
end



"""
    rem_edge!(mg::M, V1::V, V2::V) where { T, U, M <: AbstractMultiplexUGraph{T,U}, V <: MultilayerVertex{T}}


Remove edge from `src` to `dst` in `mg`.
"""
function Graphs.rem_edge!(
    mg::M, src::V, dst::V
) where {T,U,M<:AbstractMultiplexUGraph{T,U},V<:MultilayerVertex{T}}

    if !(has_vertex(mg, src) && has_vertex(mg, dst))
        println("Either vertex $src or $dst does not exist within the multiplex graph")
        return false
    elseif src.layer != dst.layer
        throw(ErrorException("Removing an edge between vertices of different layers is not allowed within a multiplex graph."))
    end

    if !has_edge(mg, src, dst)
        println("The multilayer doesn't have any edge between $src and $dst")
        return false
    end

    src_layer_cartIdxs, _src_layer = get_layer(mg, src.layer)
    dst_layer_cartIdxs, _dst_layer = get_layer(mg, dst.layer)
    # if _src_layer == _dst_layer
    layerContainingEdge_cartIdxs, layerContainingEdge = get_layer(mg, _src_layer.name)
    removed = rem_edge!(layerContainingEdge, src, dst)
    if removed
        mg.adjacency_tensor[src.node, dst.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]] = zero(
            U
        )
        mg.adjacency_tensor[dst.node, src.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]] = zero(
            U
        )
        return true
    else
        return false
    end

end


# Multilayer-specific functions
# function get_graph_of_layers end #approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
"""
    get_graph_of_layers(mg::M) where {M <: AbstractMultiplexUGraph}

Get a [`DiGraphOfGraph`](@ref) of the layers of `mg`. the weight of each edge between layers are obtained by summing all edge weights in the corresponding interlayer. See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function get_graph_of_layers(
    mg::M; normalization::String="total_edges"
) where {M<:AbstractMultiplexUGraph}
    num_layers = length(mg.layers)
    norm_factor = 1
    if cmp(normalization, "total_edges") == 0
        norm_factor = ne(mg)
    end
    n_nodes = nn(mg)
    adjacency_matrix = reshape(
        [
            i != j ? (1 / norm_factor) * n_nodes / 1 : 0.0 for
            i in 1:num_layers for j in 1:num_layers
        ],
        (num_layers, num_layers),
    )

    isapproxsymmetric(adjacency_matrix) ||
        throw(ErrorException("Adjacency / distance matrices must be symmetric"))

    symmetric_adjacency_matrix = Symmetric(adjacency_matrix)

    return GraphOfGraphs(
        getproperty.(collect(values(mg.layers)), :graph), symmetric_adjacency_matrix
    )
end
