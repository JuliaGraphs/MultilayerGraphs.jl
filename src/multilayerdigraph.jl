"""
    MultilayerDiGraph{T, U} <: AbstractMultilayerDiGraph{T,U}

A concrete type that can represent a general directed multilayer graph.

# FIELDS

- `adjacency_tensor::Array{U, 4}`: the 4-dimensional tensor that encodes all (weighted) connections within the graph. adjacency_tensor[1,2,3,4] encodes the strength of the (directed or undirected) link between node 1 in layer 3 and node 2 in layer 4.
- `layers::OrderedDict{ Tuple{Int64,Int64}, Layer{T,U,G}}`: the ordered dictionary containing all the layers of the multilayer graph. Their underlying graphs must be all directed.
- `interlayers::OrderedDict{ Tuple{Int64,Int64}, Interlayer{T,U}}`: the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all directed.
"""
mutable struct MultilayerDiGraph{T,U} <: AbstractMultilayerDiGraph{T,U}
    adjacency_tensor::Array{U,4}
    layers::OrderedDict{Tuple{Int64,Int64},Layer{T,U}}
    interlayers::OrderedDict{Tuple{Int64,Int64},Interlayer{T,U}}
end

@traitimpl IsDirected{MultilayerDiGraph}

"""
    MultilayerDiGraph(num_layers::Int64, n_nodes::Int64, min_edges::Int64, max_edges::Int64, graph_types::Vector{DataType})

Return a random `MultilayerDiGraph` with `num_layers` layers, `n_nodes` nodes and each `Layer` and `Interlayer` has a random number of edges between `min_edges` and `max_edges`. `Layers` and `Interlayers` have parametric type `graph_type`.
"""
function MultilayerDiGraph(
    num_layers::Int64,
    n_nodes::Int64,
    min_edges::Int64,
    max_edges::Int64,
    graph_types::Vector{DataType},
)

    #TODO: this method requires that all Graphs.jl extensions are such that:
    # The first type parameter of graphs types is the type of the nodes, like the T in SimpleGraph{T};
    # all graph types have a random constructor with signature (nv::Int64, ne::Int64)
    # all graph types have an empty constructor with signature ()
    all([istrait(IsDirected{graph_type}) for graph_type in graph_types]) || throw(ErrorException("Not all graph types in argument `graph_types` are directed."))
    T = get_common_type([graph_type.parameters[1] for graph_type in graph_types])
    U = get_common_type(
        eltype.([adjacency_matrix(graph_type()) for graph_type in graph_types])
    )

    layers = Layer{T,U}[]
    for i in 1:num_layers
        graph_type = rand(graph_types)
        push!(
            layers,
            Layer(Symbol("layer_$i"), graph_type(n_nodes, rand(min_edges:max_edges)); U=U),
        )
    end

    interlayers_idxs = [(i, j) for i in 1:(num_layers - 1) for j in (i + 1):num_layers]
    interlayers = [
        Interlayer(
            n_nodes,
            Symbol("interlayer_$(i)_$(j)"),
            Symbol("layer_$(i)"),
            Symbol("layer_$(j)"),
            rand(graph_types),
            rand(min_edges:max_edges);
            U=U,
        ) for (i, j) in interlayers_idxs
    ]
    return MultilayerDiGraph(layers, interlayers)::MultilayerDiGraph{T,U}
end

@traitimpl IsWeighted{MultilayerDiGraph}

"""
    MultilayerDiGraph(n_nodes::Int64, T::Type{ <: Number}, U::Type{ <: Number} )

Return a MultilayerDiGraph with `n_nodes` of type `T` nodes and an adjacency tensor eltype `U`. Use this constructor and then add Layers and Interlayers via the `add_layer!` and `add_interlayer!` methods.
"""
function MultilayerDiGraph(n_nodes::Int64, T::Type{<:Number}, U::Type{<:Number})
    adjacency_tensor = zeros(U, n_nodes, n_nodes, 0, 0)
    return MultilayerDiGraph(
        adjacency_tensor,
        OrderedDict{Tuple{Int64,Int64},Layer{T,U}}(),
        OrderedDict{Tuple{Int64,Int64},Interlayer{T,U}}(),
    )
end

"""
    MultilayerDiGraph(layers::Vector{ <: Layer{T,U}};  default_interlayer::String  = "multiplex") where {T,U} 

Construct a MultilayerDiGraph with layers `layers` and all interlayers of type `default_interlayer` (only "multiplex" is allowed).
"""
function MultilayerDiGraph(
    layers::Vector{<:Layer{T,U}}; default_interlayer::String="multiplex"
) where {T,U}
    return MultilayerDiGraph(
        layers, Interlayer{T,U}[]; default_interlayer=default_interlayer
    )
end

"""
    MultilayerDiGraph(layers::Vector{ <: Layer{T,U }}, specified_interlayers::Vector{ <: Interlayer{T,U}};  default_interlayer::String  = "multiplex") where {T, U} 

Construct a MultilayerDiGraph with layers given by `layers`. The interlayers will be constructed by default according to `default_interlayer` (only `"multiplex"` is allowed), except for those specified in `specified_interlayers`.
"""
function MultilayerDiGraph(
    layers::Vector{<:Layer{T,U}},
    specified_interlayers::Vector{<:Interlayer{T,U}};
    default_interlayer::String="multiplex",
) where {T,U}
    # Check that all layers and specified interlayers are directed
    specified_subgraphs = vcat(layers, specified_interlayers)
    all([
        istrait(IsDirected{typeof(subgraph.graph)}) for subgraph in specified_subgraphs
    ]) || throw(
        ErrorException(
            "Not all the underlying layers' and interlayers' graphs are directed."
        ),
    )

    num_nodes = nv(layers[1])
    multilayerdigraph = MultilayerDiGraph(num_nodes, T, U)

    for layer in layers
        add_layer!(
            multilayerdigraph, layer; new_default_interlayers_type=default_interlayer
        )
    end

    if !isnothing(specified_interlayers)
        for interlayer in specified_interlayers
            specify_interlayer!(multilayerdigraph, interlayer)
        end
    end

    return multilayerdigraph
end




"""
    specify_interlayer!(mg::M, layer_1::Symbol, layer_2::Symbol, graph::G; new_interlayer_name::Symbol, symmetric_interlayer_name::Symbol ,  forbidden_vertices::Tuple{Vararg{MultilayerVertex{T}}}, forbidden_edges::Tuple{Vararg{NTuple{2, MultilayerVertex{T}}}} ) where { T, U, G<: AbstractGraph{T}, M <: MultilayerDiGraph{T,U}}

Specify an interlayer which is represented by `graph` between `layer_1` and `layer_2`. The Interlayer's name can be specified via `new_interlayer_name`, as the name of its corresponding Interlayer between `layer_2` and `layer_1` (via the `symmetric_interlayer_name`) which is the same graph just with the adjacency matrix reordered.
"""
function specify_interlayer!(
    mg::M,
    layer_1::Symbol,
    layer_2::Symbol,
    graph::G;
    new_interlayer_name::String="interlayer_$(layer_1)_$(layer_2)",
    symmetric_interlayer_name::String="interlayer_$(layer_2)_$(layer_1)",
    forbidden_vertices::Tuple{Vararg{MultilayerVertex{T}}},
    forbidden_edges::Tuple{Vararg{NTuple{2,MultilayerVertex{T}}}},
) where {T,U,G<:AbstractGraph{T},M<:MultilayerDiGraph{T,U}} # 
    istrait(isDirected{G}) || throw(
        ErrorException(
            "The new interlayer's underlying graph $(graph) is undirected, so it is not compatible with a `MultilayerDiGraph`.",
        ),
    )
    interlayer = Interlayer(
        Symbol(new_interlayer_name),
        layer_1,
        layer_2,
        graph;
        forbidden_vertices=forbidden_vertices,
        forbidden_edges=forbidden_edges,
    )
    return specify_interlayer!(
        mg, interlayer; symmetric_interlayer_name=symmetric_interlayer_name
    )
end



"""
    add_edge!(mg::M, src::V, dst::V, weight::U) where { T, U <: Real, M <: MultilayerDiGraph{T,U}, V <: MultilayerVertex{T}}

Add edge from `src` to `dst` with weight `weight` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, src::V, dst::V, weight::U
) where {T,U<:Real,M<:MultilayerDiGraph{T,U},V<:MultilayerVertex{T}}
    if !(has_vertex(mg, src) && has_vertex(mg, dst))
        return false
    end

    src_layer_cartIdxs, _src_layer = get_layer(mg, src.layer)
    dst_layer_cartIdxs, _dst_layer = get_layer(mg, dst.layer)
    added = false
    if _src_layer == _dst_layer
        layerContainingEdge_cartIdxs, layerContainingEdge = get_layer(mg, _src_layer.name)
        istrait(IsWeighted{typeof(layerContainingEdge.graph)}) || throw(
            ErrorException(
                "You are trying to add a weighted edge to an unweighted layer. If you are convinced that the graph type of the layer to which the two vertices $(src) and $(dst) belong is weighted, then please submit an issue so that we may give the `IsWeighted` trait to the graph type.",
            ),
        )
        added = add_edge!(layerContainingEdge, src, dst, weight)
        if added
            mg.adjacency_tensor[
                src.node, dst.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]
            ] = weight
            return true
        else
            return false
        end
    else
        interlayerContainingEdge_cartIdxs, interlayerContainingEdge = get_interlayer(
            mg, _src_layer.name, _dst_layer.name
        )
        interlayerSymmetricContainingEdge_cartIdxs, interlayerSymmetricContainingEdge = get_interlayer(
            mg, _dst_layer.name, _src_layer.name
        )
        istrait(IsWeighted{typeof(interlayerContainingEdge.graph)}) &&
            istrait(IsWeighted{typeof(interlayerSymmetricContainingEdge.graph)}) || throw(
            ErrorException(
                "You are trying to add a weighted edge to an unweighted interlayer. If you are convinced that the graph type of the interlayer to which the two vertices $(src) and $(dst) belong is weighted, then please submit an issue so that we may give the `IsWeighted` trait to the graph type.",
            ),
        )

        added = add_edge!(interlayerContainingEdge, src, dst, weight)
        added_symmetric = add_edge!(interlayerSymmetricContainingEdge, src, dst, weight)

        added == added_symmetric || throw(
            ErrorException(
                "add_edge!; While adding and edge to an interlayer and its symmetric, one of the two already has/hasn't the said edge while the other doesn't/does. This implies that the multilayer has been corrupted. If you didn't modify its fields without relying on the provided API, please file an issue.",
            ),
        )

        if added && added_symmetric
            mg.adjacency_tensor[
                src.node, dst.node, src_layer_cartIdxs[1], dst_layer_cartIdxs[1]
            ] = weight
            return true
        else
            return false
        end
    end
    return added
end

"""
    add_edge!(mg::M, src::V, dst::V, weight::U) where { T, U <: Real, M <: MultilayerDiGraph{T,U}, V <: MultilayerVertex{T}}

Add edge from `src` to `dst` with weight `one(T)` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, src::V, dst::V
) where {T,U<:Real,M<:MultilayerDiGraph{T,U},V<:MultilayerVertex{T}}
    if !(has_vertex(mg, src) && has_vertex(mg, dst))
        return false
    end

    src_layer_cartIdxs, _src_layer = get_layer(mg, src.layer)
    dst_layer_cartIdxs, _dst_layer = get_layer(mg, dst.layer)
    added = false
    if _src_layer == _dst_layer
        layerContainingEdge_cartIdxs, layerContainingEdge = get_layer(mg, _src_layer.name)
        added = add_edge!(layerContainingEdge, src, dst)
        if added
            mg.adjacency_tensor[src.node, dst.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]] = one(
                U
            )
            return true
        else
            return false
        end
    else
        interlayerContainingEdge_cartIdxs, interlayerContainingEdge = get_interlayer(
            mg, _src_layer.name, _dst_layer.name
        )
        interlayerSymmetricContainingEdge_cartIdxs, interlayerSymmetricContainingEdge = get_interlayer(
            mg, _dst_layer.name, _src_layer.name
        )

        added = add_edge!(interlayerContainingEdge, src, dst)
        added_symmetric = add_edge!(interlayerSymmetricContainingEdge, src, dst)

        added == added_symmetric || throw(
            ErrorException(
                "add_edge!; While adding and edge to an interlayer and its symmetric, one of the two already has/hasn't the said edge while the other doesn't/does. This implies that the multilayer has been corrupted. If you didn't modify its fields without relying on the provided API, please file an issue.",
            ),
        )

        if added && added_symmetric
            mg.adjacency_tensor[src.node, dst.node, src_layer_cartIdxs[1], dst_layer_cartIdxs[1]] = one(
                U
            )
            return true
        else
            return false
        end
    end
    return added
end


"""
    rem_edge!(mg::M, src::V, dst::V) where { T, U, M <: MultilayerDiGraph{T,U}, V <: MultilayerVertex{T}}

Remove edge from `src` to `dst` in `mg`.
"""
function Graphs.rem_edge!(
    mg::M, src::V, dst::V
) where {T,U,M<:MultilayerDiGraph{T,U},V<:MultilayerVertex{T}}
    if !(has_vertex(mg, src) && has_vertex(mg, dst))
        return false
    end

    if !has_edge(mg, src, dst)
        println("The multilayer doesn't have any edge between $src and $dst")
        return false
    end

    src_layer_cartIdxs, _src_layer = get_layer(mg, src.layer)
    dst_layer_cartIdxs, _dst_layer = get_layer(mg, dst.layer)
    if _src_layer == _dst_layer
        layerContainingEdge_cartIdxs, layerContainingEdge = get_layer(mg, _src_layer.name)
        removed = rem_edge!(layerContainingEdge, src, dst)
        if removed
            mg.adjacency_tensor[src.node, dst.node, src_layer_cartIdxs[1], src_layer_cartIdxs[1]] = zero(
                U
            )
            return true
        else
            return false
        end
    else
        interlayerContainingEdge_cartIdxs, interlayerContainingEdge = get_interlayer(
            mg, _src_layer.name, _dst_layer.name
        )
        interlayerSymmetricContainingEdge_cartIdxs, interlayerSymmetricContainingEdge = get_interlayer(
            mg, _dst_layer.name, _src_layer.name
        )

        removed = rem_edge!(interlayerContainingEdge, src, dst)
        removed_symmetric = rem_edge!(interlayerSymmetricContainingEdge, src, dst)

        if removed && removed_symmetric
            mg.adjacency_tensor[src.node, dst.node, src_layer_cartIdxs[1], dst_layer_cartIdxs[1]] = zero(
                U
            )
            return true
        elseif removed == !removed_symmetric
            throw(
                ErrorException(
                    "rem_edge!; While removing and edge to an interlayer and its symmetric, one of the two already has/hasn't the said edge while the other doesn't/does. This implies that the multilayer has been corrupted. If you didn't modify its fields without relying on the provided API, please file an issue.",
                ),
            )
        elseif !removed && !removed_symmetric
            return false
        end
    end
end