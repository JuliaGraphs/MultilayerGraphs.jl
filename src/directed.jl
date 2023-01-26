# General MultilayerDiGraph Utilities
@traitfn badjlist(mg::M) where {M <: AbstractMultilayerGraph; IsDirected{M}} = mg.badjlist

# Nodes
# Vertices
"""
    _add_vertex!(mg::M, V::MultilayerVertex) where {T, U, M <: AbstractMultilayerDiGraph{T,U}}

Add MultilayerVertex `V` to multilayer graph `mg`.  If `add_node` is true and `node(mv)` is not already part of `mg`, then add `node(mv)` to `mg` before adding `mv` to `mg` instead of throwing an error. Return true if succeeds. 
"""
@traitfn function _add_vertex!(
    mg::M, V::MultilayerVertex; add_node::Bool=true
) where {T,U,M<:AbstractMultilayerGraph{T,U};IsDirected{M}}
    has_vertex(mg, V) && return false
    if add_node
        _node = node(V)
        if add_node && !has_node(mg, _node)
            add_node!(mg, _node)
        end
    else
        !has_node(mg, node(V)) && return false
    end

    n_nodes = nn(mg)

    # Re-insert the associations with the proper vertex
    v = if isempty(domain(mg.v_V_associations))
        one(T)
    else
        maximum(domain(mg.v_V_associations)) + one(T)
    end

    mg.v_V_associations[v] = get_bare_mv(V)

    mg.v_metadata_dict[v] = V.metadata

    push!(mg.fadjlist, HalfEdge{MultilayerVertex,U}[])
    push!(mg.badjlist, HalfEdge{MultilayerVertex,U}[])

    return true
end

"""
    _rem_vertex!(mg::AbstractMultilayerDiGraph, V::MultilayerVertex)

Remove [MultilayerVertex](@ref) `mv` from `mg`. Return true if succeeds, false otherwise.
"""
@traitfn function _rem_vertex!(
    mg::M, V::MultilayerVertex
) where {M <: AbstractMultilayerGraph; IsDirected{M}}
    # Check that the node exists and then that the vertex exists
    has_node(mg, V.node) || return false
    has_vertex(mg, V) || return false

    # Get the v corresponding to V, delete the association and replace it with a MissingVertex. Also substitute the metadata with an empty NamedTuple
    v = get_v(mg, V)
    forward_halfedges_to_be_deleted = mg.fadjlist[v]

    # Loop over the halfedges to be removed in  `forward_halfedges_to_be_deleted`, get their v and remove halfedges to `V` in their corresponding arrays
    for halfedge_to_be_deleted in forward_halfedges_to_be_deleted
        dst_v = get_v(mg, vertex(halfedge_to_be_deleted))
        idx_tbr = nothing

        # Don't attempt to remove twice a self loop
        if dst_v != v
            idx_tbr = findfirst(
                halfedge -> compare_multilayervertices(vertex(halfedge), V),
                mg.badjlist[dst_v],
            )
            deleteat!(mg.badjlist[dst_v], idx_tbr)
        end
    end

    backward_halfedges_to_be_deleted = mg.badjlist[v]
    # Loop over the halfedges to be removed in  `backward_halfedges_to_be_deleted`, get their v and remove halfedges to `V` in their corresponding arrays
    for halfedge_to_be_deleted in backward_halfedges_to_be_deleted
        src_v = get_v(mg, vertex(halfedge_to_be_deleted))
        idx_tbr = nothing

        # Don't attempt to remove twice a self loop
        if src_v != v
            idx_tbr = findfirst(
                halfedge -> compare_multilayervertices(vertex(halfedge), V),
                mg.fadjlist[src_v],
            )
            deleteat!(mg.fadjlist[src_v], idx_tbr)
        end
    end

    maximum_v = maximum(domain(mg.v_V_associations))

    if v != maximum_v
        mg.fadjlist[v] = mg.fadjlist[end]
        pop!(mg.fadjlist)

        mg.badjlist[v] = mg.badjlist[end]
        pop!(mg.badjlist)

        last_V = mg.v_V_associations[maximum_v]

        delete!(mg.v_V_associations, v)
        delete!(mg.v_V_associations, maximum_v)
        mg.v_V_associations[v] = last_V

        mg.v_metadata_dict[v] = mg.v_metadata_dict[maximum_v]
        delete!(mg.v_metadata_dict, maximum_v)
    else
        pop!(mg.fadjlist)
        pop!(mg.badjlist)

        delete!(mg.v_V_associations, v)
        delete!(mg.v_metadata_dict, v)
    end
    return true
end

# Edges
"""
    has_edge(mg::M, src::T, dst::T) where {T,M<:AbstractMultilayerGraph{T}; IsDirected{M}}

Return true if `mg` has edge between the `src` and `dst` (does not check edge or vertex metadata).
"""
@traitfn function Graphs.has_edge(
    mg::M, src::T, dst::T
) where {T,M<:AbstractMultilayerGraph{T};IsDirected{M}}
    # Returns true if src is a vertex of mg
    has_vertex(mg, src) || return false
    # Returns true if dst is a vertex of mg
    has_vertex(mg, dst) || return false

    halfedges_from_src = mg.fadjlist[src]
    halfedges_to_dst = mg.badjlist[dst]

    if length(halfedges_from_src) >= length(halfedges_to_dst)
        dsts_v = get_v.(Ref(mg), vertex.(halfedges_from_src))
        return dst in dsts_v
    else
        srcs_v = get_v.(Ref(mg), vertex.(halfedges_to_dst))
        return src in srcs_v
    end
end

# Overloads that make AbstractMultilayerDiGraph an extension of Graphs.jl. These are all well-inferred.
"""
    edges(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}; IsDirected{M}}

Return an list of all the edges of `mg`.
"""
@traitfn function Graphs.edges(
    mg::M
) where {T,U,M<:AbstractMultilayerGraph{T,U}; IsDirected{M}}
    edge_list = MultilayerEdge{U}[]

    for (_src_v, halfedges) in enumerate(mg.fadjlist)
        src_v = T(_src_v)
        src_V = get_rich_mv(mg, src_v)

        for halfedge in halfedges
            dst_v = get_v(mg, vertex(halfedge))

            dst_V = get_rich_mv(mg, dst_v)
            push!(edge_list, ME(src_V, dst_V, weight(halfedge), metadata(halfedge)))
        end
    end
    # Remove duplicate self loops and return
    return edge_list
end

"""
    _add_edge!(mg::M, me::E) where {T,U,E<:MultilayerEdge{<:Union{U,Nothing}},M<:AbstractMultilayerGraph{T,U}; IsDirected{M}}

Add MultilayerEdge `me` to the AbstractMultilayerDiGraph `mg`. Return true if succeeds, false otherwise.
"""
@traitfn function _add_edge!(
    mg::M, me::E
) where {
    T,U,E<:MultilayerEdge{<:Union{U,Nothing}},M<:AbstractMultilayerGraph{T,U};IsDirected{M}
}
    _src = get_bare_mv(src(me))
    _dst = get_bare_mv(dst(me))
    has_vertex(mg, _src) ||
        throw(ErrorException("Vertex $_src does not belong to the multilayer graph."))
    has_vertex(mg, _dst) ||
        throw(ErrorException("Vertex $_dst does not belong to the multilayer graph."))

    # Add edge to `edge_dict`
    src_V_idx = get_v(mg, _src)
    dst_V_idx = get_v(mg, _dst)

    _weight = isnothing(weight(me)) ? one(U) : weight(me)
    _metadata = metadata(me)

    if !has_edge(mg, _src, _dst)
        push!(mg.fadjlist[src_V_idx], HalfEdge(_dst, _weight, _metadata))
        push!(mg.badjlist[dst_V_idx], HalfEdge(_src, _weight, _metadata))
    else
        @debug "An edge between $(src(me)) and $(dst(me)) already exists"
        return false
    end

    return true
end

"""
    _rem_edge!(
    mg::M, src::MultilayerVertex, dst::MultilayerVertex
) where {M<:AbstractMultilayerGraph; IsDirected{M}}

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
@traitfn function _rem_edge!(
    mg::M, src::MultilayerVertex, dst::MultilayerVertex
) where {M <: AbstractMultilayerGraph; IsDirected{M}}

    # Perform routine checks
    has_vertex(mg, src) ||
        throw(ErrorException("Vertex $_src does not belong to the multilayer graph."))
    has_vertex(mg, dst) ||
        throw(ErrorException("Vertex $_dst does not belong to the multilayer graph."))

    has_edge(mg, src, dst) || return false

    src_V_idx = get_v(mg, src)
    dst_V_idx = get_v(mg, dst)

    _src = get_bare_mv(src)
    _dst = get_bare_mv(dst)

    src_idx_tbr = findfirst(halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx])
    deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)

    dst_idx_tbr = findfirst(halfedge -> halfedge.vertex == _src, mg.badjlist[dst_V_idx])
    deleteat!(mg.badjlist[dst_V_idx], dst_idx_tbr)

    return true
end

"""
    set_weight!(
        mg::M, src::MultilayerVertex, dst::MultilayerVertex, weight::U
    ) where {T,U,M<:AbstractMultilayerGraph{T,U}; IsDirected{M}}

Set the weight of the edge between `src` and `dst` to `weight`. Return true if succeeds (i.e. if the edge exists and the underlying graph chosen for the Layer/Interlayer where the edge lies is weighted under the `IsWeighted` trait).
"""
@traitfn function set_weight!(
    mg::M, src::MultilayerVertex, dst::MultilayerVertex, weight::U
) where {T,U,M<:AbstractMultilayerGraph{T,U};IsDirected{M}}
    # Get the subgraph descriptor for the layer containing both src and dst
    descriptor = get_subgraph_descriptor(mg, layer(src), layer(dst))
    # Check if the subgraph is weighted
    is_weighted(descriptor.null_graph) || return false
    # Check if an edge exists between src and dst
    has_edge(mg, src, dst) || return false
    # Get the halfedge from src to dst
    halfedges_from_src = mg.fadjlist[get_v(mg, src)]
    halfedge_from_src = halfedges_from_src[findfirst(
        he -> vertex(he) == dst, halfedges_from_src
    )]
    # Set the weight of the halfedge from src to dst
    halfedge_from_src.weight = weight
    # Get the halfedge from dst to src
    halfedges_to_dst = mg.badjlist[get_v(mg, dst)]
    halfedge_to_dst = halfedges_to_dst[findfirst(he -> vertex(he) == src, halfedges_to_dst)]
    # Set the weight of the halfedge from dst to src
    halfedge_to_dst.weight = weight
    return true
end

"""
    set_metadata!(
        mg::M,
        src::MultilayerVertex,
        dst::MultilayerVertex,
        metadata::Union{Tuple,NamedTuple},
    ) where {M<:AbstractMultilayerGraph; IsDirected{M}}

Set the metadata of the edge between `src` and `dst` to `metadata`. Return true if succeeds (i.e. if the edge exists and the underlying graph chosen for the Layer/Interlayer where the edge lies supports metadata at the edge level  under the `IsMeta` trait).
"""
@traitfn function set_metadata!(
    mg::M, src::MultilayerVertex, dst::MultilayerVertex, metadata::Union{Tuple,NamedTuple}
) where {M <: AbstractMultilayerGraph; IsDirected{M}}
    # Get the subgraph descriptor that corresponds to the layer of src and dst
    descriptor = get_subgraph_descriptor(mg, layer(src), layer(dst))
    # If the subgraph descriptor's null graph is true, then the edge does not exist
    is_meta(descriptor.null_graph) || return false
    # If the edge does not exist, return false
    has_edge(mg, src, dst) || return false
    # Set the halfedge from src to dst
    halfedges_from_src = mg.fadjlist[get_v(mg, src)]
    halfedge_from_src = halfedges_from_src[findfirst(
        he -> vertex(he) == dst, halfedges_from_src
    )]
    halfedge_from_src.metadata = metadata
    # Set the halfedge from dst to src
    halfedges_to_dst = mg.badjlist[get_v(mg, dst)]
    halfedge_to_dst = halfedges_to_dst[findfirst(he -> vertex(he) == src, halfedges_to_dst)]
    halfedge_to_dst.metadata = metadata
    return true
end

# Layers and Interlayers
"""
    add_layer_directedness!(
        mg::M,
        new_layer::L;
        default_interlayers_null_graph::H=SimpleGraph{T}(),
        default_interlayers_structure::String="multiplex",
    ) where {
        T,
        U,
        G<:AbstractGraph{T},
        L<:Layer{T,U,G},
        H<:AbstractGraph{T},
        M<:AbstractMultilayerGraph{T,U}; 
        IsDirected{M}
    }
    
Add layer `layer` to `mg`.

# ARGUMENTS
    
- `mg::M`: the `MultilayerDiGraph` which the new layer will be added to;
- `new_layer::L`: the new `Layer` to add to `mg`;
- `default_interlayers_null_graph::H`: upon addition of a new `Layer`, all the `Interlayer`s between the new and the existing `Layer`s are immediately created. This keyword argument specifies their `null_graph` See the `Layer` constructor for more information. Defaults to `SimpleGraph{T}()`;
- `default_interlayers_structure::String`: The structure of the `Interlayer`s created by default. May either be "multiplex" to have diagonally-coupled only interlayers, or "empty" for empty interlayers. Defaults to "multiplex".
"""
@traitfn function add_layer_directedness!(
    mg::M,
    new_layer::L;
    default_interlayers_null_graph::H=SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {
    T,
    U,
    G<:AbstractGraph{T},
    L<:Layer{T,U,G},
    H<:AbstractGraph{T},
    M<:AbstractMultilayerGraph{T,U};IsDirected{M},
}
    # Check that the layer is directed
    istrait(IsDirected{typeof(new_layer.graph)}) || throw(
        ErrorException(
            "The `new_layer`'s underlying graph $(new_layer.graph) is undirected, so it is not compatible with a `AbstractMultilayerDiGraph`.",
        ),
    )

    return _add_layer!(
        mg,
        new_layer;
        default_interlayers_null_graph=default_interlayers_null_graph,
        default_interlayers_structure=default_interlayers_structure,
    )
end

"""
    get_subgraph(
        mg::M, descriptor::LD
    ) where {T,U,LD<:LayerDescriptor{T,U}, M<:AbstractMultilayerGraph{T,U}; IsDirected{M}}

Internal function. Instantiate the Layer described by `descriptor` whose vertices and edges are contained in `mg`.
"""
@traitfn function get_subgraph(
    mg::M, descriptor::LD
) where {T,U,LD<:LayerDescriptor{T,U},M<:AbstractMultilayerGraph{T,U};IsDirected{M}}
    vs = sort([
        v for (v, mv) in collect(mg.v_V_associations) if mv.layer == descriptor.name
    ])

    _vertices = get_rich_mv.(Ref(mg), vs)

    edge_list = MultilayerEdge{U}[]

    # We may also loop over the badjlist
    for (src_v, halfedges_from_src) in zip(vs, getindex.(Ref(mg.fadjlist), vs))
        src_bare_V = mg.v_V_associations[src_v]
        for halfedge in halfedges_from_src
            dst_bare_V = vertex(halfedge)

            if dst_bare_V.layer == descriptor.name
                push!(
                    edge_list,
                    MultilayerEdge(
                        src_bare_V, dst_bare_V, weight(halfedge), metadata(halfedge)
                    ),
                )
            else
                continue
            end
        end
    end

    return Layer(descriptor, _vertices, edge_list)
end

"""
    get_subgraph(
        mg::M, descriptor::InD
    ) where {
        T,
        U,
        # G<:AbstractGraph{T},
        InD<:InterlayerDescriptor{T,U},# G},
        M<:AbstractMultilayerGraph{T,U}; 
        Is

Internal function. Instantiate the Interlayer described by `descriptor` whose vertices and edges are contained in `mg`.
"""
@traitfn function get_subgraph(
    mg::M, descriptor::InD
) where {
    T,
    U,
    # G<:AbstractGraph{T},
    InD<:InterlayerDescriptor{T,U},# G},
    M<:AbstractMultilayerGraph{T,U};IsDirected{M},
}
    layer_1_vs = T[]
    layer_2_vs = T[]

    for (v, mv) in collect(mg.v_V_associations)
        if mv.layer == descriptor.layer_1
            push!(layer_1_vs, v)

        elseif mv.layer == descriptor.layer_2
            push!(layer_2_vs, v)
        else
            continue
        end
    end

    sort!(layer_1_vs)
    sort!(layer_2_vs)

    layer_1_multilayervertices = get_rich_mv.(Ref(mg), layer_1_vs)
    layer_2_multilayervertices = get_rich_mv.(Ref(mg), layer_2_vs)

    layers_vs = vcat(layer_1_vs, layer_2_vs)

    edge_list = MultilayerEdge{U}[]

    for (src_v, halfedges_from_src) in
        zip(layer_1_vs, getindex.(Ref(mg.fadjlist), layer_1_vs))
        src_bare_V = mg.v_V_associations[src_v]

        for halfedge in halfedges_from_src
            dst_bare_V = vertex(halfedge)
            # dst_v = get_v(mg, dst_bare_V)

            # Don't take the same edge twice (except for self-loops: they will be taken twice, but are later removed using `unique`)
            if dst_bare_V.layer == descriptor.layer_2
                push!(
                    edge_list,
                    MultilayerEdge(
                        src_bare_V, dst_bare_V, weight(halfedge), metadata(halfedge)
                    ),
                )
            else
                continue
            end
        end
    end

    for (src_v, halfedges_from_src) in
        zip(layer_1_vs, getindex.(Ref(mg.fadjlist), layer_2_vs))
        src_bare_V = mg.v_V_associations[src_v]

        for halfedge in halfedges_from_src
            dst_bare_V = vertex(halfedge)
            # dst_v = get_v(mg, dst_bare_V)

            # Don't take the same edge twice (except for self-loops: they will be taken twice, but are later removed using `unique`)
            if dst_bare_V.layer == descriptor.layer_1
                push!(
                    edge_list,
                    MultilayerEdge(
                        src_bare_V, dst_bare_V, weight(halfedge), metadata(halfedge)
                    ),
                )
            else
                continue
            end
        end
    end

    return _Interlayer(
        layer_1_multilayervertices,
        layer_2_multilayervertices,
        unique(edge_list),
        descriptor,
    )
end

# Graphs.jl's internals extra overrides
"""
    degree(
        mg::M, mv::V
    ) where {M<:AbstractMultilayerGraph, V<:MultilayerVertex;  IsDirected{M}}

Return the degree of MultilayerVertex `v` within `mg`.
"""
@traitfn function Graphs.degree(
    mg::M, mv::V
) where {M<:AbstractMultilayerGraph,V<:MultilayerVertex;IsDirected{M}}
    return indegree(mg, mv) + outdegree(mg, mv)
end

"""
    inneighbors(mg::M, v::T) where {T,M<:AbstractMultilayerGraph; IsDirected{M}}


Return the list of inneighbors of `v` within `mg`.
"""
@traitfn function Graphs.inneighbors(
    mg::M, v::T
) where {T,M<:AbstractMultilayerGraph;IsDirected{M}}
    _inneighbors = T[]

    for helfedge in mg.badjlist[v]
        push!(_inneighbors, get_v(mg, vertex(helfedge)))
    end

    return _inneighbors
end

"""
    get_overlay_monoplex_graph(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}; IsDirected{M}}
    
Get overlay monoplex graph (i.e. the graph that has the same nodes as `mg` but the link between node `i` and `j` has weight equal to the sum of all edges weights between the various vertices representing `i` and `j` in `mg`, accounting for both layers and interlayers). See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
@traitfn function get_overlay_monoplex_graph(
    mg::M
) where {T,U,M<:AbstractMultilayerGraph{T,U};IsDirected{M}}
    # Convert the multilayer graph to a weight tensor
    wgt = weight_tensor(mg).array
    # Sum the weights for each edge in the multilayer graph
    projected_overlay_adjacency_matrix = sum([wgt[:, :, i, i] for i in 1:size(wgt, 3)])
    return SimpleWeightedDiGraph{T,U}(projected_overlay_adjacency_matrix)
    # Approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
end
