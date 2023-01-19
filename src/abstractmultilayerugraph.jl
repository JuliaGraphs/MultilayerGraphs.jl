#= """
    AbstractMultilayerUGraph{T,U} <: AbstractMultilayerGraph{T,U} 

Abstract type representing an undirected multilayer graph.
"""
abstract type AbstractMultilayerUGraph{T,U} <: AbstractMultilayerGraph{T,U} end =#

# Nodes

# Vertices
"""
    _add_vertex!(mg::M, V::MultilayerVertex) where {T, U, M <: AbstractMultilayerUGraph{T,U}}

Add MultilayerVertex `V` to multilayer graph `mg`.  If `add_node` is true and `node(mv)` is not already part of `mg`, then add `node(mv)` to `mg` before adding `mv` to `mg` instead of throwing an error. Return true if succeeds. 
"""
@traitfn function _add_vertex!(
    mg::M, V::MultilayerVertex;  add_node::Bool=true
) where {T,U, M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}
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

    return true
end

"""
    _rem_vertex!(mg::AbstractMultilayerUGraph, V::MultilayerVertex)

Remove [MultilayerVertex](@ref) `mv` from `mg`. Return true if succeeds, false otherwise.
"""
@traitfn function _rem_vertex!(mg::M, V::MultilayerVertex) where {M <: AbstractMultilayerGraph; !IsDirected{M}}
    # Check that the node exists and then that the vertex exists
    has_node(mg, V.node) || return false
    has_vertex(mg, V) || return false

    # Get the v corresponding to V, delete the association and replace it with a MissingVertex. Also substitute the metadata with an empty NamedTuple
    v = get_v(mg, V)

    halfedges_to_be_deleted = mg.fadjlist[v]

    # Loop over the halfedges to be removed in  `halfedges_to_be_deleted`, get their v and remove halfedges to `V` in their corresponding arrays
    for halfedge_to_be_deleted in halfedges_to_be_deleted
        dst_v = get_v(mg, vertex(halfedge_to_be_deleted))
        idx_tbr = nothing

        # Don't attempt to remove twice a self loop
        if dst_v != v
            idx_tbr = findfirst(
                halfedge -> compare_multilayervertices(vertex(halfedge), V),
                mg.fadjlist[dst_v],
            )
            deleteat!(mg.fadjlist[dst_v], idx_tbr)
        end
    end

    maximum_v = maximum(domain(mg.v_V_associations))

    if v != maximum_v
        mg.fadjlist[v] = mg.fadjlist[end]
        pop!(mg.fadjlist)

        last_V = mg.v_V_associations[maximum_v]

        delete!(mg.v_V_associations, v)
        delete!(mg.v_V_associations, maximum_v)
        mg.v_V_associations[v] = last_V

        mg.v_metadata_dict[v] = mg.v_metadata_dict[maximum_v]
        delete!(mg.v_metadata_dict, maximum_v)
    else
        pop!(mg.fadjlist)
        delete!(mg.v_V_associations, v)
        delete!(mg.v_metadata_dict, v)
    end
    return true
end

# Edges
"""
    has_edge(mg::M, src::T, dst::T) where { T, M <: AbstractMultilayerUGraph{T}}

Return true if `mg` has edge between the `src` and `dst` (does not check edge or vertex metadata).
"""
@traitfn function Graphs.has_edge(mg::M, src::T, dst::T) where {T,M<:AbstractMultilayerGraph{T}; !IsDirected{M}}
    has_vertex(mg, src) || return false
    has_vertex(mg, dst) || return false

    halfedges_from_src = mg.fadjlist[src]
    halfedges_from_dst = mg.fadjlist[dst]

    if length(halfedges_from_src) >= length(halfedges_from_dst)
        dsts_v = get_v.(Ref(mg), vertex.(halfedges_from_src))
        return dst in dsts_v
    else
        srcs_v = get_v.(Ref(mg), vertex.(halfedges_from_dst))
        return src in srcs_v
    end
end


# Overloads that make AbstractMultilayerUGraph an extension of Graphs.jl. These are all well-inferred .
"""
    edges(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}

Return an list of all the edges of `mg`.
"""
@traitfn function Graphs.edges(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}
    edge_list = MultilayerEdge{U}[]

    for (_src_v, halfedges) in enumerate(mg.fadjlist)
        src_v = T(_src_v)
        src_V = get_rich_mv(mg, src_v)

        for halfedge in halfedges
            dst_v = get_v(mg, vertex(halfedge))
            # Don't take the same edge twice, except for self loops that are later removed by the unique 
            if dst_v >= src_v
                dst_V = get_rich_mv(mg, dst_v)
                push!(edge_list, ME(src_V, dst_V, weight(halfedge), metadata(halfedge)))
            end
        end
    end

    # Remove duplicate self loops and return
    return unique(edge_list)
end


"""
    _add_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerUGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add MultilayerEdge `me` to the AbstractMultilayerUGraph `mg`. Return true if succeeds, false otherwise.
"""
@traitfn function _add_edge!(
    mg::M, me::E
) where {T,U,E<:MultilayerEdge{<:Union{U,Nothing}},M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}
    _src = get_bare_mv(src(me))
    _dst = get_bare_mv(dst(me))
    has_vertex(mg, _src) ||
        throw(ErrorException("Vertex $_src does not belong to the multilayer graph."))
    has_vertex(mg, _dst) ||
        throw(ErrorException("Vertex $_dst does not belong to the multilayer graph."))

    # Add edge to `edge_dict`
    src_v = get_v(mg, _src)
    dst_v = get_v(mg, _dst)

    _weight = isnothing(weight(me)) ? one(U) : weight(me)
    _metadata = metadata(me)

    if !has_edge(mg, _src, _dst)
        if src_v != dst_v
            push!(mg.fadjlist[src_v], HalfEdge(_dst, _weight, _metadata))
            push!(mg.fadjlist[dst_v], HalfEdge(_src, _weight, _metadata))
        else
            push!(mg.fadjlist[src_v], HalfEdge(_dst, _weight, _metadata))
        end
        return true
    else
        #=         # Should we modify weight and metadata or should we return false? This may be something to decide ecosystem-wise
                set_weight!(mg, src, dst, _weight)
                set_metadata!(mg, src, dst, _metadata) =#
        @debug "An edge between $(src(me)) and $(dst(me)) already exists"  
        
        return false
    end
end


"""
    _rem_edge!(mg::AbstractMultilayerUGraph, src::MultilayerVertex, dst::MultilayerVertex)

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
@traitfn function _rem_edge!(
    mg::M, src::MultilayerVertex, dst::MultilayerVertex
) where {M<:AbstractMultilayerGraph; !IsDirected{M}}
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

    if get_bare_mv(src) != get_bare_mv(dst)
        src_idx_tbr = findfirst(
            halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx]
        )
        deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)

        dst_idx_tbr = findfirst(halfedge -> halfedge.vertex == _src, mg.fadjlist[dst_V_idx])
        deleteat!(mg.fadjlist[dst_V_idx], dst_idx_tbr)
    else
        src_idx_tbr = findfirst(
            halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx]
        )
        deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)
    end

    return true
end


"""
    set_weight!(mg::M, src::MultilayerVertex{L1}, dst::MultilayerVertex{L2}, weight::U) where {L1 <: Symbol, L2 <: Symbol, T,U, M <: AbstractMultilayerGraph{T,U}}

Set the weight of the edge between `src` and `dst` to `weight`. Return true if succeeds (i.e. if the edge exists and the underlying graph chosen for the Layer/Interlayer where the edge lies is weighted under the `IsWeighted` trait).
"""
@traitfn function set_weight!(
    mg::M, src::MultilayerVertex, dst::MultilayerVertex, weight::U
) where {T,U,M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}
    descriptor = get_subgraph_descriptor(mg, layer(src), layer(dst))
    is_weighted(descriptor.null_graph) || return false
    has_edge(mg, src, dst) || return false

    halfedges_from_src = mg.fadjlist[get_v(mg, src)]
    halfedge_from_src = halfedges_from_src[findfirst(
        he -> vertex(he) == dst, halfedges_from_src
    )]
    halfedge_from_src.weight = weight

    halfedges_from_dst = mg.fadjlist[get_v(mg, dst)]
    halfedge_from_dst = halfedges_from_dst[findfirst(
        he -> vertex(he) == src, halfedges_from_dst
    )]
    halfedge_from_dst.weight = weight

    return true
end

"""
    set_metadata!(mg::AbstractMultilayerUGraph, src::MultilayerVertex, dst::MultilayerVertex, metadata::Union{Tuple, NamedTuple})

Set the metadata of the edge between `src` and `dst` to `metadata`. Return true if succeeds (i.e. if the edge exists and the underlying graph chosen for the Layer/Interlayer where the edge lies supports metadata at the edge level  under the `IsMeta` trait).
"""
@traitfn function set_metadata!(
    mg::M,
    src::MultilayerVertex,
    dst::MultilayerVertex,
    metadata::Union{Tuple,NamedTuple},
) where {M<:AbstractMultilayerGraph; !IsDirected{M}}
    descriptor = get_subgraph_descriptor(mg, layer(src), layer(dst))
    is_meta(descriptor.null_graph) || return false
    has_edge(mg, src, dst) || return false

    halfedges_from_src = mg.fadjlist[get_v(mg, src)]
    halfedge_from_src = halfedges_from_src[findfirst(
        he -> vertex(he) == dst, halfedges_from_src
    )]
    halfedge_from_src.metadata = metadata

    halfedges_from_dst = mg.fadjlist[get_v(mg, dst)]
    halfedge_from_dst = halfedges_from_dst[findfirst(
        he -> vertex(he) == src, halfedges_from_dst
    )]
    halfedge_from_dst.metadata = metadata

    return true
end

# Layers and Interlayers
"""
    add_layer!( mg::M,
        new_layer::L; 
        default_interlayers_null_graph::H = SimpleGraph{T}(), 
        default_interlayers_structure::String ="multiplex"
    ) where {T,U,G<:AbstractGraph{T},L<:Layer{T,U,G}, H <: AbstractGraph{T}, M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}

Add layer `layer` to `mg`.

# ARGUMENTS

- `mg::M`: the `MultilayerGraph` which the new layer will be added to;
- `new_layer::L`: the new `Layer` to add to `mg`;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: upon addition of a new `Layer`, all the `Interlayer`s between the new and the existing `Layer`s are immediately created. This keyword argument specifies their `null_graph` See the `Layer` constructor for more information. Defaults to `SimpleGraph{T}()`;
- `default_interlayers_structure::String = "multiplex"`: The structure of the `Interlayer`s created by default. May either be "multiplex" to have diagonally-coupled only interlayers, or "empty" for empty interlayers. Defaults to "multiplex".
"""
@traitfn function add_layer!(
    mg::M,
    new_layer::L;
    default_interlayers_null_graph::H=SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {
    T,
    U,
    G<:AbstractGraph{T},
    L<:Layer{T,U,G},
    H <: AbstractGraph{T},
    M<:AbstractMultilayerGraph{T,U}; 
    !IsDirected{M}
}

    @assert(eltype(default_interlayers_null_graph) == T, "The eltype of argument default_interlayers_null_graph is not $T, found $(eltype(default_interlayers_null_graph))")
    # Check that the layer is directed
    !istrait(IsDirected{typeof(new_layer.graph)}) || throw(
        ErrorException(
            "The `new_layer`'s underlying graph $(new_layer.graph) is directed, so it is not compatible with a `AbstractMultilayerUGraph`.",
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
    specify_interlayer!(
        mg::M,
        new_interlayer::In
    ) where {T,U,G<:AbstractGraph{T},In<:Interlayer{T,U,G}, M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}

Specify the interlayer `new_interlayer` as part of `mg`.
"""
@traitfn function specify_interlayer!(
    mg::M, new_interlayer::In
) where {T,U,G<:AbstractGraph{T},In<:Interlayer{T,U,G}, M<:AbstractMultilayerGraph{T,U}; !IsDirected{M} } # and(!istrait(IsDirected{M}), !istrait(IsMultiplex{M}))
    !is_directed(new_interlayer.graph) || throw( # !istrait(IsDirected{typeof(new_interlayer.graph)})
        ErrorException(
            "The `new_interlayer`'s underlying graphs $(new_interlayer.graph) is directed, so it is not compatible with a `AbstractMultilayerUGraph`.",
        ),
    )

    return _specify_interlayer!(mg, new_interlayer;)
end

"""
    get_subgraph(mg::M, descriptor::LD) where {T,U, M <: AbstractMultilayerUGraph{T,U}, LD <: LayerDescriptor{T,U}}

Internal function. Instantiate the Layer described by `descriptor` whose vertices and edges are contained in `mg`.
"""
@traitfn function get_subgraph(
    mg::M, descriptor::LD
) where {T,U,LD<:LayerDescriptor{T,U}, M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}
    vs = sort([
        v for (v, mv) in collect(mg.v_V_associations) if mv.layer == descriptor.name
    ])

    _vertices = get_rich_mv.(Ref(mg), vs)

    edge_list = MultilayerEdge{U}[]

    for (src_v, halfedges_from_src) in zip(vs, getindex.(Ref(mg.fadjlist), vs))
        src_bare_V = mg.v_V_associations[src_v]
        for halfedge in halfedges_from_src
            dst_bare_V = vertex(halfedge)
            dst_v = get_v(mg, dst_bare_V)
            # Don't take the same edge twice (except for self-loops: they will be taken twice)
            if dst_v >= src_v && dst_bare_V.layer == descriptor.name
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

    return Layer(descriptor, _vertices, unique(edge_list))
end

"""
    get_subgraph(mg::M, descriptor::InD) where {T,U, G<: AbstractGraph{T}, M <: AbstractMultilayerUGraph{T,U}, InD <: InterlayerDescriptor{T,U,G}}

Internal function. Instantiate the Interlayer described by `descriptor` whose vertices and edges are contained in `mg`.
"""
@traitfn function get_subgraph(
    mg::M, descriptor::InD
) where {
    T,
    U,
    # G<:AbstractGraph{T},
    InD<:InterlayerDescriptor{T,U}, # G},
    M<:AbstractMultilayerGraph{T,U};
    !IsDirected{M}
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

    shortest_name, longest_name = if length(layer_1_vs) >= length(layer_2_vs)
        (descriptor.layer_2, descriptor.layer_1)
    else
        (descriptor.layer_1, descriptor.layer_2)
    end
    shortest_vs = shortest_name == descriptor.layer_1 ? layer_1_vs : layer_2_vs

    edge_list = MultilayerEdge{U}[]

    for (src_v, halfedges_from_src) in
        zip(shortest_vs, getindex.(Ref(mg.fadjlist), shortest_vs))
        src_bare_V = mg.v_V_associations[src_v]

        for halfedge in halfedges_from_src
            dst_bare_V = vertex(halfedge)

            # Don't take the same edge twice (except for self-loops: they will be taken twice, but are later removed using `unique`)
            if dst_bare_V.layer == longest_name
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
    degree(mg::M, v::V) where {T,M<:AbstractMultilayerUGraph{T,<:Real},V<:MultilayerVertex}

Return the degree of MultilayerVertex `v` within `mg`.
"""
@traitfn function Graphs.degree(
    mg::M, v::V
) where {M<:AbstractMultilayerGraph, V<:MultilayerVertex; !IsDirected{M}}
    return indegree(mg, v)
end

#= """
    is_directed(mg::AbstractMultilayerUGraph)

Return `true` if `mg` is directed, `false` otherwise. 
"""
@traitfn Graphs.is_directed(mg::M) where {M<:AbstractMultilayerGraph; !IsDirected{M}} = false

"""
    is_directed(m::M) where { M <: Type{ <: AbstractMultilayerUGraph}}

Return `true` if `mg` is directed, `false` otherwise. 
"""
@traitfn Graphs.is_directed(mg::M) where {M<:Type{<:AbstractMultilayerGraph}; !IsDirected{M} }  = false =#

"""
    inneighbors(mg::M, v::T) where {T,M<:AbstractMultilayerUGraph{T,<:Real}}

Return the list of inneighbors of `v` within `mg`.
"""
@traitfn function Graphs.inneighbors(mg::M, v::T) where {T,M<:AbstractMultilayerGraph; !IsDirected{M}}
    return outneighbors(mg, v)
end

# Multilayer-specific functions
# function get_overlay_monoplex_graph end #approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
"""
    get_overlay_monoplex_graph(mg::M) where {M<: AbstractMultilayerUGraph}
Get overlay monoplex graph (i.e. the graph that has the same nodes as `mg` but the link between node `i` and `j` has weight equal to the sum of all edges weights between the various vertices representing `i` and `j` in `mg`, accounting for only within-layer edges). See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
@traitfn function get_overlay_monoplex_graph(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}
    wgt = weight_tensor(mg).array
    projected_overlay_adjacency_matrix = sum([wgt[:, :, i, i] for i in 1:size(wgt, 3)])
    return SimpleWeightedGraph{T,U}(projected_overlay_adjacency_matrix)
end

"""
    von_neumann_entropy(mg::M) where {T,U,  M <: AbstractMultilayerUGraph{T, U}}

Compute the Von Neumann entropy of `mg`, according to [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Only for undirected multilayer graphs.
"""
@traitfn function von_neumann_entropy(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}; !IsDirected{M}}
    wgt = weight_tensor(mg).array

    num_nodes = length(nodes(mg))
    num_layers = length(mg.layers)

    # Multistrength tensor
    Δ = ein"ijkm,ik,nj,om -> njom"(
        wgt, ones(T, size(wgt)[[1, 3]]), δk(num_nodes), δk(num_layers)
    )

    # Trace of Δ
    tr_Δ = ein"iikk->"(Δ)[]

    # Multilayer laplacian tensor
    L = Δ .- wgt

    # Multilayer density tensor
    ρ = (1.0 / tr_Δ) .* L

    eigvals, eigvects = tensoreig(ρ, [2, 4], [1, 3])

    #=  # Check that we are calculating the right eigenvalues
        lhs = ein"ijkm,ik -> jm"(ρ,eigvects[:,:,1])
        rhs = eigvals[1].*eigvects[:,:,1]
        @assert all(lhs .≈ rhs)
        # Indeed we are =#

    Λ = get_diagonal_weight_tensor(eigvals, size(wgt))

    # Correct for machine precision
    Λ[Λ .< eps()] .= 0.0

    log2Λ = log2.(Λ)

    log2Λ[isinf.(log2Λ)] .= 0

    return -ein"ijkm,jimk ->"(Λ, log2Λ)[]
end
