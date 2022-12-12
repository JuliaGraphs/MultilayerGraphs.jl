# ne
# has_edge
# edges

"""
    AbstractMultilayerUGraph{T,U} <: AbstractMultilayerGraph{T,U} 

Abstract type representing an undirected multilayer graph.
"""
abstract type AbstractMultilayerUGraph{T,U} <: AbstractMultilayerGraph{T,U} end

# Nodes

# Vertices
"""
    add_vertex!(mg::M, V::MultilayerVertex) where {T, U, M <: AbstractMultilayerUGraph{T,U}}

Add MultilayerVertex `V` to multilayer graph `mg`, provided that `node(V)` is a `Node` of `mg`. Return true if succeeds.
"""
function Graphs.add_vertex!(mg::M, V::MultilayerVertex) where {T, U, M <: AbstractMultilayerUGraph{T,U}}
    !has_node(mg, node(V)) && return false
    has_vertex(mg, V) && return false

    n_nodes = nn(mg)

    # Re-insert the associations with the proper vertex
    v = isempty(domain(mg.v_V_associations)) ? one(T) : maximum(domain(mg.v_V_associations)) + one(T)

    mg.v_V_associations[v] = get_bare_mv(V)

    mg.v_metadata_dict[v] = V.metadata

    push!(mg.fadjlist, HalfEdge{MultilayerVertex,U}[])

    return true
end

"""
    rem_vertex!(mg::AbstractMultilayerUGraph, V::MultilayerVertex)

Remove [MultilayerVertex](@ref) `mv` from `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.rem_vertex!(mg::AbstractMultilayerUGraph, V::MultilayerVertex)
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
            idx_tbr = findfirst(halfedge -> compare_multilayervertices(vertex(halfedge), V), mg.fadjlist[dst_v])
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
function Graphs.has_edge(mg::M, src::T, dst::T) where { T, M <: AbstractMultilayerUGraph{T}}

    has_vertex(mg,src) || return false
    has_vertex(mg,dst) || return false


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


"""
    set_weight!(mg::M, src::MultilayerVertex{L1}, dst::MultilayerVertex{L2}, weight::U) where {L1 <: Symbol, L2 <: Symbol, T,U, M <: AbstractMultilayerGraph{T,U}}

Set the weight of the edge between `src` and `dst` to `weight`. Return true if succeeds (i.e. if the edge exists and the underlying graph chosen for the Layer/Interlayer where the edge lies is weighted under the `IsWeighted` trait).
"""
function set_weight!(mg::M, src::MultilayerVertex, dst::MultilayerVertex, weight::U) where { T,U, M <: AbstractMultilayerUGraph{T,U}}
    descriptor = get_subgraph_descriptor(mg, layer(src), layer(dst))
    is_weighted(descriptor.null_graph) || return false
    has_edge(mg, src, dst) ||  return false

    halfedges_from_src = mg.fadjlist[get_v(mg, src)]
    halfedge_from_src = halfedges_from_src[findfirst(he -> vertex(he) == dst, halfedges_from_src)]
    halfedge_from_src.weight = weight

    halfedges_from_dst = mg.fadjlist[get_v(mg, dst)]
    halfedge_from_dst = halfedges_from_dst[findfirst(he -> vertex(he) == src, halfedges_from_dst)]
    halfedge_from_dst.weight = weight

    return true

end



"""
    set_metadata!(mg::M, src::MultilayerVertex{L1}, dst::MultilayerVertex{L2}, metadata::U) where {L1 <: Symbol, L2 <: Symbol, T,U, M <: AbstractMultilayerGraph{T,U}}

Set the metadata of the edge between `src` and `dst` to `metadata`. Return true if succeeds (i.e. if the edge exists and the underlying graph chosen for the Layer/Interlayer where the edge lies supports metadata at the edge level  under the `IsMeta` trait).
"""
function set_metadata!(mg::M, src::MultilayerVertex, dst::MultilayerVertex, metadata::Union{Tuple, NamedTuple}) where M <: AbstractMultilayerUGraph
    descriptor = get_subgraph_descriptor(mg, layer(src), layer(dst))
    is_meta(descriptor.null_graph) || return false
    has_edge(mg, src, dst) ||  return false

    halfedges_from_src = mg.fadjlist[get_v(mg, src)]
    halfedge_from_src = halfedges_from_src[findfirst(he -> vertex(he) == dst, halfedges_from_src)]
    halfedge_from_src.metadata = metadata

    halfedges_from_dst = mg.fadjlist[get_v(mg, dst)]
    halfedge_from_dst = halfedges_from_dst[findfirst(he -> vertex(he) == src, halfedges_from_dst)]
    halfedge_from_dst.metadata = metadata

    return true

end


# Overloads that make AbstractMultilayerUGraph an extension of Graphs.jl. These are all well-inferred .
"""
    edges(mg::M) where {T,U,M<:AbstractMultilayerUGraph{T,U}}

Return an list of all the edges of `mg`.
"""
function Graphs.edges(mg::M) where {T,U,M<:AbstractMultilayerUGraph{T,U}}

    edge_list = MultilayerEdge{U}[]

    for (_src_v,halfedges) in enumerate(mg.fadjlist)
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

# Layers and Interlayers
"""
    add_layer!(mg::M,layer::L; interlayers_type = "multiplex") where { T, U, G<: AbstractGraph{T}, M <: AbstractMultilayerUGraph{T, U}, L <: Layer{T,U,G}}

Add layer `layer` to `mg`.

# ARGUMENTS

- `mg::M`: the `MultilayerGraph` which the new layer will be added to;
- `new_layer::L`: the new `Layer` to add to `mg`;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: upon addition of a new `Layer`, all the `Interlayer`s between the new and the existing `Layer`s are immediately created. This keyword argument specifies their `null_graph` See the `Layer` constructor for more information. Defaults to `SimpleGraph{T}()`;
- `default_interlayers_structure::String = "multiplex"`: The structure of the `Interlayer`s created by deafault. May either be "multiplex" to have diagonally-coupled only interlayers, or "empty" for empty interlayers. Defaults to "multiplex".
"""
function add_layer!(
    mg::M, new_layer::L; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String ="multiplex"
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerUGraph{T,U},L<:Layer{T,U,G}, H <: AbstractGraph{T}}
    # Check that the layer is directed
    !istrait(IsDirected{typeof(new_layer.graph)}) || throw(
        ErrorException(
            "The `new_layer`'s underlying graph $(new_layer.graph) is directed, so it is not compatible with a `AbstractMultilayerUGraph`.",
        ),
    )

    _add_layer!(mg, new_layer; default_interlayers_null_graph = default_interlayers_null_graph, default_interlayers_structure = default_interlayers_structure)
end

"""
    specify_interlayer!(mg::M, new_interlayer::In; symmetric_interlayer_name::String) where { T, U, G<: AbstractGraph{T}, M <: AbstractMultilayerUGraph{T, U}, In <: Interlayer{T,U,G}}

Specify the interlayer `new_interlayer` as part of `mg`.
"""
function specify_interlayer!(
    mg::M,
    new_interlayer::In
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerUGraph{T,U},In<:Interlayer{T,U,G}}

    !istrait(IsDirected{typeof(new_interlayer.graph)}) || throw(
        ErrorException(
            "The `new_interlayer`'s underlying graphs $(new_interlayer.graph) is directed, so it is not compatible with a `AbstractMultilayerUGraph`.",
        ),
    )

    return _specify_interlayer!(
        mg, new_interlayer;
    )
end

"""
    get_subgraph(mg::M, descriptor::LD) where {T,U, M <: AbstractMultilayerUGraph{T,U}, LD <: LayerDescriptor{T,U}}

Internal function. Instantiate the Layer described by `descriptor` whose vertices and edges are contained in `mg`.
"""
function get_subgraph(mg::M, descriptor::LD) where {T,U, M <: AbstractMultilayerUGraph{T,U}, LD <: LayerDescriptor{T,U}}
    vs = sort([v for (v,mv) in collect(mg.v_V_associations) if mv.layer == descriptor.name])

    _vertices = get_rich_mv.(Ref(mg), vs)#MultilayerVertex{descriptor.name}[]

#=     for (v,mv) in zip(vs, getindex.(Ref(mg.v_V_associations),vs))
        push!(_vertices, get_rich_mv(mg, v)) #MV(mv.node, mv.layer, mg.v_metadata_dict[v] )
    end
 =#
    edge_list = MultilayerEdge{U}[]#MultilayerEdge{U}[MultilayerEdge(src, vertex(halfedge),  weight(halfedge), metadata(halfedge)) for (src,halfedge) in zip(_vertices,mg.fadjlist[vs]) if vertex(halfedge).layer == descriptor.name && get_v(mg, vertex(halfedge)) >= get_v(mg,src) ]
 
    for (src_v,halfedges_from_src) in zip(vs, getindex.(Ref(mg.fadjlist), vs) )
       
        src_bare_V = mg.v_V_associations[src_v]
        for halfedge in halfedges_from_src
            dst_bare_V = vertex(halfedge)
            dst_v = get_v(mg, dst_bare_V)
            # Don't take the same edge twice (except for self-loops: they will be taken twice)
            if dst_v >= src_v && dst_bare_V.layer == descriptor.name
                push!(edge_list, MultilayerEdge(src_bare_V, dst_bare_V,  weight(halfedge), metadata(halfedge)))
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
function get_subgraph(mg::M, descriptor::InD) where {T,U, G<: AbstractGraph{T}, M <: AbstractMultilayerUGraph{T,U}, InD <: InterlayerDescriptor{T,U,G}}
    layer_1_vs = T[]
    layer_2_vs = T[]

    for (v,mv) in collect(mg.v_V_associations)
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

    shortest_name, longest_name = length(layer_1_vs) >= length(layer_2_vs) ?  (descriptor.layer_2, descriptor.layer_1) : (descriptor.layer_1, descriptor.layer_2)
    shortest_vs = shortest_name == descriptor.layer_1 ? layer_1_vs : layer_2_vs

    edge_list = MultilayerEdge{U}[]

    for (src_v,halfedges_from_src) in zip(shortest_vs, getindex.(Ref(mg.fadjlist), shortest_vs) )
        src_bare_V = mg.v_V_associations[src_v]

        for halfedge in halfedges_from_src
            dst_bare_V = vertex(halfedge)
            
            # Don't take the same edge twice (except for self-loops: they will be taken twice, but are later removed using `unique`)
            if dst_bare_V.layer == longest_name
                push!(edge_list, MultilayerEdge(src_bare_V, dst_bare_V,  weight(halfedge), metadata(halfedge)))
            else
                continue
            end
        end
    end

    return _Interlayer(layer_1_multilayervertices, layer_2_multilayervertices, unique(edge_list), descriptor)
end

# Graphs.jl's internals extra overrides
"""
    degree(mg::M, v::V) where {T,M<:AbstractMultilayerUGraph{T,<:Real},V<:MultilayerVertex}

Return the degree of MultilayerVertex `v` within `mg`.
"""
Graphs.degree(mg::M, v::V ) where {T,M<:AbstractMultilayerUGraph{T,<:Real},V<:MultilayerVertex} = indegree(mg, v)

"""
    is_directed(mg::AbstractMultilayerUGraph)

Return `true` if `mg` is directed, `false` otherwise. 
"""
Graphs.is_directed(mg::AbstractMultilayerUGraph) = false

"""
    is_directed(m::M) where { M <: Type{ <: AbstractMultilayerUGraph}}

Return `true` if `mg` is directed, `false` otherwise. 
"""
Graphs.is_directed(mg::M) where {M<:Type{<:AbstractMultilayerUGraph}} = false

"""
    inneighbors(mg::M, v::T) where {M <: AbstractMultilayerUGraph{T} } where { T <: Integer}

Return the list of inneighbors of `v` within `mg`.
"""
Graphs.inneighbors(mg::M, v::T) where {T,M<:AbstractMultilayerUGraph{T,<:Real}} = outneighbors(mg, v)


# Multilayer-specific functions
# function get_overlay_monoplex_graph end #approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
"""
    get_overlay_monoplex_graph(mg::M) where {M<: AbstractMultilayerUGraph}
Get overlay monoplex graph (i.e. the graph that has the same nodes as `mg` but the link between node `i` and `j` has weight equal to the sum of all edges weights between the various vertices representing `i` and `j` in `mg`, accounting for only within-layer edges). See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function get_overlay_monoplex_graph(mg::M) where {T,U,M<:AbstractMultilayerUGraph{T,U}}
    wgt = weight_tensor(mg).array
    projected_overlay_adjacency_matrix = sum([
        wgt[:, :, i, i] for i in 1:size(wgt, 3)
    ])
    return SimpleWeightedGraph{T,U}(
        projected_overlay_adjacency_matrix
    )
end

"""
    von_neumann_entropy(mg::M) where {T,U,  M <: AbstractMultilayerUGraph{T, U}}

Compute the Von Neumann entropy of `mg`, according to [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Only for undirected multilayer graphs.
"""
function von_neumann_entropy(mg::M) where {T,U,M<:AbstractMultilayerUGraph{T,U}}

    wgt = weight_tensor(mg).array

    num_nodes = length(nodes(mg))
    num_layers = length(mg.layers)

    # Multistrength tensor
    Δ = ein"ijkm,ik,nj,om -> njom"(
        wgt,
        ones(T, size(wgt)[[1, 3]]),
        δk(num_nodes),
        δk(num_layers),
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