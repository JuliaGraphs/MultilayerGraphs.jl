# ==
# nl
# nIn
# has_layer
# _add_layer!
# rem_layer!
# _specify_interlayer!
# get_interlayer
# get_layer_idx

# Base.eltype
# weighttype
# edgetype

# has_vertex
# multilayer_vertices
# mv_inneighbors
# inneighbors
# mv_outneighbors
# outneighbors
# nv
# nv_withmissing
# vertices
# nodes

"""
    AbstractMultilayerGraph{T <: Integer, U <: Real} <: AbstractGraph{T}

An abstract type for multilayer graphs. It is a subtype of AbstractGraph and its concrete subtypes may extend Graphs.jl.

Its concrete subtypes must have the following fields:

- `idx_N_associations::Bijection{Int64,Node}`:;
- `v_V_associations::Bijection{T,<:MultilayerVertex}`:;
- `v_metadata_dict::Dict{T,<:Union{<:Tuple,<:NamedTuple}}`:;
- `layers`: an indexable collection of `Layer`s.
- `interlayers`:a collection of `Interlayer`s;
- `layers_names`: a collection of the names of the layers.
- `subgraphs_names`: a collection of the names of all the subgraphs.
- `fadjlist::Vector{Vector{HalfEdge{<:MultilayerVertex,<:Union{Nothing,U}}}}`: the forward adjacency list.
"""
abstract type AbstractMultilayerGraph{T<:Integer,U<:Real} <: AbstractGraph{T} end

# General MultilayerGraph Utilities
fadjlist(mg::AbstractMultilayerGraph) = mg.fadjlist

# Nodes
"""
    nodes(mg::AbstractMultilayerGraph

Return the nodes of the AbstractMultilayerGraph `mg`, in order of addition.
"""
function nodes(mg::AbstractMultilayerGraph)
    return Node[node for (_, node) in sort!(collect(mg.idx_N_associations); by=first)]
end

"""
    nn(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of nodes in `mg`.
"""
nn(mg::AbstractMultilayerGraph) = length(nodes(mg))

"""
    has_node(mg::AbstractMultilayerGraph, n::Node)

Return true if `n` is a node of `mg`.
"""
has_node(mg::AbstractMultilayerGraph, n::Node) = n in image(mg.idx_N_associations)

"""
    _add_node!(mg::AbstractMultilayerGraph, n::Node; add_vertex_to_layers::Union{Vector{Symbol}, Symbol} = Symbol[])

Add node `n` to `mg`. Return true if succeeds. Additionally, add a corresponding vertex to all layers whose name is listed in `add_vertex_to_layers`. If `add_vertex_to_layers == :all`, then a corresponding vertex is added to all layers.
"""
function _add_node!(
    mg::AbstractMultilayerGraph,
    n::Node;
    add_vertex_to_layers::Union{Vector{Symbol},Symbol}=Symbol[],
)
    !has_node(mg, n) || return false

    maximum_idx =
        isempty(domain(mg.idx_N_associations)) ? 0 : maximum(domain(mg.idx_N_associations))

    mg.idx_N_associations[maximum_idx + 1] = n

    if add_vertex_to_layers == :all
        for layer_name in mg.layers_names
            _add_vertex!(mg, MV(n, layer_name))
        end
    elseif add_vertex_to_layers isa Vector{Symbol}
        for layer_name in add_vertex_to_layers
            _add_vertex!(mg, MV(n, layer_name))
        end
    end

    return true
end

"""
    _rem_node!(mg::AbstractMultilayerGraph, n::Node)

Remove node `n` to `mg`. Return true if succeeds.
"""
function _rem_node!(mg::AbstractMultilayerGraph, n::Node)
    has_node(mg, n) || return false

    Vs_tbr = MultilayerVertex[]
    for V in image(mg.v_V_associations)
        if V.node == n
            push!(Vs_tbr, V)
        end
    end

    for mv in Vs_tbr
        _rem_vertex!(mg, mv)
    end

    idx_tbr = mg.idx_N_associations(n)

    delete!(mg.idx_N_associations, idx_tbr)

    return true
end

# Vertices
"""
    eltype(::M) where {T,M<:AbstractMultilayerGraph{T}}

Return the vertex type of `mg`.
"""
Base.eltype(::M) where {T,M<:AbstractMultilayerGraph{T}} = T

"""
    has_vertex(mg::M, v::T) where {T,M <: AbstractMultilayerGraph{T}}

Return true if `v` is in mg, else false.
"""
function Graphs.has_vertex(mg::M, v::T) where {T,M<:AbstractMultilayerGraph{T}}
    return v in domain(mg.v_V_associations)
end

"""
    has_vertex(mg::AbstractMultilayerGraph, mv::MultilayerVertex) 

Return true if `mv` is in `mg`, else false.
"""
function Graphs.has_vertex(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return get_bare_mv(mv) in image(mg.v_V_associations)
end

"""
    mv_vertices(mg::AbstractMultilayerGraph)

Return a list of the `MultilayerVertex`s contained in `mg`.
"""
mv_vertices(mg::AbstractMultilayerGraph) = [get_rich_mv(mg, v) for v in vertices(mg)]

"""
    nv(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of vertices in `mg`, excluding the missing vertices.
"""
Graphs.nv(mg::AbstractMultilayerGraph) = length(mg.v_V_associations)

"""
    vertices(mg::M) where {M<:AbstractMultilayerGraph}

Return the collection of the vertices of `mg`.
"""
function Graphs.vertices(mg::M) where {M<:AbstractMultilayerGraph}
    return sort(collect(domain(mg.v_V_associations)))
end

"""
    get_metadata(mg::AbstractMultilayerGraph, bare_mv::MultilayerVertex)

Return the metadata associated to `MultilayerVertex` mv (regardless of metadata assigned to `bare_mv`).
"""
function get_metadata(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return mg.v_metadata_dict[get_v(mg, mv)]
end

"""
    set_metadata!(mg::AbstractMultilayerGraph, mv::MultilayerVertex, metadata::Union{Tuple, NamedTuple}) 

Set the metadata of vertex `mv` to `metadata`. Return true if succeeds
"""
function set_metadata!(
    mg::AbstractMultilayerGraph, mv::MultilayerVertex, metadata::Union{Tuple,NamedTuple}
)
    descriptor = mg.layers[get_layer_idx(mg, layer(mv))]
    is_meta(descriptor.null_graph) || return false
    has_vertex(mg, mv) || return false
    mg.v_metadata_dict[get_v(mg, mv)] = metadata

    return true
end

# Edges
"""
    edgetype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}

Return the edge type for `mg`.
"""
Graphs.edgetype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} = MultilayerEdge{U}

"""
    ne(mg::AbstractMultilayerGraph)

Return the number of edges in `mg`.
"""
Graphs.ne(mg::AbstractMultilayerGraph) = length(edges(mg))

"""
    has_edge(mg::AbstractMultilayerGraph, edge::MultilayerEdge) 

Return true if `mg` has an edge between the source and the destination of `edge` (does not check edge or vertex metadata).
"""
function Graphs.has_edge(mg::AbstractMultilayerGraph, edge::MultilayerEdge)
    return has_edge(mg, get_v(mg, src(edge)), get_v(mg, dst(edge)))
end

"""
    has_edge(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)

Return true if `mg` has edge between the `src` and `dst` (does not check edge or vertex metadata).
"""
function Graphs.has_edge(
    mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex
)
    return has_edge(mg, get_v(mg, src), get_v(mg, dst))
end

"""
    add_edge!(mg::M, src::T, dst::T; weight::Union{Nothing, U} = one(U), metadata::Union{Tuple,NamedTuple} = NamedTuple() ) where {T,U, M <: AbstractMultilayerGraph{T,U}} 

Internal method. Add a MultilayerEdge between `src` and `dst` with weight `weight` and metadata `metadata`. Return true if succeeds, false otherwise.
"""
function Graphs.add_edge!(
    mg::M,
    src::T,
    dst::T;
    weight::Union{Nothing,U}=one(U),
    metadata::Union{Tuple,NamedTuple}=NamedTuple(),
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    return add_edge!(
        mg, ME(mg.v_V_associations[src], mg.v_V_associations[dst], weight, metadata)
    )
end

"""
    add_edge!(mg::M, src::V, dst::V; weight::Union{Nothing, U} = one(U), metadata::Union{Tuple,NamedTuple} = NamedTuple() ) where {T,U, M <: AbstractMultilayerGraph{T,U}, V <: MultilayerVertex}

Add a MultilayerEdge between `src` and `dst` with weight `weight` and metadata `metadata`. Return true if succeeds, false otherwise.
"""
function Graphs.add_edge!(
    mg::M,
    src::V,
    dst::V;
    weight::Union{Nothing,U}=one(U),
    metadata::Union{Tuple,NamedTuple}=NamedTuple(),
) where {T,U,M<:AbstractMultilayerGraph{T,U},V<:MultilayerVertex}
    return add_edge!(mg, ME(src, dst, weight, metadata))
end

"""
    rem_edge!(mg::M, src::T, dst::T) where {T, M <: AbstractMultilayerGraph{T}}

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.rem_edge!(mg::M, src::T, dst::T) where {T,M<:AbstractMultilayerGraph{T}}
    return rem_edge!(mg, mg.v_V_associations[src], mg.v_V_associations[dst])
end

"""
    rem_edge!(mg::AbstractMultilayerGraph, me::MultilayerEdge)

Remove edge from `src(me)` to `dst(me)` from `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.rem_edge!(mg::AbstractMultilayerGraph, me::MultilayerEdge)
    return rem_edge!(mg, src(me), dst(me))
end

"""
    get_halfegde(mg::M, src::MultilayerVertex, dst::MultilayerVertex) where M <: AbstractMultilayerGraph

Internal function. Return the `HalfEdge`, if it exists, between `src` and `dst`. Error if there is no `HalfEdge`.
"""
function get_halfegde(
    mg::M, src::MultilayerVertex, dst::MultilayerVertex
) where {M<:AbstractMultilayerGraph}
    has_edge(mg, src, dst) ||
        throw(ErrorException("There is no `HalfEdge` from `src` to `dst`"))
    halfedges_from_src = fadjlist(mg)[get_v(mg, src)]
    return halfedges_from_src[findfirst(
        halfedge -> vertex(halfedge) == dst, halfedges_from_src
    )]
end

"""
    get_metadata(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)

Return the metadata associated to the `MultilayerEdge` from `src` to `dst`.
"""
function get_metadata(
    mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex
)
    return get_halfegde(mg, src, dst).metadata
end

"""
    get_weight(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)h

Return the weight associated to the `MultilayerEdge` from `src` to `dst`.
"""
function SimpleWeightedGraphs.get_weight(
    mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex
)
    return get_halfegde(mg, src, dst).weight
end

# Layers and Interlayers

"""
    nl(mg::AbstractMultilayerGraph)

Return the number of layers in `mg`.
"""
nl(mg::AbstractMultilayerGraph) = length(mg.layers)

"""
    nIn(mg::AbstractMultilayerGraph)

Return the number of interlayers in `mg`.
"""
nIn(mg::AbstractMultilayerGraph) = length(mg.interlayers)

"""
    has_layer(mg::AbstractMultilayerGraph, layer_name::Symbol)

Return true in `layer_name` is a name of a `[Layer](@ref)` of `mg`.
"""
has_layer(mg::AbstractMultilayerGraph, layer_name::Symbol) = layer_name in mg.layers_names

"""
    _add_layer!(mg::M,new_layer::L; new_default_interlayers_type::H) where { T, U, G <: AbstractGraph{T}, H <: AbstractGraph{T}, M <: AbstractMultilayerGraph{T, U}, L <: Layer{T,U,G}

Internal function. It is called by the `add_layer!` API functions, which needs to specify the default interlayer graph type.
"""
function _add_layer!(
    mg::M,
    new_layer::L;
    default_interlayers_null_graph::H,
    default_interlayers_structure::String="multiplex",
) where {
    T,
    U,
    M<:AbstractMultilayerGraph{T,U},
    G<:AbstractGraph{T},
    L<:Layer{T,U,G},
    H<:AbstractGraph{T},
}

    # Check that the new layer has a name different from all the existing ones
    new_layer.name ∉ mg.subgraphs_names || throw(
        ErrorException(
            "The new layer has the same name as an existing layer within the multilayer graph. Layers' names must be unique. Existing layers names are $(mg.layers_names).",
        ),
    )

    # Check that the default_interlayers_null_graph argument is indeed empty
    if !(nv(default_interlayers_null_graph) == ne(default_interlayers_null_graph) == 0)
        throw(
            ErrorException(
                "The `default_interlayer_empty_graph` has not been assigned to an empty graph. Expected 0 vertices and 0 edges, found $(nv(default_interlayers_null_graph)) vertices and $(ne(default_interlayers_null_graph)) edges.",
            ),
        )
    end

    # Add the new layer
    push!(mg.layers, new_layer.descriptor)

    # Add the new nodes that the new layer represents
    new_nodes = setdiff(nodes(new_layer), nodes(mg))
    for new_node in new_nodes
        add_node!(mg, new_node)
    end

    # Add vertices
    for vertex in mv_vertices(new_layer)
        _add_vertex!(mg, vertex)
    end

    # Add edges
    for edge in edges(new_layer)
        add_edge!(mg, edge)
    end

    # Add default interlayers
    if default_interlayers_structure == "multiplex"
        for layer_descriptor in mg.layers
            if layer_descriptor.name != new_layer.name
                _specify_interlayer!(
                    mg,
                    multiplex_interlayer(
                        new_layer,
                        getproperty(mg, layer_descriptor.name),
                        default_interlayers_null_graph,
                    ),
                )
            end
        end
    elseif default_interlayers_structure == "empty"
        for layer_descriptor in mg.layers
            if layer_descriptor.name != new_layer.name
                _specify_interlayer!(
                    mg,
                    empty_interlayer(
                        new_layer,
                        getproperty(mg, layer_descriptor.name),
                        default_interlayers_null_graph,
                    ),
                )
            end
        end
    else
        throw(
            ErrorException(
                "Default interlayer structured as '$default_interlayers_structure' not yet implemented. Only 'multiplex' and 'null' are available.",
            ),
        )
    end

    return true
end

"""
    rem_layer!(mg::AbstractMultilayerGraph, layer_name::Symbol; remove_nodes::Bool = false)

Remove layer `layer_name` from multilayer graph `mg`. If `remove_nodes` is true, also remove from the multilayer graph all the nodes associated with the layer. Warning: this action has multilayer-wide consequences, amd may inadvertently remove vertices and edges that were meant to be kept.
"""
function rem_layer!(
    mg::AbstractMultilayerGraph, layer_name::Symbol; remove_nodes::Bool=false
)
    layer_idx = get_layer_idx(mg, layer_name)
    isnothing(layer_idx) && return false

    layer_tbr = getproperty(mg, layer_name)

    if remove_nodes
        for node in nodes(layer_tbr)
            rem_node!(mg, node)
        end
    else
        for mv in mv_vertices(layer_tbr)
            rem_vertex!(mg, mv)
        end
    end

    deleteat!(mg.layers, layer_idx)

    keys_tbr = Set{Symbol}[]
    for connected_layers_set in keys(mg.interlayers)
        if layer_name ∈ connected_layers_set
            push!(keys_tbr, connected_layers_set)
        end
    end

    delete!.(Ref(mg.interlayers), keys_tbr)
    return true
end

"""
    function _specify_interlayer!(
        mg::M, new_interlayer::In
    ) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerGraph{T,U},In<:Interlayer{T,U,G}}

Internal function. It is called by the `specify_interlayer!` API functions.
"""
function _specify_interlayer!(
    mg::M, new_interlayer::In
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerGraph{T,U},In<:Interlayer{T,U,G}}
    all(in.([new_interlayer.layer_1, new_interlayer.layer_2], Ref(mg.layers_names))) ||
    #  throw(
    #      ErrorException(
        throw(
            ErrorException(
                "The new interlayer connects two layers that are not (one or both) part of the multilayer graph. Make sure you spelled the `layer_1` and `layer_2` arguments of the `Interlayer` correctly. Available layers are $(mg.layers_names), found $(new_interlayer.layer_1) and $(new_interlayer.layer_2).",
            ),
        )
    #      ),
    #  )

    # Check that it has the correct number of nodes on both layers
    (
        isempty(
            setdiff(
                Set(new_interlayer.layer_1_nodes),
                Set(nodes(Base.getproperty(mg, new_interlayer.layer_1))),
            ),
        ) && isempty(
            setdiff(
                Set(new_interlayer.layer_2_nodes),
                Set(nodes(Base.getproperty(mg, new_interlayer.layer_2))),
            ),
        )
    ) || throw(
        ErrorException(
            "The nodes in the interlayer $(new_interlayer.name) do not correspond to the nodes in the respective layers $(new_interlayer.layer_1) and $(new_interlayer.layer_2). Found $( setdiff(Set(new_interlayer.layer_1_nodes), Set(nodes(Base.getproperty(mg, new_interlayer.layer_1)))) ) and $(setdiff(Set(new_interlayer.layer_2_nodes), Set(nodes(Base.getproperty(mg, new_interlayer.layer_2)))))",
        ),
    )

    # A rem_interlayer! function may not exist since there always must be all interlayers. We then proceed to effectively remove the interlayer here
    key = Set([new_interlayer.layer_1, new_interlayer.layer_2])
    if haskey(mg.interlayers, key)
        existing_interlayer = getproperty(mg, mg.interlayers[key].name)

        for edge in edges(existing_interlayer)
            rem_edge!(mg, edge)
        end
    end

    mg.interlayers[Set([new_interlayer.layer_1, new_interlayer.layer_2])] =
        new_interlayer.descriptor

    for edge in edges(new_interlayer)
        @assert add_edge!(mg, edge)
        # assert success
    end

    return true
end

"""
    get_interlayer(
        mg::AbstractMultilayerGraph, layer_1_name::Symbol, 
        layer_2_name::Symbol
    )

Return the `Interlayer` between `layer_1` and `layer_2`.
"""
function get_interlayer(
    mg::AbstractMultilayerGraph, layer_1_name::Symbol, layer_2_name::Symbol
)
    layer_1_name ∈ mg.layers_names || throw(
        ErrorException(
            "$layer_1_name doesn't belong to the multilayer graph. Available layers are $(mg.layers_names).",
        ),
    )
    layer_2_name ∈ mg.layers_names || throw(
        ErrorException(
            "$layer_2_name doesn't belong to the multilayer graph. Available layers are $(mg.layers_names).",
        ),
    )
    layer_1_name != layer_2_name || throw(
        ErrorException(
            "`layer_1` argument is the same as `layer_2`. There is no interlayer between a layer and itself.",
        ),
    )

    names = [layer_1_name, layer_2_name]
    for interlayer_descriptor in values(mg.interlayers)
        if all(interlayer_descriptor.layers_names .== names) #issetequal(interlayer_descriptor.layers_names, names_set)
            return get_subgraph(mg, interlayer_descriptor)
        elseif all(interlayer_descriptor.layers_names .== reverse(names))
            interlayer = get_subgraph(mg, interlayer_descriptor)
            return get_symmetric_interlayer(
                interlayer; symmetric_interlayer_name=String(interlayer.name) * "_rev"
            )
        end
    end
end

"""
    get_layer_idx(mg::M, layer_name::Symbol) where {T, U, M <: AbstractMultilayerGraph{T, U}}

Return the index of the `Layer` whose name is `layer_name` within `mg.layers`.
"""
function get_layer_idx(
    mg::M, layer_name::Symbol
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    idx = findfirst(descriptor -> descriptor.name == layer_name, mg.layers)
    if !isnothing(idx)
        return idx
    else
        return nothing
    end
end

"""
    get_subgraph_descriptor(mg::M, layer_1_name::Symbol, layer_2_name::Symbol) where {T,U,M<:AbstractMultilayerGraph{T,U}}  

Return the descriptor associated to the interlayer connecting `layer_1` to `layer_2` (or to the Layer named `layer_1` if `layer_1` == `layer_2`)
"""
function get_subgraph_descriptor(
    mg::M, layer_1_name::Symbol, layer_2_name::Symbol
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    if layer_1_name == layer_2_name
        idx = get_layer_idx(mg, layer_1_name)
        if !isnothing(idx)
            return mg.layers[idx]
        else
            throw(
                ErrorException(
                    "The multilayer graph does not contain any Layer named $(layer_1_name). Available layers are $(mg.layers_names).",
                ),
            )
        end
    else
        layer_1_name ∈ mg.layers_names || throw(
            ErrorException(
                "$layer_1_name does nto belong to the multilayer graph. Available layers are $(mg.layers_names).",
            ),
        )
        layer_2_name ∈ mg.layers_names || throw(
            ErrorException(
                "$layer_2_name does nto belong to the multilayer graph. Available layers are $(mg.layers_names).",
            ),
        )

        return mg.interlayers[Set([layer_1_name, layer_2_name])]
    end
end

# Graphs.jl's internals and ecosystem extra overrides
"""
    indegree( mg::AbstractMultilayerGraph, v::MultilayerVertex)      

Get the indegree of vertex `v` in `mg`.
"""
function Graphs.indegree(mg::AbstractMultilayerGraph, mv::V) where {V<:MultilayerVertex}
    return length(inneighbors(mg, mv))
end

"""
    indegree( mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}

Get the vector of indegrees of vertices `vs` in `mg`.
"""
function Graphs.indegree(
    mg::AbstractMultilayerGraph, vs::AbstractVector{<:MultilayerVertex}=mv_vertices(mg)
)
    return [indegree(mg, x) for x in vs]
end

"""
    outdegree(mg::AbstractMultilayerGraph, mv::MultilayerVertex)

Get the outdegree of vertex `v` in `mg`.
"""
function Graphs.outdegree(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return length(outneighbors(mg, mv))
end

"""
    outdegree(mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} 

Get the vector of outdegrees of vertices `vs` in `mg`.
"""
function Graphs.outdegree(
    mg::AbstractMultilayerGraph, vs::AbstractVector{<:MultilayerVertex}=mv_vertices(mg)
)
    return [outdegree(mg, x) for x in vs]
end

"""
    degree(mg::AbstractMultilayerGraph, vs::AbstractVector{<:MultilayerVertex}=vertices(mg)) 

Get the degree of vertices `vs` in `mg`.
"""
function Graphs.degree(
    mg::AbstractMultilayerGraph, vs::AbstractVector{<:MultilayerVertex}=mv_vertices(mg)
)
    return degree.(Ref(mg), vs)
end

"""
    inneighbors( mg::AbstractMultilayerGraph, mv::MultilayerVertex ) 

Return the list of inneighbors of `mv` within `mg`.
"""
function Graphs.inneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return inneighbors(mg, get_v(mg, mv))
end

"""
    outneighbors( mg::AbstractMultilayerGraph, mv::MultilayerVertex)

Return the list of outneighbors of `v` within `mg`.
"""
function Graphs.outneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return outneighbors(mg, get_v(mg, mv))
end

"""
    outneighbors(mg::M, v::T) where {T, M<:AbstractMultilayerGraph{T}}

Return the list of outneighbors of `v` within `mg`.
"""
function Graphs.outneighbors(mg::M, v::T) where {T,M<:AbstractMultilayerGraph{T}}
    _outneighbors = T[]

    for helfedge in mg.fadjlist[v]
        push!(_outneighbors, get_v(mg, vertex(helfedge)))
    end

    return _outneighbors
end

"""
    neighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex) 

Get the neighbors of vertex `mv` in `mg`. Reduces to `outneighbors` for both directed and undirected multilayer graphs.
"""
Graphs.neighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex) = outneighbors(mg, mv)

"""
    weighttype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} 

Return the weight type of `mg` (i.e. the eltype of the weight tensor or the supra-adjacency matrix).
"""
weighttype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} = U

# Multilayer-specific methods
"""
    mv_inneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    
Return the list of `MultilayerVertex` inneighbors of `mv` within `mg`.
"""
function mv_inneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return getindex.(Ref(mg.v_V_associations), inneighbors(mg, mv))
end

"""
    mv_outneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)

Return the list of `MultilayerVertex` outneighbors of `mv` within `mg`.
"""
function mv_outneighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return getindex.(Ref(mg.v_V_associations), outneighbors(mg, mv))
end

"""
    mv_neighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)

Return the list of `MultilayerVertex` neighbors of `mv` within `mg`.
"""
function mv_neighbors(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return getindex.(Ref(mg.v_V_associations), neighbors(mg, mv))
end

"""
    get_supra_weight_matrix_from_weight_tensor(weight_tensor::Array{U, 4}) where { U <: Real}

Internal method. Convert a weight tensor into the corresponding supra-adjacency matrix.
"""
function get_supra_weight_matrix_from_weight_tensor(
    weight_tensor::Array{U,4}
) where {U<:Real}
    N = size(weight_tensor, 1)
    L = size(weight_tensor, 3)
    supra_weight_matrix = Array{U}(undef, N * L, N * L)
    for i in 1:L
        for j in 1:L
            supra_weight_matrix[(N * (i - 1) + 1):(N * i), (N * (j - 1) + 1):(N * j)] .= weight_tensor[
                1:N, 1:N, i, j
            ]
        end
    end
    return supra_weight_matrix
end

"""
    get_weight_tensor_from_supra_weight_matrix(mg::M, supra_weight_matrix::S) where {T, U, S <: Array{U, 2}, M <: AbstractMultilayerGraph{T,U} }

Internal method. Convert a supra-adjacency matrix into the corresponding weight tensor.
"""
function get_weight_tensor_from_supra_weight_matrix(
    mg::M, supra_weight_matrix::S
) where {T,U,S<:Array{U,2},M<:AbstractMultilayerGraph{T,U}}
    N = nn(mg)
    L = nl(mg)
    weight_tensor = zeros(U, N, N, L, L)
    if N != 0
        for i in 1:L
            for j in 1:L
                weight_tensor[1:N, 1:N, i, j] .= supra_weight_matrix[
                    (N * (i - 1) + 1):(N * i), (N * (j - 1) + 1):(N * j)
                ]
            end
        end
        return weight_tensor
    else
        return weight_tensor
    end
end

"""
    weight_tensor(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Compute the weight tensor of `mg`. Return an object of type `WeightTensor`.
"""
function weight_tensor(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    N = nn(mg)
    L = nl(mg)

    _size = (N, N, L, L)

    _weight_tensor = zeros(U, _size...)

    v_V_associations = Bijection{Int,Union{MissingVertex,MultilayerVertex}}()

    for (_src_v, halfedges_from_src) in enumerate(mg.fadjlist)
        if !isempty(halfedges_from_src)
            src_v = T(_src_v)

            src_bare_V = mg.v_V_associations[src_v]
            src_n_idx = mg.idx_N_associations(src_bare_V.node)
            src_layer_idx = get_layer_idx(mg, src_bare_V.layer)

            for halfedge in halfedges_from_src
                dst_bare_V = vertex(halfedge)

                dst_n_idx = mg.idx_N_associations(dst_bare_V.node)
                dst_layer_idx = get_layer_idx(mg, dst_bare_V.layer)
                _weight_tensor[src_n_idx, dst_n_idx, src_layer_idx, dst_layer_idx] = weight(
                    halfedge
                )
            end
        end
    end

    return WeightTensor(_weight_tensor, mg.layers_names, mg.idx_N_associations)
end

"""
    get_v_V_associations_withmissings(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Internal function. Return the `v_V_associations` for `mg` taking into account missing vertices (the order of the vertices of each layer is induced by the order of the nodes in `mg.idx_N_associations`).
"""
function get_v_V_associations_withmissings(
    mg::M
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    v_V_associations = Bijection{T,Union{MissingVertex,MultilayerVertex}}()

    n_nodes = nn(mg)
    n_layers = nl(mg)

    for i in 1:(n_nodes * n_layers)
        v_V_associations[i] = MissingVertex()
    end

    for V in mv_vertices(mg)
        if !(V isa MissingVertex)
            layer_idxs = get_layer_idx(mg, V.layer)
            v = T(mg.idx_N_associations(V.node) + (layer_idxs[1] - 1) * n_nodes)
            delete!(v_V_associations, v)
            v_V_associations[v] = get_bare_mv(V)
        end
    end

    return v_V_associations
end

"""
    supra_weight_matrix(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Compute the supra weight matrix of `mg`. Return an object of type `SupraWeightMatrix`
"""
function supra_weight_matrix(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    n_nodes = nn(mg)
    n_layers = nl(mg)

    v_V_associations = get_v_V_associations_withmissings(mg)

    _supra_weight_matrix = zeros(U, n_nodes * n_layers, n_nodes * n_layers)
    for (_src_v, halfedges) in enumerate(mg.fadjlist)
        src_v = v_V_associations(mg.v_V_associations[_src_v])

        for halfedge in halfedges
            _weight = isnothing(weight(halfedge)) ? one(U) : weight(halfedge)
            dst_v = v_V_associations(vertex(halfedge))
            # The loop goes through the entire fadjlist, so [dst_v, src_v] will be obtained at some point
            _supra_weight_matrix[src_v, dst_v] = _weight
        end
    end

    return SupraWeightMatrix(_supra_weight_matrix, v_V_associations)
end

"""
    metadata_tensor(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Compute the weight tensor of `mg`. Return an object of type `WeightTensor`.
"""
function metadata_tensor(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    N = nn(mg)
    L = nl(mg)

    _metadata_tensor = Array{Union{Nothing,Tuple,NamedTuple}}(nothing, N, N, L, L)

    for (_src_v, halfedges_from_src) in enumerate(mg.fadjlist)
        if !isempty(halfedges_from_src)
            src_v = T(_src_v)

            src_bare_V = mg.v_V_associations[src_v]
            src_n_idx = mg.idx_N_associations(src_bare_V.node)
            src_layer_idx = get_layer_idx(mg, src_bare_V.layer)

            for halfedge in halfedges_from_src
                dst_bare_V = vertex(halfedge)

                dst_n_idx = mg.idx_N_associations(dst_bare_V.node)
                dst_layer_idx = get_layer_idx(mg, dst_bare_V.layer)
                _metadata_tensor[src_n_idx, dst_n_idx, src_layer_idx, dst_layer_idx] = metadata(
                    halfedge
                )
            end
        end
    end

    return MetadataTensor(_metadata_tensor, mg.layers_names, mg.idx_N_associations)
end

"""
    mean_degree(mg::AbstractMultilayerGraph)

Return the mean of the degree sequence of `mg`.
"""
mean_degree(mg::AbstractMultilayerGraph) = mean(degree(mg))

"""
    degree_second_moment(mg::AbstractMultilayerGraph) 

Calculate the second moment of the degree sequence of `mg`.
"""
degree_second_moment(mg::AbstractMultilayerGraph) = mean(degree(mg) .^ 2)

"""
    degree_variance(mg::AbstractMultilayerGraph)

Return the variance of the degree sequence of `mg`.
"""
degree_variance(mg::AbstractMultilayerGraph) = var(degree(mg))

"""
    multilayer_global_clustering_coefficient(
        mg::AbstractMultilayerGraph, 
        norm_factor::Union{Float64,Symbol}=:max
    )

Return the complete multilayer global clustering coefficient, equal to the ratio of realized triplets over all possible triplets, including those whose every or some edges belong to interlayers, normalized by `norm_factor`. If `norm_factor == :max`, then the ratio is normalized by `maximum(array(weight_tensor(mg)))`, else it is not normalized. This function does not override Graphs.jl's `global_clustering_coefficient`, since the latter does not consider cliques where two nodes are the same node but in different layers/interlayers. See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function multilayer_global_clustering_coefficient(
    mg::AbstractMultilayerGraph, norm_factor::Union{Float64,Symbol}=:max
)
    wgt = weight_tensor(mg).array

    _normalization_inverse = 1.0
    if norm_factor == :max
        _normalization_inverse = 1.0 / maximum(wgt)
    end

    A_right = wgt .- get_diagonal_elements(wgt)

    # DeDomenico2013 numerator implementation. Inconsistent with both Wikipedia and Graphs.jl's implementation
    num = ein"ijkm,jnmo,niok ->"(A_right, A_right, A_right)[]

    # Wikipedia-informed denominator implementation (consistent with Graphs.jl's global_clustering_coefficient)
    # ntriangles = 0 
    # for vertex in vertices(mg)
    #     k = degree(mg, vertex)
    #     ntriangles += k * (k - 1)
    # end

    # DeDomenico2013 denominator implementation
    F = ones(size(wgt)...) .- multilayer_kronecker_delta(size(wgt))
    den = ein"ijkm,jnmo,niok ->"(A_right, F, A_right)[]

    return _normalization_inverse * (num / den)
end

"""
    multilayer_weighted_global_clustering_coefficient(mg::M, norm_factor::Union{Float64, Symbol} = :max) where {M <: AbstractMultilayerGraph}

Return the complete multilayer global clustering coefficient, equal to the ratio of realized triplets over all possible triplets, including those whose every or some edges belong to interlayers, normalized by `norm_factor`. Each triplets contributes for `w[1]` if all of its vertices are in one layer, `w[2]` if its vertices span two layers, and `w[3]` if they span 3 layers. If `norm_factor == :max`, then the ratio is normalized by `maximum(array(weight_tensor(mg)))`, else it is not normalized. This function does not override Graphs.jl's `global_clustering_coefficient`, since the latter does not consider cliques where two nodes are the same node but in different layers/interlayers. See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function multilayer_weighted_global_clustering_coefficient(
    mg::M, w::Vector{Float64}, norm_factor::Union{Float64,Symbol}=:max
) where {M<:AbstractMultilayerGraph} #This is well defined for both weighted and unweighted multilayer graphs
    wgt = weight_tensor(mg).array

    sum(w) == 1 ||
        throw(ErrorException("Weight vector `w` does not sum to 1. Found $(sum(w))."))
    _normalization_inverse = 1.0
    if norm_factor == :max
        _normalization_inverse = 1.0 / maximum(wgt)
    end

    A_right = wgt .- get_diagonal_elements(wgt)
    num_layers = size(wgt, 3)

    num = ein"ijkm,jnmo,niok,skmo,s ->"(A_right, A_right, A_right, δ_Ω(num_layers), w)[]

    F = ones(size(wgt)...) .- multilayer_kronecker_delta(size(wgt))

    den = ein"ijkm,jnmo,niok,skmo,s ->"(A_right, F, A_right, δ_Ω(num_layers), w)[]

    return _normalization_inverse * (num / den)
end

"""
    overlay_clustering_coefficient(
        mg::AbstractMultilayerGraph, 
        norm_factor::Union{Float64,Symbol}=:max
    )

Return the overlay clustering coefficient as calculated in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). If `norm_factor == :max`, then the ratio is normalized by `maximum(array(weight_tensor(mg)))`, else it is not normalized. 
"""
function overlay_clustering_coefficient(
    mg::AbstractMultilayerGraph, norm_factor::Union{Float64,Symbol}=:max
)
    wgt = weight_tensor(mg).array

    _normalization_inverse = 1.0
    if norm_factor == :max
        _normalization_inverse = 1.0 / (maximum(ein"ijkl->ij"(wgt)) / length(mg.layers))
        # Check that we are using OMEinsum correctly.
        # @assert all(ein"ijkl->ij"(mg.array) .== dropdims(sum(mg.array, dims = (3,4)), dims = (3,4)))
    end

    num = ein"ij,jm,mi ->"(ein"ijkm->ij"(wgt), ein"ijkm->ij"(wgt), ein"ijkm->ij"(wgt))[]

    F = ones(size(wgt)...) - multilayer_kronecker_delta(size(wgt))

    den = ein"ij,jm,mi ->"(ein"ijkm->ij"(wgt), ein"ijkm->ij"(F), ein"ijkm->ij"(wgt))[]

    return _normalization_inverse * (num / den)
end

"""
    eigenvector_centrality(
        mg::M;
        norm::String = "1",
        tol::Float64 = 1e-6,
        maxiter::Int64 = 2000
        ) where {T, U, M <: AbstractMultilayerGraph{T, U}}

Calculate the eigenvector centrality of `mg` via an iterative algorithm. The `norm` parameter may be `"1"` or `"n"`,  and respectively the eigenvector centrality will be normalized to 1 or further divided by the number of nodes of `mg`. The `tol` parameter terminates the approximation when two consecutive iteration differ by no more than  `tol`. The `maxiters` parameter terminates the algorithm when it goes beyond `maxiters` iterations.

The returned values are: the eigenvector centrality and the relative error at each algorithm iteration, that is, the summed absolute values of the componentwise differences between the centrality computed at the current iteration minus the centrality computed at the previous iteration.

Note: in the limit case of a monoplex graph, this function outputs a eigenvector centrality vector that coincides the one outputted by Graphs.jl's `eigenvector_centrality`.
"""
function Graphs.eigenvector_centrality(
    mg::M; weighted::Bool=true, norm::String="1", tol::Float64=1e-6, maxiter::Int64=2000
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    swm_m = nothing
    v_V_associations_withmissings = nothing
    if weighted
        swm = supra_weight_matrix(mg)
        swm_m = swm.array
        v_V_associations_withmissings = swm.v_V_associations
    else
        swm = supra_weight_matrix(mg)
        swm_m = swm.array .!= zero(U)
        v_V_associations_withmissings = swm.v_V_associations
    end

    num_nodes = length(nodes(mg))
    X = ones(Float64, num_nodes * length(mg.layers))

    err = 1.0
    errs = Float64[]
    iter = 0
    while err > tol && iter < maxiter
        new_X = ein"ij,i -> j"(swm_m, X)
        new_X = new_X ./ sqrt(sum(abs2, new_X))
        err = sum(abs.(X .- new_X))
        push!(errs, err)
        X .= new_X
        iter += 1
    end

    if norm == "1"
        X = X ./ sum(X)
    elseif norm == "n"
        X = X ./ (sum(X) / num_nodes)
    end

    # Reorder centrality values ot march the order multilayer vertices are given in mv_vertices(mg)
    layers_mvs_ordered = [
        couple[2] for couple in sort(collect(v_V_associations_withmissings); by=first) if
        !(couple[2] isa MissingVertex)
    ]

    _sortperm = sortperm(mg.v_V_associations.(layers_mvs_ordered))

    X = vec(X)[[
        couple[1] for couple in sort(collect(v_V_associations_withmissings); by=first) if
        !(couple[2] isa MissingVertex)
    ]]

    X = X[_sortperm]

    return X, errs
end

"""
    modularity(
        mg::M, 
        c::Matrix{Int64}; 
        null_model::Union{String,Array{U,4}} = "degree"
    ) where {T, U, M <: AbstractMultilayerGraph{T,U}}

Calculate the modularity of `mg`, as shown in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function Graphs.modularity(
    mg::M, c::Matrix{Int64}; null_model::Union{String,Array{U,4}}="degree"
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    wgt = weight_tensor(mg).array

    # Check that c has the correct size
    n_nodes = length(nodes(mg))
    n_layers = length(mg.layers)
    size(c) == (n_nodes, n_layers) || throw(
        ErrorException(
            "The size of the community matrix does not match (nn(mg),length(mg.layers)), found $(size(c)) and $((nn(mg),length(mg.layers))).",
        ),
    )

    # Build S
    n_communities = length(unique(c))
    S = Array{Bool}(undef, n_nodes, n_layers, n_communities)
    for (i, community) in enumerate(unique(c))
        S[:, :, i] .= (c .== community)
    end

    # Build P
    P = Array{Float64}(undef, size(wgt))
    tot_links = length(edges(mg))
    if typeof(null_model) != String && size(null_model) == size(P)
        P .= null_model
    elseif typeof(null_model) != String && size(null_model) != size(P)
        throw(
            ErrorException(
                "size of `null_model` does not match the size of the adjacency tensor. Got $(size(null_model)) and $(size(P)) respectively.",
            ),
        )
    elseif null_model == "degree"
        for cart_idx in CartesianIndices(P)
            layer_1_idx = cart_idx[3]
            layer_2_idx = cart_idx[4]

            mv1 = MultilayerVertex(
                mg.idx_N_associations[cart_idx[1]], mg.layers[layer_1_idx].name
            )
            mv2 = MultilayerVertex(
                mg.idx_N_associations[cart_idx[2]], mg.layers[layer_2_idx].name
            )
            if has_vertex(mg, mv1) && has_vertex(mg, mv2) # haskey(mg.v_V_associations.finv, mv1 ) && haskey(mg.v_V_associations.finv, mv2 )
                P[cart_idx] = (degree(mg, mv1) * degree(mg, mv2)) / (2 * tot_links - 1)
            else
                P[cart_idx] = 0.0
            end
        end
    else
        throw(ErrorException("Null model '$null_model' not implemented."))
    end

    # Build B
    B = wgt .- P

    # Build K
    K = ein"ijkm,jimk -> "(wgt, ones(T, size(wgt)...))[]

    return (1 / K) * ein"ija,ikjm,kma->"(S, B, S)[]
end

# Base overloads
"""
    Base.(==)(x::AbstractMultilayerGraph, y::AbstractMultilayerGraph)

Overload equality for `AbstractMultilayerGraph`s.
"""
function Base.:(==)(x::AbstractMultilayerGraph, y::AbstractMultilayerGraph)
    typeof(x) == typeof(y) || false

    for field in fieldnames(typeof(x))
        if @eval $x.$field != $y.$field
            return false
        end
    end
    return true
end

# Utilities
"""
    get_v(mg::AbstractMultilayerGraph, V::MultilayerVertex)

Internal method. Get the Integer label associated to `mv` within `mg.v_V_associations`.
"""
function get_v(mg::AbstractMultilayerGraph, mv::MultilayerVertex)
    return mg.v_V_associations(get_bare_mv(mv))
end

"""
    get_rich_mv(mg::M, i::T) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Return `V` together with its metadata.
"""
function get_rich_mv(
    mg::M, i::T; perform_checks::Bool=false
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    if perform_checks
        haskey(mg.v_V_associations, i) ||
            throw(ErrorException("$i is not a vertex of the multilayer graph."))
    end

    bare_V = mg.v_V_associations[i]
    return MV(bare_V.node, bare_V.layer, mg.v_metadata_dict[i])
end

# Console print utilities
function to_string(x::AbstractMultilayerGraph)
    unionall_type = typeof(x).name.wrapper
    parameters = typeof(x).parameters

    layers_names = name.(x.layers)
    layers_underlying_graphs = typeof.(graph.(x.layers))
    layers_table = pretty_table(
        hcat(layers_names, layers_underlying_graphs);
        title="### LAYERS",
        column_labels=(["NAME", "UNDERLYING GRAPH"]),
        alignment=:c,
        column_label_alignment=:c,
        style=TextTableStyle(; first_line_column_label=crayon"yellow bold"),
        table_format=TextTableFormat(; @text__all_horizontal_lines()),
    )

    interlayers_names = name.(values(x.interlayers))
    interlayers_underlying_graphs = typeof.(graph.(values(x.interlayers)))
    interlayer_layer_1s = getproperty.(values(x.interlayers), Ref(:layer_1))
    interlayer_layer_2s = getproperty.(values(x.interlayers), Ref(:layer_2))
    interlayer_tranfers = getproperty.(values(x.interlayers), Ref(:transfer_vertex_metadata))

    interlayers_table = pretty_table(
        hcat(
            interlayers_names,
            interlayer_layer_1s,
            interlayer_layer_2s,
            interlayers_underlying_graphs,
            interlayer_tranfers,
        );
        title="### INTERLAYERS",
        column_labels=([
            "NAME", "LAYER 1", "LAYER 2", "UNDERLYING GRAPH", "TRANSFER VERTEX METADATA"
        ]),
        alignment=:c,
        column_label_alignment=:c,
        style=TextTableStyle(; first_line_column_label=crayon"yellow bold"),
        table_format=TextTableFormat(; @text__all_horizontal_lines()),
    )

    return """
           `$unionall_type` with vertex type `$(parameters[1])` and weight type `$(parameters[2])`.

           $layers_table

           $interlayers_table
           """
end
Base.show(io::IO, x::AbstractMultilayerGraph) = print(io, to_string(x))

"""
    getproperty(mg::AbstractMultilayerGraph, f::Symbol)
"""
function Base.getproperty(mg::AbstractMultilayerGraph, f::Symbol)
    if f in (
        :v_V_associations,
        :fadjlist,
        :idx_N_associations,
        :layers,
        :interlayers,
        :v_metadata_dict,
    ) # :weight_tensor, :supra_weight_matrix, 
        Base.getfield(mg, f)
    elseif f == :badjlist && is_directed(mg)
        Base.getfield(mg, f)
    elseif f == :edge_list
        return edges(mg)
    elseif f == :subgraphs
        return merge(mg.layers, mg.interlayers)
    elseif f == :layers_names
        return [layer.name for layer in mg.layers]
    elseif f == :interlayers_names
        return [interlayer.name for interlayer in values(mg.interlayers)]
    elseif f == :subgraphs_names
        return vcat(mg.layers_names, mg.interlayers_names)
    else
        for descriptor in mg.layers
            if descriptor.name == f
                return get_subgraph(mg, descriptor)
            end
        end

        for descriptor in values(mg.interlayers)
            if descriptor.name == f
                return get_subgraph(mg, descriptor)
            end
        end
    end
end
