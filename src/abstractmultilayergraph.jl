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
# Graphs.edgetype

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
"""
abstract type AbstractMultilayerGraph{T<:Integer,U<:Real} <: AbstractGraph{T} end


# Nodes
"""
    nodes(mg::M) where {M <: AbstractMultilayerGraph}

Return the nodes of the AbstractMultilayerGraph `mg`, in order of addition.
"""
nodes(mg::M) where {M<:AbstractMultilayerGraph} = [couple[2] for couple in sort(collect(mg.idx_N_associations), by = first)]

"""
    nn(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of nodes in `mg`.
"""
nn(mg::M) where {M<:AbstractMultilayerGraph} = length(nodes(mg))

"""
    has_node(mg::M, n::Node) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Return true if `n` is a node of `mg`.
"""
has_node(mg::M, n::Node) where {T,U, M <: AbstractMultilayerGraph{T,U}} = n in image(mg.idx_N_associations)

"""
    add_node!(mg::M, n::Node)  where {T,U, M <: AbstractMultilayerGraph{T,U}}

Add node `n` to `mg`. Return true if succeeds.
"""
function add_node!(mg::M, n::Node)  where {T,U, M <: AbstractMultilayerGraph{T,U}}
    !has_node(mg, n) || return false

    maximum_idx = isempty(domain(mg.idx_N_associations)) ?  0 : maximum(domain(mg.idx_N_associations))

    mg.idx_N_associations[maximum_idx + 1] = n

    return true

end

"""
    rem_node!(mg::M, n::Node)  where {T,U, M <: AbstractMultilayerGraph{T,U}}

Remove node `n` to `mg`. Return true if succeeds.
"""
function rem_node!(mg::M, n::Node)  where {T,U, M <: AbstractMultilayerGraph{T,U}}
    has_node(mg, n) || return false
    
    idx_tbr = mg.idx_N_associations(n)

    delete!(mg.idx_N_associations, idx_tbr)

    Vs_tbr = MultilayerVertex[]
    for V in image(mg.v_V_associations)
        if V.node == n
            push!(Vs_tbr,V)
        end
    end

    rem_vertex!.(Ref(mg), Vs_tbr)

    return true
end

# Vertices
"""
    Base.eltype(mg::M) where {M <: AbstractMultilayerGraph}

Return the vertex type of `mg`.
"""
Base.eltype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} = T


"""
    has_vertex(mg::M, v::T) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Return true if `v` is in mg, else false.
"""
Graphs.has_vertex(mg::M, v::T ) where {T,U, M <: AbstractMultilayerGraph{T,U}} = v in domain(mg.v_V_associations) # && !(mg.v_V_associations[v] isa MissingVertex)

"""
    has_vertex(mg::M, mv::MultilayerVertex) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Return true if `mv` is in `mg`, else false.
"""
Graphs.has_vertex(mg::M, mv::MultilayerVertex) where {T,U, M <: AbstractMultilayerGraph{T,U}} = get_bare_mv(mv) in image(mg.v_V_associations)
    
"""
    mv_vertices(mg::AbstractMultilayerGraph)

Return a list of the `MultilayerVertex`s contained in `mg`.
"""
mv_vertices(mg::AbstractMultilayerGraph)  =  [get_rich_mv(mg, v) for v in vertices(mg)]

"""
    nv(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of vertices in `mg`, excluding the missing vertices.
"""
Graphs.nv(mg::M) where {M<:AbstractMultilayerGraph} = length(mg.v_V_associations) #length([mv for mv in image(mg.v_V_associations) if !(mv isa MissingVertex)])

#= """
    nv_withmissing(mg::M) where {M<:AbstractMultilayerGraph}

Return the number of vertices of `mg`, including the missing vertices.
"""
nv_withmissing(mg::M) where {M<:AbstractMultilayerGraph} = length(mg.v_V_associations) =#

"""
    vertices(mg::M) where {M<:AbstractMultilayerGraph}

Return the collection of the vertices of `mg`.
"""
Graphs.vertices(mg::M) where {M<:AbstractMultilayerGraph} = sort(collect(domain(mg.v_V_associations)))

"""
    get_metadata(mg::M, mv::MultilayerVertex) where M <: AbstractMultilayerGraph

Return the metadata associated to `MultilayerVertex` mv.
"""
get_metadata(mg::M, mv::MultilayerVertex) where M <: AbstractMultilayerGraph = mg.v_metadata_dict[get_v(mg, mv)]

"""
    set_metadata!(mg::M, mv::MultilayerVertex, metadata::Union{Tuple, NamedTuple}) where M <: AbstractMultilayerGraph

Set the metadata of vertex `mv` to `metadata`. Return true if succeeds
"""
function set_metadata!(mg::M, mv::MultilayerVertex, metadata::Union{Tuple, NamedTuple}) where M <: AbstractMultilayerGraph
    descriptor = mg.layers[get_layer_idx(mg, layer(mv))]#get_subgraph_descriptor(mg, layer(mv))
    is_meta(descriptor.null_graph) || return false
    has_vertex(mg, mv) || return false
    mg.v_metadata_dict[get_v(mg, mv)] = metadata

    return true

end


# Edges
"""
    edgetype(mg::M) where {M <: AbstractMultilayerGraph}

Return the edge type for `mg`.
"""
Graphs.edgetype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} = MultilayerEdge{U}

"""
    ne(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of edges in `mg`.
"""
Graphs.ne(mg::M) where {M<:AbstractMultilayerGraph} = length(edges(mg))

"""
    has_edge(mg::AbstractMultilayerGraph, edge::MultilayerEdge) 

Return true if `mg` has an edge between the source and the destination of `edge` (does not check edge or vertex metadata).
"""
Graphs.has_edge(mg::AbstractMultilayerGraph, edge::MultilayerEdge)  = has_edge(mg,get_v(mg,src(edge)), get_v(mg,dst(edge)) )

"""
    has_edge(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)

Return true if `mg` has edge between the `src` and `dst` (does not check edge or vertex metadata).
"""
Graphs.has_edge(mg::AbstractMultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)  = has_edge(mg,get_v(mg,src), get_v(mg,dst) )

"""
    add_edge!(mg::M, src::T, dst::T; weight::Union{Nothing, U} = one(U), metadata::Union{Tuple,NamedTuple} = NamedTuple() ) where {T,U, M <: AbstractMultilayerGraph{T,U}} 

Internal method. Add a MultilayerEdge between `src` and `dst` with weight `weight` and metadata `metadata`. Return true if succeeds, false otherwise.
"""
Graphs.add_edge!(mg::M, src::T, dst::T; weight::Union{Nothing, U} = one(U), metadata::Union{Tuple,NamedTuple} = NamedTuple() ) where {T,U, M <: AbstractMultilayerGraph{T,U}} = add_edge!(mg, ME(mg.v_V_associations[src], mg.v_V_associations[dst], weight, metadata))

"""
    add_edge!(mg::M, src::V, dst::V; weight::Union{Nothing, U} = one(U), metadata::Union{Tuple,NamedTuple} = NamedTuple() ) where {T,U, M <: AbstractMultilayerGraph{T,U}, V <: MultilayerVertex}

Add a MultilayerEdge between `src` and `dst` with weight `weight` and metadata `metadata`. Return true if succeeds, false otherwise.
"""
Graphs.add_edge!(mg::M, src::V, dst::V; weight::Union{Nothing, U} = one(U), metadata::Union{Tuple,NamedTuple} = NamedTuple() ) where {T,U, M <: AbstractMultilayerGraph{T,U}, V <: MultilayerVertex} = add_edge!(mg, ME(src, dst, weight, metadata))

"""
    rem_edge!(mg::M, src::T, dst::T) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
Graphs.rem_edge!(mg::M, src::T, dst::T) where {T,U, M <: AbstractMultilayerGraph{T,U}} = rem_edge!(mg, mg.v_V_associations[src], mg.v_V_associations[dst])

"""
    rem_edge!(mg::M, src::T, dst::T) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Remove edge from `src(me)` to `dst(me)` from `mg`. Return true if succeeds, false otherwise.
"""
Graphs.rem_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerGraph{T,U}, E <: MultilayerEdge} = rem_edge!(mg, src(me), dst(me))


"""
    get_halfegde(mg::M, src::MultilayerVertex, dst::MultilayerVertex) where M <: AbstractMultilayerGraph

Internal function. Return the `HalfEdge`, if it exists, between `src` and `dst`. Error if there is no `HalfEdge`.
"""
function get_halfegde(mg::M, src::MultilayerVertex, dst::MultilayerVertex) where M <: AbstractMultilayerGraph
    has_edge(mg, src, dst) || throw(ErrorException("There is no `HalfEdge` from `src` to `dst`"))
    halfedges_from_src = fadjlist(mg)[get_v(mg,src)]
    return halfedges_from_src[findfirst(halfedge -> vertex(halfedge) == dst,halfedges_from_src)]
end


"""
    get_metadata(mg::M, mv::MultilayerVertex) where M <: AbstractMultilayerGraph

Return the metadata associated to the `MultilayerEdge` from `src` to `dst`.
"""
get_metadata(mg::M, src::MultilayerVertex, dst::MultilayerVertex) where M <: AbstractMultilayerGraph  = get_halfegde(mg, src, dst).metadata

"""
    get_weight(mg::M, src::MultilayerVertex, dst::MultilayerVertex) where M <: AbstractMultilayerGraph

Return the weight associated to the `MultilayerEdge` from `src` to `dst`.
"""
SimpleWeightedGraphs.get_weight(mg::M, src::MultilayerVertex, dst::MultilayerVertex) where M <: AbstractMultilayerGraph = get_halfegde(mg, src, dst).weight


# Layers and Interlayers

"""
    nl(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of layers in `mg`.
"""
nl(mg::M) where {M<:AbstractMultilayerGraph} = length(mg.layers)

"""
    nIn(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of interlayers in `mg`.
"""
nIn(mg::M) where {M<:AbstractMultilayerGraph} = length(mg.interlayers)

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
    mg::M, new_layer::L; default_interlayers_null_graph::H, default_interlayers_structure::String ="multiplex"
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
        throw(ErrorException("The `default_interlayer_empty_graph` has not been assigned to an empty graph. Expected 0 vertices and 0 edges, found $(nv(default_interlayers_null_graph)) vertices and $(ne(default_interlayers_null_graph)) edges."))
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
        add_vertex!(mg, vertex)
    end

    # Add edges
    for edge in edges(new_layer)
        add_edge!(mg, edge)
    end

    # Add default interlayers
    if default_interlayers_structure == "multiplex"
        for layer_descriptor in mg.layers
            if layer_descriptor.name != new_layer.name
                _specify_interlayer!( mg, multiplex_interlayer(new_layer, getproperty(mg, layer_descriptor.name), default_interlayers_null_graph) )
            end
        end
    elseif default_interlayers_structure == "empty"
        for layer_descriptor in mg.layers
            if layer_descriptor.name != new_layer.name
                _specify_interlayer!( mg, empty_interlayer(new_layer, getproperty(mg, layer_descriptor.name), default_interlayers_null_graph) )
            end
        end
    else
        throw(ErrorException("Default interlayer structured as '$default_interlayers_structure' not yet implemented. Only 'multiplex' and 'null' are available."))
    end

    return true

end

"""
    rem_layer!(mg::AbstractMultilayerGraph, layer_name::Symbol; remove_nodes::Bool = false)

Remove layer `layer_name` from multilayer graph `mg`. If `remove_nodes` is true, also remove from the multilayer graph all the nodes associated with the layer. Warning: this action has multilayer-wide consequences, amd may inadvertently remove vertices and edges that were meant to be kept.
"""
function rem_layer!(mg::AbstractMultilayerGraph, layer_name::Symbol; remove_nodes::Bool = false)
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

    delete!.(Ref(mg.interlayers),  keys_tbr)
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
         throw(ErrorException("The new interlayer connects two layers that are not (one or both) part of the multilayer graph. Make sure you spelled the `layer_1` and `layer_2` arguments of the `Interlayer` correctly. Available layers are $(mg.layers_names), found $(new_interlayer.layer_1) and $(new_interlayer.layer_2)."))
#      ),
#  )

 # Check that it has the correct number of nodes on both layers
 (isempty(setdiff(Set(new_interlayer.layer_1_nodes), Set(nodes(Base.getproperty(mg, new_interlayer.layer_1)))) ) &&  isempty(setdiff(Set(new_interlayer.layer_2_nodes), Set(nodes(Base.getproperty(mg, new_interlayer.layer_2)))) )) || throw(ErrorException("The nodes in the interlayer $(new_interlayer.name) do not correspond to the nodes in the respective layers $(new_interlayer.layer_1) and $(new_interlayer.layer_2). Found $( setdiff(Set(new_interlayer.layer_1_nodes), Set(nodes(Base.getproperty(mg, new_interlayer.layer_1)))) ) and $(setdiff(Set(new_interlayer.layer_2_nodes), Set(nodes(Base.getproperty(mg, new_interlayer.layer_2)))))"))

# A rem_interlayer! function may not exist since there always must be all interlayers. We then proceed to effectively remove the interlayer here
key = Set([new_interlayer.layer_1, new_interlayer.layer_2])
if haskey(mg.interlayers, key)
    existing_interlayer = getproperty(mg, mg.interlayers[key].name)

    for edge in edges(existing_interlayer)
        rem_edge!(mg, edge)
    end
end

mg.interlayers[Set([new_interlayer.layer_1, new_interlayer.layer_2])] = new_interlayer.descriptor

for edge in edges(new_interlayer)
    success = add_edge!(mg, edge)
    @assert success
end

 return true
end



"""
    get_interlayer(mg::M, layer_1::Symbol, layer_2::Symbol) where {M <: AbstractMultilayerGraph}

Return the `Interlayer` between `layer_1` and `layer_2`.
"""
function get_interlayer(
    mg::M, layer_1_name::Symbol, layer_2_name::Symbol
) where {M<:AbstractMultilayerGraph}

    layer_1_name ∈ mg.layers_names || throw(ErrorException("$layer_1_name doesn't belong to the multilayer graph. Available layers are $(mg.layers_names)."))
    layer_2_name ∈ mg.layers_names || throw(ErrorException("$layer_2_name doesn't belong to the multilayer graph. Available layers are $(mg.layers_names)."))
    layer_1_name != layer_2_name || throw(ErrorException("`layer_1` argument is the same as `layer_2`. There is no interlayer between a layer and itself."))

    names = [layer_1_name, layer_2_name]
    for interlayer_descriptor in values(mg.interlayers)
        if all(interlayer_descriptor.layers_names .== names) #issetequal(interlayer_descriptor.layers_names, names_set)
            return get_subgraph(mg, interlayer_descriptor)
        elseif all(interlayer_descriptor.layers_names .== reverse(names))
            interlayer = get_subgraph(mg, interlayer_descriptor)
            return get_symmetric_interlayer(interlayer; symmetric_interlayer_name = String(interlayer.name)*"_rev" )
        end
    end
end

"""
    get_layer_idx(mg::M, layer_name::Symbol) where {T, U, M <: AbstractMultilayerGraph{T, U}}

Return the index of the `Layer` whose name is `layer_name` within `mg.layers`.
"""
function get_layer_idx(mg::M, layer_name::Symbol) where {T,U,M<:AbstractMultilayerGraph{T,U}}

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
function get_subgraph_descriptor(mg::M, layer_1_name::Symbol, layer_2_name::Symbol) where {T,U,M<:AbstractMultilayerGraph{T,U}}

    if layer_1_name == layer_2_name
        idx = get_layer_idx(mg, layer_1_name)
        if !isnothing(idx)
            return mg.layers[idx] 
        else
            throw(ErrorException("The multilayer graph does not contain any Layer named $(layer_1_name). Available layers are $(mg.layers_names)."))
        end

    else
        layer_1_name ∈ mg.layers_names || throw(ErrorException("$layer_1_name does nto belong to the multilayer graph. Available layers are $(mg.layers_names)."))
        layer_2_name ∈ mg.layers_names ||throw(ErrorException("$layer_2_name does nto belong to the multilayer graph. Available layers are $(mg.layers_names)."))

        return mg.interlayers[Set([layer_1_name, layer_2_name])]
    end
end

# Graphs.jl's internals and ecosystem extra overrides
"""
    indegree( mg::M, v::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}       

Get the indegree of vertex `v` in `mg`.
"""
Graphs.indegree( mg::M, v::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = length(inneighbors(mg, v))

"""
    indegree( mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}

Get the vector of indegrees of vertices `vs` in `mg`.
"""
Graphs.indegree( mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = [indegree(mg, x) for x in vs]

"""
    outdegree(mg::M, v::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}

Get the outdegree of vertex `v` in `mg`.
"""
Graphs.outdegree(mg::M, v::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = length(outneighbors(mg, v))
     
"""
    outdegree(mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} 

Get the vector of outdegrees of vertices `vs` in `mg`.
"""
Graphs.outdegree(mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = [outdegree(mg, x) for x in vs]

"""
    degree(mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}

Get the degree of vertices `vs` in `mg`.
"""
Graphs.degree(mg::M, vs::AbstractVector{V}=vertices(mg)) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = [degree(mg, x) for x in vs]

"""
    inneighbors( mg::M, mv::V ) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}

Return the list of inneighbors of `mv` within `mg`.
"""
Graphs.inneighbors( mg::M, mv::V ) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = inneighbors(mg, get_v(mg, mv))

"""
    outneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}

Return the list of outneighbors of `v` within `mg`.
"""
Graphs.outneighbors( mg::M, mv::V ) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = outneighbors(mg, get_v(mg, mv))


"""
    outneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}

Return the list of outneighbors of `v` within `mg`.
"""
function Graphs.outneighbors(
    mg::M, v::T
) where {M<:AbstractMultilayerGraph{T,<:Real}} where {T}

    _outneighbors = T[]

    for helfedge in mg.fadjlist[v]
        push!(_outneighbors, get_v(mg, vertex(helfedge)))
    end

    return _outneighbors
end

"""
    neighbors(mg::M, v::V) where {T, M <: AbstractMultilayerGraph{T, <: Real}, V <: MultilayerVertex}

Get the neighbors of vertices `vs` in `mg`. Reduces to `outneighbors` for both directed and undirected multilayer graphs.
"""
Graphs.neighbors(mg::M, v::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = outneighbors(mg, v)

"""
    weighttype(mg::M) where {M <: AbstractMultilayerGraph}

Return the weight type of `mg` (i.e. the eltype of the weight tensor or the supra-adjacency matrix).
"""
weighttype(::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} = U

# Multilayer-specific methods
"""
    mv_inneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}
    
Return the list of `MultilayerVertex` inneighbors of `v` within `mg`.
"""
mv_inneighbors(mg::M, mv::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = getindex.(Ref(mg.v_V_associations), inneighbors(mg, mv))

"""
    mv_outneighbors( mg::M, mv::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex}

Return the list of `MultilayerVertex` outneighbors of `v` within `mg`.
"""
mv_outneighbors( mg::M, mv::V) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex} = getindex.(Ref(mg.v_V_associations), outneighbors(mg, mv))

"""
    get_supra_weight_matrix_from_weight_tensor(weight_tensor::Array{U, 4}) where { U <: Real}

Internal method. Convert a weight tensor into the corresponding supra-adjacency matrix.
"""
function get_supra_weight_matrix_from_weight_tensor(weight_tensor::Array{U, 4}) where { U <: Real}
    N = size(weight_tensor,1)
    L = size(weight_tensor,3)
    supra_weight_matrix  = Array{U}(undef, N*L, N*L)
    for i in 1:L
        for j in 1:L
            supra_weight_matrix[(N*(i - 1) + 1):N*i, (N*(j - 1)+1):N*j ] .= weight_tensor[1:N, 1:N, i, j]
        end
    end
    return supra_weight_matrix
end

"""
    get_weight_tensor_from_supra_weight_matrix(mg::M, supra_weight_matrix::S) where {T, U, S <: Array{U, 2}, M <: AbstractMultilayerGraph{T,U} }

Internal method. Convert a supra-adjacency matrix into the corresponding weight tensor.
"""
function get_weight_tensor_from_supra_weight_matrix(mg::M, supra_weight_matrix::S) where {T, U, S <: Array{U, 2}, M <: AbstractMultilayerGraph{T,U} }
    N = nn(mg)
    L = nl(mg)
    weight_tensor  = zeros(U, N, N, L, L)
    if N != 0
        for i in 1:L
            for j in 1:L
                weight_tensor[1:N, 1:N, i, j]  .= supra_weight_matrix[(N*(i - 1) + 1):N*i, (N*(j - 1) + 1):N*j] 
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
function weight_tensor(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}
    N = nn(mg)
    L = nl(mg)

    _size = (N,N,L,L)

    _weight_tensor = zeros(U, _size...)

    v_V_associations = Bijection{Int, Union{MissingVertex, MultilayerVertex}}()
    

    for (_src_v, halfedges_from_src) in enumerate(mg.fadjlist)
        if !isempty(halfedges_from_src)
            src_v = T(_src_v)

            src_bare_V    = mg.v_V_associations[src_v]
            src_n_idx     = mg.idx_N_associations(src_bare_V.node)
            src_layer_idx = get_layer_idx(mg, src_bare_V.layer)
            

            for halfedge in halfedges_from_src
                dst_bare_V = vertex(halfedge)

                dst_n_idx = mg.idx_N_associations(dst_bare_V.node)
                dst_layer_idx = get_layer_idx(mg, dst_bare_V.layer)
#=                 vec_idx = cartIndexTovecIndex( (src_n_idx,dst_n_idx,src_layer_idx ,dst_layer_idx ) ,_size)
                v_V_associations[vec_idx] =  =#
                _weight_tensor[src_n_idx,dst_n_idx,src_layer_idx ,dst_layer_idx] = weight(halfedge)

            end
        end
    end

    return WeightTensor(_weight_tensor, mg.layers_names, mg.idx_N_associations)
end

"""
    get_v_V_associations_withmissings(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Internal function. Return the `v_V_associations` for `mg` taking into account missing vertices (the order of the vertices of each layer is induced by the order of the nodes in `mg.idx_N_associations`).
"""
function get_v_V_associations_withmissings(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}

    v_V_associations = Bijection{T,  Union{MissingVertex, MultilayerVertex} }()

    n_nodes   = nn(mg)
    n_layers  = nl(mg)



    for i in 1:(n_nodes*n_layers)
        v_V_associations[i] = MissingVertex()
    end

    for V in mv_vertices(mg)
        if !(V isa MissingVertex)
            layer_idxs = get_layer_idx(mg, V.layer)
            v = T(mg.idx_N_associations(V.node) + (layer_idxs[1]-1)*n_nodes)
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
function supra_weight_matrix(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}

    n_nodes   = nn(mg)
    n_layers  = nl(mg)

    v_V_associations = get_v_V_associations_withmissings(mg)

    _supra_weight_matrix = zeros(U,n_nodes*n_layers, n_nodes*n_layers)
    for (_src_v,halfedges) in enumerate(mg.fadjlist)
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
function metadata_tensor(mg::M) where {T,U, M <: AbstractMultilayerGraph{T,U}}
    N = nn(mg)
    L = nl(mg)

    _metadata_tensor = Array{Union{Nothing,Tuple,NamedTuple}}( nothing, N, N, L, L)

    for (_src_v, halfedges_from_src) in enumerate(mg.fadjlist)
        if !isempty(halfedges_from_src)
            src_v = T(_src_v)

            src_bare_V    = mg.v_V_associations[src_v]
            src_n_idx     = mg.idx_N_associations(src_bare_V.node)
            src_layer_idx = get_layer_idx(mg, src_bare_V.layer)
            

            for halfedge in halfedges_from_src
                dst_bare_V = vertex(halfedge)

                dst_n_idx = mg.idx_N_associations(dst_bare_V.node)
                dst_layer_idx = get_layer_idx(mg, dst_bare_V.layer)
                _metadata_tensor[src_n_idx, dst_n_idx, src_layer_idx, dst_layer_idx] = metadata(halfedge)

            end
        end
    end

    return MetadataTensor(_metadata_tensor, mg.layers_names, mg.idx_N_associations)
end

"""
    mean_degree(mg::M) where { M <: AbstractMultilayerGraph}

Return the mean of the degree sequence of `mg`.
"""
mean_degree(mg::M) where {M<:AbstractMultilayerGraph} = mean(degree(mg))

"""
    degree_second_moment(mg::M) where { M <: AbstractMultilayerGraph}

Calculate the second moment of the degree sequence of `mg`.
"""
degree_second_moment(mg::M) where {M<:AbstractMultilayerGraph} = mean(degree(mg) .^ 2)

"""
    degree_variance(mg::M) where { M <: AbstractMultilayerGraph}

Return the variance of the degree sequence of `mg`.
"""
degree_variance(mg::M) where {M<:AbstractMultilayerGraph} = var(degree(mg))

"""
    multilayer_clustering_coefficient(mg::M, norm_factor::Union{Float64, Symbol} = :max) where {M <: AbstractMultilayerGraph}

Return the complete multilayer global clustering coefficient, equal to the ratio of realized triplets over all possible triplets, including those whose every or some edges belong to interlayers, normalized by `norm_factor`. If `norm_factor == :max`, then the ratio is normalized by `maximum(mg.array)`. This function does not override Graphs.jl's `global_clustering_coefficient`, since the latter does not consider cliques where two nodes are the same node but in different layers/interlayers. See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function multilayer_global_clustering_coefficient(
    mg::M, norm_factor::Union{Float64,Symbol}=:max
) where {M<:AbstractMultilayerGraph}
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
    F =
        ones(size(wgt)...) .-
        multilayer_kronecker_delta(size(wgt))
    den = ein"ijkm,jnmo,niok ->"(A_right, F, A_right)[]

    return _normalization_inverse * (num / den)
end

"""
    multilayer_weighted_global_clustering_coefficient(mg::M, norm_factor::Union{Float64, Symbol} = :max) where {M <: AbstractMultilayerGraph}

Return the complete multilayer global clustering coefficient, equal to the ratio of realized triplets over all possible triplets, including those whose every or some edges belong to interlayers, normalized by `norm_factor`. Each triplets contributes for `w[1]` if all of its vertices are in one layer, `w[2]` if its vertices span two layers, and `w[3]` if they span 3 layers. If `norm_factor == :max`, then the ratio is normalized by `maximum(mg.array)`. This function does not override Graphs.jl's `global_clustering_coefficient`, since the latter does not consider cliques where two nodes are the same node but in different layers/interlayers. See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
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

    F =
        ones(size(wgt)...) .-
        multilayer_kronecker_delta(size(wgt))

    den = ein"ijkm,jnmo,niok,skmo,s ->"(A_right, F, A_right, δ_Ω(num_layers), w)[]

    return _normalization_inverse * (num / den)
end

"""
    overlay_clustering_coefficient(mg::M, norm_factor::Union{Float64, Symbol} = :max) where {M <: AbstractMultilayerGraph}

Return the overlay clustering coefficient as calculated in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function overlay_clustering_coefficient(
    mg::M, norm_factor::Union{Float64,Symbol}=:max
) where {M<:AbstractMultilayerGraph}

    wgt = weight_tensor(mg).array

    _normalization_inverse = 1.0
    if norm_factor == :max
        _normalization_inverse =
            1.0 / (maximum(ein"ijkl->ij"(wgt)) / length(mg.layers))
        # Check that we are using OMEinsum correctly.
        # @assert all(ein"ijkl->ij"(mg.array) .== dropdims(sum(mg.array, dims = (3,4)), dims = (3,4)))
    end

    num = ein"ij,jm,mi ->"(
        ein"ijkm->ij"(wgt),
        ein"ijkm->ij"(wgt),
        ein"ijkm->ij"(wgt),
    )[]

    F =
        ones(size(wgt)...) -
        multilayer_kronecker_delta(size(wgt))

    den = ein"ij,jm,mi ->"(
        ein"ijkm->ij"(wgt),
        ein"ijkm->ij"(F),
        ein"ijkm->ij"(wgt),
    )[]

    return _normalization_inverse * (num / den)
end

"""
    eigenvector_centrality(mg::M; norm::String = "1", tol::Float64 = 1e-6, maxiter::Int64 = 2000) where {T, U, M <: AbstractMultilayerGraph{T, U}}

Calculate the eigenvector centrality of `mg` via an iterative algorithm. The `norm` parameter may be `"1"` or `"n"`,  and respectively the eigenvector centrality will be normalized to 1 or further divided by the number of nodes of `mg`. The `tol` parameter terminates the approximation when two consecutive iteration differ by no more than  `tol`. The `maxiters` parameter terminates the algorithm when it goes beyond `maxiters` iterations.

The returned values are: the eigenvector centrality and the relative error at each algorithm iteration, that is, the summed absolute values of the componentwise differences between the centrality computed at the current iteration minus the centrality computed at the previous iteration.

Note: in the limit case of a monoplex graph, this function outputs a eigenvector centrality vector that coincides the one outputted by Graphs.jl's `eigenvector_centrality`.
"""
function Graphs.eigenvector_centrality(
    mg::M; weighted::Bool = true,  norm::String="1", tol::Float64=1e-6, maxiter::Int64=2000
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
    layers_mvs_ordered = [couple[2] for couple in sort(collect(v_V_associations_withmissings), by = first) if !(couple[2] isa MissingVertex)]

    _sortperm = sortperm(mg.v_V_associations.(layers_mvs_ordered))

    X = vec(X)[[couple[1] for couple in sort(collect(v_V_associations_withmissings), by = first) if !(couple[2] isa MissingVertex)]]

    X = X[_sortperm]

    return X, errs
end

"""
    modularity(mg::M, c::Matrix{Int64}; null_model::Union{String,Array{U,4}} = "degree") where {T, U, M <: AbstractMultilayerGraph{T,U}}

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
                    P[cart_idx] =
                    (
                        degree(
                            mg,
                            mv1,
                        ) * degree(
                            mg,
                            mv2,
                        )
                    ) / (2 * tot_links - 1)
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
get_v(mg::AbstractMultilayerGraph, mv::MultilayerVertex) = mg.v_V_associations(get_bare_mv(mv))


"""
    get_rich_mv(mg::M, i::T) where {T,U, M <: AbstractMultilayerGraph{T,U}}

Return `V` together with its metadata.
"""
function get_rich_mv(mg::M, i::T; perform_checks::Bool = false) where {T,U, M <: AbstractMultilayerGraph{T,U}}
    if perform_checks
        haskey(mg.v_V_associations,i) || throw(ErrorException("$i is not a vertex of the multilayer graph"))
    end

    bare_V = mg.v_V_associations[i]
    return MV(bare_V.node, bare_V.layer, mg.v_metadata_dict[i])
end