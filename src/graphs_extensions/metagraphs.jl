function __add_vertex!(
    g::AbstractMetaGraph{T}; metadata::Union{Tuple,NamedTuple}=NamedTuple()
) where {T<:Integer}
    return add_vertex!(
        g, Dict(Symbol(pair.first) => pair.second for pair in pairs(metadata))
    )
end

function _get_vertex_metadata(g::AbstractMetaGraph{T}, vertex::T) where {T}
    return NamedTuple(props(g, vertex))
end

# MetaGraphs.jl extra overrides. They are necessary since one may not modify an edge's metadata via add_edge!, for this kind of graphs
"""
    set_prop!(subgraph::S, prop, val) where {S <:AbstractSubGraph}
"""
function MetaGraphs.set_prop!(subgraph::S, prop, val) where {S<:AbstractSubGraph}
    return set_prop!(subgraph.graph, prop, val)
end
"""
    set_prop!(layer::L, v::MultilayerVertex, prop, val) where {L <: Layer} 
"""
function MetaGraphs.set_prop!(layer::L, v::MultilayerVertex, prop, val) where {L<:Layer}
    return set_prop!(layer.graph, get_v(layer, v), prop, val)
end
"""
set_prop!(subgraph::S, s::MultilayerVertex, d::MultilayerVertex,  prop, val) where {S <:AbstractSubGraph}
"""
function MetaGraphs.set_prop!(
    subgraph::S, s::MultilayerVertex, d::MultilayerVertex, prop, val
) where {S<:AbstractSubGraph}
    return set_prop!(subgraph.graph, get_v(subgraph, s), get_v(subgraph, d), prop, val)
end

"""
    get_prop(subgraph::S, prop) where {S <:AbstractSubGraph}
"""
function MetaGraphs.get_prop(subgraph::S, prop) where {S<:AbstractSubGraph}
    return get_prop(subgraph.graph, prop)
end
"""
"""
function MetaGraphs.get_prop(
    subgraph::S, v::MultilayerVertex, prop
) where {S<:AbstractSubGraph}
    return get_prop(subgraph.graph, get_v(subgraph, v), prop)
end
"""
    get_prop(subgraph::S, v::MultilayerVertex, prop) where {S <:AbstractSubGraph}
"""
function MetaGraphs.get_prop(
    subgraph::S, s::MultilayerVertex, d::MultilayerVertex, prop
) where {S<:AbstractSubGraph}
    return get_prop(subgraph.graph, get_v(subgraph, s), get_v(subgraph, d), prop)
end

function _add_edge!(
    g::AbstractMetaGraph{T},
    src::T,
    dst::T;
    weight::W=nothing,
    metadata::Union{Tuple,NamedTuple}=NamedTuple(),
) where {T<:Integer,W<:Union{<:Real,Nothing}}
    (isnothing(weight) || weight == 1.0) ||
        @warn "Trying to add a weighted edge to an unweighted graph of type $(typeof(g)). Weight $weight will be ignored."
    return add_edge!(
        g, src, dst, Dict(Symbol(pair.first) => pair.second for pair in pairs(metadata))
    ) #Symbol(string(key)) => value for (key, value) in pairs(metadata)
end

function _get_edge_weight(
    g::AbstractMetaGraph{T}, src::T, dst::T, weighttype::Type{U}
) where {T,U<:Real}
    return one(U)
end

function _get_edge_metadata(g::AbstractMetaGraph{T}, src::T, dst::T) where {T}
    return NamedTuple(props(g, src, dst))
end

function _set_metadata!(
    g::AbstractMetaGraph{T}, src::T, dst::T, metadata::NamedTuple
) where {T}
    return set_props!(g, src, dst, Dict(key => value for (key, value) in pairs(metadata)))
end

weights(g::AbstractMetaGraph{T}) where {T} = adjacency_matrix(g)

# is_directed(g::G) where {G <: AbstractMetaGraph} = is_directed(G)
