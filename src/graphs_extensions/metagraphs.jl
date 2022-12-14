#= """
    MetaGraph{T,U}()

Empty MetaGraph with vertex type `T` and adjacency matrix eltype `U`. The underlying graph is a SimpleGraph. 
"""
MetaGraphs.MetaGraph{T,U}() where {T,U} =  MetaGraph{T,U}(SimpleGraph{T}())

"""
    MetaDiGraph{T,U}()

Empty MetaDiGraph with vertex type `T` and adjacency matrix eltype `U`. The underlying graph is a SimpleDiGraph. 
"""
MetaGraphs.MetaDiGraph{T,U}() where {T,U} =  MetaDiGraph{T,U}(SimpleDiGraph{T}())

"""
    MetaGraph{T,U}(n_vertices::Integer, n_edges::Integer)

Random MetaGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. the underlying graph is a SimpleGraph. 
"""
MetaGraphs.MetaGraph{T,U}(n_vertices::Integer, n_edges::Integer) where {T,U} = MetaGraph{T,U}(SimpleGraph(n_vertices, n_edges))

"""
    MetaGraph(n_vertices::Integer, n_edges::Integer; T = Int64, U = Float64)

Random MetaGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. the underlying graph is a SimpleGraph. 
"""
MetaGraphs.MetaGraph(n_vertices::Integer, n_edges::Integer; T = Int64, U = Float64) = MetaGraph{T,U}(SimpleGraph(n_vertices, n_edges))


"""
    MetaDiGraph{T,U}(n_vertices::Integer, n_edges::Integer)

Randoms MetaDiGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. the underlying graph is a SimpleDiGraph. 
"""
MetaGraphs.MetaDiGraph{T,U}(n_vertices::Integer, n_edges::Integer) where {T,U} =  MetaDiGraph{T,U}(SimpleDiGraph(n_vertices, n_edges))

"""
    MetaGraph(n_vertices::Integer, n_edges::Integer; T = Int64, U = Float64)

Randoms MetaGraph with `n_vertices` vertices and `n_edges` edges, vertex type `T` and adjacency matrix eltype `U`. the underlying graph is a SimpleDiGraph. 
"""
MetaGraphs.MetaDiGraph(n_vertices::Integer, n_edges::Integer; T = Int64, U = Float64) =  MetaDiGraph{T,U}(SimpleDiGraph(n_vertices, n_edges))


"""
    MetaGraph{T,U}(adjm::Matrix) where {T,U}

MetaGraph with adjacency matrix `adjm`, vertex type `T` and adjacency matrix eltype `U`. The underlying graph is a SimpleGraph. 
"""
MetaGraphs.MetaGraph{T,U}(adjm::Union{Matrix,SparseMatrixCSC}) where {T,U} =  MetaGraph{T,U}(SimpleGraph(adjm))

"""
    MetaGraph{T,U}(n_vertices::Integer, n_edges::Integer)

MetaDiGraph with adjacency matrix `adjm`, vertex type `T` and adjacency matrix eltype `U`. The underlying graph is a SimpleDiGraph. 
"""
MetaGraphs.MetaDiGraph{T,U}(adjm::Union{Matrix,SparseMatrixCSC}) where {T,U} =  MetaDiGraph{T,U}(SimpleDiGraph(adjm)) =#

__add_vertex!(g::AbstractMetaGraph{T}; metadata::Union{Tuple,NamedTuple} = NamedTuple()) where {T <: Integer} = add_vertex!(g, Dict(Symbol(pair.first) => pair.second for pair in pairs(metadata)))

_get_vertex_metadata(g::AbstractMetaGraph{T}, vertex::T) where T = NamedTuple(props(g,vertex))

# MetaGraphs.jl extra overrides. They are necessary since one may not modify an edge's metadata via add_edge!, for this kind of graphs
"""
    set_prop!(subgraph::S, prop, val) where {S <:AbstractSubGraph}
"""
MetaGraphs.set_prop!(subgraph::S, prop, val) where {S <:AbstractSubGraph} = set_prop!(subgraph.graph, prop, val)
"""
    set_prop!(layer::L, v::MultilayerVertex, prop, val) where {L <: Layer} 
"""
MetaGraphs.set_prop!(layer::L, v::MultilayerVertex, prop, val) where {L <: Layer} = set_prop!(layer.graph, get_v(layer,v), prop, val)
"""
set_prop!(subgraph::S, s::MultilayerVertex, d::MultilayerVertex,  prop, val) where {S <:AbstractSubGraph}
"""
MetaGraphs.set_prop!(subgraph::S, s::MultilayerVertex, d::MultilayerVertex,  prop, val) where {S <:AbstractSubGraph} = set_prop!(subgraph.graph, get_v(subgraph,s), get_v(subgraph,d), prop, val)

"""
    get_prop(subgraph::S, prop) where {S <:AbstractSubGraph}
"""
MetaGraphs.get_prop(subgraph::S, prop) where {S <:AbstractSubGraph}  = get_prop(subgraph.graph, prop)
"""
"""
MetaGraphs.get_prop(subgraph::S, v::MultilayerVertex, prop) where {S <:AbstractSubGraph} = get_prop(subgraph.graph, get_v(subgraph,v), prop)
"""
    get_prop(subgraph::S, v::MultilayerVertex, prop) where {S <:AbstractSubGraph}
"""
MetaGraphs.get_prop(subgraph::S, s::MultilayerVertex, d::MultilayerVertex, prop) where {S <:AbstractSubGraph} = get_prop(subgraph.graph, get_v(subgraph,s), get_v(subgraph,d), prop,)


function _add_edge!(g::AbstractMetaGraph{T}, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple,NamedTuple} = NamedTuple()) where {T <: Integer, W<: Union{<: Real, Nothing}}
    (isnothing(weight) || weight == 1.0) || @warn "Trying to add a weighted edge to an unweighted graph of type $(typeof(g)). Weight $weight will be ignored."
    add_edge!(g, src, dst, Dict(key => value for (key,value) in pairs(metadata) ))
end

_get_edge_weight(g::AbstractMetaGraph{T}, src::T, dst::T, weighttype::Type{U} ) where {T, U <: Real} = one(U) 

_get_edge_metadata(g::AbstractMetaGraph{T}, src::T, dst::T ) where T = NamedTuple(props(g,src,dst)) 

_set_metadata!(g::AbstractMetaGraph{T}, src::T, dst::T, metadata::NamedTuple) where T = set_props!(g, src, dst, Dict(key => value for (key,value) in pairs(metadata)))

Graphs.weights(g::AbstractMetaGraph{T}) where T = adjacency_matrix(g)