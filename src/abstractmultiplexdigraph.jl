"""
"""
abstract type AbstractMultiplexDiGraph{T,U} <: AbstractMultilayerDiGraph{T,U} end


# Edges
"""
    add_edge_specialized!(mg::M, me::E) where {T,U, M <: AbstractMultiplexDiGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add MultilayerEdge `me` to the directed multiplex graph `mg`. Return true if succeeds, false otherwise.
"""
function add_edge_specialized!(
    mg::M, me::E
) where {T,U,M<:AbstractMultiplexDiGraph{T,U},E<:MultilayerEdge{<:Union{U,Nothing}}}
    layer(_src) == layer(_dst) || throw(ArgumentError("Cannot add an edge between two `MultilayerVertex`s that belong to two different layers in a multiplex graph"))
    add_edge_directed!(mg, me)
end