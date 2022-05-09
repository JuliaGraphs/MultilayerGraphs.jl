# CUSTOM TYPE/CLASS
# This file ddefines the custom type `Layer` and straightforwardly makes it compatible with the Graphs.jl ecosystem. The reason ro have thsi custom type is to have a way not to break other people's code when modification to layers are required.


"""
    mutable struct Layer{G <: AbstractGraph}

Represents a layer in a multilayer graph. 
"""
mutable struct Layer{G <: AbstractGraph}
    name::Symbol
    graph::G
end

# Exted all Graphs.jl required methods (https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/)
Graphs.edges(l::L) where { L <: Layer} = edges(l.graph)

Base.eltype(l::L) where { L <: Layer}  = Base.eltype(l.graph)

Graphs.edgetype(l::L) where { L <: Layer} = edgetype(l.graph)

Graphs.has_edge(l::L, s::Integer, d::Integer) where { L <: Layer} = has_edge(l.graph, s, d)

Graphs.has_vertex(l::L, v::Integer) where { L <: Layer} = has_vertex(l.graph, v)

Graphs.inneighbors(l::L, v::Integer) where { L <: Layer} = inneighbors(l.graph, v)

Graphs.ne(l::L) where { L <: Layer} = ne(l.graph)

Graphs.nv(l::L) where { L <: Layer} = nv(l.graph)

Graphs.outneighbors(l::L, v::Integer) where { L <: Layer} = outneighbors(l.graph)

Graphs.vertices(l::L) where { L <: Layer} = vertices(l.graph)

Graphs.is_directed(l::L) where { L <: Layer} = is_directed(l.graph)

Graphs.is_directed(::Type{L}) where {L <: Layer} = is_directed(typeof(L).parameters[1])