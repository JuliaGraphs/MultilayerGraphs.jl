# CUSTOM TYPE/CLASS
# This file defines the custom type `InterLayer` and straightforwardly makes it compatible with the Graphs.jl ecosystem. The reason to have this custom type is to have a way not to break other people's code when modification to interlayers are required.

"""
    mutable struct InterLayer{G <: AbstractGraph}

Represents an interlayer in a multilayer graph. 
"""
mutable struct InterLayer{G <: AbstractGraph}
    name::Symbol
    layer_1::Symbol
    layer_2::Symbol
    graph::G
    
    function InterLayer(name::Symbol,layer_1::Symbol, layer_2::Symbol, graph::G ) where {G <: AbstractGraph}
        # Check that the linke dlayers are different, esle error
        @assert layer_1 != layer_2

        return new{G}(name, layer_1, layer_2, graph)
    end

end

# Exted all Graphs.jl required methods (https://juliagraphs.org/Graphs.jl/dev/ecosystem/interface/)
Graphs.edges(interlayer::I) where { I <: InterLayer} = edges(I.graph)

Base.eltype(interlayer::I) where { I <: InterLayer}  = Base.eltype(I.graph)

Graphs.edgetype(interlayer::I) where { I <: InterLayer} = edgetype(I.graph)

Graphs.has_edge(interlayer::I, s::Integer, d::Integer) where { I <: InterLayer} = has_edge(I.graph, s, d)

Graphs.has_vertex(interlayer::I, v::Integer) where { I <: InterLayer} = has_vertex(I.graph, v)

Graphs.inneighbors(interlayer::I, v::Integer) where { I <: InterLayer} = inneighbors(I.graph, v)

Graphs.ne(interlayer::I) where { I <: InterLayer} = ne(I.graph)

Graphs.nv(interlayer::I) where { I <: InterLayer} = nv(interlayer.graph)

Graphs.outneighbors(interlayer::I, v::Integer) where { I <: InterLayer} = outneighbors(interlayer.graph)

Graphs.vertices(interlayer::I) where { I <: InterLayer} = vertices(interlayer.graph)

Graphs.is_directed(interlayer::I) where { I <: InterLayer} = is_directed(interlayer.graph)

Graphs.is_directed(interlayer_type::Type{I}) where {I <: InterLayer} = is_directed(interlayer_type.parameters[1])