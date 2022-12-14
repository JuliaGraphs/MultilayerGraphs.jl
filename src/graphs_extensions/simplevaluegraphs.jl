function __add_vertex!(g::SimpleValueGraphs.AbstractValGraph{T}; metadata::Union{Tuple,NamedTuple} = NamedTuple()) where {T <: Integer} 
    if isempty(metadata) 
        add_vertex!(g)
    else
        add_vertex!(g, metadata)
    end
end

_get_vertex_metadata(g::SimpleValueGraphs.AbstractValGraph{T}, vertex::T) where T = get_vertexval(g, vertex, :)

# SimpleValueGraphs.jl extra overrides. They are necessary since add_vertex!, conversely to add_edge!, returns false immediately if the node already exists in the layer (instead of proceeding to possibly modify the (weight and) metadata as add_edge! would do), so we need to extend the functions that allow for modifying vertices
"""
    set_vertexval!(layer::L, n::MultilayerVertex, prop::Symbol, val) where {T,U,G <: SimpleValueGraphs.AbstractValGraph{T}, L <: Layer{T,U,G}} 
"""
SimpleValueGraphs.set_vertexval!(layer::L, n::MultilayerVertex, prop::Symbol, val) where {T,U,G <: SimpleValueGraphs.AbstractValGraph{T}, L <: Layer{T,U,G}} = set_vertexval!(layer.graph, layer.v_V_associations(n), prop, val) 

function _add_edge!(g::SimpleValueGraphs.AbstractValGraph{T}, src::T, dst::T; weight::W = nothing, metadata::Union{Tuple,NamedTuple} = NamedTuple()) where {T <: Integer, W<: Union{<: Real, Nothing}}

    (isnothing(weight) || weight == 1.0) || @warn "Trying to add a weighted edge to an unweighted graph of type $(typeof(g)). Weight $weight will be ignored."

    if has_edge(g, src, dst)
        if metadata isa NamedTuple
            for (key,value) in pairs(metadata)
                set_edgeval!(g,src,dst,key,value)
            end
            return true
        elseif metadata isa Tuple
            for (i,value) in enumerate(metadata)
                set_edgeval!(g,src,dst,i,value)
            end
            return true
        end
    else
        add_edge!(g, src, dst, metadata)
    end
end

_get_edge_weight(g::SimpleValueGraphs.AbstractValGraph{T}, src::T, dst::T, weighttype::Type{U} ) where {T, U <: Real} = one(U)

_get_edge_metadata(g::SimpleValueGraphs.AbstractValGraph{T}, src::T, dst::T ) where T   = get_edgeval(g, src, dst, :) 

function _set_metadata!(g::SimpleValueGraphs.AbstractValGraph{T}, src::T, dst::T, metadata::Union{Tuple, NamedTuple}) where T
    
    if metadata isa NamedTuple
        for (key,value) in pairs(metadata)
            set_edgeval!(g,src,dst,key,value)
        end
        return true
    elseif metadata isa Tuple
        for (i,value) in enumerate(metadata)
            set_edgeval!(g,src,dst,i,value)
        end
        return true
    end
end

Graphs.weights(g::SimpleValueGraphs.AbstractValGraph{T}) where T = adjacency_matrix(g)