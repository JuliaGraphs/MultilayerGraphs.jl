# This file is the entry point for scripts each dedicated to make one Graphs.jl extension compatible with MultilayerGraphs.jl
include("graphs.jl")
include("simpleweightedgraphs.jl")
include("metagraphs.jl")
include("simplevaluegraphs.jl")

"""
    _edges(g::Union{Graphs.SimpleGraphs.AbstractSimpleGraph{T}, AbstractSimpleWeightedGraph{T}, AbstractMetaGraph{T}, SimpleValueGraphs.AbstractValGraph{T} }, weighttype::Type{U}) where {T,U}

Internal function. It serves as a unified interface between MultilayerGraphs.jl `edges` method and the homonymous methods of the other packages.
"""
function _edges(
    g::Union{
        Graphs.SimpleGraphs.AbstractSimpleGraph{T},
        AbstractSimpleWeightedGraph{T},
        AbstractMetaGraph{T},
        SimpleValueGraphs.AbstractValGraph{T},
    },
    weighttype::Type{U},
) where {T,U}
    return [
        (
            src(edge),
            dst(edge),
            _get_edge_weight(g, src(edge), dst(edge), weighttype),
            _get_edge_metadata(g, src(edge), dst(edge)),
        ) for edge in edges(g)
    ]
end
