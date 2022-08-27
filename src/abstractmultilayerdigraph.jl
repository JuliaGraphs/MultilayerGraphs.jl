"""
    AbstractMultilayerDiGraph{T,U} <: AbstractMultilayerGraph{T,U} 

Abstract type representing an undirected multilayer graph.
"""
abstract type AbstractMultilayerDiGraph{T,U} <: AbstractMultilayerGraph{T,U} end

"""
    add_layer!(mg::M,layer::L; interlayers_type = "multiplex") where { T, undef, M <: AbstractMultilayerDiGraph{T, G, U}, L <: Layer{T,U,G}}

Add layer `layer` to `mg`. Also add `Interlayer`s of type `interlayers_type` (can only be `"multiplex"`) between the new layer and all the other ones. 
"""
function add_layer!(
    mg::M, new_layer::L; new_default_interlayers_type="multiplex"
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerDiGraph{T,U},L<:Layer{T,U,G}}
    #Check that the layer is directed
    istrait(IsDirected{typeof(new_layer.graph)}) || throw(
        ErrorException(
            "The `new_layer`'s underlying graph $(new_layer.graph) is undirected, so it is not compatible with a `AbstractMultilayerDiGraph`.",
        ),
    )

    if new_default_interlayers_type == "multiplex"
        _add_layer!(mg, new_layer; new_default_interlayers_type=SimpleDiGraph{T})
    end
end

"""
    specify_interlayer!(mg::M, new_interlayer::In; symmetric_interlayer_name::Symbol) where { T, U, G<: AbstractGraph{T}, M <: AbstractMultilayerDiGraph{T, U}, In <: Interlayer{T,G}; IsDirected{M}}

Specify the interlayer `new_interlayer` as part of `mg`. The underlying graph of `new_interlayer` must be directed.
"""
function specify_interlayer!(
    mg::M,
    new_interlayer::In;
    symmetric_interlayer_name::String="interlayer_$(new_interlayer.layer_2)_$(new_interlayer.layer_1)",
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerDiGraph{T,U},In<:Interlayer{T,U,G}}
    istrait(IsDirected{typeof(new_interlayer.graph)}) || throw(
        ErrorException(
            "The `new_interlayer`'s underlying graph $(new_interlayer.graph) is undirected, so it is not compatible with a `AbstractMultilayerDiGraph`.",
        ),
    )
    return _specify_interlayer!(
        mg, new_interlayer; symmetric_interlayer_name=symmetric_interlayer_name
    )
end

# Graphs.jl's internals extra overrides
function Graphs.degree(
    mg::M, v::V
) where {T,M<:AbstractMultilayerDiGraph{T,<:Real},V<:MultilayerVertex{T}}
    return indegree(mg, v) + outdegree(mg, v)
end

"""
    add_edge!(mg::M, me::E) where { T, U <: Real, M <: AbstractMultilayerDiGraph{T,U}, E <: MultilayerEdge{MultilayerVertex{T},U}}

Add weighted edge `me` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, me::E
) where {T,U<:Real,M<:AbstractMultilayerDiGraph{T,U},E<:MultilayerEdge{MultilayerVertex{T},U}}
    return add_edge!(mg, src(me), dst(me), weight(me))
end

"""
    add_edge!(mg::M, me::E) where { T, U, M <: AbstractMultilayerDiGraph{T,U}, E <: MultilayerEdge{MultilayerVertex{T},Nothing}}

Add unweighted edge `me` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, me::E
) where {T,U,M<:AbstractMultilayerDiGraph{T,U},E<:MultilayerEdge{MultilayerVertex{T},Nothing}}
    return add_edge!(mg, src(me), dst(me))
end


"""
    is_directed(mg::M) where { M <: AbstractMultilayerDiGraph}

Returns `true` if `mg` is directed, `false` otherwise. 
"""
Graphs.is_directed(mg::M) where {M<:AbstractMultilayerDiGraph} = true

"""
    is_directed(mg::M) where { M <: Type{ <: AbstractMultilayerDiGraph}}

Returns `true` if `mg` is directed, `false` otherwise. 
"""
Graphs.is_directed(mg::M) where {M<:Type{<:AbstractMultilayerDiGraph}} = true



# TODO:
# it may be not well inferred
# function get_projected_monoplex_graph end #approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
"""
    get_projected_monoplex_graph(mg::M) where {M<: AbstractMultilayerDiGraph}

Get projected monoplex graph (i.e. the graph that as the same nodes as `mg` but the link between node `i` and `j` has weight equal to the sum of all edges weights between the various vertices representing `i` and `j` in `mg`, accounting for both layers and interlayers). See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function get_projected_monoplex_graph(mg::M) where {M<:AbstractMultilayerDiGraph}
    projected_monoplex_adjacency_matrix = dropdims(
        sum(mg.adjacency_tensor; dims=(3, 4)); dims=(3, 4)
    ) #::Matrix{eltype(mg.adjacency_tensor)}
    return SimpleWeightedDiGraph{M.parameters[1],M.parameters[2]}(
        projected_monoplex_adjacency_matrix
    )
end

# function get_overlay_monoplex_graph end #approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
"""
    get_overlay_monoplex_graph(mg::M) where {M<: AbstractMultilayerDiGraph}

Get overlay monoplex graph (i.e. the graph that has the same nodes as `mg` but the link between node `i` and `j` has weight equal to the sum of all edges weights between the various vertices representing `i` and `j` in `mg`, accounting for both layers and interlayers). See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function get_overlay_monoplex_graph(mg::M) where {M<:AbstractMultilayerDiGraph}
    projected_overlay_adjacency_matrix = sum([
        mg.adjacency_tensor[:, :, i, i] for i in 1:size(mg.adjacency_tensor, 3)
    ])
    return SimpleWeightedDiGraph{M.parameters[1],M.parameters[2]}(
        projected_overlay_adjacency_matrix
    )
end
