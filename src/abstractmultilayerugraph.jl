
"""
    AbstractMultilayerUGraph{T,U} <: AbstractAbstractMultilayerUGraph{T,U} 

Abstract type representing an undirected multilayer graph.
"""
abstract type AbstractMultilayerUGraph{T,U} <: AbstractAbstractMultilayerUGraph{T,U} end



"""
    add_layer!(mg::M,layer::L; interlayers_type = "multiplex") where { T, U, G<: AbstractGraph{T}, M <: AbstractMultilayerUGraph{T, U}, L <: Layer{T,U,G}}

Add layer `layer` to `mg`. Also add `Interlayer`s of type `interlayers_type` (can only be `"multiplex"`) between the new layer and all the other ones. 
"""
function add_layer!(
    mg::M, new_layer::L; new_default_interlayers_type="multiplex"
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerUGraph{T,U},L<:Layer{T,U,G}}
    # Check that the layer is directed
    !istrait(IsDirected{typeof(new_layer.graph)}) || throw(
        ErrorException(
            "The `new_layer`'s underlying graph $(new_layer.graph) is directed, so it is not compatible with a `AbstractMultilayerUGraph`.",
        ),
    )

    if new_default_interlayers_type == "multiplex"
        _add_layer!(mg, new_layer; new_default_interlayers_type=SimpleGraph{T})
    end
end


# Graphs.jl's internals extra overrides
function Graphs.degree(
    mg::M, v::V
) where {T,M<:AbstractMultilayerUGraph{T,<:Real},V<:MultilayerVertex{T}}
    return indegree(mg, v)
end


"""
    add_edge!(mg::M, me::E) where { T, U <: Real, M <: AbstractMultilayerUGraph{T,U}, E <: MultilayerEdge{MultilayerVertex{T},U}}

Add weighted edge `me` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, me::E
) where {T,U<:Real,M<:AbstractMultilayerUGraph{T,U},E<:MultilayerEdge{MultilayerVertex{T},U}}
    return add_edge!(mg, src(me), dst(me), weight(me))
end

"""
    add_edge!(mg::M, me::E) where { T, U, M <: AbstractMultilayerUGraph{T,U}, E <: MultilayerEdge{MultilayerVertex{T},Nothing}}

Add unweighted edge `me` to `mg`.
"""
function Graphs.add_edge!(
    mg::M, me::E
) where {T,U,M<:AbstractMultilayerUGraph{T,U},E<:MultilayerEdge{MultilayerVertex{T},Nothing}}
    return add_edge!(mg, src(me), dst(me))
end


"""
    is_directed(m::M) where { M <: AbstractMultilayerUGraph}

Returns `true` if `m` is directed, `false` otherwise. 
"""
Graphs.is_directed(m::M) where {M<:AbstractMultilayerUGraph} = false

"""
    is_directed(m::M) where { M <: Type{ <: AbstractMultilayerUGraph}}

Returns `true` if `m` is directed, `false` otherwise. 
"""
Graphs.is_directed(m::M) where {M<:Type{<:AbstractMultilayerUGraph}} = false


# Multilayer-specific functions
# TODO:
# it may be not well inferred
# function get_projected_monoplex_graph end #approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
"""
    get_projected_monoplex_graph(mg::M) where {M<: AbstractMultilayerUGraph}

Get projected monoplex graph (i.e. that graph that as the same nodes as `mg` but the link between node `i` and `j` has weight equal to the sum of all edges weights between the various vertices representing `i` and `j` in `mg`, accounting for both layers and interlayers). See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function get_projected_monoplex_graph(mg::M) where {M<:AbstractMultilayerUGraph}
    projected_monoplex_adjacency_matrix = dropdims(
        sum(mg.adjacency_tensor; dims=(3, 4)); dims=(3, 4)
    )

    isapproxsymmetric(projected_monoplex_adjacency_matrix) ||
        throw(ErrorException("Adjacency / distance matrices must be symmetric"))

    symmetric_projected_monoplex_adjacency_matrix = Symmetric(
        projected_monoplex_adjacency_matrix
    )

    return SimpleWeightedGraph{M.parameters[1],M.parameters[2]}(
        symmetric_projected_monoplex_adjacency_matrix
    )
end

# function get_overlay_monoplex_graph end #approach taken from https://github.com/JuliaGraphs/Graphs.jl/blob/7152d540631219fd51c43ab761ec96f12c27680e/src/core.jl#L124
"""
    get_overlay_monoplex_graph(mg::M) where {M<: AbstractMultilayerUGraph}


Get overlay monoplex graph (i.e. the graph that has the same nodes as `mg` but the link between node `i` and `j` has weight equal to the sum of all edges weights between the various vertices representing `i` and `j` in `mg`, accounting for only within-layer edges). See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function get_overlay_monoplex_graph(mg::M) where {M<:AbstractMultilayerUGraph}
    projected_overlay_adjacency_matrix = sum([
        mg.adjacency_tensor[:, :, i, i] for i in 1:size(mg.adjacency_tensor, 3)
    ])
    return SimpleWeightedGraph{M.parameters[1],M.parameters[2]}(
        projected_overlay_adjacency_matrix
    )
end




"""
    von_neumann_entropy(mg::M) where {T,U,  M <: AbstractMultilayerUGraph{T, U}}

Compute the Von Neumann entropy of `mg`, according to [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Only for undirected multilayer graphs.
"""
function von_neumann_entropy(mg::M) where {T,U,M<:AbstractMultilayerUGraph{T,U}}
    num_nodes = length(nodes(mg))
    num_layers = length(mg.layers)

    # multistrength tensor
    Δ = ein"ijkm,ik,nj,om -> njom"(
        mg.adjacency_tensor,
        ones(T, size(mg.adjacency_tensor)[[1, 3]]),
        δ(num_nodes),
        δ(num_layers),
    )

    # trace of Δ
    tr_Δ = ein"iikk->"(Δ)[]

    # multilayer laplacian tensor
    L = Δ .- mg.adjacency_tensor

    # multilayer density tensor
    ρ = (1.0 / tr_Δ) .* L

    eigvals, eigvects = tensoreig(ρ, [2, 4], [1, 3])

    #=     # Check that we are calculating the right eigenvalues
        lhs = ein"ijkm,ik -> jm"(ρ,eigvects[:,:,1])
        rhs = eigvals[1].*eigvects[:,:,1]
        @assert all(lhs .≈ rhs)
        # Indeed we are =#

    Λ = get_diagonal_adjacency_tensor(eigvals, size(mg.adjacency_tensor))

    # Correct for machine precision
    Λ[Λ .< eps()] .= 0.0

    log2Λ = log2.(Λ)

    log2Λ[isinf.(log2Λ)] .= 0

    return -ein"ijkm,jimk ->"(Λ, log2Λ)[]
end



