"""
    AbstractMultilayerGraph{T <: Integer, U <: Real} <: AbstractGraph{T}

An abstract type for multilayer graphs, it must contain the fields: `adjacency_tensor`, `layers`, `Interlayers`. It is a subtype of AbstractGraph and its concrete subtypes may extend Graphs.jl.
"""
abstract type AbstractMultilayerGraph{T<:Integer,U<:Real} <: AbstractGraph{T} end

"""
    Base.getproperty(mg::M, f::Symbol) where { M <: AbstractMultilayerGraph }
"""
function Base.getproperty(mg::M, f::Symbol) where {M<:AbstractMultilayerGraph}
    if f == :adjacency_tensor
        Base.getfield(mg, :adjacency_tensor)
    elseif f == :layers
        Base.getfield(mg, :layers)
    elseif f == :interlayers
        Base.getfield(mg, :interlayers)
    elseif f == :graphs
        return merge(mg.layers, mg.interlayers)
    elseif f == :layers_names
        return [layer.name for layer in values(mg.layers)]
    elseif f == :interlayers_names
        return [interlayer.name for layer in values(mg.interlayers)]
    else
        try
            collect(values(mg.layers))[findfirst(
                layer -> layer.name == f, collect(values(mg.layers))
            )]
        catch
            collect(values(mg.interlayers))[findfirst(
                interlayer -> interlayer.name == f, collect(values(mg.interlayers))
            )]
        end
    end
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

# TODO:
# write performant loops and adjacency_tensor structures following https://docs.huihoo.com/julia/0.3/manual/performance-tips/index.html#access-arrays-in-memory-order-along-columns

# ADD AND REMOVE (INTER)LAYERS
"""
    _add_layer!(mg::M,new_layer::L; new_default_interlayers_type::H) where { T, U, G <: AbstractGraph{T}, H <: AbstractGraph{T}, M <: AbstractMultilayerGraph{T, U}, L <: Layer{T,U,G}

Internal function. It is called by the `add_layer!` API functions, which needs to specify the default interlayer graph type (SimpleGraph for undirected multilayer graphs, SimpleDiGraph for directed multilayer graphs).
"""
function _add_layer!(
    mg::M, new_layer::L; new_default_interlayers_type::H
) where {
    T,
    U,
    M<:AbstractMultilayerGraph{T,U},
    G<:AbstractGraph{T},
    L<:Layer{T,U,G},
    H<:Type{<:AbstractGraph{T}},
}

    # Check that the new layer has a name different from all the existing ones
    new_layer.name ∉ mg.layers_names || throw(
        ErrorException(
            "The new layer has the same name as an existing layer within the multilayer graph. Layers' names must be unique.",
        ),
    )

    # Get number of nodes and of (already existing) layers
    n_nodes = size(mg.adjacency_tensor, 1)
    n_layers = size(mg.adjacency_tensor, 3)
    # Get index tuple for the new layer
    last_layer_tuple_idx =
        length(collect(keys(mg.layers))) == 0 ? (0, 0) : collect(keys(mg.layers))[end]
    new_layer_tuple_idx = last_layer_tuple_idx .+ 1
    # Add the new layer to mg.layers
    push!(mg.layers, new_layer_tuple_idx => new_layer)
    # Add components to the mg.adjacency_tensor
    new_adjacency_tensor = Array{U}(undef, n_nodes, n_nodes, n_layers + 1, n_layers + 1)
    new_adjacency_tensor[:, :, 1:n_layers, 1:n_layers] .= mg.adjacency_tensor
    new_adjacency_tensor[:, :, end, end] = adjacency_matrix(new_layer.graph)
    mg.adjacency_tensor = new_adjacency_tensor

    # Specify all the interlayers
    for layer in values(mg.layers)
        if layer.name != new_layer.name
            specify_interlayer!(
                mg,
                multiplex_interlayer(
                    nv(new_layer),
                    Symbol("interlayer_$(new_layer.name)_$(layer.name)"),
                    new_layer.name,
                    layer.name,
                    new_default_interlayers_type;
                    U=U,
                    forbidden_vertices=vcat(
                        new_layer.forbidden_vertices, layer.forbidden_vertices
                    ),
                    forbidden_edges=NTuple{2,MultilayerVertex{T}}[],
                ),
            )
        end
    end
end

"""
    _specify_interlayer!(mg::M, new_interlayer::In; symmetric_interlayer_name::String) where { T, U, G <: AbstractGraph{T}, M <: AbstractMultilayerGraph{T, U}, In <: Interlayer{T,U,G}}

Internal function. It is called by the `specify_interlayer!` API functions.
"""
function _specify_interlayer!(
    mg::M, new_interlayer::In; symmetric_interlayer_name::String
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerGraph{T,U},In<:Interlayer{T,U,G}}
    num_vertices = length(new_interlayer.layer_1_vertices)
    # Check that the `new_interlayer` has the same number of vertices as the `mg`
    num_vertices == size(mg.adjacency_tensor, 1) || throw(
        ErrorException(
            "The new interlayer does not have the same number of vertices as the other interlayers.",
        ),
    )
    # Check that the specified layers are part of `mg`
    all(in.([new_interlayer.layer_1, new_interlayer.layer_2], Ref(mg.layers_names))) ||
        throw(
            ErrorException(
                "The new interlayer connects two layers that are not (one or both) part of the multilayer graph. Make sure you spelled the `layer_1` and `layer_2` arguments of the `Interlayer` correctly.",
            ),
        )

    idxs_cart_1, layer_1 = get_layer(mg, new_interlayer.layer_1)
    idxs_cart_2, layer_2 = get_layer(mg, new_interlayer.layer_2)
    idxs_tup = (idxs_cart_1[1], idxs_cart_2[1])

    mg.interlayers[idxs_tup] = new_interlayer
    new_interlayer_adjm = adjacency_matrix(new_interlayer.graph)
    mg.adjacency_tensor[:, :, idxs_tup...] .= @views new_interlayer_adjm[
        1:num_vertices, (num_vertices + 1):end
    ]

    new_interlayer_symmetric = get_symmetric_interlayer(
        new_interlayer; symmetric_interlayer_name=symmetric_interlayer_name
    )
    mg.interlayers[reverse(idxs_tup)] = new_interlayer_symmetric
    new_interlayer_symmmetric_adjm = adjacency_matrix(new_interlayer_symmetric.graph)
    return mg.adjacency_tensor[:, :, reverse(idxs_tup)...] .= @views new_interlayer_symmmetric_adjm[
        1:num_vertices, (num_vertices + 1):end
    ] # (num_vertices+1):end, 1:num_vertices
end

"""
    get_layer(mg::M, layer_name::Symbol) where {T, U, M <: AbstractMultilayerGraph{T, U}}

Return `(idxs_cart,layer)` where `layer_name` is the name of the `Layer` being queried and `idxs_cart` is the `CartesianIndex` such that `mg.adjacency_tensor[:,:,idxs_cart]` is the adjacency matrix of such layer.
"""
function get_layer(mg::M, layer_name::Symbol) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    idxs_cart = nothing
    output_layer = nothing
    for (idxs_tup, layer) in collect(mg.layers)
        if layer.name == layer_name
            output_layer = layer
            idxs_cart = CartesianIndex(idxs_tup...)
            break
        end
    end

    if isnothing(output_layer)
        throw(
            ErrorException(
                "Layer $(layer_name) does not correspond to a Layer in the AbstractMultilayerGraph `mg`. Existing layers are $(tuple(layer.name for layer in values(mg.layers) ) )",
            ),
        )
    else
        return idxs_cart::CartesianIndex{2}, output_layer
    end
end

"""
    get_interlayer(mg::M, layer_1::Symbol, layer_2::Symbol) where {M <: AbstractMultilayerGraph}

Return `(idxs_cart,interlayer)` where `interlayer` is the `Interlayer` between `layer_1` and `layer_2` and `idxs_cart` is the `CartesianIndex` such that `mg.adjacency_tensor[:,:,idxs_cart]` is the adjacency matrix of `interlayer`.
"""
function get_interlayer(
    mg::M, layer_1::Symbol, layer_2::Symbol
) where {M<:AbstractMultilayerGraph}
    idxs_cart = nothing
    output_interlayer = nothing
    for (idxs_tup, interlayer) in collect(mg.interlayers)
        if all(interlayer.layers .== [layer_1, layer_2])
            output_interlayer = interlayer
            idxs_cart = CartesianIndex(idxs_tup...)
            break
        end
    end

    if isnothing(idxs_cart)
        throw(
            ErrorException(
                "At least one of the specified layers $((layer_1,layer_2)) does not correspond to a Layer in the AbstractMultilayerGraph `mg`. Existing layers are $(Tuple(layer.name for layer in values(mg.layers) ) )",
            ),
        )
    end

    return idxs_cart::CartesianIndex{2}, output_interlayer::Interlayer{M.parameters[1]}
end

"""
    get_subgraph(mg::M, names...) where {M <: AbstractMultilayerGraph}

If `length(names) == 1`, then this is equivalent to calling `get_layer(mg,names[1])` else if `length(names) == 2`, this is equivalent to calling `get_interlayer(mg,names[1],names[2])`. Throws an `ErrorException` in all other cases. 
"""
function get_subgraph(mg::M, names::Symbol...) where {M<:AbstractMultilayerGraph}
    names_unique = unique(names)
    if length(names_unique) == 1
        return get_layer(mg, names_unique...)
    elseif length(names_unique) == 2
        return get_interlayer(mg, names_unique...)
    else
        throw(
            ErrorException(
                "Function `get_subgraph` can be called with either 1 symbol (to get a layer and its adjacency matrix's `CartesianIndex` inside `mg.adjacency_tensor`) or 2 symbols (to get an interlayer and its adjacency matrix's `CartesianIndex` inside `mg.adjacency_tensor`). Called with $(length(names)) symbols.",
            ),
        )
    end
end

# Overloads that make AbstractMultilayerGraph an extension of Graphs.jl. These are all well-inferred .
"""
    edges(mg::M) where {M <: AbstractMultilayerDiGraph}

Return an iterator over all the edges of `mg`. The iterators first loops over all layers' edges (in the order they are given in `mg.layers`), then over all interlayers' edges (in the order they are given in `mg.interlayers`).
"""
function Graphs.edges(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    output = MultilayerEdge[]
    already_seen_interlayers_idxs = Set{Int64}[]
    for (idxs, graph) in collect(mg.graphs)
        if Set(idxs) ∉ already_seen_interlayers_idxs
            push!(output, collect(edges(graph))...)
            push!(already_seen_interlayers_idxs, Set(idxs))
        end
    end

    return output
end

"""
    Base.eltype(mg::M) where {M <: AbstractMultilayerGraph}

Return the vertex type of `mg`.
"""
Base.eltype(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} = T

"""
    adjm_eltype(mg::M) where {M <: AbstractMultilayerGraph}

Return the vertex type of `mg`.
"""
adjm_eltype(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}} = U

"""
    edgetype(mg::M) where {M <: AbstractMultilayerGraph}

Return the edge type for `mg`.
"""
function Graphs.edgetype(mg::M) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    return MultilayerEdge{MultilayerVertex{T},U}
end

"""
    has_edge(mg::M,s::T,d::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}

Return `true` if there is an edge between `s` and `d`, `false` otherwise. In case `s` and `d` belong to two different layers, but the link is found only in one of the two corresponding interlayers between them, throws an error.
"""
function Graphs.has_edge(
    mg::M, s::V, d::V
)::Bool where {T,U,M<:AbstractMultilayerGraph{T,U},V<:MultilayerVertex{T}}
    if s.layer == d.layer
        cart_index, layer = get_layer(mg, s.layer)
        return has_edge(layer, s, d)
    else
        cart_index, interlayer = get_interlayer(mg, s.layer, d.layer)
        cart_index_symmetric, interlayer_symmetric = get_interlayer(mg, d.layer, s.layer)
        hasEdge = has_edge(interlayer, s, d)
        hasEdge_symmetric = has_edge(interlayer_symmetric, s, d)

        hasEdge == hasEdge_symmetric || throw(
            ErrorException(
                "Found an edge in an interlayer which is not reflected in its symmetric. An error must have been occurred in past executions. If you didn't modify the multilayer's fields without using the provided API, please file a bug report.",
            ),
        )

        return hasEdge
    end
end

"""
    has_vertex(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}

Return `true` if `v` is a vertex of `mg`.
"""
function Graphs.has_vertex(
    mg::M, v::V
)::Bool where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}} # where { T <: Integer}
    graph = v.layer
    _has_vertex = eval(:(!($v in $mg.$graph.forbidden_vertices)))
    return _has_vertex
end

"""
    inneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}

Return the list of inneighbors of `v` within `mg`, looping first over all layers (in the order they are given in `mg.layers`), then over all interlayers (in the order they are given in `mg.interlayers`).
"""
function Graphs.inneighbors(
    mg::M, v::V
) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}}
    _inneighbors = MultilayerVertex{M.parameters[1]}[]
    layer_name = v.layer
    push!(_inneighbors, inneighbors(eval(:($mg.$layer_name)), v)...)

    for (idx_tup, interlayer) in collect(mg.interlayers)
        if interlayer.layer_1 == v.layer
            interlayer_name = interlayer.name
            interlayer_inneighbors = inneighbors(eval(:($mg.$interlayer_name)), v)
            push!(_inneighbors, interlayer_inneighbors...)
        end
    end
    return _inneighbors
end

"""
    ne(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of edges in `mg`.
"""
function Graphs.ne(mg::M) where {M<:AbstractMultilayerGraph}
    _ne = zero(Int64)
    for layer in values(mg.layers)
        _ne += ne(layer)::Int64
    end

    already_seen_interlayers_idxs = Set{Int64}[]
    for (idxs, interlayer) in collect(mg.interlayers)
        idxs_set = Set(idxs)
        if idxs_set ∉ already_seen_interlayers_idxs
            _ne += ne(interlayer)::Int64
            push!(already_seen_interlayers_idxs, idxs_set)
        else
            continue
        end
    end

    return _ne
end

"""
    nv(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of vertices in `mg`.
"""
function Graphs.nv(mg::M)::Int64 where {M<:AbstractMultilayerGraph}
    _nv = zero(Int64)
    for layer in values(mg.layers)
        _nv += nv(layer)
    end

    return _nv
end

"""
    nn(mg::M) where {M <: AbstractMultilayerGraph }

Return the number of nodes in `mg`.
"""
nn(mg::M) where {M<:AbstractMultilayerGraph} = length(nodes(mg))

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
    outneighbors(mg::M, v::T) where {M <: AbstractMultilayerGraph{T} } where { T <: Integer}

Return the list of outneighbors of `v` within `mg`, looping first over all layers (in the order they are given in `mg.layers`), then over all interlayers (in the order they are given in `mg.interlayers`).
"""
function Graphs.outneighbors(
    mg::M, v::V
) where {M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}} where {T}
    _outneighbors = MultilayerVertex{M.parameters[1]}[]
    layer_name = v.layer
    push!(_outneighbors, outneighbors(eval(:($mg.$layer_name)), v)...)

    for (idx_tup, interlayer) in collect(mg.interlayers)
        if interlayer.layer_1 == v.layer
            interlayer_name = interlayer.name
            interlayer_outneighbors = outneighbors(eval(:($mg.$interlayer_name)), v)
            push!(_outneighbors, interlayer_outneighbors...)
        end
    end

    return _outneighbors
end

"""
    vertices(mg::M) where {M <: AbstractMultilayerGraph{ <: Integer, <: AbstractSimpleGraph}}

Return the collection of the vertices of `mg`.
"""
function Graphs.vertices(mg::M) where {M<:AbstractMultilayerGraph}
    output = MultilayerVertex{M.parameters[1]}[]
    for layer in values(mg.layers)
        push!(output, vertices(layer)...)
    end

    return output
end

"""
    nodes(mg::M) where {M <: AbstractMultilayerGraph}

Return the nodes of the AbstractMultilayerGraph `mg`.
"""
function nodes(mg::M) where {M<:AbstractMultilayerGraph}
    return [vertex.node for vertex in vertices(mg.layers[(1, 1)])]
end

# Graphs.jl's internals extra overrides
"""
    indegree(mg::M, v::V) where {T, M <: AbstractMultilayerGraph{T, <: Real}, V <: MultilayerVertex{T}}

Get the indegree of vertex `v` in `mg`.
"""
function Graphs.indegree(
    mg::M, v::V
) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}}
    return length(inneighbors(mg, v))
end
"""
    indegree(mg::M, vs::AbstractVector{V} = vertices(mg)) where {T,M <: AbstractMultilayerGraph{T, <: Real}, V <: MultilayerVertex{T}} 

Get the vector of indegrees of vertices `vs` in `mg`.
"""
function Graphs.indegree(
    mg::M, vs::AbstractVector{V}=vertices(mg)
) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}}
    return [indegree(mg, x) for x in vs]
end

"""
    outdegree(mg::M, v::V) where {T, M <: AbstractMultilayerGraph{T, <: Real}, V <: MultilayerVertex{T}}

Get the outdegree of vertex `v` in `mg`.
"""
function Graphs.outdegree(
    mg::M, v::V
) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}}
    return length(outneighbors(mg, v))
end
"""
outdegree(mg::M, vs::AbstractVector{V} = vertices(mg)) where {T,M <: AbstractMultilayerGraph{T, <: Real}, V <: MultilayerVertex{T}} 

Get the vector of outdegrees of vertices `vs` in `mg`.
"""
function Graphs.outdegree(
    mg::M, vs::AbstractVector{V}=vertices(mg)
) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}}
    return [outdegree(mg, x) for x in vs]
end

"""
    degree(mg::M, v::V) where {T, M <: AbstractMultilayerGraph{T, <: Real}, V <: MultilayerVertex{T}}

Get the degree of vertices `vs` in `mg`.
"""
function Graphs.degree(
    mg::M, vs::AbstractVector{V}=vertices(mg)
) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}}
    return [degree(mg, x) for x in vs]
end

"""
    neighbors(mg::M, v::V) where {T, M <: AbstractMultilayerGraph{T, <: Real}, V <: MultilayerVertex{T}}

Get the neighbors of vertices `vs` in `mg`. Reduces to `outneighbors` for both directed and undirected multilayer graphs.
"""
function Graphs.neighbors(
    mg::M, v::V
) where {T,M<:AbstractMultilayerGraph{T,<:Real},V<:MultilayerVertex{T}}
    return outneighbors(mg, v)
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

Return the complete multilayer global clustering coefficient, equal to the ratio of realized triplets over all possible triplets, including those whose every or some edges belong to interlayers, normalized by `norm_factor`. If `norm_factor == :max`, then the ratio is normalized by `maximum(mg.adjacency_tensor)`. This function overrides Graphs.jl's `global_clustering_coefficient`, since the latter does not consider cliques where two nodes are the same node but in different layers/interlayers. See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function multilayer_global_clustering_coefficient(
    mg::M, norm_factor::Union{Float64,Symbol}=:max
) where {M<:AbstractMultilayerGraph}
    _normalization_inverse = 1.0
    if norm_factor == :max
        _normalization_inverse = 1.0 / maximum(mg.adjacency_tensor)
    end

    A_right = mg.adjacency_tensor .- get_diagonal_elements(mg.adjacency_tensor)

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
        ones(size(mg.adjacency_tensor)...) .-
        multilayer_kronecker_delta(size(mg.adjacency_tensor))
    den = ein"ijkm,jnmo,niok ->"(A_right, F, A_right)[]

    return _normalization_inverse * (num / den)
end

"""
    multilayer_weighted_global_clustering_coefficient(mg::M, norm_factor::Union{Float64, Symbol} = :max) where {M <: AbstractMultilayerGraph}

Return the complete multilayer global clustering coefficient, equal to the ratio of realized triplets over all possible triplets, including those whose every or some edges belong to interlayers, normalized by `norm_factor`. Each triplets contributes for `w[1]` if all of its vertices are in one layer, `w[2]` if its vertices span two layers, and `w[3]` if they span 3 layers. If `norm_factor == :max`, then the ratio is normalized by `maximum(mg.adjacency_tensor)`. This function overrides Graphs.jl's `global_clustering_coefficient`, since the latter does not consider cliques where two nodes are the same node but in different layers/interlayers. See [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function multilayer_weighted_global_clustering_coefficient(
    mg::M, w::Vector{Float64}, norm_factor::Union{Float64,Symbol}=:max
) where {M<:AbstractMultilayerGraph} #This is well defined for both weighted and unweighted multilayer graphs
    sum(w) == 1 ||
        throw(ErrorException("Weight vector `w` does not sum to 1. Found $(sum(w))."))
    _normalization_inverse = 1.0
    if norm_factor == :max
        _normalization_inverse = 1.0 / maximum(mg.adjacency_tensor)
    end

    A_right = mg.adjacency_tensor .- get_diagonal_elements(mg.adjacency_tensor)
    num_layers = size(mg.adjacency_tensor, 3)

    num = ein"ijkm,jnmo,niok,skmo,s ->"(A_right, A_right, A_right, δ_Ω(num_layers), w)[]

    F =
        ones(size(mg.adjacency_tensor)...) .-
        multilayer_kronecker_delta(size(mg.adjacency_tensor))

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
    _normalization_inverse = 1.0
    if norm_factor == :max
        _normalization_inverse =
            1.0 / (maximum(ein"ijkl->ij"(mg.adjacency_tensor)) / length(mg.layers))
        # Check that we are using OMEinsum correctly.
        # @assert all(ein"ijkl->ij"(mg.adjacency_tensor) .== dropdims(sum(mg.adjacency_tensor, dims = (3,4)), dims = (3,4)))
    end

    num = ein"ij,jm,mi ->"(
        ein"ijkm->ij"(mg.adjacency_tensor),
        ein"ijkm->ij"(mg.adjacency_tensor),
        ein"ijkm->ij"(mg.adjacency_tensor),
    )[]

    F =
        ones(size(mg.adjacency_tensor)...) -
        multilayer_kronecker_delta(size(mg.adjacency_tensor))

    den = ein"ij,jm,mi ->"(
        ein"ijkm->ij"(mg.adjacency_tensor),
        ein"ijkm->ij"(F),
        ein"ijkm->ij"(mg.adjacency_tensor),
    )[]

    return _normalization_inverse * (num / den)
end

"""
    eigenvector_centrality(mg::M; norm::String = "1", tol::Float64 = 1e-6, maxiter::Int64 = 2000) where {T, U, M <: AbstractMultilayerGraph{T, U}}

Calculate the eigenvector centrality of `mg` via an iterative algorithm. The `norm` parameter may be `"1"` or `"n"`,  and respectively the eigenvector centrality will be normalized to 1 or further divided by the number of nodes of `mg`. The `tol` parameter terminates the approximation when two consecutive iteration differ by no more than  `tol`. The `maxiters` parameter terminates the algorithm when it goes beyond `maxiters` iterations.

The returned values are: the eigenvector centrality and the relative error at each algorithm iteration, that is, the summed absolute values of the componentwise differences between the centrality computed at the current iteration minus the centrality computed at the previous iteration.


Note: in the limit case of a monoplex graph, this function outputs a eigenvector centrality vector that is a multiple of the one outputted by Graphs.jl's `eigenvector_centrality`. It is currently under investigation.
"""
function Graphs.eigenvector_centrality(
    mg::M; norm::String="1", tol::Float64=1e-6, maxiter::Int64=2000
) where {T,U,M<:AbstractMultilayerGraph{T,U}}
    num_nodes = length(nodes(mg))
    X = ones(Float64, num_nodes, length(mg.layers))

    err = 1.0
    errs = Float64[]
    iter = 0
    while err > tol && iter < maxiter
        new_X = ein"ijkm,ik -> jm"(mg.adjacency_tensor, X)
        if norm == "1"
            new_X = new_X ./ sum(new_X)
        elseif norm == "n"
            new_X = new_X ./ (sum(new_X) / num_nodes)
        end
        err = sum(abs.(X .- new_X))
        push!(errs, err)
        X .= new_X
        iter += 1
    end
    return X, errs
end

#= function eigenvector_centrality_2(mg::M; norm::String = "1", tol::Float64 = 1e-6) where {T, U, M <: AbstractMultilayerGraph{T, U}}

    eigvals, eigvects = tensoreig(mg.adjacency_tensor,[2,4],[1,3])

    eigvects = [eigvects[:,:,i] for i in 1:size(eigvects,3)]

    println("eigvects = ")
    display(eigvects)

    eigvals_eigvects_sorted = sort(collect(zip(abs.(eigvals),eigvects)); by=first)

    println("eigvals_sorted = $([couple[1] for couple in sort(collect(zip(abs.(eigvals),eigvects)); by=first)])" )

    println("eigvects_sorted = ")
    display(eigvects_sorted)

end =#

"""
    modularity(mg::M, c::Matrix{Int64}; null_model::Union{String,Array{U,4}} = "degree") where {T, U, M <: AbstractMultilayerGraph{T,U}}

Calculate the modularity of `mg`, as shown in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).
"""
function Graphs.modularity(
    mg::M, c::Matrix{Int64}; null_model::Union{String,Array{U,4}}="degree"
) where {T,U,M<:AbstractMultilayerGraph{T,U}}

    # Check that c has the correct size
    n_nodes = length(nodes(mg))
    n_layers = length(mg.layers)
    size(c) == (n_nodes, n_layers) || throw(
        ErrorException(
            "The size of the community matrix does not match (nn(mg),length(mg.layers)), found $(size(c)) and $((nn(mg),length(mg.layers)))",
        ),
    )

    # Build S
    n_communities = length(unique(c))
    S = Array{Bool}(undef, n_nodes, n_layers, n_communities)
    for (i, community) in enumerate(unique(c))
        S[:, :, i] .== (c .== community)
    end

    # Build P
    P = Array{Float64}(undef, size(mg.adjacency_tensor)) # similar(mg.adjacency_tensor, element_type = Float64)
    tot_links = length(edges(mg))
    if typeof(null_model) != String && size(null_model) == size(P)
        P .= null_model
    elseif typeof(null_model) != String && size(null_model) != size(P)
        throw(
            ErrorException(
                "size of `null_model` does not match the size of the adjacency tensor. Got $(size(null_model)) and $(size(P)) respectively .",
            ),
        )
    elseif null_model == "degree"
        for cart_idx in CartesianIndices(P)
            layer_1_idx = cart_idx[3]
            layer_2_idx = cart_idx[4]
            P[cart_idx] =
                (
                    degree(
                        mg,
                        MultilayerVertex(
                            cart_idx[1], mg.layers[(layer_1_idx, layer_1_idx)].name
                        ),
                    ) * degree(
                        mg,
                        MultilayerVertex(
                            cart_idx[2], mg.layers[(layer_2_idx, layer_2_idx)].name
                        ),
                    )
                ) / (2 * tot_links - 1)
        end
    else
        throw(ErrorException("$null_model not implemented."))
    end

    # Build B
    B = mg.adjacency_tensor .- P

    # Build K
    K = ein"ijkm,jimk -> "(mg.adjacency_tensor, ones(T, size(mg.adjacency_tensor)...))[]

    return (1 / K) * ein"ija,ikjm,kma->"(S, B, S)[]
end
