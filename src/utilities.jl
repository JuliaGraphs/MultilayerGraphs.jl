# Base extensions
#= """
    getindex(od::O, key::Int64) where {O <: OrderedDict}

Return the `key`th pair of `od`, following `od`'s order.
"""
Base.getindex(od::O, key::Int64) where {O<:OrderedDict} = collect(values(od))[key] =#

# Non package-specific utilities
"""
    get_common_type(types::Vector{<: DataType})

Return the minimum common supertype of `types`.
"""
function get_common_type(types::Vector{<:DataType})
    promoted_types = DataType[]
    if length(types) > 1
        for i in 1:(length(types) - 1)
            push!(promoted_types, promote_type(types[i], types[i + 1]))
        end
    else
        return types[1]
    end

    return length(promoted_types) == 1 ? promoted_types[1] : get_common_type(promoted_types)
end

"""
    check_unique(vec::Union{Missing, <: Vector, <: Tuple})

Return `true` if all elements in `vec` are unique or if `ismissing(vec)`, else return `false`.
"""
function check_unique(vec::Union{Missing,<:Vector,<:Tuple})
    return ismissing(vec) ? true : length(setdiff(unique(vec), vec)) == 0
end

"""
    multilayer_kronecker_delta(dims::NTuple{4,Int64})

Return a 4-dimensional Kronecker delta with size equal to `dims`.
"""
function multilayer_kronecker_delta(dims::NTuple{4,Int64})
    output = Array{Int64}(undef, dims...)
    ranges = [1:dim for dim in dims]
    idxs = Iterators.product(ranges...)

    for idx in idxs
        vertices_indexes = idx[1:2]
        layers_indexes = idx[3:4]
        output[idx...] =
            if all(length(unique(vertices_indexes)) == length(unique(layers_indexes)) .== 1)
                1
            else
                0
            end
    end

    return output
end

"""
    get_diagonal_weight_tensor(arr,dims)

Return a tensor whose size is `dims` and has `arr` on the diagonal.
"""
function get_diagonal_weight_tensor(arr::Vector{T}, dims) where {T}
    output = Array{T}(undef, dims...)
    ranges = [1:dim for dim in dims]
    idxs = Iterators.product(ranges...)

    i = 1
    for idx in idxs
        vertices_indexes = idx[1:2]
        layers_indexes = idx[3:4]
        on_diag = all(
            length(unique(vertices_indexes)) == length(unique(layers_indexes)) .== 1
        )
        if on_diag
            output[idx...] = arr[i]
            i += 1
        else
            output[idx...] = zero(T)
        end
    end

    return output
end

"""
    mutable struct δk{T} <: AbstractVector{T}

The Kronecker delta.

# FIELDS

- `N::Int64`: the number of dimensions;
- `representation::Matrix{Int64}`: the matrix representing the Kronecker delta;
- `T`: the return type when called `δk[i,j]`.

### CONSTRUCTORS

    δk{T}(N::Int64) where {T <: Number}

Inner constructor that only requires `N` and the `eltype`.
"""
mutable struct δk{T} <: AbstractVector{T}
    N::Int64
    representation::Matrix{Int64}

    # Inner constructor that only requires N and the eltype.

    function δk{T}(N::Int64) where {T<:Number}
        out = new{T}(N)
        representation = [out[h, k] for h in 1:N, k in 1:N]
        out.representation = representation
        return out
    end
end

"""
    δk(N::Int64)

Outer constructor that only requires `N`.
"""
δk(N::Int64) = δk{Int64}(N)

"""
    getindex(d::δk{T}, h::Int64, k::Int64) where T

`getindex` dispatch that allows to easily construct the `representation` field of `δ_Ω` inside its inner constructor.
"""

Base.getindex(d::δk{T}, h::Int64, k::Int64) where {T} = I[h, k] ? one(T) : zero(T)

"""
    getindex(d::δk{T}, i::Int) where T

The getindex called by OMEinsum.jl.
"""
Base.getindex(d::δk{T}, i::Int) where {T} = d.representation[i]

"""
    size(d::δk)

Override required by OMEinsum.jl.
"""
Base.size(d::δk) = (d.N, d.N)

"""
    struct δ_1{T<: Number}

The `δ_1` from [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Evaluate it via the notation `[i,j]`.

# FIELDS

- `N:Int64`: the dimensionality of `δ_1`;
- `T`: the return type.

# CONSTRUCTORS

    δ_1{T<: Number}(N::Int64)
"""
struct δ_1{T<:Number}
    N::Int64
end
function Base.getindex(d::δ_1{T}, h::Int64, k::Int64, l::Int64) where {T}
    return Bool(I[h, l] * I[h, k] * I[k, l]) ? one(T) : zero(T)
end
Base.size(d::δ_1) = (d.N, d.N)

"""
    struct δ_2{T<: Number}

The `δ_2` from [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Evaluate it via the notation `[i,j]`.

# FIELDS

- `N:Int64`: the dimensionality of `δ_2`;
- `T`: the return type.

# CONSTRUCTORS

    δ_2{T<: Number}(N::Int64)
"""
struct δ_2{T<:Number}
    N::Int64
end
function Base.getindex(d::δ_2{T}, h::Int64, k::Int64, l::Int64) where {T}
    return if Bool(
        (1 - I[h, l]) * I[h, k] + (1 - I[k, l]) * I[h, l] + (1 - I[h, k]) * I[k, l]
    )
        one(T)
    else
        zero(T)
    end
end
Base.size(d::δ_2) = (d.N, d.N)

"""
    struct δ_3{T<: Number}

The `δ_3` from [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Evaluate it via the notation `[i,j]`.

# FIELDS

- `N:Int64`: the dimensionality of `δ_3`;
- `T`: the return type.

# CONSTRUCTORS

    δ_3{T<: Number}(N::Int64)
"""
struct δ_3{T<:Number}
    N::Int64
end
function Base.getindex(d::δ_3{T}, h::Int64, k::Int64, l::Int64) where {T}
    return Bool((1 - I[h, l]) * (1 - I[k, l]) * (1 - I[h, k])) ? one(T) : zero(T)
end
Base.size(d::δ_3) = (d.N, d.N)

"""
    δ_Ω{T} <: AbstractVector{T}

Struct that represents the `δ_Ω` defined in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).

# FIELDS

- `δ_1::δ_1{T}`: Instance of `δ_1`;
- `δ_2::δ_2{T}`: Instance of `δ_2`;
- `δ_3::δ_3{T}`: Instance of `δ_3`;
- `N::Int64 `  : Maximum index (number of layers);
- `representation::Array{Int64,4}`: Multidimensional-array representation of `δ_Ω`.
"""
mutable struct δ_Ω{T} <: AbstractVector{T}
    δ_1::δ_1{T}
    δ_2::δ_2{T}
    δ_3::δ_3{T}
    N::Int64
    representation::Array{Int64,4}

    # Inner constructor that only requires N and the eltype.
    function δ_Ω{T}(N::Int64) where {T<:Number}
        out = new{T}(δ_1{T}(N), δ_2{T}(N), δ_3{T}(N), N)
        representation = [out[Ω, h, k, l] for Ω in 1:3, h in 1:N, k in 1:N, l in 1:N]
        out.representation = representation
        return out
    end
end

# Outer constructor that only requires N
δ_Ω(N::Int64) = δ_Ω{Int64}(N)

# getindex dispatch that allows to easily construct the `representation` field of δ_Ω inside its inner constructor
function Base.getindex(d::δ_Ω{T}, Ω::Int64, h::Int64, k::Int64, l::Int64) where {T}
    if Ω == 1
        return d.δ_1[h, k, l]
    elseif Ω == 2
        return d.δ_2[h, k, l]
    elseif Ω == 3
        return d.δ_3[h, k, l]
    end
end

# The getindex called by OMEinsum.jl
Base.getindex(d::δ_Ω{T}, i::Int) where {T} = d.representation[i]

# Override required by OMEinsum.jl
Base.size(d::δ_Ω{T}) where {T} = (d.N >= 3 ? 3 : d.N, d.N, d.N, d.N)

"""
    get_diagonal_elements(arr::Array{T,N}) where {T,N}
"""
function get_diagonal_elements(arr::Array{T,N}) where {T,N}
    output = similar(arr)

    for cart_idx in CartesianIndices(arr)
        vertices_indexes = Tuple(cart_idx)[1:2]
        layers_indexes = Tuple(cart_idx)[3:4]
        output[cart_idx] =
            if length(unique(vertices_indexes)) == length(unique(layers_indexes)) == 1
                arr[cart_idx]
            else
                0.0
            end
    end

    return output
end

#Taken from https://gist.github.com/Datseris/1b1aa1287041cab1b2dff306ddc4f899
"""
    get_concrete_subtypes(type::Type)

Return an array of all concrete subtypes of `type` (which should be abstract).
"""
function get_concrete_subtypes(type::Type)
    out = DataType[]
    return _subtypes!(out, type)
end

function _subtypes!(out, type::Type)
    if isconcretetype(type)
        push!(out, type)
    else
        foreach(T -> _subtypes!(out, T), InteractiveUtils.subtypes(type))
    end
    return out
end

# Check if struct is completely initialized
"""
    iscompletelyinitialized(obj::Any)

Check whether `obj` is completely initialized.
"""
function iscompletelyinitialized(obj::Any)
    completely_initialized = true
    for field_name in fieldnames(typeof(obj))
        try
            @eval $obj.$field_name
        catch
            completely_initialized = false
            break
        end
    end
    return completely_initialized
end

"""
    get_oom(x::Number)

Get the order of magnitude of `x`.
"""
function get_oom(x::Number)
    float_oom = floor(log10(x))
    # if float_oom == -Inf
    #     return  0
    if !isinf(float_oom)
        return Int64(float_oom)
    else
        return float_oom
    end
end

"""
    isapproxsymmetric(A::Matrix{T}) where {T <: Real }

"Check whether `A` is approximately symmetric (within `eps(T)`).
"""
function isapproxsymmetric(A::Matrix{T}) where {T<:Real}
    return all(get_oom.(abs.(A .- A')) .<= get_oom(eps(T)))
end

"""
    isapproxsymmetric(A::Matrix{T}) where {T <: Real }

Check whether `A` is symmetric (within `zero(T)`).
"""
isapproxsymmetric(A::Matrix{T}) where {T<:Integer} = all(abs.(A .- A') .<= zero(T))

"""
    isdigraphical(indegree_sequence::AbstractVector{<:Integer}, outdegree_sequence::AbstractVector{<:Integer},)
Check whether the given indegree sequence and outdegree sequence are digraphical, that is whether they can be the indegree and outdegree sequence of a digraph.

### Implementation Notes
According to Fulkerson-Chen-Anstee theorem, a sequence ``\\{(a_1, b_1), ...,(a_n, b_n)\\}`` (sorted in descending order of a) is graphic iff the sum of vertex degrees is even and  ``\\sum_{i = 1}^{n} a_i = \\sum_{i = 1}^{n} b_i\\}`` and the sequence obeys the property -
```math
\\sum_{i=1}^{r} a_i \\leq \\sum_{i=r+1}^n min(r-1,b_i) + \\sum_{i=r+1}^n min(r,b_i)
```
for each integer r <= n-1. 
"""
function isdigraphical(
    indegree_sequence::AbstractVector{<:Integer},
    outdegree_sequence::AbstractVector{<:Integer},
)
    # Check whether the degree sequences have the same length 
    length(indegree_sequence) == length(outdegree_sequence) || throw(
        ArgumentError("The indegree and outdegree sequences must have the same length.")
    )
    # Check whether the degree sequence is empty
    !(isempty(indegree_sequence) && isempty(outdegree_sequence)) || return true
    # Check whether the degree sequences have only non-negative values
    all(indegree_sequence .>= 0) || throw(
        ArgumentError("The indegree sequence must contain non-negative integers only.")
    )
    all(outdegree_sequence .>= 0) || throw(
        ArgumentError("The outdegree sequence must contain non-negative integers only.")
    )

    n = length(indegree_sequence)
    n == length(outdegree_sequence) || return false

    sum(indegree_sequence) == sum(outdegree_sequence) || return false

    _sortperm = sortperm(indegree_sequence; rev=true)

    sorted_indegree_sequence = indegree_sequence[_sortperm]
    sorted_outdegree_sequence = outdegree_sequence[_sortperm]

    indegree_sum = zero(Int64)
    outdegree_min_sum = zero(Int64)

    cum_min = zero(Int64)

    # The following approach, which requires substituting the line
    # cum_min = sum([min(sorted_outdegree_sequence[i], r) for i in (1+r):n])
    # with the line
    # cum_min -= mindeg[r]
    # inside the for loop below, work as well, but the values of `cum_min` at each iteration differ. To be on the safe side we implemented it as in https://en.wikipedia.org/wiki/Fulkerson%E2%80%93Chen%E2%80%93Anstee_theorem
    #=     mindeg = Vector{Int64}(undef, n)
        @inbounds for i = 1:n
            mindeg[i] = min(i, sorted_outdegree_sequence[i])
        end
        cum_min = sum(mindeg) =#
    # Similarly for `outdegree_min_sum`.

    @inbounds for r in 1:(n - 1)
        indegree_sum += sorted_indegree_sequence[r]
        outdegree_min_sum = sum([min(sorted_outdegree_sequence[i], r - 1) for i in 1:r])
        cum_min = sum([min(sorted_outdegree_sequence[i], r) for i in (1 + r):n])
        cond = indegree_sum <= (outdegree_min_sum + cum_min)
        cond || return false
    end

    return true
end

#= """
    _random_undirected_configuration(empty_mg::M, degree_sequence::Vector{ <: Integer}) where {T,U,M <: MultilayerGraph{T,U}}

Internal function. Returns a `MultilayerEdge` list compatible with `empty_mg`, using a relatively inefficient algorithm.
"""
function _random_undirected_configuration(empty_mg::M, degree_sequence::Vector{ <: Integer}, allow_self_loops::Bool) where {T,U,M <: MultilayerGraph{T,U}}

    # Get all MultilayerVertexs
    mvs = mv_vertices(empty_mg)

    edge_list = MultilayerEdge{U}[]
    # Boolean that states if the wiring was successful
    success = false
    @info "Looping through wirings to find one that works..."
    # Loop until a successful wiring is found
    while !success

        mvs_degree_dict = Dict(mv => deg for (mv,deg) in zip(mvs,degree_sequence) if deg != 0)
        edge_list = MultilayerEdge[]
        for src in mvs

            if src in keys(mvs_degree_dict)
               try
                    dsts = nothing
                    if allow_self_loops
                        dsts = sample(collect(keys(mvs_degree_dict)), mvs_degree_dict[src]; replace = false) # This would be correct  but we have no "isgraphical" function that takes into account self loops. This section of the code is thus disabled.
                    else
                        dsts = sample(collect(setdiff(keys(mvs_degree_dict), [src])), mvs_degree_dict[src]; replace = false)
                    end

                    for dst in dsts

                        mvs_degree_dict[dst] = mvs_degree_dict[dst] - 1

                        if mvs_degree_dict[dst] == 0
                            delete!(mvs_degree_dict, dst)
                        end

                        descriptor = get_subgraph_descriptor(empty_mg, src.layer, dst.layer)

                        push!(edge_list, MultilayerEdge(src, dst, descriptor.default_edge_weight(src,dst), descriptor.default_edge_metadata(src,dst)))
                    end

                    delete!(mvs_degree_dict, src)
                catch e
                    if cmp(e.msg,"Cannot draw more samples without replacement.") == 0
                        break
                    else
                        throw(e)
                    end
                end
            else
                continue
            end
        end
        success = length(mvs_degree_dict) == 0
    end

    return edge_list

end

"""
    _random_directed_configuration(empty_mg::M, indegree_sequence::Vector{ <: Integer}, outdegree_sequence::Vector{ <: Integer}, allow_self_loops::Bool) where {T,U,M <: MultilayerDiGraph{T,U}}

Internal function. Returns a `MultilayerEdge` list compatible with `empty_mg`, using a relatively inefficient algorithm.
"""
function _random_directed_configuration(empty_mg::M, indegree_sequence::Vector{ <: Integer}, outdegree_sequence::Vector{<:Integer}, allow_self_loops::Bool) where {T,U,M <: MultilayerDiGraph{T,U}}

        # Get all MultilayerVertexs
        mvs = mv_vertices(empty_mg)

        edge_list = MultilayerEdge{U}[]
        # Boolean that states if the wiring was successful
        success = false
        @info "Looping through wirings to find one that works..."
        # Loop until a successful wiring is found
        while !success

            mvs_indegree_dict = Dict(mv => indeg for (mv,indeg) in zip(mvs,indegree_sequence) if indeg != 0)
            mvs_outdegree_dict = Dict(mv => outdeg for (mv,outdeg) in zip(mvs,outdegree_sequence) if outdeg != 0)
            edge_list = MultilayerEdge[]
            for src in mvs

                if src in keys(mvs_outdegree_dict)
                   try
                        dsts = nothing
                        if allow_self_loops
                            dsts = sample(collect(keys(mvs_indegree_dict)), mvs_outdegree_dict[src]; replace = false) # This would be correct  but we have no "isgraphical" function that takes into account self loops. This section of the code is thus disabled.
                        else
                            dsts = sample(collect(setdiff(keys(mvs_indegree_dict), [src])), mvs_outdegree_dict[src]; replace = false)
                        end

                        for dst in dsts

                            mvs_indegree_dict[dst] = mvs_indegree_dict[dst] - 1

                            if mvs_indegree_dict[dst] == 0
                                delete!(mvs_indegree_dict, dst)
                            end

                            descriptor = get_subgraph_descriptor(empty_mg, src.layer, dst.layer)

                            push!(edge_list, MultilayerEdge(src, dst, descriptor.default_edge_weight(src,dst), descriptor.default_edge_metadata(src,dst)))
                        end

                        delete!(mvs_outdegree_dict, src)
                    catch e
                        if cmp(e.msg,"Cannot draw more samples without replacement.") == 0
                            break
                        else
                            throw(e)
                        end
                    end
                else
                    continue
                end
            end
            success = length(mvs_indegree_dict) == 0 && length(mvs_outdegree_dict) == 0
        end
    return edge_list
end =#

"""
    cartIndexTovecIndex(cart_index::CartesianIndex, tensor_size::NTuple{N, <: Integer} ) where N

Internal function. Converts `cart_index` to an integer index such that it corresponds to the same element under flattening of the tensor whose size is `tensor_size`.
"""
function cartIndexTovecIndex(
    cart_index::Union{NTuple{N,Integer},CartesianIndex}, tensor_size::NTuple{N,<:Integer}
) where {N}
    return cart_index[1] +
           sum(collect(Tuple(cart_index)[2:end] .- 1) .* cumprod(tensor_size[1:(end - 1)]))
end

"""
    havel_hakimi_graph_generator(degree_sequence::AbstractVector{<:Integer})

Returns a simple graph with a given finite degree sequence of non-negative integers generated via the Havel-Hakimi algorithm which works as follows:
1. successively connect the node of highest degree to other nodes of highest degree;
2. sort the remaining nodes by degree in decreasing order;
3. repeat the procedure.

## References
1. [Hakimi (1962)](https://doi.org/10.1137/0110037);
2. [Wikipedia](https://en.wikipedia.org/wiki/Havel%E2%80%93Hakimi_algorithm).
"""
function havel_hakimi_graph_generator(degree_sequence::AbstractVector{<:Integer})
    # Check whether the degree sequence has only non-negative values
    all(degree_sequence .>= 0) ||
        throw(ArgumentError("The degree sequence must contain non-negative integers only."))
    # Instantiate an empty simple graph 
    graph = SimpleGraph(length(degree_sequence))
    # Create a (vertex, degree) ordered dictionary
    vertices_degrees_dict = OrderedDict(
        vertex => degree for (vertex, degree) in enumerate(degree_sequence)
    )
    # Havel-Hakimi algorithm  
    while (any(values(vertices_degrees_dict) .!= 0))
        # Sort the new sequence in non-increasing order
        vertices_degrees_dict = OrderedDict(
            sort(collect(vertices_degrees_dict); by=last, rev=true)
        )
        # Remove the first vertex and distribute its stabs
        max_vertex, max_degree = popfirst!(vertices_degrees_dict)
        # Check whether the new sequence has only positive values
        all(collect(values(vertices_degrees_dict))[1:max_degree] .> 0) ||
            throw(ErrorException("The degree sequence is not graphical."))
        # Connect the node of highest degree to other nodes of highest degree 
        for vertex in collect(keys(vertices_degrees_dict))[1:max_degree]
            add_edge!(graph, max_vertex, vertex)
            vertices_degrees_dict[vertex] -= 1
        end
    end
    # Return the simple graph
    return graph
end

"""
    lexicographical_order_lt(A::Vector{T}, B::Vector{T}) where T

The less than (lt) function that implements lexicographical order.

See [Wikipedia](https://en.wikipedia.org/wiki/Lexicographic_order).
"""
function lexicographical_order_lt(
    A::Union{Vector{T},NTuple{N,T}}, B::Union{Vector{T},NTuple{M,T}}
) where {N,M,T}
    A_dc = deepcopy(A)
    B_dc = deepcopy(B)
    diff_length = length(A) - length(B)
    if diff_length >= 0
        A_dc = vcat(A_dc, repeat([-Inf], diff_length))
    else
        B_dc = vcat(B_dc, repeat([-Inf], abs(diff_length)))
    end
    for (a, b) in zip(A_dc, B_dc)
        if a != b
            a < b
        end
    end
end

"""
lexicographical_order_ntuple(A::NTuple{N,T}, B::NTuple{M,T}) where {N,T}

The less than (lt) function that implements lexicographical order for `NTuple` of equal length.

See [Wikipedia](https://en.wikipedia.org/wiki/Lexicographic_order).
"""
function lexicographical_order_ntuple(A::NTuple{N,T}, B::NTuple{N,T}) where {N,T}
    for (a, b) in zip(A, B)
        if a != b
            return a < b
        end
    end

    return false
end

"""
    kleitman_wang_graph_generator(indegree_sequence::AbstractVector{<:Integer},outdegree_sequence::AbstractVector{<:Integer})

Returns a simple directed graph with given finite in-degree and out-degree sequences of non-negative integers generated via the Kleitman-Wang algorithm, that works like follows:
1. Sort the indegree-outdegree pairs in lexicographical order;
2. Select a pair that has strictly positive outdegree, say the i-th pairs that has outdegree = b_i;
3. Subtract 1 to the first b_i highest indegrees (the i-th being excluded), and set b_i to 0;
4. Repeat from 1. until all indegree-outdegree pairs are of the form (0.0).

## References
- [Wikipedia](https://en.wikipedia.org/wiki/Kleitman%E2%80%93Wang_algorithms)
- [Kleitman and Wang (1973)](https://doi.org/10.1016/0012-365X(73)90037-X)
"""
function kleitman_wang_graph_generator(
    indegree_sequence::AbstractVector{<:Integer},
    outdegree_sequence::AbstractVector{<:Integer},
)
    length(indegree_sequence) == length(outdegree_sequence) || throw(
        ArgumentError(
            "The provided `indegree_sequence` and `outdegree_sequence` must be of the dame length.",
        ),
    )
    # Check whether the indegree_sequence and outdegree_sequence have only non-negative values
    all(indegree_sequence .>= 0) || throw(
        ArgumentError(
            "The `indegree_sequence` sequence must contain non-negative integers only."
        ),
    )
    all(outdegree_sequence .>= 0) || throw(
        ArgumentError(
            "The `outdegree_sequence` sequence must contain non-negative integers only."
        ),
    )

    # Instantiate an empty simple graph 
    graph = SimpleDiGraph(length(indegree_sequence))
    # Create a (vertex, degree) ordered dictionary
    S = zip(deepcopy(indegree_sequence), deepcopy(outdegree_sequence))
    vertices_degrees_dict = OrderedDict(i => tup for (i, tup) in enumerate(S))
    # Kleitman-Wang algorithm  
    while (any(Iterators.flatten(values(vertices_degrees_dict)) .!= 0))
        # Sort the new sequence in non-increasing lexicographical order
        vertices_degrees_dict = OrderedDict(
            sort(
                collect(vertices_degrees_dict);
                by=last,
                lt=lexicographical_order_ntuple,
                rev=true,
            ),
        )
        # Find a vertex with positive outdegree,a nd temporarily remove it from `vertices_degrees_dict`
        i, (a_i, b_i) = 0, (0, 0)
        for (_i, (_a_i, _b_i)) in collect(deepcopy(vertices_degrees_dict))
            if _b_i != 0
                i, a_i, b_i = (_i, _a_i, _b_i)
                delete!(vertices_degrees_dict, _i)
                break
            end
        end
        # Connect the vertex found above to other nodes of highest degree 
        for (v, degs) in collect(vertices_degrees_dict)[1:b_i]
            add_edge!(graph, i, v)
            vertices_degrees_dict[v] = (degs[1] - 1, degs[2])
        end
        # Check whether the new sequence has only positive values
        all(
            collect(Iterators.flatten(collect(values(vertices_degrees_dict))))[1:b_i] .>= 0
        ) || throw(
            ErrorException("The in-degree and out-degree sequences are not digraphical."),
        )
        # Reinsert the vertex, with zero outdegree
        vertices_degrees_dict[i] = (a_i, 0)
    end
    return graph
end


"""
    sample_graphical_degree_sequence(degree_distribution::UnivariateDistribution, n::Integer)

Sample a graphical degree sequence for a graph with `n` vertices from `degree_distribution`
"""
function sample_graphical_degree_sequence(degree_distribution::UnivariateDistribution, n::Integer)

    minimum(support(degree_distribution)) >= 0 || throw(
        ErrorException(
            "The `degree_distribution` must have positive support. Found $(support(degree_distribution)).",
        ),
    )


    @info "Trying to sample a graphical sequence from the provided distribution $degree_distribution..."

    degree_sequence = nothing
    acceptable = false

    while !acceptable
        degree_sequence = round.(Ref(Int), rand(degree_distribution, n))

        acceptable = all(0 .<= degree_sequence .< n)

        if acceptable
            acceptable = isgraphical(degree_sequence)
        end
    end

    return degree_sequence

end


"""
    sample_digraphical_degree_sequences(indegree_distribution::UnivariateDistribution, outdegree_distribution::UnivariateDistribution, n::Integer)

Sample a digraphical couple `(indegree_sequence, outdegree_sequence)` for a graph with `n` vertices from respectively `indegree_distribution` and `outdegree_distribution`.
"""
function sample_digraphical_degree_sequences(indegree_distribution::UnivariateDistribution, outdegree_distribution::UnivariateDistribution, n::Integer)

    minimum(support(indegree_distribution)) >= 0 && minimum(support(outdegree_distribution)) >= 0  || throw(
        ErrorException(
            "Both the `indegree_distribution` and the `outdegree_distribution` must have positive support. Found $(support(indegree_distribution)) and $(support(outdegree_distribution)).",
        ),
    )

    indegree_sequence = Vector{Int64}(undef, n)
    outdegree_sequence = Vector{Int64}(undef, n)
    acceptable = false

    @info "Trying to sample a digraphical sequence from the two provided distributions $indegree_distribution and $outdegree_distribution..."
    while !acceptable
        indegree_sequence .=
            round.(Ref(Int), rand(indegree_distribution, n))
        outdegree_sequence .=
            round.(Ref(Int), rand(outdegree_distribution, n))

        acceptable = all(0 .<= vcat(indegree_sequence, outdegree_sequence) .< n)

        if acceptable
            acceptable = isdigraphical(indegree_sequence, outdegree_sequence)
        end
    end

    return (indegree_sequence,outdegree_sequence)

end


"""
    get_valtypes(init_function::Function)

Return the `vertexval_types` (or `edgeval_types`), deduced from from the `init_function` (which must be the `vertexval_init` or the `edgeval_init` function) for SimpleValueGraphs.jl's structs.
"""
function get_valtypes(init_function::Function)

    returned_type = Base.return_types(init_function)[1]

    if returned_type <: Tuple
        tuple(returned_type.parameters...)
    elseif returned_type <: NamedTuple
        returned_namedtuple_type = returned_type.parameters

        type_pars = tuple(returned_namedtuple_type[end].parameters...)

        NamedTuple{returned_namedtuple_type[1]}(type_pars)
    end

end
