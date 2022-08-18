# Base extensions
"""
    getindex(od::O, key::Int64) where {O <: OrderedDict}

Return the `key`th pair of `od`, following `od`'s order.
"""
Base.getindex(od::O, key::Int64) where {O<:OrderedDict} = collect(values(od))[key]

# Non package-specific utilities
"""
    get_common_type(types::Vector{<: DataType})

Returns the minimum common supertype of `types`.
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
    is_weighted(g::G) where { G <: AbstractGraph}

Check whether g is weighted.
"""
is_weighted(g::G) where {G<:AbstractGraph} = is_weighted(typeof(g))
function is_weighted(g::G) where {G<:Type{<:AbstractGraph}}
    return g <: AbstractSimpleWeightedGraph ||
           g <: MultilayerGraph{<:Integer,<:AbstractSimpleWeightedGraph,<:Real}
end

"""
    multilayer_kronecker_delta(dims...)

Returns a 4 dimensional kronecker delta with size equal to `dims`.
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
    get_diagonal_adjacency_tensor(arr,dims)

Returns a tensor whose size is `dims` and has `arr` on the diagonal.
"""
function get_diagonal_adjacency_tensor(arr::Vector{T}, dims) where {T}
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
    mutable struct δ{T} <: AbstractVector{T}

The kronecker delta.

# FIELDS

- `N::Int64`: the number of dimensions;
- `representation::Matrix{Int64}`: the matrix representing the kronecker delta;
- `T`: the return type when called ad δ[i,j].

# CONSTRUCTORS

    δ{T}(N::Int64) where {T <: Number}

Inner constructor that only requires N and the eltype.
"""
mutable struct δ{T} <: AbstractVector{T}
    N::Int64
    representation::Matrix{Int64}

    # Inner constructor that only requires N and the eltype.
    function δ{T}(N::Int64) where {T<:Number}
        out = new{T}(N)
        representation = [out[h, k] for h in 1:N, k in 1:N]
        out.representation = representation
        return out
    end
end

"""
    δ(N::Int64)
Outer constructor that only requires N
"""
δ(N::Int64) = δ{Int64}(N)

"""
    getindex(d::δ{T}, h::Int64, k::Int64) where T 

getindex dispatch that allows to easily construct the `representation` field of δ_Ω inside its inner constructor
"""

Base.getindex(d::δ{T}, h::Int64, k::Int64) where {T} = I[h, k] ? one(T) : zero(T)

"""
    getindex(d::δ{T}, i::Int) where T

The getindex called by OMEinsum
"""
Base.getindex(d::δ{T}, i::Int) where {T} = d.representation[i]

"""
    size(d::δ)

Override required by OMEinsum
"""
Base.size(d::δ) = (d.N, d.N)

"""
    struct δ_1{T<: Number}

The δ_1 from [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Evaluate it via the notation [i,j].

# FIELDS

- `N:Int64`: the dimensionality of δ_1;
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

The δ_2 from [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Evaluate it via the notation [i,j].

# FIELDS

- `N:Int64`: the dimensionality of δ_2;
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

The δ_3 from [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022). Evaluate it via the notation [i,j].

# FIELDS

- `N:Int64`: the dimensionality of δ_3;
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

Struct that represents the δ_Ω defined in [De Domenico et al. (2013)](https://doi.org/10.1103/PhysRevX.3.041022).

# FIELDS

- `δ_1::δ_1{T}`: Instance of δ_1;
- `δ_2::δ_2{T}`: Instance of δ_2;
- `δ_3::δ_3{T}`: Instance of δ_3;
- `N::Int64 `  : Maximum index (number of layers);
- `representation::Array{Int64,4}`: Multidimensional-array representation of δ_Ω.
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

# The getindex called by OMEinsum
Base.getindex(d::δ_Ω{T}, i::Int) where {T} = d.representation[i]

# Override required by OMEinsum
Base.size(d::δ_Ω{T}) where {T} = (d.N >= 3 ? 3 : d.N, d.N, d.N, d.N)

"""
    get_diagonal_elements(arr::Array{T,N}) where {T,N}
"""
function get_diagonal_elements(arr::Array{T,N}) where {T,N}
    output = similar(arr)
    # ranges = [1:dim for dim in dims]
    # idxs = Iterators.product(ranges...)

    # for idx in idxs
    #     output[idx...] = length(unique(idx)) == 1 ? 1 : 0
    # end

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
