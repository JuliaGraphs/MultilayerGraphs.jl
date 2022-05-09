# Base extensions

"""
"""
Base.getindex(od::O, key::Int64) where {O <: OrderedDict} = collect(values(od))[key]

# Non package-specific utilities

"""
"""
function get_common_type(types::Vector{<: DataType})
    promoted_types = DataType[]
    for i in 1:(length(types)-1)
        push!(promoted_types, promote_type(types[i], types[i+1]))
    end

    length(promoted_types) == 1 ? promoted_types[1] : get_common_type(promoted_types)
end

"""
"""
check_unique(vec::Union{Missing, <: Vector, <: Tuple}) = ismissing(vec) ? true : length(setdiff(unique(vec), vec)) == 0 

