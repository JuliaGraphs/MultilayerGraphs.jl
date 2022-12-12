"""
    AbstractMatrixRepresentation{T,U}

An abstract type encoding a generic matrix representation of the links and metadata of a multilayer graph. 

Concrete subtypes must have an `array` field (a matrix of eltype U) and a `v_V_associations` (a `Bijection{T, Union{MissingVertex, MultilayerVertex}}`).

# PARAMETRIC TYPES

- `T`: the type of the internal representation of vertices;
- `U`: the weight type of the multilayer graph.
"""
abstract type AbstractMatrixRepresentation end

"""
    getindex(amr::AbstractMatrixRepresentation, src_mv::MultilayerVertex, dst_mv::MultilayerVertex)
"""
function Base.getindex(amr::AbstractMatrixRepresentation, src_mv::MultilayerVertex, dst_mv::MultilayerVertex)
    src_idx = amr.v_V_associations(get_bare_mv(src_mv))
    dst_idx = amr.v_V_associations(get_bare_mv(dst_mv))

    return amr.array[src_idx, dst_idx]
end

"""
    getindex(amr::WeightTensor, src_tup::Tuple{String, Symbol}, dst_tup::Tuple{String, Symbol})
"""
function Base.getindex(amr::AbstractMatrixRepresentation, src_tup::Tuple{String, Symbol}, dst_tup::Tuple{String, Symbol})
    src_idx = amr.v_V_associations(MV(Node(src_tup[1]), src_tup[2]))
    dst_idx = amr.v_V_associations(MV(Node(dst_tup[1]), dst_tup[2]))

    return amr.array[src_idx, dst_idx]
end

"""
    array(amr::AbstractMatrixRepresentation) 

Return the array representation of `amr`.
"""
array(amr::AbstractMatrixRepresentation) = amr.array