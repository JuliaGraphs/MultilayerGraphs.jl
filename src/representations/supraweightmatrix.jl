"""
    SupraWeightMatrix{T,U}

A concrete type representing the (supra) weight matrix of a multilayer graph. It takes into account missing vertices by default. Look at the EXAMPLES section to learn how to use it.

# EXAMPLES

```julia
# Assuming a MultilayerGraph named mg is defined, and that mv1 and mv2 are two of its `MultilayerVertex`ss
swm = SupraWeightMatrix(mg)
# One may access te corresponding SupraWeightMatrix's entry via:
swm[mv1, mv2]
```
"""
struct SupraWeightMatrix{T,U} <: AbstractMatrixRepresentation
    array::Array{U,2}
    v_V_associations::Bijection{T,<:Union{<:MissingVertex,<:MultilayerVertex}}
end
