"""
    WeightTensor{U}

Concrete type representing the weight tensor of a multilayer graph.  Look at the EXAMPLES section below to learn how to use it.

# EXAMPLES

```julia
# Assuming a MultilayerGraph named mg is defined, and that mv1 and mv2 are two of its `MultilayerVertex`ss
wt = WeightTensor(mg)
# One may access te corresponding WeightTensor's entry via:
wt[mv1, mv2]
```
"""
struct WeightTensor{U<:Real} <: AbstractTensorRepresentation
    array::Array{U,4}
    layers_names::Vector{Symbol}
    idx_N_associations::Bijection{Int64,Node}
end
