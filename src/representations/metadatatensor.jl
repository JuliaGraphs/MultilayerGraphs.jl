"""
    MetadataTensor{U}

Concrete type representing the metadata tensor of a multilayer graph.  Look at the EXAMPLES section below to learn how to use it.

# EXAMPLES

```julia
# Assuming a MultilayerGraph named mg is defined, and that mv1 and mv2 are two of its `MultilayerVertex`s
mt = MetadataTensor(mg)
# One may access te corresponding MetadataTensor's entry via:
mt[mv1, mv2]
```
"""
struct MetadataTensor <: AbstractTensorRepresentation
    array::Array{ <: Union{Nothing,Tuple, NamedTuple},4}
    layers_names::Vector{Symbol}
    idx_N_associations::Bijection{Int64, Node}
end