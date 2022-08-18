"""
    AbstractVertex{T}

An abstract type representing a graph vertex.
"""
abstract type AbstractVertex{T} end

"""
    AbstractMultilayerVertex{T} <: AbstractVertex{T}

An abstract type representing an abstract MultilayerGraph vertex.
"""
abstract type AbstractMultilayerVertex{T} <: AbstractVertex{T} end

"""
    MultilayerVertex{T <: Integer} <: AbstractMultilayerVertex{T}

A struct representing a vertex of a MultilayerGraph.

# FIELDS
- `vertex::T` : the node that the MultilayerVertex represents;
- `layer::Symbol`:  the layer the MultilayerVertex belongs to.
"""
struct MultilayerVertex{T<:Integer} <: AbstractMultilayerVertex{T}
    node::T
    layer::Symbol
end
