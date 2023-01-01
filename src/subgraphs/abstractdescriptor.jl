"""
    AbstractDescriptor{T<:Integer,U<:Real,G<:AbstractGraph{T}}

An abstract tyep representing an abstract notion of subgraph descriptor
"""
abstract type AbstractDescriptor{T<:Integer,U<:Real,G<:AbstractGraph{T}} end

name(x::AbstractDescriptor)  = x.name
graph(x::AbstractDescriptor) = x.null_graph