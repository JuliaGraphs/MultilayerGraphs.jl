"""
    abstract type AbstractGraphOfGraphs{T <: AbstractGraph} <: AbstractGraph{T} end

An abstract type representing a graph of graphs.
"""
abstract type AbstractGraphOfGraphs{T<:AbstractGraph} <: AbstractGraph{T} end