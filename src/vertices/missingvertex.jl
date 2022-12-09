"""
    MissingVertex

A mutable struct that acts as a placeholder for a vertex that is missing in a Layer. It is mutable so that it may be added more than once to `Bijections` from Bijections.jl.s
"""
mutable struct MissingVertex <: AbstractVertex end