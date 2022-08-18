"""
AbstractMultiplexGraph{T}
"""
abstract type AbstractMultiplexGraph{T,U} <: AbstractMultilayerGraph{T,U} end

"""
MultiplexGraph(layers::Vector{ <: Layer })
"""
MultiplexGraph(layers::Vector{<:Layer}) = MultilayerGraph(layers, InterLayer[])

"""
MultiplexGraph(layers::Vector{ <: Layer })
"""
MultiplexDiGraph(layers::Vector{<:Layer}) = MultilayerDiGraph(layers, InterLayer[])
