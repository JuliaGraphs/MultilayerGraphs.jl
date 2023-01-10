"""
    add_node!(mg::M, n::Node) where {M <: AbstractMultilayerGraph; IsUncoupled{M}}

Add Node `n` to the 
"""
@traits add_node!(mg::M, n::Node) where {M <: AbstractMultilayerGraph, istrait(IsUncoupled{M})} = _add_node!(mg, n; add_vertex_to_layers = :all)