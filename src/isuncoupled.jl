"""
    add_node!(mg::M, n::Node) where {M <: AbstractMultilayerGraph; IsUncoupled{M}}

Add Node `n` to the 
"""
@traitfn add_node!(mg::M, n::Node) where {M <: AbstractMultilayerGraph; IsUncoupled{M}} = _add_node!(mg, n; add_vertex_to_layers = :all)