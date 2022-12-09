# Define traits from SimpleTraits.jl. Right now, no trait except `IsDirected` is used by MultilayerGraphs.jl
"""
    IsWeighted{X}

Trait that discerns between weighted and unweighted graphs. A graph type should take the `IsWeighted` trait IF AND ONLY IF it implements the signature add_edge!(src,dst,weight). Otherwise it should not.
"""
@traitdef IsWeighted{X}

@traitimpl IsWeighted{SimpleWeightedGraphs.AbstractSimpleWeightedGraph}

"""
    is_weighted(g::G) where { G <: AbstractGraph}

Check whether `g` is weighted.
"""
is_weighted(g::G) where {G<:AbstractGraph} = is_weighted(G)

"""
    is_weighted(g::G) where {G<:Type{<:AbstractGraph}} 

Check whether `g` is weighted.
"""
is_weighted(g::G) where {G<:Type{<:AbstractGraph}} = istrait(IsWeighted{g})



"""
    IsMeta{X}

Trait that discerns between graphs that sport edge and vertex  metadata.
"""
@traitdef IsMeta{X}

@traitimpl IsMeta{MetaGraphs.AbstractMetaGraph}
@traitimpl IsMeta{SimpleValueGraphs.AbstractValGraph}

"""
    is_meta(g::G) where {G <: AbstractGraph}

Check whether `g` supports edge AND vertex metadata.
"""
is_meta(g::G) where {G<:AbstractGraph} = is_meta(typeof(g))

"""
    is_meta(g::G) where {G<:Type{<:AbstractGraph}}

Check whether `g` supports edge AND vertex metadata.
"""
is_meta(g::G) where {G<:Type{<:AbstractGraph}} = istrait(IsMeta{g})