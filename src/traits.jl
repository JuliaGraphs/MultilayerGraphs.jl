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

Trait that discerns between graphs that sport edge and vertex metadata.
"""
@traitdef IsMeta{X}

@traitimpl IsMeta{MetaGraphs.AbstractMetaGraph}
@traitimpl IsMeta{SimpleValueGraphs.AbstractValGraph}

"""
    is_meta(g::G) where {G <: AbstractGraph}

Check whether `g` supports edge AND vertex metadata.
"""
is_meta(g::G) where {G<:AbstractGraph} = is_meta(G)

"""
    is_meta(g::G) where {G<:Type{<:AbstractGraph}}

Check whether `g` supports edge AND vertex metadata.
"""
is_meta(g::G) where {G<:Type{<:AbstractGraph}} = istrait(IsMeta{g})



# IsMultiplex{X}

#= """
    is_multiplex(mg::M) where {M<:Type{<:AbstractGraph}}

Check whether `mg` is a multiplex graph.
"""
is_multiplex(mg::M) where {M<:Type{<:AbstractGraph}} = istrait(IsMultiplex{mg}) =#

"""
    IsMultiplex{X}

Trait that characterizes multilayer graphs that have multiplex-like behavior (i.e. diagonal and immutable interlayers).
"""
@traitdef IsMultiplex{X}

# IsMultiplexDirected
@traitdef IsMultiplexDirected{X}

is_multiplex_directed(X) = !istrait(IsMultiplex{X}) && istrait(IsDirected{X})
@traitimpl IsMultiplexDirected{X} <- is_multiplex_directed(X)

# IsMultiplexNotDirected
@traitdef IsMultiplexNotDirected{X}

is_multiplex_not_directed(X) = istrait(IsMultiplex{X}) && !istrait(IsDirected{X})
@traitimpl IsMultiplexNotDirected{X} <- is_multiplex_not_directed(X)
# IsNotMultiplexDirected
@traitdef IsNotMultiplexDirected{X}

is_not_multiplex_directed(X) =  !istrait(IsMultiplex{X}) && istrait(IsDirected{X})
@traitimpl IsNotMultiplexDirected{X} <- is_not_multiplex_directed(X)

# IsNotMultiplexNotDirected
@traitdef IsNotMultiplexNotDirected{X}

is_not_multiplex_not_directed(X) = !istrait(IsMultiplex{X}) && !istrait(IsDirected{X})
@traitimpl IsNotMultiplexNotDirected{X} <- is_not_multiplex_not_directed(X)






# IsUncoupled


"""
    IsUncoupled{X}

Trait that characterizes multilayer graphs that have uncoupled layers (i.e. no inter-layer edges).
"""
@traitdef IsUncoupled{X}

