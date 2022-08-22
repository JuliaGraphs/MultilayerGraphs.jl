# Define traits from SimpleTraits.jl
"""
    IsWeighted{X}

Trait that discerns between weighted and unweighted graphs.
"""
@traitdef IsWeighted{X}

@traitimpl IsWeighted{SimpleWeightedGraph}
@traitimpl IsWeighted{SimpleWeightedDiGraph}

@traitimpl IsWeighted{MetaGraph}
@traitimpl IsWeighted{MetaDiGraph}

@traitimpl IsWeighted{ValGraph}
@traitimpl IsWeighted{ValOutDiGraph}
@traitimpl IsWeighted{ValDiGraph}
