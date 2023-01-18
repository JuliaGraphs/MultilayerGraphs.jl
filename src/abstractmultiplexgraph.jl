# Due to the impossibility to perform multiple inheritance in Julia, if AbstractMultiplex(U/Di)Graph wish to inherit all methods that apply to AbstractMultilayer(U/Di)Graph then they may not subtype any AbstractMultiplexGraph (since they have to subtype AbstractMultilayer(U/Di)Graph). Thus,  AbstractMultiplex(U/Di)Grap will not be defined, and, in exchange, the methods that need a dispatch for general multiplex graphs (of whatever kind) will be trait-ed with the IsMultiplex trait.

# Anyway, both SimpleTraits.jl and WhereTraits.jl do not allow for implementing a trait-based sub ecosystem right now. 

# SimpleTraits.jl cannot dispatch on multple traits in a single method and cannot create dispatched of the same function over different traits (only dispatching over ONE trait and its negation is allowed). Unfortunately SimpleTraots.jl is no longer maintained

# WhereTraits.jl's @traits is not compatible with the current Graphs.jl interface

# Thus we temporarily resort to using a mix of type hierarchy and SimpleTraits.jl's basic IsDirected traits.

"""
"""
abstract type AbstractMultiplexGraph{T,U} <: AbstractMultilayerGraph{T,U} end



# Layers and Interlayers
function specify_interlayer!( mg::M, new_interlayer::In
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultiplexGraph{T,U},In<:Interlayer{T,U,G} } #istrait(IsMultiplex{M})
     throw(ErrorException("Interlayers are automatically specified in multiplex graphs and they should not be modified in any way. If you need to modify interlayers, perhaps a standard multilayer graph will better suit your use case."))
end


