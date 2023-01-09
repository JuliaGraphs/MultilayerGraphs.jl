# Due to the impossibility to perform multiple inheritance in Julia, if AbstractMultiplex(U/Di)Graph wish to inherit all methods that apply to AbstractMultilayer(U/Di)Graph then they may not subtype any AbstractMultiplexGraph (since they have to subtype AbstractMultilayer(U/Di)Graph). Thus,  AbstractMultiplex(U/Di)Grap will not be defined, and, in exchange, the methods that need a dispatch for general multiplex graphs (of whatever kind) will be trait-ed with the IsMultiplex trait.



# Layers and Interlayers
@traitfn function specify_interlayer!(
    mg::M, new_interlayer::In
) where {T,U,G<:AbstractGraph{T},M<:AbstractMultilayerGraph{T,U},In<:Interlayer{T,U,G}; IsMultiplex{M}}
     throw(ErrorException("Interlayers are automatically specified in multiplex graphs and they should not be modified in any way. If you need to modify interlayers, perhaps a standard multilayer graph will better suit your use case."))
end
