"""
    abstract type AbstractInterlayerDescriptor{ T <: Integer, U <: Real, G <: AbstractGraph{T}} end

An abstract type representing an `Interlayer` descriptor object.
"""
abstract type AbstractInterlayerDescriptor{ T <: Integer, U <: Real, G <: AbstractGraph{T}} end

"""
    struct InterlayerDescriptor{T,U,G} <: AbstractInterlayerDescriptor{T,U,G}

Custom concrete type representing an `Interlayer` descriptor object.

# FIELDS

-`name::Symbol`:.
-`layer_1::Symbol`:.
-`layer_2::Symbol`:.
-`null_graph::G`:.
-`default_edge_weight::Function`:.
-`default_edge_metadata::Function`:.
-`transfer_vertex_metadata::Bool`:.
"""
struct InterlayerDescriptor{T,U,G} <: AbstractInterlayerDescriptor{T,U,G}
    name::Symbol
    layer_1::Symbol
    layer_2::Symbol
    null_graph::G
    default_edge_weight::Function
    default_edge_metadata::Function
    transfer_vertex_metadata::Bool

    # Override inner constructor to check that the empty graph is indeed empty
    function InterlayerDescriptor( name::Symbol, layer_1::Symbol, layer_2::Symbol, null_graph::G, default_edge_weight::Function, default_edge_metadata::Function, transfer_vertex_metadata::Bool, weighttype::Type{U}) where {T,U, G <: AbstractGraph{T}} #  edge_metadata_function::Function, edge_weight_function::Function, 
        (layer_1 != layer_2) || throw(ErrorException("`layer_1` must be different from `layer_2`. Found $layer_1 and $layer_2 ."))
        nv(null_graph) == ne(null_graph) == 0 || throw(ErrorException("The provided graph is not empty. It contains $(nv(null_graph)) vertices and $(ne(null_graph)). Expected 0 and 0."))
        return new{T,weighttype,G}(name, layer_1, layer_2, null_graph, default_edge_weight, default_edge_metadata, transfer_vertex_metadata) # edge_weight_function, edge_metadata_function,
    end
end

# Outer constructor that does not require name, edge_metadata_function and transfer_vertex_metadata
InterlayerDescriptor(layer_1::Symbol, layer_2::Symbol, null_graph::G, weighttype::Type{U}; default_edge_weight::Function = (src, dst) -> nothing, default_edge_metadata::Function = (src, dst) -> NamedTuple(),  transfer_vertex_metadata::Bool = false, name::Symbol = Symbol("interlayer_$(layer_1)_$(layer_2)") ) where {T,U, G <: AbstractGraph{T}}= InterlayerDescriptor(name, layer_1, layer_2, null_graph, default_edge_weight, default_edge_metadata, transfer_vertex_metadata, weighttype) # edge_weight_function::Function = (src,dst) -> one(U), edge_metadata_function::Function = (src,dst) -> NamedTuple(),

# This function is called whenever one tries to access a field of an object of type InterlayerDescriptor
function Base.getproperty(descriptor::InterlayerDescriptor, f::Symbol)
    # If the field is one of the ones listed below, return its value
    if f ∈ (:name, :layer_1, :layer_2, :null_graph, :transfer_vertex_metadata, :default_edge_metadata,:default_edge_weight )
        Base.getfield(descriptor, f)
    # Otherwise, if the field is layers_names, return the names of the two layers
    elseif f == :layers_names
        [descriptor.layer_1, descriptor.layer_2]
    end
end