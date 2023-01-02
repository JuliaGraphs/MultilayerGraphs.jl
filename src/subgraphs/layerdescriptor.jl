"""
    abstract type AbstractLayerDescriptor{ T <: Integer, U <: Real, G <: AbstractGraph{T}} end

Abstract type representing a `Layer` descriptor object.
"""
abstract type AbstractLayerDescriptor{T<:Integer,U<:Real,G<:AbstractGraph{T}} <:
              AbstractDescriptor{T,U,G} end

"""
    struct LayerDescriptor{T,U,G} <: AbstractLayerDescriptor{T,U,G}

Custom concrete type representing a `Layer` descriptor object.

# FIELDS

-`name::Symbol`.
-`null_graph::G`.
-`default_vertex_metadata::Function`.
-`default_edge_weight::Function`.
-`default_edge_metadata::Function`.
"""
struct LayerDescriptor{T,U,G} <: AbstractLayerDescriptor{T,U,G}
    name::Symbol
    null_graph::G
    default_vertex_metadata::Function
    default_edge_weight::Function
    default_edge_metadata::Function

    # Override inner constructor to check that the `null_graph` is indeed empty and to assign the parametric type
    function LayerDescriptor(
        name::Symbol,
        null_graph::G,
        weighttype::Type{U};
        default_vertex_metadata::Function=mv -> NamedTuple(),
        default_edge_weight::Function=(src, dst) -> nothing,
        default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    ) where {T,U,G<:AbstractGraph{T}}
        nv(null_graph) == ne(null_graph) == 0 || throw(
            ErrorException(
                "The provided graph is not empty. It contains $(nv(null_graph)) vertices and $(ne(null_graph)). Expected 0 and 0.",
            ),
        )
        return new{T,weighttype,G}(
            name,
            null_graph,
            default_vertex_metadata,
            default_edge_weight,
            default_edge_metadata,
        )
    end
end

# Console print utilities
function to_string(x::LayerDescriptor)
    parameters = typeof(x).parameters
    return """
           Layer\t$(name(x))
           underlying_graph: $(typeof(graph(x)))
           vertex_type: $(parameters[1])
           weight_type: $(parameters[2]) 
           """
end
Base.show(io::IO, x::LayerDescriptor) = print(io, to_string(x))
