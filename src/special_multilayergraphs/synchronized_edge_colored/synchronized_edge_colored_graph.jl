"""
    SynchronizedEdgeColoredGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultilayerGraph{T,U}

A concrete type that can represent a general edge colored graph, that is synchronized i.e. that it represents every node in each layer. Thus:
-  `add_node!` will always add the corresponding vertex in all layers;
-  `add_vertex!` and `rem_vertex!` are not available for this type;
-  All Interlayers automatically added by `add_layer!` are empty simple graphs.
- `specify_interlayer!` is not available.

Its internal fields aren't meant to be modified by the user. Please prefer the provided API.
"""
mutable struct SynchronizedEdgeColoredGraph{T,U} <: AbstractSynchronizedEdgeColoredGraph{T,U}
    layers::Vector{LayerDescriptor{T,U}} # vector containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
    interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}} # the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
    v_V_associations::Bijection{T,<:MultilayerVertex} # A Bijection from Bijections.jl that associates numeric vertices to `MultilayerVertex`s.
    idx_N_associations::Bijection{Int64,Node} # A Bijection from Bijections.jl that associates Int64 to `Node`s.
    fadjlist::Vector{Vector{HalfEdge{<:MultilayerVertex,<:Union{Nothing,U}}}} # the forward adjacency list of the SynchronizedEdgeColoredGraph. It is a vector of vectors of `HalfEdge`s. Its i-th element are the `HalfEdge`s that originate from `v_V_associations[i]`.
    v_metadata_dict::Dict{T,<:Union{<:Tuple,<:NamedTuple}} #  A Dictionary that associates numeric vertices to their metadata.
end

# Traits
@traitimpl IsWeighted{SynchronizedEdgeColoredGraph}
@traitimpl IsMeta{SynchronizedEdgeColoredGraph}


"""
    is_directed(m::M) where { M <: Type{ <: SynchronizedEdgeColoredGraph}}

Return `false`
"""
Graphs.is_directed(mg::M) where {M<:Type{<:SynchronizedEdgeColoredGraph}}  = false


# Constructors
"""
    SynchronizedEdgeColoredGraph(
        layers::Vector{<:Layer{T,U}}
    ) where {T,U, H <: AbstractGraph{T}}   

Construct a SynchronizedEdgeColoredGraph with layers given by `layers`. The interlayers will be constructed by default as empty.

# ARGUMENTS

- `layers::Vector{<:Layer{T,U}}`: The (ordered) list of layers the multilayer graph will have;
"""
function SynchronizedEdgeColoredGraph(
    layers::Vector{<:Layer{T,U}},
) where {T,U}

    multilayergraph = SynchronizedEdgeColoredGraph(T, U)

    for layer in deepcopy(layers)
        add_layer!(
            multilayergraph,
            layer
        )
    end


    return multilayergraph
end


"""
    SynchronizedEdgeColoredGraph(T::Type{<:Number}, U::Type{<:Number})

Return a null SynchronizedEdgeColoredGraph with with vertex type `T` weighttype `U`. Use this constructor and then add Layers via the `add_layer!`method.
"""
function SynchronizedEdgeColoredGraph(T::Type{<:Number}, U::Type{<:Number})
    return SynchronizedEdgeColoredGraph{T,U}(
        LayerDescriptor{T,U}[],
        OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}}(),
        Bijection{T,MultilayerVertex}(),
        Bijection{Int64,Node}(),
        Vector{HalfEdge{MultilayerVertex,<:Union{Nothing,U}}}[],
        Dict{T,Union{Tuple,NamedTuple}}(),
    )
end