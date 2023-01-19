"""
    MultiplexGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultilayerGraph{T,U}

A concrete type that can represent a general (undirected) multiplex graph. Differently from `MultilayerGraph`, this type will error if:
1. The user tries to specify one of its interlayers;
2. The user tries to add and inter-layer edge which is not diagonal;
3. The user tries to remove an inter-layer edge.

Its internal fields aren't meant to be modified by the user. Please prefer the provided API.
"""
mutable struct MultiplexGraph{T,U} <: AbstractMultiplexGraph{T,U}
    layers::Vector{LayerDescriptor{T,U}} # vector containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
    interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}} # the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
    v_V_associations::Bijection{T,<:MultilayerVertex} # A Bijection from Bijections.jl that associates numeric vertices to `MultilayerVertex`s.
    idx_N_associations::Bijection{Int64,Node} # A Bijection from Bijections.jl that associates Int64 to `Node`s.
    fadjlist::Vector{Vector{HalfEdge{<:MultilayerVertex,<:Union{Nothing,U}}}} # the forward adjacency list of the MultiplexGraph. It is a vector of vectors of `HalfEdge`s. Its i-th element are the `HalfEdge`s that originate from `v_V_associations[i]`.
    v_metadata_dict::Dict{T,<:Union{<:Tuple,<:NamedTuple}} #  A Dictionary that associates numeric vertices to their metadata.
end

# Traits
@traitimpl IsWeighted{MultiplexGraph}
@traitimpl IsMeta{MultiplexGraph}

# Constructors
"""
    MultiplexGraph(
        layers::Vector{<:Layer{T,U}};
        default_interlayers_null_graph::H = SimpleGraph{T}().
    ) where {T,U, H <: AbstractGraph{T}}   

Construct an undirected multiplex graph. The interlayers will have underlying graph type `default_interlayers_null_graph` and of course have only diagonal couplings.

# ARGUMENTS

- `layers::Vector{<:Layer{T,U}}`: The (ordered) list of layers the multilayer graph will have;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: Sets the underlying graph for the interlayers that are to be automatically specified. Defaults to `SimpleGraph{T}()`. See the `Interlayer` constructors for more information;
"""
function MultiplexGraph(
    layers::Vector{<:Layer{T,U}};
    default_interlayers_null_graph::H=SimpleGraph{T}()
) where {T,U,H<:AbstractGraph{T}}
    multilplexgraph = MultiplexGraph(T, U)

    for layer in deepcopy(layers)
        add_layer!(
            multilplexgraph,
            layer;
            default_interlayers_null_graph=default_interlayers_null_graph,
            default_interlayers_structure="multiplex",
        )
    end

    return multilplexgraph
end


"""
    MultiplexGraph(T::Type{<:Number}, U::Type{<:Number})

Return a null MultiplexGraph with with vertex type `T` weighttype `U`. Use this constructor and then add Layers via the `add_layer!`.
"""
function MultiplexGraph(T::Type{<:Number}, U::Type{<:Number})
    return MultiplexGraph{T,U}(
        LayerDescriptor{T,U}[],
        OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}}(),
        Bijection{T,MultilayerVertex}(),
        Bijection{Int64,Node}(),
        Vector{HalfEdge{MultilayerVertex,<:Union{Nothing,U}}}[],
        Dict{T,Union{Tuple,NamedTuple}}(),
    )
end


# General MultiplexGraph Utilities
fadjlist(mg::MultiplexGraph) = mg.fadjlist

# Nodes

# Vertices

# Edges

"""
    add_edge_specialized!(mg::M, me::E) where {T,U, M <: MultiplexGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add MultilayerEdge `me` to the multiplex graph `mg`. Return true if succeeds, false otherwise.
"""
function add_edge_specialized!(
    mg::M, me::E
) where {T,U,M<:MultiplexGraph{T,U},E<:MultilayerEdge{<:Union{U,Nothing}}}
    layer(_src) == layer(_dst) || throw(ArgumentError("Cannot add an edge between two `MultilayerVertex`s that belong to two different layers in a multiplex graph"))
    _add_edge!(mg, me)
end



"""
    rem_edge_specialized!(mg::MultiplexGraph, src::MultilayerVertex, dst::MultilayerVertex)

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
function rem_edge_specialized!(
    mg::MultiplexGraph, src::MultilayerVertex, dst::MultilayerVertex)
    (layer(src) != layer(dst) && node(src) == node(dst)) && throw(ArgumentError("Cannot remove an edge between two `MultilayerVertex`s that belong to two different layers in a multiplex graph"))
    _rem_edge!(mg, src, dst)

end

# Layers and Interlayers


# Graphs.jl's extensions

# Multilayer-specific methodss