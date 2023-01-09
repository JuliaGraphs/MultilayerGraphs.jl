
"""
MultiplexDiGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultilayerGraph{T,U}

A concrete type that can represent a general multilayer graph. Its internal fields aren't meant to be modified by the user. Please prefer the provided API.
"""
mutable struct MultiplexDiGraph{T,U} <: AbstractMultilayerDiGraph{T,U}
layers::Vector{LayerDescriptor{T,U}} # vector containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}} # the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
v_V_associations::Bijection{T,<:MultilayerVertex} # A Bijection from Bijections.jl that associates numeric vertices to `MultilayerVertex`s.
idx_N_associations::Bijection{Int64,Node} # A Bijection from Bijections.jl that associates Int64 to `Node`s.
fadjlist::Vector{Vector{HalfEdge{<:MultilayerVertex,<:Union{Nothing,U}}}} # the forward adjacency list of the MultiplexDiGraph. It is a vector of vectors of `HalfEdge`s. Its i-th element are the `HalfEdge`s that originate from `v_V_associations[i]`.
v_metadata_dict::Dict{T,<:Union{<:Tuple,<:NamedTuple}} #  A Dictionary that associates numeric vertices to their metadata.
end

# Traits
@traitimpl IsWeighted{MultiplexDiGraph}
@traitimpl IsMeta{MultiplexDiGraph}

# Constructors
"""
MultiplexDiGraph(
    layers::Vector{<:Layer{T,U}};
    default_interlayers_null_graph::H = SimpleGraph{T}().
) where {T,U, H <: AbstractGraph{T}}   

Construct an undirected multiplex graph. The interlayers will have underlying graph type `default_interlayers_null_graph` and of course have only diagonal couplings.

# ARGUMENTS

- `layers::Vector{<:Layer{T,U}}`: The (ordered) list of layers the multilayer graph will have;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: Sets the underlying graph for the interlayers that are to be automatically specified. Defaults to `SimpleGraph{T}()`. See the `Interlayer` constructors for more information;
"""
function MultiplexDiGraph(
layers::Vector{<:Layer{T,U}};
default_interlayers_null_graph::H=SimpleGraph{T}()
) where {T,U,H<:AbstractGraph{T}}
multilplexgraph = MultiplexDiGraph(T, U)

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

@traitimpl IsMultiplex{MultiplexDiGraph}


"""
MultiplexDiGraph(T::Type{<:Number}, U::Type{<:Number})

Return a null MultiplexDiGraph with with vertex type `T` weighttype `U`. Use this constructor and then add Layers via the `add_layer!`.
"""
function MultiplexDiGraph(T::Type{<:Number}, U::Type{<:Number})
return MultiplexDiGraph{T,U}(
    LayerDescriptor{T,U}[],
    OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}}(),
    Bijection{T,MultilayerVertex}(),
    Bijection{Int64,Node}(),
    Vector{HalfEdge{MultilayerVertex,<:Union{Nothing,U}}}[],
    Dict{T,Union{Tuple,NamedTuple}}(),
)
end

#= """
MultiplexDiGraph(layers::Vector{<:Layer{T,U}}; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String="multiplex") where {T,U, H <: AbstractGraph{T}}

Construct a MultiplexDiGraph with layers `layers` and all interlayers with structure `default_interlayers_structure` (only "multiplex" and "empty" are allowed) and type `default_interlayers_null_graph`.
"""
function MultiplexDiGraph(
layers::Vector{<:Layer{T,U}};
default_interlayers_null_graph::H=SimpleGraph{T}(),
default_interlayers_structure::String="multiplex",
) where {T,U,H<:AbstractGraph{T}}
return MultiplexDiGraph(
    layers,
    Interlayer{T,U}[];
    default_interlayers_null_graph=default_interlayers_null_graph,
    default_interlayers_structure=default_interlayers_structure,
)
end
=#
# General MultiplexDiGraph Utilities
fadjlist(mg::MultiplexDiGraph) = mg.fadjlist

# Nodes

# Vertices

# Edges

"""
add_edge_specialized!(mg::M, me::E) where {T,U, M <: MultiplexDiGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add MultilayerEdge `me` to the multiplex graph `mg`. Return true if succeeds, false otherwise.
"""
function add_edge_specialized!(
mg::M, me::E
) where {T,U,M<:MultiplexDiGraph{T,U},E<:MultilayerEdge{<:Union{U,Nothing}}}
layer(_src) == layer(_dst) || throw(ArgumentError("Cannot add an edge between two `MultilayerVertex`s that belong to two different layers in a multiplex graph"))
add_edge_directed!(mg, me)
end



"""
rem_edge_specialized!(mg::MultiplexDiGraph, src::MultilayerVertex, dst::MultilayerVertex)

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
function rem_edge_specialized!(
mg::MultiplexDiGraph, src::MultilayerVertex, dst::MultilayerVertex)
(layer(src) != layer(dst) && node(src) == node(dst)) && throw(ArgumentError("Cannot remove an edge between two `MultilayerVertex`s that belong to two different layers in a multiplex graph"))
rem_edge_directed!(mg, src, dst)

end

# Layers and Interlayers


# Graphs.jl's extensions

# Multilayer-specific methods
# "empty graph" could be the correct way of calling a graph with no edges: https://math.stackexchange.com/questions/320859/what-is-the-term-for-a-graph-on-n-vertices-with-no-edges

#= # Base overloads
"""
Base.getproperty(mg::M, f::Symbol) where { M <: MultiplexDiGraph }
"""
function Base.getproperty(mg::MultiplexDiGraph, f::Symbol)
if f in (
    :v_V_associations,
    :fadjlist,
    :idx_N_associations,
    :layers,
    :interlayers,
    :v_metadata_dict,
) # :weight_tensor, :supra_weight_matrix, 
    Base.getfield(mg, f)
elseif f == :edge_list
    return edges(mg)
elseif f == :subgraphs
    return merge(mg.layers, mg.interlayers)
elseif f == :layers_names
    return [layer.name for layer in mg.layers]
elseif f == :interlayers_names
    return [interlayer.name for interlayer in values(mg.interlayers)]
elseif f == :subgraphs_names
    return vcat(mg.layers_names, mg.interlayers_names)
else
    for descriptor in mg.layers
        if descriptor.name == f
            return get_subgraph(mg, descriptor)
        end
    end

    for descriptor in values(mg.interlayers)
        if descriptor.name == f
            return get_subgraph(mg, descriptor)
        end
    end
end
end =#
