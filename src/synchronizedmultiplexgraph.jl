"""
    SynchronizedMultiplexGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultilayerGraph{T,U}

A concrete type that can represent an undirected multiplex graph with the following special functionalities: 

1. When adding a `Node` via `add_node!, a corresponding `MultilayerVertex` is created in all layers. `add_vertex!` and `remove_vertex!` are thus disabled.




Its internal fields aren't meant to be modified by the user. Please prefer the provided API.
"""
mutable struct SynchronizedMultiplexGraph{T,U} <: AbstractMultiplexUGraph{T,U}
    layers::Vector{LayerDescriptor{T,U}} # vector containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
    interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U, SimpleGraph{T}}} # the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
    v_V_associations::Bijection{T,<:MultilayerVertex} # A Bijection from Bijections.jl that associates numeric vertices to `MultilayerVertex`s.
    idx_N_associations::Bijection{Int64,Node} # A Bijection from Bijections.jl that associates Int64 to `Node`s.
    fadjlist::Vector{Vector{HalfEdge{<:MultilayerVertex,<:Union{Nothing,U}}}} # the forward adjacency list of the SynchronizedMultiplexGraph. It is a vector of vectors of `HalfEdge`s. Its i-th element are the `HalfEdge`s that originate from `v_V_associations[i]`.
    v_metadata_dict::Dict{T,<:Union{<:Tuple,<:NamedTuple}} #  A Dictionary that associates numeric vertices to their metadata.
end

# Traits
@traitimpl IsWeighted{SynchronizedMultiplexGraph}
@traitimpl IsMeta{SynchronizedMultiplexGraph}

# Constructors
"""
    SynchronizedMultiplexGraph(
        layers::Vector{<:Layer{T,U}};
        default_interlayers_null_graph::H = SimpleGraph{T}()
    ) where {T,U, H <: AbstractGraph{T}}   

Construct a `SynchronizedMultiplexGraph` with layers given by `layers`. The interlayers will have underlying graph goven by `default_interlayer` and will of course have only diagonal couplings.

# ARGUMENTS

- `layers::Vector{<:Layer{T,U}}`: The (ordered) list of layers the multilayer graph will have;
- `specified_interlayers::Vector{<:Interlayer{T,U}}`: The list of interlayers specified by the user. Note that the user does not need to specify all interlayers, as the unspecified ones will be automatically constructed using the indications given by the `default_interlayers_null_graph` and `default_interlayers_structure` keywords;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: Sets the underlying graph for the interlayers that are to be automatically specified. Defaults to `SimpleGraph{T}()`. See the `Layer` constructors for more information;
- `default_interlayers_structure::String = "multiplex"`: Sets the structure of the interlayers that are to be automatically specified. May be "multiplex" for diagonally coupled interlayers, or "empty" for empty interlayers (no edges).  "multiplex". See the `Interlayer` constructors for more information.
"""
function SynchronizedMultiplexGraph(
    layers::Vector{<:Layer{T,U}},
    specified_interlayers::Vector{<:Interlayer{T,U}};
    default_interlayers_null_graph::H=SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U,H<:AbstractGraph{T}}
    multilayergraph = SynchronizedMultiplexGraph(T, U)

    for layer in deepcopy(layers)
        add_layer!(
            multilayergraph,
            layer;
            default_interlayers_null_graph=default_interlayers_null_graph,
            default_interlayers_structure=default_interlayers_structure,
        )
    end

    if !isnothing(specified_interlayers)
        for interlayer in deepcopy(specified_interlayers)
            specify_interlayer!(multilayergraph, interlayer)
        end
    end

    return multilayergraph
end

"""
    SynchronizedMultiplexGraph(
        empty_layers::Vector{<:Layer{T,U}},
        empty_interlayers::Vector{<:Interlayer{T,U}},
        degree_distribution::UnivariateDistribution;
        allow_self_loops::Bool = false,
        default_interlayers_null_graph::H = SimpleGraph{T}(),
    ) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

Return a random SynchronizedMultiplexGraph that has `empty_layers` as layers and `empty_interlayers` as specified interlayers. `empty_layers` and `empty_interlayers` must respectively be `Layer`s and `Interlayer`s with whatever number of vertices but no edges (if any edge is found, an error is thrown). The  degree distribution of the returned random `SynchronizedMultiplexGraph` is given by `degree_distribution`, which must have a support that only contains positive numbers for obvious reasons. `allow_self_loops = true` allows for self loops t be present in the final random SynchronizedMultiplexGraph. `default_interlayers_null_graph` controls the `null_graph` argument passed to automatically-generated interlayers. 
"""
function SynchronizedMultiplexGraph(
    empty_layers::Vector{<:Layer{T,U}},
    empty_interlayers::Vector{<:Interlayer{T,U}},
    degree_distribution::UnivariateDistribution;
    allow_self_loops::Bool=false,
    default_interlayers_null_graph::H=SimpleGraph{T}(),
) where {T<:Integer,U<:Real,H<:AbstractGraph{T}}
    !allow_self_loops || throw(
        ErrorException(
            "`allow_self_loops` must currently be set to `false`. The configuration model algorithm does not support self-loops yet.",
        ),
    )

    empty_multilayergraph = SynchronizedMultiplexGraph(
        empty_layers,
        empty_interlayers;
        default_interlayers_null_graph=default_interlayers_null_graph,
        default_interlayers_structure="empty",
    )

    n = nv(empty_multilayergraph)

    degree_sequence = sample_graphical_degree_sequence(degree_distribution, n)

    return SynchronizedMultiplexGraph(
        empty_multilayergraph, degree_sequence; allow_self_loops=false, perform_checks=false
    )
end

"""
    SynchronizedMultiplexGraph(empty_multilayergraph::SynchronizedMultiplexGraph{T,U},
    degree_sequence::Vector{<:Integer}; 
    allow_self_loops::Bool = false,
    perform_checks::Bool = false) where {T,U}

Return a random `SynchronizedMultiplexGraph` with degree sequence `degree_sequence`. `allow_self_loops` controls the presence of self-loops, while if `perform_checks` is true, the `degree_sequence` os checked to be graphical.
"""
function SynchronizedMultiplexGraph(
    empty_multilayergraph::SynchronizedMultiplexGraph{T,U},
    degree_sequence::Vector{<:Integer};
    allow_self_loops::Bool=false,
    perform_checks::Bool=true,
) where {T,U}
    (allow_self_loops && perform_checks) &&
        @warn "Checks for graphicality and coherence with the provided `empty_multilayergraph` are currently performed without taking into account self-loops. Thus said checks may fail even though the provided `degree_sequence` may be graphical when one allows for self-loops within the multilayer graph to be present. If you are sure that the provided `degree_sequence` is indeed graphical under those circumstances, you may want to disable checks by setting `perform_checks = false`. We apologize for the inconvenient."

    _multilayergraph = deepcopy(empty_multilayergraph)

    ne(_multilayergraph) == 0 || throw(
        ErrorException(
            "The `empty_multilayergraph` argument should be an empty SynchronizedMultiplexGraph. Found $(ne(_multilayergraph)) edges.",
        ),
    )

    if perform_checks
        n = nv(_multilayergraph)
        n == length(degree_sequence) || throw(
            ErrorException(
                "The number of vertices of the provided empty SynchronizedMultiplexGraph does not match the length of the degree sequence. Found $(nv(_multilayergraph)) and $(length(degree_sequence)).",
            ),
        )

        isgraphical(degree_sequence) ||
            throw(ArgumentError("degree_sequence must be graphical."))
    end

    # edge_list = _random_undirected_configuration(_multilayergraph, degree_sequence, allow_self_loops)
    equivalent_graph = havel_hakimi_graph_generator(degree_sequence)

    edge_list = [
        ME(
            empty_multilayergraph.v_V_associations[src(edge)],
            empty_multilayergraph.v_V_associations[dst(edge)],
        ) for edge in edges(equivalent_graph)
    ]

    for edge in edge_list
        add_edge!(_multilayergraph, edge)
    end

    return _multilayergraph
end

"""
    SynchronizedMultiplexGraph(T::Type{<:Number}, U::Type{<:Number})

Return a null SynchronizedMultiplexGraph with with vertex type `T` weighttype `U`. Use this constructor and then add Layers and Interlayers via the `add_layer!` and `specify_interlayer!` methods.
"""
function SynchronizedMultiplexGraph(T::Type{<:Number}, U::Type{<:Number})
    return SynchronizedMultiplexGraph{T,U}(
        LayerDescriptor{T,U}[],
        OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}}(),
        Bijection{T,MultilayerVertex}(),
        Bijection{Int64,Node}(),
        Vector{HalfEdge{MultilayerVertex,<:Union{Nothing,U}}}[],
        Dict{T,Union{Tuple,NamedTuple}}(),
    )
end

"""
    SynchronizedMultiplexGraph(layers::Vector{<:Layer{T,U}}; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String="multiplex") where {T,U, H <: AbstractGraph{T}}

Construct a SynchronizedMultiplexGraph with layers `layers` and all interlayers with structure `default_interlayers_structure` (only "multiplex" and "empty" are allowed) and type `default_interlayers_null_graph`.
"""
function SynchronizedMultiplexGraph(
    layers::Vector{<:Layer{T,U}};
    default_interlayers_null_graph::H=SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U,H<:AbstractGraph{T}}
    return SynchronizedMultiplexGraph(
        layers,
        Interlayer{T,U}[];
        default_interlayers_null_graph=default_interlayers_null_graph,
        default_interlayers_structure=default_interlayers_structure,
    )
end

# General SynchronizedMultiplexGraph Utilities
fadjlist(mg::SynchronizedMultiplexGraph) = mg.fadjlist

# Nodes

# Vertices

# Edges

"""
    add_edge_specialized!(mg::M, me::E) where {T,U, M <: MultiplexGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add MultilayerEdge `me` to the MultiplexGraph `mg`. Return true if succeeds, false otherwise.
"""
function add_edge_specialized!(
    mg::M, me::E
) where {T,U,M<:MultiplexGraph{T,U},E<:MultilayerEdge{<:Union{U,Nothing}}}
    layer(_src) == layer(_dst) || throw(ArgumentError("Cannot add an edge between two `MultilayerVertex`s that belong to two different layers in a multiplex graph"))
    add_edge_undirected!(mg, me)
end



"""
    rem_edge!(mg::AbstractMultiplexUGraph, src::MultilayerVertex, dst::MultilayerVertex)

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.rem_edge!(
    mg::AbstractMultiplexUGraph, src::MultilayerVertex, dst::MultilayerVertex
)
    # Perform routine checks
    has_vertex(mg, src) ||
        throw(ErrorException("Vertex $_src does not belong to the multilayer graph."))
    has_vertex(mg, dst) ||
        throw(ErrorException("Vertex $_dst does not belong to the multilayer graph."))

    has_edge(mg, src, dst) || return false

    src_V_idx = get_v(mg, src)
    dst_V_idx = get_v(mg, dst)

    _src = get_bare_mv(src)
    _dst = get_bare_mv(dst)

    if get_bare_mv(src) != get_bare_mv(dst)
        src_idx_tbr = findfirst(
            halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx]
        )
        deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)

        dst_idx_tbr = findfirst(halfedge -> halfedge.vertex == _src, mg.fadjlist[dst_V_idx])
        deleteat!(mg.fadjlist[dst_V_idx], dst_idx_tbr)
    else
        src_idx_tbr = findfirst(
            halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx]
        )
        deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)
    end

    return true
end

# Layers and Interlayers

# Graphs.jl's extensions

# Multilayer-specific methods
# "empty graph" could be the correct way of calling a graph with no edges: https://math.stackexchange.com/questions/320859/what-is-the-term-for-a-graph-on-n-vertices-with-no-edges

# Base overloads
"""
    Base.getproperty(mg::M, f::Symbol) where { M <: SynchronizedMultiplexGraph }
"""
function Base.getproperty(mg::SynchronizedMultiplexGraph, f::Symbol)#  where {T,U,M<:AbstractMultiplexUGraph{T,U}}
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
end
