"""
    MultilayerGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultilayerGraph{T,U}

A concrete type that can represent a general multilayer graph. Its internal fields aren't meant to be modified by the user. Please prefer the provided API.
"""
mutable struct MultilayerGraph{T,U} <: AbstractMultilayerGraph{T,U}
    layers::Vector{LayerDescriptor{T,U}} # vector containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
    interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}} # the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
    v_V_associations::Bijection{T,<:MultilayerVertex} # A Bijection from Bijections.jl that associates numeric vertices to `MultilayerVertex`s.
    idx_N_associations::Bijection{Int64,Node,Dict{Int64,Node},Dict{Node,Int64}} # A Bijection from Bijections.jl that associates Int64 to `Node`s.
    fadjlist::Vector{Vector{HalfEdge{<:MultilayerVertex,<:Union{Nothing,U}}}} # the forward adjacency list of the MultilayerGraph. It is a vector of vectors of `HalfEdge`s. Its i-th element are the `HalfEdge`s that originate from `v_V_associations[i]`.
    v_metadata_dict::Dict{T,<:Union{<:Tuple,<:NamedTuple}} #  A Dictionary that associates numeric vertices to their metadata.
end

# Traits
@traitimpl IsWeighted{MultilayerGraph}
@traitimpl IsMeta{MultilayerGraph}

"""
    is_directed(m::M) where { M <: Type{ <: MultilayerGraph}}

Return `false`
"""
Graphs.is_directed(mg::M) where {M<:Type{<:MultilayerGraph}} = false

# Constructors
"""
    MultilayerGraph(
        layers::Vector{<:Layer{T,U}},
        specified_interlayers::Vector{<:Interlayer{T,U}};
        default_interlayers_null_graph::H = SimpleGraph{T}(),
        default_interlayers_structure::String="multiplex",
    ) where {T,U, H <: AbstractGraph{T}}   

Construct a MultilayerGraph with layers given by `layers`. The interlayers will be constructed by default according to `default_interlayer` where only `"multiplex"` and `"empty"` are allowed, except for those specified in `specified_interlayers`. `default_interlayer = "multiplex"` will imply that unspecified interlayers will have only diagonal couplings, while  `default_interlayer = "multiplex"` will produced interlayers that have no couplings.

# ARGUMENTS

- `layers::Vector{<:Layer{T,U}}`: The (ordered) list of layers the multilayer graph will have;
- `specified_interlayers::Vector{<:Interlayer{T,U}}`: The list of interlayers specified by the user. Note that the user does not need to specify all interlayers, as the unspecified ones will be automatically constructed using the indications given by the `default_interlayers_null_graph` and `default_interlayers_structure` keywords;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: Sets the underlying graph for the interlayers that are to be automatically specified. Defaults to `SimpleGraph{T}()`. See the `Interlayer` constructors for more information;
- `default_interlayers_structure::String = "multiplex"`: Sets the structure of the interlayers that are to be automatically specified. May be "multiplex" for diagonally coupled interlayers, or "empty" for empty interlayers (no edges).  "multiplex". See the `Interlayer` constructors for more information.
"""
function MultilayerGraph(
    layers::Vector{<:Layer{T,U}},
    specified_interlayers::Vector{<:Interlayer{T,U}};
    default_interlayers_null_graph::H=SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U,H<:AbstractGraph{T}}
    multilayergraph = MultilayerGraph(T, U)

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
    MultilayerGraph(
        empty_layers::Vector{<:Layer{T,U}},
        empty_interlayers::Vector{<:Interlayer{T,U}},
        degree_distribution::UnivariateDistribution;
        allow_self_loops::Bool = false,
        default_interlayers_null_graph::H = SimpleGraph{T}(),
    ) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

Return a random MultilayerGraph that has `empty_layers` as layers and `empty_interlayers` as specified interlayers. `empty_layers` and `empty_interlayers` must respectively be `Layer`s and `Interlayer`s with whatever number of vertices but no edges (if any edge is found, an error is thrown). The  degree distribution of the returned random `MultilayerGraph` is given by `degree_distribution`, which must have a support that only contains positive numbers for obvious reasons. `allow_self_loops = true` allows for self loops t be present in the final random MultilayerGraph. `default_interlayers_null_graph` controls the `null_graph` argument passed to automatically-generated interlayers. 
"""
function MultilayerGraph(
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

    empty_multilayergraph = MultilayerGraph(
        empty_layers,
        empty_interlayers;
        default_interlayers_null_graph=default_interlayers_null_graph,
        default_interlayers_structure="empty",
    )

    n = nv(empty_multilayergraph)

    degree_sequence = sample_graphical_degree_sequence(degree_distribution, n)

    return MultilayerGraph(
        empty_multilayergraph, degree_sequence; allow_self_loops=false, perform_checks=false
    )
end

"""
    MultilayerGraph(empty_multilayergraph::MultilayerGraph{T,U},
    degree_sequence::Vector{<:Integer}; 
    allow_self_loops::Bool = false,
    perform_checks::Bool = false) where {T,U}

Return a random `MultilayerGraph` with degree sequence `degree_sequence`. `allow_self_loops` controls the presence of self-loops, while if `perform_checks` is true, the `degree_sequence` os checked to be graphical.
"""
function MultilayerGraph(
    empty_multilayergraph::MultilayerGraph{T,U},
    degree_sequence::Vector{<:Integer};
    allow_self_loops::Bool=false,
    perform_checks::Bool=true,
) where {T,U}
    (allow_self_loops && perform_checks) &&
        @warn "Checks for graphicality and coherence with the provided `empty_multilayergraph` are currently performed without taking into account self-loops. Thus said checks may fail even though the provided `degree_sequence` may be graphical when one allows for self-loops within the multilayer graph to be present. If you are sure that the provided `degree_sequence` is indeed graphical under those circumstances, you may want to disable checks by setting `perform_checks = false`. We apologize for the inconvenient."

    _multilayergraph = deepcopy(empty_multilayergraph)

    ne(_multilayergraph) == 0 || throw(
        ErrorException(
            "The `empty_multilayergraph` argument should be an empty MultilayerGraph. Found $(ne(_multilayergraph)) edges.",
        ),
    )

    if perform_checks
        n = nv(_multilayergraph)
        n == length(degree_sequence) || throw(
            ErrorException(
                "The number of vertices of the provided empty MultilayerGraph does not match the length of the degree sequence. Found $(nv(_multilayergraph)) and $(length(degree_sequence)).",
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
    MultilayerGraph(T::Type{<:Number}, U::Type{<:Number})

Return a null MultilayerGraph with with vertex type `T` weighttype `U`. Use this constructor and then add Layers and Interlayers via the `add_layer!` and `specify_interlayer!` methods.
"""
function MultilayerGraph(T::Type{<:Number}, U::Type{<:Number})
    return MultilayerGraph{T,U}(
        LayerDescriptor{T,U}[],
        OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}}(),
        Bijection{T,MultilayerVertex}(),
        Bijection{Int64,Node}(),
        Vector{HalfEdge{MultilayerVertex,<:Union{Nothing,U}}}[],
        Dict{T,Union{Tuple,NamedTuple}}(),
    )
end

"""
    MultilayerGraph(layers::Vector{<:Layer{T,U}}; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String="multiplex") where {T,U, H <: AbstractGraph{T}}

Construct a MultilayerGraph with layers `layers` and all interlayers with structure `default_interlayers_structure` (only "multiplex" and "empty" are allowed) and type `default_interlayers_null_graph`.
"""
function MultilayerGraph(
    layers::Vector{<:Layer{T,U}};
    default_interlayers_null_graph::H=SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U,H<:AbstractGraph{T}}
    return MultilayerGraph(
        layers,
        Interlayer{T,U}[];
        default_interlayers_null_graph=default_interlayers_null_graph,
        default_interlayers_structure=default_interlayers_structure,
    )
end

# Nodes

"""
    add_node!(mg::MultilayerGraph, n::Node; add_vertex_to_layers::Union{Vector{Symbol}, Symbol} = Symbol[])

Add node `n` to `mg`. Return true if succeeds. Additionally, add a corresponding vertex to all layers whose name is listed in `add_vertex_to_layers`. If `add_vertex_to_layers == :all`, then a corresponding vertex is added to all layers.
"""
function add_node!(
    mg::MultilayerGraph,
    n::Node;
    add_vertex_to_layers::Union{Vector{Symbol},Symbol}=Symbol[],
)
    return _add_node!(mg, n; add_vertex_to_layers=add_vertex_to_layers)
end

"""
    rem_node!(mg::MultilayerGraph, n::Node)

Remove node `n` to `mg`. Return true if succeeds.
"""
rem_node!(mg::MultilayerGraph, n::Node) = _rem_node!(mg, n)

# Vertices

"""
    add_vertex!(mg::MultilayerGraph, mv::MultilayerVertex; add_node::Bool = true)

Add MultilayerVertex `mv` to multilayer graph `mg`. If `add_node` is true and `node(mv)` is not already part of `mg`, then add `node(mv)` to `mg` before adding `mv` to `mg` instead of throwing an error.
"""
function Graphs.add_vertex!(mg::MultilayerGraph, mv::MultilayerVertex; add_node::Bool=true)
    return _add_vertex!(mg, mv; add_node=add_node)
end

"""
    rem_vertex!(mg::MultilayerGraph, V::MultilayerVertex)

Remove [MultilayerVertex](@ref) `mv` from `mg`. Return true if succeeds, false otherwise.
"""
Graphs.rem_vertex!(mg::MultilayerGraph, V::MultilayerVertex) = _rem_vertex!(mg, V)

# Edges

"""
    add_edge!(mg::M, me::E) where {T,U, M <: MultilayerGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add a MultilayerEdge between `src` and `dst` with weight `weight` and metadata `metadata`. Return true if succeeds, false otherwise.
"""
function Graphs.add_edge!(
    mg::M, me::E
) where {T,U,M<:MultilayerGraph{T,U},E<:MultilayerEdge{<:Union{U,Nothing}}}
    return _add_edge!(mg, me)
end

"""
    rem_edge!(mg::MultilayerGraph, me::MultilayerEdge)

Remove edge from `src(me)` to `dst(me)` from `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.rem_edge!(mg::MultilayerGraph, src::MultilayerVertex, dst::MultilayerVertex)
    return _rem_edge!(mg, src, dst)
end

# Layers and Interlayers
"""
    add_layer!(
        mg::M,
        new_layer::L;
        default_interlayers_null_graph::H=SimpleGraph{T}(),
        default_interlayers_structure::String="multiplex",
    ) where {
        T,
        U,
        G<:AbstractGraph{T},
        L<:Layer{T,U,G},
        H<:AbstractGraph{T},
        M<:MultilayerGraph{T,U}
    }

Add layer `layer` to `mg`.

# ARGUMENTS

- `mg::M`: the `MultilayerGraph` which the new layer will be added to;
- `new_layer::L`: the new `Layer` to add to `mg`;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: upon addition of a new `Layer`, all the `Interlayer`s between the new and the existing `Layer`s are immediately created. This keyword argument specifies their `null_graph` See the `Layer` constructor for more information. Defaults to `SimpleGraph{T}()`;
- `default_interlayers_structure::String = "multiplex"`: The structure of the `Interlayer`s created by default. May either be "multiplex" to have diagonally-coupled only interlayers, or "empty" for empty interlayers. Defaults to "multiplex".
"""
function add_layer!(
    mg::M,
    new_layer::L;
    default_interlayers_null_graph::H=SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {
    T,U,G<:AbstractGraph{T},L<:Layer{T,U,G},H<:AbstractGraph{T},M<:MultilayerGraph{T,U}
}
    return add_layer_directedness!(
        mg,
        new_layer;
        default_interlayers_null_graph=default_interlayers_null_graph,
        default_interlayers_structure=default_interlayers_structure,
    )
end

"""
    specify_interlayer!(
        mg::M,
        new_interlayer::In
    ) where {T,U,G<:AbstractGraph{T},In<:Interlayer{T,U,G}, M<:MultilayerGraph{T,U}; !IsDirected{M}}

Specify the interlayer `new_interlayer` as part of `mg`.
"""
@traitfn function specify_interlayer!(
    mg::M, new_interlayer::In
) where {
    T,U,G<:AbstractGraph{T},In<:Interlayer{T,U,G},M<:MultilayerGraph{T,U};!IsDirected{M}
} # and(!istrait(IsDirected{M}), !istrait(IsMultiplex{M}))
    !is_directed(new_interlayer.graph) || throw( # !istrait(IsDirected{typeof(new_interlayer.graph)})
        ErrorException(
            "The `new_interlayer`'s underlying graphs $(new_interlayer.graph) is directed, so it is not compatible with a `MultilayerGraph`.",
        ),
    )

    return _specify_interlayer!(mg, new_interlayer;)
end
