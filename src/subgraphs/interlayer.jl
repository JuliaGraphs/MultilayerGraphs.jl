# CUSTOM TYPE/CLASS
# This file defines the custom type `Interlayer` and straightforwardly makes it compatible with the Graphs.jl ecosystem. The reason to have this custom type is to have a way not to break other people's code when modification to interlayers are required.

"""
    AbstractInterlayer{T,U,G}

An abstract type representing a generic Interlayer.

# PARAMETRIC TYPES

- `T`: the node type;
- `U`: the adjacency matrix/tensor eltype;
- `G`: the underlying graph type.
"""
abstract type AbstractInterlayer{T,U,G} <: AbstractSubGraph{T,U,G} end

"""
    Interlayer{T<:Integer,U<:Real,G<:AbstractGraph{T}} <: AbstractInterlayer{T,U,G}

Represents an interlayer in a `Multilayer(Di)Graph`. Its type hierarchy is: Interlayer <: AbstractInterlayer <: AbstractSubGraph .
"""
mutable struct Interlayer{T<:Integer,U<:Real,G<:AbstractGraph{T}} <:
               AbstractInterlayer{T,U,G}
    descriptor::InterlayerDescriptor{T,U,G}
    graph::G
    v_V_associations::Bijection{T,<:MultilayerVertex}
end

# Old inner constructor that has been removed in favor of the current inner constructor since not all graph packages implement == for the concrete graphs they define
"""
    _Interlayer(
        layer_1_multilayervertices::Vector{<: MultilayerVertex},
        layer_2_multilayervertices::Vector{<: MultilayerVertex},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}
        descriptor::InterlayerDescriptor{T,U,G}
        
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Internal constructor used with InterlayerDescriptor.
"""
function _Interlayer(
    layer_1_multilayervertices::Vector{<:MultilayerVertex},
    layer_2_multilayervertices::Vector{<:MultilayerVertex},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    },
    descriptor::InterlayerDescriptor{T,U,G},
) where {T<:Integer,U<:Real,G<:AbstractGraph{T}}
    graph, v_V_associations = _compute_interlayer_graph(
        layer_1_multilayervertices, layer_2_multilayervertices, edge_list, descriptor
    )

    return Interlayer(descriptor, graph, v_V_associations)
end

"""
    _Interlayer(
        layer_1_multilayervertices::Vector{<: MultilayerVertex},
        layer_2_multilayervertices::Vector{<: MultilayerVertex},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}
        descriptor::InterlayerDescriptor{T,U,G}
        
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Internal constructor used with InterlayerDescriptor.
"""
function _Interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    },
    descriptor::InterlayerDescriptor{T,U,G},
) where {T<:Integer,U<:Real,G<:AbstractGraph{T}}
    return _Interlayer(mv_vertices(layer_1), mv_vertices(layer_2), edge_list, descriptor)
end

"""
    Interlayer(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        null_graph::G,
        edge_list::Vector{ <: MultilayerEdge{<: Union{U, Nothing}}};
        default_edge_weight::Function = (x,y) -> nothing,
        default_edge_metadata::Function = (x,y) -> NamedTuple(),
        transfer_vertex_metadata::Bool = false,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Vector{ <: MultilayerEdge{<: Union{U, Nothing}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;
- `null_graph::G`: the Interlayer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors.

"""
function Interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G,
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_weight::Function=(x, y) -> nothing,
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real,G<:AbstractGraph{T}}
    descriptor = InterlayerDescriptor(
        name(layer_1),
        name(layer_2),
        null_graph,
        U;
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        name=interlayer_name,
    )

    return _Interlayer(layer_1, layer_2, edge_list, descriptor)
end

"""
    Interlayer(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64,
        null_graph::G;
        default_edge_weight::Function = (x,y) -> nothing,
        default_edge_metadata::Function = (x,y) -> NamedTuple(),
        interlayer_name::Symbol,
        transfer_vertex_metadata::Bool = false
    ) where {T<:Integer, U <: Union{Nothing, <: Real}, G<:AbstractGraph{T}}

Return a random `Interlayer`.

# ARGUMENTS
- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;
- `null_graph::G`: the Interlayer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors.
"""
function Interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64,
    null_graph::G;
    default_edge_weight::Function=(x, y) -> nothing,
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
    transfer_vertex_metadata::Bool=false,
) where {T<:Integer,U<:Union{Nothing,<:Real},G<:AbstractGraph{T}}
    nv_1 = nv(layer_1)
    nv_2 = nv(layer_2)

    directed = is_directed(null_graph)
    maxe = directed ? 2 * nv_1 * nv_2 : nv_1 * nv_2

    @assert(
        ne <= maxe,
        "The number of required edges, $ne, is greater than the number of edges the provided graph supports i.e. $maxe"
    )

    layer_1_name = name(layer_1)
    layer_2_name = name(layer_2)

    layer_1_multilayervertices = collect(mv_vertices(layer_1))
    layer_2_multilayervertices = collect(mv_vertices(layer_2))
    all_vertices = vcat(layer_1_multilayervertices, layer_2_multilayervertices)

    edge_list = NTuple{2,MultilayerVertex}[]
    fadjlist = Dict{MultilayerVertex,Vector{MultilayerVertex}}()

    for i in 1:ne
        rand_vertex_1 = rand(
            setdiff(
                all_vertices,
                [
                    mv for (mv, _list) in collect(fadjlist) if
                    (layer(mv) == layer_1_name && length(_list) == nv_2) ||
                    (layer(mv) == layer_2_name && length(_list) == nv_1)
                ],
            ),
        )

        if !haskey(fadjlist, rand_vertex_1)
            fadjlist[rand_vertex_1] = MultilayerVertex{layer_2_name}[]
        end

        rand_vertex_2 = if layer(rand_vertex_1) == layer_1_name
            rand(setdiff(layer_2_multilayervertices, fadjlist[rand_vertex_1])) # rand(layer_2_multilayervertices)
        else
            rand(setdiff(layer_1_multilayervertices, fadjlist[rand_vertex_1])) # rand(layer_2_multilayervertices)
        end # rand(layer_2_multilayervertices)

        push!(fadjlist[rand_vertex_1], rand_vertex_2)

        if !is_directed(null_graph)
            if !haskey(fadjlist, rand_vertex_2)
                fadjlist[rand_vertex_2] = MultilayerVertex{layer_1_name}[]
            end
            push!(fadjlist[rand_vertex_2], rand_vertex_1)
        end

        push!(edge_list, (rand_vertex_1, rand_vertex_2))
    end

    edge_list = if istrait(IsDirected{G})
        [rand() < 0.5 ? tup : reverse(tup) for tup in edge_list]#MultilayerEdge
    else
        edge_list#MultilayerEdge
    end#MultilayerEdge

    descriptor = InterlayerDescriptor(
        name(layer_1),
        name(layer_2),
        null_graph,
        U;
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        name=interlayer_name,
    )

    return _Interlayer(
        layer_1_multilayervertices, layer_2_multilayervertices, edge_list, descriptor
    )
end

"""
    multiplex_interlayer(
        nv::Int64,
        interlayer_name::Symbol,
        layer_1::Symbol, 
        layer_2::Symbol, 
        graph_type::Type{G}; 
        forbidden_vertices::Vector{MultilayerVertex}, forbidden_edges::Vector{NTuple{2, MultilayerVertex}}
    ) where {T <: Union{ <: Integer, AbstractVertex}, G <: AbstractGraph{T}}

Return an `Interlayer{T,U,G}` that has edges only between vertices that represent the same node.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `null_graph::G`: the Interlayer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors;
"""
function multiplex_interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G;
    default_edge_weight::Function=(x, y) -> nothing,
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real,G<:AbstractGraph{T}}
    return _multiplex_interlayer(
        mv_vertices(layer_1),
        mv_vertices(layer_2),
        null_graph,
        U;
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    _multiplex_interlayer(
        layer_1_multilayervertices::Vector{MultilayerVertex{L1}},
        layer_2_multilayervertices::Vector{MultilayerVertex{L2}},
        null_graph::G,
        weighttype::Type{U};
        default_edge_weight::Function=(x, y) -> nothing,
        default_edge_metadata::Function=(x, y) -> NamedTuple(),
        transfer_vertex_metadata::Bool=false,
        interlayer_name::Symbol,
    ) where {L1,L2,T<:Integer,U<:Real,G<:AbstractGraph{T}}

Internal method for multiplex_interlayer.
"""
function _multiplex_interlayer(
    layer_1_multilayervertices::Vector{MultilayerVertex{L1}},
    layer_2_multilayervertices::Vector{MultilayerVertex{L2}},
    null_graph::G,
    weighttype::Type{U};
    default_edge_weight::Function=(x, y) -> nothing,
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {L1,L2,T<:Integer,U<:Real,G<:AbstractGraph{T}}
    common_nodes = intersect(
        [mv.node for mv in layer_1_multilayervertices],
        [mv.node for mv in layer_2_multilayervertices],
    )
    edge_list = nothing
    if istrait(IsDirected{typeof(null_graph)})
        _edge_list = [
            MultilayerEdge(
                MultilayerVertex(node, L1),
                MultilayerVertex(node, L2),
                default_edge_weight(MV(node, L1), MV(node, L2)),
                default_edge_metadata(MV(node, L1), MV(node, L2)),
            ) for node in common_nodes
        ]
        edge_list = vcat(_edge_list..., reverse.(_edge_list)...)
    else
        edge_list = [
            MultilayerEdge(
                MultilayerVertex(node, L1),
                MultilayerVertex(node, L2),
                default_edge_weight(MV(node, L1), MV(node, L2)),
                default_edge_metadata(MV(node, L1), MV(node, L2)),
            ) for node in common_nodes
        ]
    end

    descriptor = InterlayerDescriptor(
        L1,
        L2,
        null_graph,
        weighttype;
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        name=interlayer_name,
    )

    return _Interlayer(
        layer_1_multilayervertices, layer_2_multilayervertices, edge_list, descriptor
    )
end

"""
    empty_interlayer(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        null_graph::G;
        default_edge_weight::Function = (x,y) -> nothing,
        default_edge_metadata::Function = (x,y) -> NamedTuple(),
        interlayer_name::Symbol),
        transfer_vertex_metadata::Bool = false
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Construct an empty interlayer (i.e. an interlayer with no edges).

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `null_graph::G`: the Interlayer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors;
"""
function empty_interlayer(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    null_graph::G;
    default_edge_weight::Function=(x, y) -> nothing,
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
    transfer_vertex_metadata::Bool=false,
) where {T<:Integer,U<:Real,G<:AbstractGraph{T}}
    return _empty_interlayer(
        collect(mv_vertices(layer_1)),
        collect(mv_vertices(layer_2)),
        null_graph,
        U;
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    _empty_interlayer(
        layer_1_multilayervertices::Vector{MultilayerVertex{L1}},
        layer_2_multilayervertices::Vector{MultilayerVertex{L2}},
        null_graph::G,
        weighttype::Type{U};
        interlayer_name::Symbol = Symbol("interlayer_(layer_1.name)_(layer_2.name)"),
        transfer_vertex_metadata::Bool = false
    ) where {L1, L2, T<:Integer, U <: Real, G<:AbstractGraph{T}}

Internal method for `empty_interlayer`. Return an `Interlayer{T,U,G}` with vertices  `layer_1_multilayervertices` and `layer_2_multilayervertices` with no edges.
"""
function _empty_interlayer(
    layer_1_multilayervertices::Vector{MultilayerVertex{L1}},
    layer_2_multilayervertices::Vector{MultilayerVertex{L2}},
    null_graph::G,
    weighttype::Type{U};
    default_edge_weight::Function=(x, y) -> nothing,
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {L1,L2,T<:Integer,U<:Real,G<:AbstractGraph{T}}
    edge_list = MultilayerEdge{U}[]
    descriptor = InterlayerDescriptor(
        L1,
        L2,
        null_graph,
        weighttype;
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        name=interlayer_name,
    )

    return _Interlayer(
        layer_1_multilayervertices, layer_2_multilayervertices, edge_list, descriptor
    )
end

"""
    interlayer_simplegraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{<:Union{U, Nothing}}}, Vector{ <:Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `SimpleGraph` from `Graphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simplegraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1, layer_2, SimpleGraph{T}(), edge_list; interlayer_name=interlayer_name
    )
end

"""
    interlayer_simplegraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `SimpleGraph` from `Graphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simplegraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1, layer_2, ne, SimpleGraph{T}(); interlayer_name=interlayer_name
    )
end

"""
    interlayer_simpledigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `SimpleDiGraph` from `Graphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simpledigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1, layer_2, SimpleDiGraph{T}(), edge_list; interlayer_name=interlayer_name
    )
end

"""
    interlayer_simpledigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `SimpleDiGraph` from `Graphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simpledigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1, layer_2, ne, SimpleDiGraph{T}(); interlayer_name=interlayer_name
    )
end

"""
    interlayer_simpleweightedgraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        default_edge_weight::Function=(x, y) -> nothing,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `SimpleWeightedGraph` from `SimpleWeightedGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simpleweightedgraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_weight::Function=(x, y) -> nothing,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        SimpleWeightedGraph{T,U}(),
        edge_list;
        default_edge_weight=default_edge_weight,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_simpleweightedgraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        default_edge_weight::Function=(x, y) -> nothing,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `SimpleWeightedGraph` from `SimpleWeightedGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simpleweightedgraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    default_edge_weight::Function=(x, y) -> nothing,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        ne,
        SimpleWeightedGraph{T,U}();
        default_edge_weight=default_edge_weight,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_simpleweighteddigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        default_edge_weight::Function=(x, y) -> nothing,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `SimpleWeightedDiGraph` from `SimpleWeightedGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simpleweighteddigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_weight::Function=(x, y) -> nothing,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        SimpleWeightedGraph{T,U}(),
        edge_list;
        default_edge_weight=default_edge_weight,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_simpleweighteddigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        default_edge_weight::Function=(x, y) -> nothing,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges underlying graph is a `SimpleWeightedDiGraph` from `SimpleWeightedGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_simpleweighteddigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    default_edge_weight::Function=(x, y) -> nothing,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        ne,
        SimpleWeightedDiGraph{T,U}();
        default_edge_weight=default_edge_weight,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_metagraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        default_edge_metadata::Function=(x, y) -> NamedTuple(),
        transfer_vertex_metadata::Bool=false,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `MetaGraph` from `MetaGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_metagraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        MetaGraph{T,U}(),
        edge_list;
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_metagraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        default_edge_metadata::Function=(x, y) -> NamedTuple(),
        transfer_vertex_metadata::Bool=false,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `MetaGraph` from `MetaGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_metagraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        ne,
        MetaGraph{T,U}();
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_metadigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        default_edge_metadata::Function=(x, y) -> NamedTuple(),
        transfer_vertex_metadata::Bool=false,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `MetaDiGraph` from `MetaGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_metadigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        MetaDiGraph{T,U}(),
        edge_list;
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_metadigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        default_edge_metadata::Function=(x, y) -> NamedTuple(),
        transfer_vertex_metadata::Bool=false,
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `MetaDiGraph` from `MetaGraphs.jl`, with vertex type `T` and weight type `U`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `transfer_vertex_metadata::Bool`:if true, vertex metadata found in both connected layers are carried over to the vertices of the Interlayer. NB: not all choice of underlying graph may support this feature. Graphs types that don't support metadata or that pose limitations to it may result in errors;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_metadigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    transfer_vertex_metadata::Bool=false,
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    return Interlayer(
        layer_1,
        layer_2,
        ne,
        MetaDiGraph{T,U}();
        default_edge_metadata=default_edge_metadata,
        transfer_vertex_metadata=transfer_vertex_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_valgraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        default_edge_metadata::Function=(x, y) -> NamedTuple()
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `ValGraph` from `SimpleValueGraphs.jl`, with vertex type `T`. By default, `transfer_vertex_metadata` is set to `false`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_valgraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValGraph{T}(
        SimpleGraph{T}(); edgeval_types=edgeval_types, edgeval_init=default_edge_metadata
    )

    return Interlayer(
        layer_1,
        layer_2,
        graph,
        edge_list;
        default_edge_metadata=default_edge_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_valgraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        default_edge_metadata::Function=(x, y) -> NamedTuple()
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `ValGraph` from `SimpleValueGraphs.jl`, with vertex type `T`. By default, `transfer_vertex_metadata` is set to `false`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_valgraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValGraph{T}(
        SimpleGraph{T}(); edgeval_types=edgeval_types, edgeval_init=default_edge_metadata
    )

    return Interlayer(
        layer_1,
        layer_2,
        ne,
        graph;
        default_edge_metadata=default_edge_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_valoutdigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        default_edge_metadata::Function=(x, y) -> NamedTuple()
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `ValOutDiGraph` from `SimpleValueGraphs.jl`, with vertex type `T`. By default, `transfer_vertex_metadata` is set to `false`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_valoutdigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValOutDiGraph{T}(
        SimpleDiGraph{T}(); edgeval_types=edgeval_types, edgeval_init=default_edge_metadata
    )

    return Interlayer(
        layer_1,
        layer_2,
        graph,
        edge_list;
        default_edge_metadata=default_edge_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_valoutdigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        default_edge_metadata::Function=(x, y) -> NamedTuple()
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `ValOutDiGraph` from `SimpleValueGraphs.jl`, with vertex type `T`. By default, `transfer_vertex_metadata` is set to `false`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_valoutdigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValOutDiGraph{T}(
        SimpleDiGraph{T}(); edgeval_types=edgeval_types, edgeval_init=default_edge_metadata
    )

    return Interlayer(
        layer_1,
        layer_2,
        ne,
        graph;
        default_edge_metadata=default_edge_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_valdigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}};
        default_edge_metadata::Function=(x, y) -> NamedTuple()
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Constructor for Interlayer whose underlying graph is a `ValDiGraph` from `SimpleValueGraphs.jl`, with vertex type `T`. By default, `transfer_vertex_metadata` is set to `false`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `edge_list::Union{<:Vector{<:MultilayerEdge{ <: Union{U,Nothing}}}, <:Vector{ <: Tuple{<:MultilayerVertex, <:MultilayerVertex}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_valdigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    };
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValDiGraph{T}(
        SimpleDiGraph{T}(); edgeval_types=edgeval_types, edgeval_init=default_edge_metadata
    )

    return Interlayer(
        layer_1,
        layer_2,
        graph,
        edge_list;
        default_edge_metadata=default_edge_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    interlayer_valdigraph(
        layer_1::Layer{T,U},
        layer_2::Layer{T,U},
        ne::Int64;
        default_edge_metadata::Function=(x, y) -> NamedTuple()
        interlayer_name::Symbol
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}}

Return a random Interlayer with `ne` edges whose underlying graph is a `ValDiGraph` from `SimpleValueGraphs.jl`, with vertex type `T`. By default, `transfer_vertex_metadata` is set to `false`.

# ARGUMENTS

- `layer_1::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `layer_2::Layer{T,U}`: one of the two layers connected by the Interlayer;
- `ne::Int64`: The number of edges of the Interlayer;

# KWARGS

- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `interlayer_name::Symbol`: The name of the Interlayer. Defaults to Symbol("interlayer_(layer_1.name)_(layer_2.name)");

"""
function interlayer_valdigraph(
    layer_1::Layer{T,U},
    layer_2::Layer{T,U},
    ne::Int64;
    default_edge_metadata::Function=(x, y) -> NamedTuple(),
    interlayer_name::Symbol=Symbol("interlayer_$(layer_1.name)_$(layer_2.name)"),
) where {T<:Integer,U<:Real}
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValDiGraph{T}(
        SimpleDiGraph{T}(); edgeval_types=edgeval_types, edgeval_init=default_edge_metadata
    )

    return Interlayer(
        layer_1,
        layer_2,
        ne,
        graph;
        default_edge_metadata=default_edge_metadata,
        interlayer_name=interlayer_name,
    )
end

"""
    is_interlayer_bipartite(graph::G, v_V_associations::Bijection{T, MultilayerVertex}) where {T, G <: AbstractGraph{T}}

Check if `graph`, whose vertices are interpreted via `v_V_associations` can be the underlying bipartite graph of an Interlayer graph.
"""
function is_interlayer_bipartite(
    graph::G, v_V_associations::Bijection{T,MultilayerVertex}
) where {T,G<:AbstractGraph{T}}
    layers = unique([mv.layer for mv in image(v_V_associations)])
    # If it is an empty interlayer, we assume it is bipartite
    length(layers) == 0 && return true
    length(layers) == 2 || throw(
        ErrorException(
            "The interlayer cannot be bipartite since more than two layers are involved. Found $layers in v_V_associations.",
        ),
    )
    # This loop iterates through all edges of the graph
    for edge in edges(graph)
        # This if statement checks if the source and destination of the edge are on the same layer
        if !(v_V_associations[src(edge)].layer != v_V_associations[dst(edge)].layer)
            # If they are, the graph is not a valid layering
            return false
        end
    end
    return true
end

"""
    is_multiplex_interlayer(interlayer::In) where {In <: Interlayer}

Check that Interlayer `interlayer` is a multiplex-type Interlayer.
"""
function is_multiplex_interlayer(interlayer::Interlayer)
    if is_directed(interlayer)
        for node in intersect(interlayer.layer_1_nodes, interlayer.layer_2_nodes)
            (
                has_edge(
                    interlayer, MV(node, interlayer.layer_1), MV(node, interlayer.layer_2)
                ) && has_edge(
                    interlayer, MV(node, interlayer.layer_2), MV(node, interlayer.layer_1)
                )
            ) || return false
        end
        return true
    else
        # Loop through all nodes in both layers.
        for node in intersect(interlayer.layer_1_nodes, interlayer.layer_2_nodes)
            # Check if there is an edge between the nodes in both layers.
            has_edge(
                interlayer, MV(node, interlayer.layer_1), MV(node, interlayer.layer_2)
            ) || return false
        end
        return true
    end
end

"""
    has_node( interlayer::Interlayer, n::Node )

Return `true` if `n` is a `Node` of `interlayer`.
"""
function has_node(interlayer::Interlayer, n::Node)
    return n in nodes(interlayer.layer_1) || n in nodes(interlayer.layer_2)
end

"""
    has_vertex(interlayer::Interlayer, v::MultilayerVertex)

Return `true` if `v` is a vertex of `interlayer`.
"""
function Graphs.has_vertex(interlayer::Interlayer, mv::MultilayerVertex)
    return get_bare_mv(mv) ∈ collect(image(interlayer.v_V_associations))
end

"""
    add_edge!(interlayer::Interlayer, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...)

Add edge from vertex `src` to vertex `dst` to Interlayer `interlayer`.Returns true if succeeds, but will fail (return false) if `src` and `dst` belong to the same layer, since interlayers are always bipartite. This method supports the uniform and transparent interfaces. See the [Edges](@ref edges_tut_subg) section of the Tutorial.
"""
function Graphs.add_edge!(
    interlayer::Interlayer, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...
)
    src_bare = get_bare_mv(src)
    dst_bare = get_bare_mv(dst)
    !has_vertex(interlayer, src) &&
        throw(ErrorException("Vertex $(src) does not belong to the interlayer."))
    !has_vertex(interlayer, dst) &&
        throw(ErrorException("Vertex $(dst) does not belong to the interlayer."))
    src_bare.layer != dst_bare.layer || return false

    if !has_edge(interlayer, src_bare, dst_bare)
        success = false
        if isempty(args) &&
            length(kwargs) == 2 &&
            issetequal(Set([:weight, :metadata]), Set(keys(kwargs)))
            success = add_edge_standard!(
                interlayer,
                src_bare,
                dst_bare;
                weight=values(kwargs).weight,
                metadata=values(kwargs).metadata,
            )
        elseif length(args) == length(kwargs) == 0
            success = add_edge_standard!(
                interlayer,
                src_bare,
                dst_bare;
                weight=interlayer.default_edge_weight(src_bare, dst_bare),
                metadata=interlayer.default_edge_metadata(src_bare, dst_bare),
            )
        else
            success = add_edge!(
                interlayer.graph,
                interlayer.v_V_associations(src_bare),
                interlayer.v_V_associations(dst_bare),
                args...;
                kwargs...,
            )
        end
        return success
    else
        return false
    end
end

# Base overloads
"""
    Base.(==)(x::Interlayer, y::Interlayer)

Overload equality for `Interlayer`s.
"""
function Base.:(==)(x::Interlayer, y::Interlayer)
    for field in fieldnames(Interlayer)
        if @eval $x.$field != $y.$field
            return false
        end
    end
    return true
end

"""
    Base.getproperty(interlayer::In, f::Symbol) where { In <: Interlayer}
"""
function Base.getproperty(interlayer::In, f::Symbol) where {In<:Interlayer}
    if f ∈ (:descriptor, :graph, :v_V_associations)
        Base.getfield(interlayer, f)
    elseif f ∈ (
        :name,
        :layer_1,
        :layer_2,
        :null_graph,
        :default_edge_weight,
        :default_edge_metadata,
        :transfer_vertex_metadata,
    )
        Base.getfield(interlayer.descriptor, f)
    elseif f == :edge_list
        edges(interlayer)
    elseif f == :layers_names
        [interlayer.layer_1, interlayer.layer_2]
    elseif f == :layer_1_nodes
        [mv.node for mv in interlayer.layer_1_bare_multilayer_vertices]
    elseif f == :layer_2_nodes
        [mv.node for mv in interlayer.layer_2_bare_multilayer_vertices]
    elseif f == :layer_1_bare_multilayer_vertices
        pairs_sorted = sort(collect(interlayer.v_V_associations); by=first)
        [mv for (v, mv) in pairs_sorted if mv.layer == interlayer.layer_1]
    elseif f == :layer_2_bare_multilayer_vertices
        pairs_sorted = sort(collect(interlayer.v_V_associations); by=first)
        [mv for (v, mv) in pairs_sorted if mv.layer == interlayer.layer_2]
    elseif f == :multilayer_vertices
        vcat(
            collect(interlayer.layer_1_bare_multilayer_vertices),
            collect(interlayer.layer_2_bare_multilayer_vertices),
        )
    elseif f == :V_v_associations
        active_inv(interlayer.v_V_associations)
    elseif f == :v_V_associations
        Base.getfield(interlayer, :v_V_associations)
    else
        Base.getfield(interlayer, f)
    end
end

"""
    get_symmetric_interlayer(interlayer::In; symmetric_interlayer_name::String) where{T,U,G, In <: Interlayer{T,U,G}}

Return the `Interlayer` corresponding to `interlayer` where `layer_1` and `layer_2` are swapped. Its interlayer_name will be `symmetric_interlayer_name` (defaults to `interlayer_(interlayer.layer_2)_(interlayer.layer_1)`).
"""
function get_symmetric_interlayer(
    interlayer::In; symmetric_interlayer_name::String=String(interlayer.name) * "_rev"
) where {T,U,G,In<:Interlayer{T,U,G}}
    # Create a symmetric interlayer descriptor
    symmetric_descriptor = InterlayerDescriptor(
        Symbol(symmetric_interlayer_name),
        interlayer.layer_2,
        interlayer.layer_1,
        interlayer.null_graph,
        interlayer.default_edge_weight,
        interlayer.default_edge_metadata,
        interlayer.transfer_vertex_metadata,
        U,
    )
    # Create a symmetrized interlayer with the symmetric descriptor
    symmetrized_interlayer = Interlayer(
        symmetric_descriptor, interlayer.graph, interlayer.v_V_associations
    )
    # Recompute the interlayer, using the layer 1 and layer 2 vertices from the original interlayer
    recompute_interlayer!(
        symmetrized_interlayer,
        symmetrized_interlayer.layer_1_bare_multilayer_vertices,
        symmetrized_interlayer.layer_2_bare_multilayer_vertices,
    )
    # Return the symmetrized interlayer
    return symmetrized_interlayer
end

"""
    _compute_interlayer_graph(mvs_layer_1::Vector{MultilayerVertex{L1}}, mvs_layer_2::Vector{MultilayerVertex{L2}}, edge_list::Vector{ <: MultilayerEdge}, descriptor::InterlayerDescriptor{T,U,G}) where {L1, L2, T, U,  G <: AbstractGraph{T}} 

Compute the interlayer between the layer's vertex array `mvs_layer_1` and the layer's vertex array `mvs_layer_2`.  
"""
function _compute_interlayer_graph(
    mvs_layer_1::Vector{MultilayerVertex{L1}},
    mvs_layer_2::Vector{MultilayerVertex{L2}},
    edge_list::Union{
        <:Vector{<:MultilayerEdge{<:Union{U,Nothing}}},
        <:Vector{<:Tuple{<:MultilayerVertex,<:MultilayerVertex}},
    },
    descriptor::InterlayerDescriptor{T,U,G},
) where {L1,L2,T,U,G<:AbstractGraph{T}}
    (L1 != L2) || throw(
        ErrorException(
            "Expected all vertices in mvs_layer_1 (and mvs_layer_2) to belong to two different layers.",
        ),
    )
    L1 == descriptor.layer_1 || throw(
        ErrorException(
            "Expected vertices in mvs_layer_1 to belong to layer $(descriptor.layer_1). Found $L1.",
        ),
    )
    L2 == descriptor.layer_2 || throw(
        ErrorException(
            "Expected vertices in mvs_layer_2 to belong to layer $(descriptor.layer_2). Found $L2.",
        ),
    )

    graph = deepcopy(descriptor.null_graph)
    v_V_associations = Bijection{T,MultilayerVertex}()

    for (v, mv) in enumerate(vcat(mvs_layer_1, mvs_layer_2))
        if descriptor.transfer_vertex_metadata
            __add_vertex!(graph; metadata=mv.metadata)
        else
            __add_vertex!(graph)
        end

        v_V_associations[v] = MultilayerVertex(mv.node, mv.layer)
    end

    if typeof(edge_list) <: Vector{<:MultilayerEdge}
        for edge in edge_list
            _src = v_V_associations(get_bare_mv(src(edge)))
            _dst = v_V_associations(get_bare_mv(dst(edge)))
            if !has_edge(graph, _src, _dst)
                _add_edge!(graph, _src, _dst; weight=weight(edge), metadata=metadata(edge))
            end
        end
    else
        for tup in edge_list
            src_mv, dst_mv = tup
            _src = v_V_associations(get_bare_mv(src_mv))
            _dst = v_V_associations(get_bare_mv(dst_mv))
            _add_edge!(
                graph,
                _src,
                _dst;
                weight=descriptor.default_edge_weight(src_mv, dst_mv),
                metadata=descriptor.default_edge_metadata(src_mv, dst_mv),
            )
        end
    end
    return graph, v_V_associations
end

"""
    recompute_interlayer!(interlayer::In, mvs_layer_1::Vector{<: MultilayerVertex}, mvs_layer_2::Vector{<: MultilayerVertex} ) where {T,U, In <: Interlayer{T,U}}

Recompute the interlayer `interlayer` between the layer's vertex array `mvs_layer_1` and the layer's vertex array `mvs_layer_2`.      
"""
function recompute_interlayer!(
    interlayer::In,
    mvs_layer_1::Vector{<:MultilayerVertex},
    mvs_layer_2::Vector{<:MultilayerVertex},
) where {T,U,In<:Interlayer{T,U}}
    graph, v_V_associations = _compute_interlayer_graph(
        mvs_layer_1,
        mvs_layer_2,
        Vector{MultilayerEdge{U}}(collect(edges(interlayer))),
        interlayer.descriptor,
    )
    interlayer.graph = graph
    return interlayer.v_V_associations = v_V_associations
end

# Console print utilities
function to_string(x::Interlayer)
    parameters = typeof(x).parameters
    return """
           Interlayer\t$(name(x))
           layer_1: $(x.layer_1)
           layer_2: $(x.layer_2)
           underlying graph: $(typeof(x.graph))
           vertex type : $(parameters[1]) 
           weight type : $(parameters[2]) 
           nv : $(nv(x))
           ne : $(ne(x))
           """
end
Base.show(io::IO, x::Interlayer) = print(io, to_string(x))
