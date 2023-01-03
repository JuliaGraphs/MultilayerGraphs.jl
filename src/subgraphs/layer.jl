# CUSTOM TYPE/CLASS
# This file defines the custom type `Layer` and straightforwardly makes it compatible with the Graphs.jl ecosystem. The reason to have this custom type is to have a way not to break other people's code when modification to layers are required.

"""
    AbstractLayer{T,U,G}

An abstract type representing a generic Layer.

# FIELDS

- `T`: the node type;
- `U`: the `MultilayerEdge` weight eltype;
- `G`: the underlying graph type.
"""
abstract type AbstractLayer{T,U,G} <: AbstractSubGraph{T,U,G} end

"""
    mutable struct Layer{T <: Integer, U <: Real, G <: AbstractGraph{T}} <: AbstractLayer{T,U,G}

Represents a layer in a `Multilayer(Di)Graph`. Its type hierarchy is: Layer <: AbstractLayer <: AbstractSubGraph .
"""
mutable struct Layer{T<:Integer,U<:Real,G<:AbstractGraph{T}} <: AbstractLayer{T,U,G}
    descriptor::LayerDescriptor{T,U,G}
    graph::G
    v_V_associations::Bijection{T,<:MultilayerVertex}
    # Inner constructor that performs checks on request
    function Layer(
        descriptor::LayerDescriptor{T,U,G},
        graph::G,
        v_V_associations::Bijection{T,<:MultilayerVertex};
        check_consistency=true,
    ) where {T,U,G}
        if check_consistency
            # Check that the graph type is the same as the one in the descriptor
            typeof(descriptor.null_graph) == typeof(graph) || throw(
                ErrorException(
                    "Graph types between the provided `descriptor` and `graph` cannot differ. Found $(typeof(descriptor.null_graph)) and $(typeof(graph)).",
                ),
            )
            # Check that the graph vertices are the same as the ones in the association
            all(vertices(graph) .== sort(domain(v_V_associations))) || throw(
                ErrorException(
                    "The graph has a different set of vertices w.r.t. the domain of `v_V_associations`. Found $(vertices(graph)) and $(sort(domain(v_V_associations))).",
                ),
            )
        end
        return new{T,U,G}(descriptor, graph, v_V_associations)
    end
end

"""
    Layer(
        descriptor::LayerDescriptor{T}, 
        vertices::Union{<:Vector{<:MultilayerVertex}, Vector{Node}}, 
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}) where {T <: Integer}

Constructor for `Layer`.

# ARGUMENTS

- `descriptor::LayerDescriptor{T}`;
- `vertices::Union{Vector{<:MultilayerVertex}, Vector{<:Node}}`;
- `edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2, MultilayerVertex{nothing}}}}`;
"""
function Layer(
    descriptor::LayerDescriptor{T},
    vertices::Union{<:Vector{<:MultilayerVertex},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}},
) where {T<:Integer}
    # First check that the vertices are of the correct type
    if typeof(vertices) <: Vector{<:MultilayerVertex}
        par = eltype(vertices).parameters[1]
        (isnothing(par) || par == descriptor.name) || throw(
            ErrorException(
                "`vertices` should be a `Vector{MultilayerVertex{:$(descriptor.name)}}` or a `Vector{MultilayerVertex{nothing}}` or a `Vector{Node}}`. Found $(typeof(vertices))",
            ),
        )
    end
    # Create the layer
    layer = Layer(
        descriptor,
        deepcopy(descriptor.null_graph),
        Bijection{T,MultilayerVertex{descriptor.name}}();
        check_consistency=false,
    )
    # Add the vertices one by one
    if typeof(vertices) <: Vector{<:MultilayerVertex}
        for mv in vertices
            add_vertex!(layer, mv)
        end
    else
        for node in vertices
            add_vertex!(layer, MV(node, descriptor.default_vertex_metadata(MV(node))))
        end
    end
    # Add the edges
    if typeof(edge_list) <: Vector{<:MultilayerEdge}
        for edge in edge_list
            add_edge!(layer, edge)
        end
    else
        for tup in edge_list
            src_mv, dst_mv = tup
            add_edge!(
                layer,
                ME(
                    src_mv,
                    dst_mv,
                    descriptor.default_edge_weight(src_mv, dst_mv),
                    descriptor.default_edge_metadata(src_mv, dst_mv),
                ),
            )
        end
    end
    return layer
end

"""
    Layer(
        name::Symbol, 
        vertices::Vector{<: MultilayerVertex}, 
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}, 
        null_graph::G, 
        weighttype::Type{U};  
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_weight::Function = (src, dst) -> one(U), d
        default_edge_metadata::Function = (src, dst) -> NamedTuple()) where {T <: Integer, U <: Real,  G <: AbstractGraph{T}}

Constructor for `Layer`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;
- `null_graph::G`: the Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown;
- `weighttype::Type{U}`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function Layer(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}},
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_weight::Function=(src, dst) -> one(U),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
) where {T<:Integer,U<:Real,G<:AbstractGraph{T}}
    descriptor = LayerDescriptor(
        name,
        null_graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
    )

    return Layer(descriptor, vertices, edge_list)
end

"""
    Layer(
        name::Symbol,
        vertices::Union{V, N},
        ne::Int64,
        null_graph::G,
        weighttype::Type{U};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_weight::Function = (src, dst) -> nothing,
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        allow_self_loops::Bool = false
    ) where {T<:Integer,U<:Real,G<:AbstractGraph{T}, V <: Vector{MultilayerVertex{nothing}}, N <: Vector{Node}}

Return a random `Layer`.

# ARGUMENTS

- `vertices::Union{V, N}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `name::Symbol`: The name of the Layer
- `ne::Int64`: The number of edges of the Layer
- `null_graph::G`: the Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
- `weighttype::Type{U}`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted);

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
- `allow_self_loops::Bool`: whether to allow self loops to be generated or not. Defaults to `false`.
"""
function Layer(
    name::Symbol,
    vertices::Vector{Union{V,N}},#Union{V, N},
    ne::Int64,
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_weight::Function=(src, dst) -> nothing,
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    allow_self_loops::Bool=false,
) where {T<:Integer,U<:Real,G<:AbstractGraph{T},V<:MultilayerVertex{nothing},N<:Node} #V <: Vector{MultilayerVertex{nothing}}, N <: Vector{Node}}
    _nv = length(vertices)
    @assert(
        length(unique(vertices)) == _nv, "The argument `vertices` must be a unique list"
    )

    directed = is_directed(null_graph)
    maxe = directed ? _nv * (_nv - 1) : _nv * (_nv - 1) ÷ 2

    @assert(
        ne <= maxe,
        "The number of required edges, $ne, is greater than the number of edges the provided graph supports i.e. $maxe"
    )

    vertex_type = @isdefined(V) ? MultilayerVertex : Node

    edge_list = NTuple{2,MultilayerVertex{nothing}}[]  #MultilayerEdge
    fadjlist = Dict{vertex_type,Vector{vertex_type}}()

    max_links_per_vertex = _nv - 1 # directed ? _nv - 1 : _nv-1
    for i in 1:ne
        # Generate a random vertex
        rand_vertex_1 = rand(
            setdiff(
                vertices,
                [
                    mv for
                    mv in keys(fadjlist) if length(fadjlist[mv]) == max_links_per_vertex
                ],
            ),
        )

        if !haskey(fadjlist, rand_vertex_1)
            fadjlist[rand_vertex_1] = vertex_type[]
        end

        # Generate another random vertex
        # If we don't allow self loops, keep generating until we get a vertex that isn't the same as the previous one.
        rand_vertex_2 = if !allow_self_loops
            rand(setdiff(vertices, vcat([rand_vertex_1], fadjlist[rand_vertex_1])))
        else
            rand(setdiff(vertices, vcat(fadjlist[rand_vertex_1])))
        end

        push!(fadjlist[rand_vertex_1], rand_vertex_2)
        if !is_directed(null_graph)
            if !haskey(fadjlist, rand_vertex_2)
                fadjlist[rand_vertex_2] = vertex_type[]
            end
            push!(fadjlist[rand_vertex_2], rand_vertex_1)
        end

        # Add the edge to the edge list. Convert it to a MultilayerVertex if a list of Nodes was given as `vertices`
        if @isdefined(V)
            push!(edge_list, (rand_vertex_1, rand_vertex_2))
        else
            push!(edge_list, (MV(rand_vertex_1), MV(rand_vertex_2)))
        end
    end

    descriptor = LayerDescriptor(
        name,
        null_graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
    )

    layer = Layer(descriptor, vertices, edge_list)

    return layer
end

"""
    Layer(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        degree_distribution::UnivariateDistribution,
        null_graph::G,
        weighttype::Type{U};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_weight::Function = (src, dst) -> nothing,
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}; !IsDirected{G}}

Returns an undirected `Layer` whose degree sequence is sampled from `degree_distribution`. Graph realization is performed using the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `degree_distribution::UnivariateDistribution`: The distribution to sample degrees from.
- `null_graph::G`: the Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
- `weighttype::Type{U}`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted);

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
"""
@traitfn function Layer(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    degree_distribution::UnivariateDistribution,
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_weight::Function=(src, dst) -> nothing,
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
) where {T<:Integer,U<:Real,G<:AbstractGraph{T};!IsDirected{G}}
    degree_sequence = sample_graphical_degree_sequence(
        degree_distribution, length(vertices)
    )

    return Layer(
        name,
        vertices,
        degree_sequence,
        null_graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    Layer(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        degree_sequence::Vector{<:Integer},
        null_graph::G,
        weighttype::Type{U};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_weight::Function = (src, dst) -> nothing,
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}; !IsDirected{G}}

Returns an undirected `Layer` with given `degree_sequence` realized using the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `degree_sequence::Vector{<:Integer}`: The degree sequence of the vertices, in the order they are given in `vertices`.
- `null_graph::G`: the Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
- `weighttype::Type{U}`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted);

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
"""
@traitfn function Layer(
    name::Symbol,
    vertices::Union{V,N},#Vector{MultilayerVertex{nothing}},
    degree_sequence::Vector{<:Integer},
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_weight::Function=(src, dst) -> nothing,
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
) where {
    T<:Integer,
    U<:Real,
    G<:AbstractGraph{T},
    V<:Vector{MultilayerVertex{nothing}},
    N<:Vector{Node};!IsDirected{G},
}
    descriptor = LayerDescriptor(
        name,
        null_graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
    )

    equivalent_graph = havel_hakimi_graph_generator(degree_sequence)

    edge_list = NTuple{2,MultilayerVertex{nothing}}[]

    if @isdefined(V)
        for edge in edges(equivalent_graph)
            src_mv = vertices[src(edge)]
            dst_mv = vertices[dst(edge)]
            push!(edge_list, (src_mv, dst_mv))
        end
    else
        for edge in edges(equivalent_graph)
            src_mv = MV(vertices[src(edge)])
            dst_mv = MV(vertices[dst(edge)])
            push!(edge_list, (src_mv, dst_mv))
        end
    end

    layer = Layer(descriptor, vertices, edge_list)

    return layer
end

"""
    Layer(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution,
        null_graph::G,
        weighttype::Type{U};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_weight::Function = (src, dst) -> nothing,
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}; IsDirected{G}}

Returns a directed `Layer` whose indegree and oudegree sequences are sampled from `indegree_distribution` and `outdegree_distribution`. Graph realization is performed using the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `indegree_distribution::UnivariateDistribution`: The distribution to sample indegrees from.
- `outdegree_distribution::UnivariateDistribution`: The distribution to sample outdegrees from.
- `null_graph::G`: the Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
- `weighttype::Type{U}`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted);

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
"""
@traitfn function Layer(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},#Vector{MultilayerVertex{nothing}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution,
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_weight::Function=(src, dst) -> nothing,
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
) where {T<:Integer,U<:Real,G<:AbstractGraph{T};IsDirected{G}}
    indegree_sequence, outdegree_sequence = sample_digraphical_degree_sequences(
        indegree_distribution, outdegree_distribution, length(vertices)
    )

    return Layer(
        name,
        vertices,
        indegree_sequence,
        outdegree_sequence,
        null_graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    Layer(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        indegree_sequence::Vector{<:Integer},
        outdegree_sequence::Vector{<:Integer},
        null_graph::G,
        weighttype::Type{U};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_weight::Function = (src, dst) -> nothing,
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
    ) where {T<:Integer, U <: Real, G<:AbstractGraph{T}; IsDirected{G}}

Returns an directed `Layer` with given `indegree_sequence` and `outdegree_sequence` realized using the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `indegree_sequence::Vector{<:Integer}`: The indegree sequence of the vertices, in the order they are given in `vertices`.
- `outdegree_sequence::Vector{<:Integer}`: The outdegree sequence of the vertices, in the order they are given in `vertices`.
- `null_graph::G`: the Layer's underlying graph type, which must be passed as a null graph. If it is not, an error will be thrown.
- `weighttype::Type{U}`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted);

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`;
"""
@traitfn function Layer(
    name::Symbol,
    vertices::Union{V,N}, #Vector{MultilayerVertex{nothing}},
    indegree_sequence::Vector{<:Integer},
    outdegree_sequence::Vector{<:Integer},
    null_graph::G,
    weighttype::Type{U};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_weight::Function=(src, dst) -> nothing,
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
) where {
    T<:Integer,
    U<:Real,
    G<:AbstractGraph{T},
    V<:Vector{MultilayerVertex{nothing}},
    N<:Vector{Node};IsDirected{G},
}
    descriptor = LayerDescriptor(
        name,
        null_graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_weight=default_edge_weight,
        default_edge_metadata=default_edge_metadata,
    )

    equivalent_graph = kleitman_wang_graph_generator(indegree_sequence, outdegree_sequence)

    edge_list = NTuple{2,MultilayerVertex{nothing}}[]

    if @isdefined(V)
        for edge in edges(equivalent_graph)
            src_mv = vertices[src(edge)]
            dst_mv = vertices[dst(edge)]
            push!(edge_list, (src_mv, dst_mv))
        end
    else
        for edge in edges(equivalent_graph)
            src_mv = MV(vertices[src(edge)])
            dst_mv = MV(vertices[dst(edge)])
            push!(edge_list, (src_mv, dst_mv))
        end
    end

    layer = Layer(descriptor, vertices, edge_list)

    return layer
end

# Quick constructors
#####
"""
    layer_simplegraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimpleGraph` from `Graphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simplegraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(name, vertices, edge_list, SimpleGraph{vertextype}(), weighttype)
end

"""
    layer_simplegraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        degree_distribution::UnivariateDistribution;
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimpleGraph` from `Graphs.jl` with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `degree_distribution::UnivariateDistribution`: The degree distribution from which the degree sequence is sampled ;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simplegraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    degree_distribution::UnivariateDistribution;
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(name, vertices, degree_distribution, SimpleGraph{vertextype}(), weighttype)
end

"""
    layer_simplegraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `SimpleGraph` from `Graphs.jl` with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simplegraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(name, vertices, ne, SimpleGraph{vertextype}(), weighttype)
end

## SimpleDiGraph
"""
    layer_simpledigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimpleDiGraph` from `Graphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpledigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(name, vertices, edge_list, SimpleDiGraph{vertextype}(), weighttype)
end

"""
    layer_simpledigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution;
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimplDiGraph{vertextype}` from `Graphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `indegree_distribution::UnivariateDistribution`: The degree distribution from which the indegree sequence is sampled ;
- `outdegree_distribution::UnivariateDistribution`: The degree distribution from which the outdegree sequence is sampled ;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpledigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution;
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        indegree_distribution,
        outdegree_distribution,
        SimpleDiGraph{vertextype}(),
        weighttype,
    )
end

"""
    layer_simpledigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `SimplDiGraph{vertextype}` from `Graphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpledigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(name, vertices, ne, SimpleDiGraph{vertextype}(), weighttype)
end

## SimpleWeightedGraph
"""
    layer_simpleweightedgraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        default_edge_weight::Function = (src,dst) -> nothing,
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimpleWeightedGraph` from `SimpleWeightedGraphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpleweightedgraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    default_edge_weight::Function=(src, dst) -> nothing,
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        edge_list,
        SimpleWeightedGraph{vertextype,weighttype}(),
        weighttype;
        default_edge_weight=default_edge_weight,
    )
end

"""
    layer_simpleweightedgraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        degree_distribution::UnivariateDistribution;
        default_edge_weight::Function = (src,dst) -> nothing,
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimpleWeightedGraph` from `SimpleWeightedGraphs.jl`. with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `degree_distribution::UnivariateDistribution`: The degree distribution from which the degree sequence is sampled ;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpleweightedgraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    degree_distribution::UnivariateDistribution;
    default_edge_weight::Function=(src, dst) -> nothing,
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        degree_distribution,
        SimpleWeightedGraph{vertextype,weighttype}(),
        weighttype;
        default_edge_weight=default_edge_weight,
    )
end

"""
    layer_simpleweightedgraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        default_edge_weight::Function = (src,dst) -> nothing,
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `SimpleWeightedGraph` from `SimpleWeightedGraphs.jl`. with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpleweightedgraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    default_edge_weight::Function=(src, dst) -> nothing,
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        ne,
        SimpleWeightedGraph{vertextype,weighttype}(),
        weighttype;
        default_edge_weight=default_edge_weight,
    )
end

## SimpleWeightedDiGraph
"""
    layer_simpleweighteddigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        default_edge_weight::Function = (src,dst) -> nothing,
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimpleWeightedDiGraph` from `SimpleWeightedGraphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpleweighteddigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    default_edge_weight::Function=(src, dst) -> nothing,
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        edge_list,
        SimpleWeightedDiGraph{vertextype,weighttype}(),
        weighttype;
        default_edge_weight=default_edge_weight,
    )
end

"""
    layer_simpleweighteddigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution;
        default_edge_weight::Function = (src,dst) -> nothing,
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `SimpleWeightedDiGraph` from `SimpleWeightedGraphs.jl`, with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `indegree_distribution::UnivariateDistribution`: The degree distribution from which the indegree sequence is sampled ;
- `outdegree_distribution::UnivariateDistribution`: The degree distribution from which the outdegree sequence is sampled ;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpleweighteddigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution;
    default_edge_weight::Function=(src, dst) -> nothing,
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        indegree_distribution,
        outdegree_distribution,
        SimpleWeightedDiGraph{vertextype,weighttype}(),
        weighttype;
        default_edge_weight=default_edge_weight,
    )
end

"""
    layer_simpleweighteddigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        default_edge_weight::Function = (src,dst) -> nothing,
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}
Return a random `Layer` with `ne` edges whose underlying graph is a `SimpleWeightedDiGraph` from `SimpleWeightedGraphs.jl`, with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;
- `default_edge_weight::Function`: Function that takes a pair of `MultilayerVertex`s and returns an edge weight of type `weighttype` or `nothing` (which is compatible with unweighted underlying graphs and corresponds to `one(weighttype)` for weighted underlying graphs). Defaults to `(src, dst) -> nothing`;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_simpleweighteddigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    default_edge_weight::Function=(src, dst) -> nothing,
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        ne,
        SimpleWeightedDiGraph{vertextype,weighttype}(),
        weighttype;
        default_edge_weight=default_edge_weight,
    )
end

## MetaGraph
"""
    layer_metagraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `MetaGraph` from `MetaGraphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_metagraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        edge_list,
        MetaGraph{vertextype,weighttype}(),
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_metagraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        degree_distribution::UnivariateDistribution;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `MetaGraph` from `MetaGraphs.jl` with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `degree_distribution::UnivariateDistribution`: The degree distribution from which the degree sequence is sampled ;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_metagraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    degree_distribution::UnivariateDistribution;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        degree_distribution,
        MetaGraph{vertextype,weighttype}(),
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_metagraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `MetaGraph` from `MetaGraphs.jl` with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_metagraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        ne,
        MetaGraph{vertextype,weighttype}(),
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

## MetaDiGraph
"""
    layer_metadigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `MetaDiGraph` from `MetaGraphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_metadigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        edge_list,
        MetaDiGraph{vertextype,weighttype}(),
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_metadigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `MetaDiGraph` from `MetaDiGraphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `indegree_distribution::UnivariateDistribution`: The degree distribution from which the indegree sequence is sampled ;
- `outdegree_distribution::UnivariateDistribution`: The degree distribution from which the outdegree sequence is sampled ;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_metadigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution,
        MetaDiGraph{vertextype,weighttype}(),
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_metadigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `MetaDiGraph` from `MetaDiGraphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_metadigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    return Layer(
        name,
        vertices,
        ne,
        MetaDiGraph{vertextype,weighttype}(),
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

## ValGraph
"""
    layer_valgraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `ValGraph` from `SimpleValueGraphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_valgraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_returned_type = Base.return_types(default_vertex_metadata)[1]
    edgeval_returned_type = Base.return_types(default_edge_metadata)[1]

    vertexval_types = if vertexval_returned_type <: Tuple
        tuple(vertexval_returned_type.parameters...)
    elseif vertexval_returned_type <: NamedTuple
        returned_namedtuple_type = vertexval_returned_type.parameters

        type_pars = tuple(returned_namedtuple_type[end].parameters...)

        NamedTuple{returned_namedtuple_type[1]}(type_pars)
    end

    edgeval_types = if edgeval_returned_type <: Tuple
        tuple(edgeval_returned_type.parameters...)
    elseif edgeval_returned_type <: NamedTuple
        returned_namedtuple_type = edgeval_returned_type.parameters

        type_pars = tuple(returned_namedtuple_type[end].parameters...)

        NamedTuple{returned_namedtuple_type[1]}(type_pars)
    end

    graph = ValGraph{vertextype}(
        SimpleGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        edge_list,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_valgraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        degree_distribution::UnivariateDistribution;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `ValGraph` from `SimpleValueGraphs.jl`. with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `degree_distribution::UnivariateDistribution`: The degree distribution from which the degree sequence is sampled ;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_valgraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    degree_distribution::UnivariateDistribution;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValGraph{vertextype}(
        SimpleGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        degree_distribution,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_valgraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `ValGraph` from `SimpleValueGraphs.jl`. with a degree sequence sampled from `degree_distribution`. Realization is performed via the Havel-Hakimi algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_valgraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValGraph{vertextype}(
        SimpleGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        ne,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

## ValOutDiGraph
"""
    layer_valoutdigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `ValOutDiGraph` from `SimpleValueGraphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_valoutdigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValOutDiGraph{vertextype}(
        SimpleDiGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        edge_list,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_valoutdigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `ValOutDiGraph` from `SimpleValueGraphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `indegree_distribution::UnivariateDistribution`: The degree distribution from which the indegree sequence is sampled ;
- `outdegree_distribution::UnivariateDistribution`: The degree distribution from which the outdegree sequence is sampled ;

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`. Do not type this function's arguments;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`,  Do not type this function's arguments;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_valoutdigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValOutDiGraph{vertextype}(
        SimpleDiGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        indegree_distribution,
        outdegree_distribution,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_valoutdigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `ValOutDiGraph` from `SimpleValueGraphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`. Do not type this function's arguments;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`,  Do not type this function's arguments;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_valoutdigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValOutDiGraph{vertextype}(
        SimpleDiGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        ne,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

## ValDiGraph
"""
    layer_valdigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}};
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `ValDiGraph` from `SimpleValueGraphs.jl`.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `edge_list::Union{Vector{<:MultilayerEdge}, Vector{NTuple{2, MultilayerVertex{nothing}}}}`: The list of `MultilayerEdge`s. It may be a vector of `MultilayerEdge`s or a Vector of 2-tuples of `MultilayerVertex`s. In the latter case, the weight and the metadata of the `MultilayerEdge` to be added are computed respectively via the `default_edge_weight` and `default_edge_metadata` functions;


# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)
"""
function layer_valdigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    edge_list::Union{Vector{<:MultilayerEdge},Vector{NTuple{2,MultilayerVertex{nothing}}}};
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValDiGraph{vertextype}(
        SimpleDiGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        edge_list,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_valdigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Constructor for a `Layer` whose underlying graph is a `ValDiGraph` from `SimpleValueGraphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `indegree_distribution::UnivariateDistribution`: The degree distribution from which the indegree sequence is sampled ;
- `outdegree_distribution::UnivariateDistribution`: The degree distribution from which the outdegree sequence is sampled ;

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`. Do not type this function's arguments;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`,  Do not type this function's arguments;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_valdigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValDiGraph{vertextype}(
        SimpleDiGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        indegree_distribution,
        outdegree_distribution,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    layer_valdigraph(
        name::Symbol,
        vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
        ne::Int64;
        default_vertex_metadata::Function = mv -> NamedTuple(),
        default_edge_metadata::Function = (src, dst) -> NamedTuple(),
        vertextype::Type{T} = Int64,
        weighttype::Type{U} = Float64
    ) where {T<:Integer,U<:Real}

Return a random `Layer` with `ne` edges whose underlying graph is a `ValDiGraph` from `SimpleValueGraphs.jl` with a indegree and outdegree sequences respectively sampled from `indegree_distribution` and `outdegree_distribution`. Realization is performed via the Kleitman-Wang algorithm.

# ARGUMENTS

- `name::Symbol`: The name of the Layer;
- `vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}}`: The `MultilayerVertex`s of the Layer. May be a vector of `MultilayerVertex{nothing}`s or a vector of `Node`s. In the latter case, the metadata of the `MultilayerVertex` to be added are computed via the `default_vertex_metadata` before the vertex is added (the function will act on each element of `MV.(vertices)`);
- `ne::Int64`: The number of edges of the Layer;

# KWARGS
-` default_vertex_metadata::Function`: Function that takes a `MultilayerVertex` and returns a `Tuple` or a `NamedTuple` containing the vertex metadata. defaults to `mv -> NamedTuple()`. Do not type this function's arguments;
- `default_edge_metadata::Function`: Function that takes a pair of `MultilayerVertex`s and  returns a `Tuple` or a `NamedTuple` containing the edge metadata, that will be called when `add_edge!(mg,src,dst, args...; kwargs...)` is called without the `metadata` keyword argument, and when generating the edges in this constructor. Defaults to  `(src, dst) -> NamedTuple()`,  Do not type this function's arguments;
- `vertextype::Type{T} = Int64`: The type of the underlying integer labels associated to vertices.
- `weighttype::Type{U} = Float64`: The type of the `MultilayerEdge` weights (even when the underlying Layer's graph is unweighted, we need to specify a weight type since the `MultilayerGraph`s will always be weighted)

"""
function layer_valdigraph(
    name::Symbol,
    vertices::Union{Vector{MultilayerVertex{nothing}},Vector{Node}},
    ne::Int64;
    default_vertex_metadata::Function=mv -> NamedTuple(),
    default_edge_metadata::Function=(src, dst) -> NamedTuple(),
    vertextype::Type{T}=Int64,
    weighttype::Type{U}=Float64,
) where {T<:Integer,U<:Real}
    vertexval_types = get_valtypes(default_vertex_metadata)
    edgeval_types = get_valtypes(default_edge_metadata)

    graph = ValDiGraph{vertextype}(
        SimpleDiGraph{vertextype}();
        vertexval_types=vertexval_types,
        edgeval_types=edgeval_types,
        vertexval_init=default_vertex_metadata,
        edgeval_init=default_edge_metadata,
    )

    return Layer(
        name,
        vertices,
        ne,
        graph,
        weighttype;
        default_vertex_metadata=default_vertex_metadata,
        default_edge_metadata=default_edge_metadata,
    )
end

"""
    has_node(layer::Layer, n::Node)

Return `true` if `n` is a node of `layer`.
"""
has_node(layer::Layer, n::Node) = MV(n, layer.name) ∈ image(layer.v_V_associations)

"""
    has_vertex(layer::Layer, mv::MultilayerVertex)

Return `true` if `v` is a vertex of `layer`.
"""
function Graphs.has_vertex(layer::Layer, mv::MultilayerVertex)
    return MV(node(mv), name(layer)) ∈ collect(image(layer.v_V_associations))
end

# TODO:
# Implement a MultilayerVertex constructor that leaves the .layer field unspecified, for ease of use of the following function
"""
    add_vertex!(layer::Layer, mv::MultilayerVertex) 

Add vertex to layer `layer`. 
"""
function Graphs.add_vertex!(layer::Layer, mv::MultilayerVertex)
    (isnothing(mv.layer) || mv.layer == layer.name) || throw(
        ErrorException("The multilayer vertex $mv cannot belong to layer $(layer.name)."),
    )
    return add_vertex!(layer, mv.node; metadata=mv.metadata)
end

"""
    add_vertex!(layer::L, n::Node, args...; kwargs...) where {T, L <: Layer{T}}      

Add vertex associated with node `n` to layer `layer`. This method supports the uniform and transparent interfaces. See the [Vertices](@ref) section of the Tutorial.
"""
function Graphs.add_vertex!(layer::L, n::Node, args...; kwargs...) where {T,L<:Layer{T}}
    # Check if the vertex is already in the layer
    has_node(layer, n) && return false

    success = false
    if isempty(args) &&
        length(kwargs) == 1 &&
        issetequal(Set([:metadata]), Set(keys(kwargs)))
        # If only the metadata is provided, call the standard add_vertex method
        success = add_vertex_standard!(layer; metadata=values(kwargs).metadata)
    elseif length(args) == length(kwargs) == 0
        # If no arguments or keyword arguments are provided, use the layer's default vertex metadata
        success = add_vertex_standard!(
            layer; metadata=layer.default_vertex_metadata(MV(n, layer.name))
        )
    else
        # Otherwise, call the generic add_vertex method
        success = add_vertex!(layer.graph, args...; kwargs...)
    end

    # Check if the vertex is already in the layer
    if success
        # Get the last vertex in the layer
        last_vertex = if length(layer.v_V_associations) == 0
            zero(T)
        else
            maximum(domain(layer.v_V_associations))
        end
        # Add the vertex to the layer
        layer.v_V_associations[last_vertex + one(T)] = MV(n, layer.name)
        return true
    else
        return false
    end
end

"""
    add_vertex_standard!(layer::Layer; metadata::Union{Tuple, NamedTuple}= NamedTuple())

Add vertex with metadata to layer `layer`.
"""
function add_vertex_standard!(layer::Layer; metadata::Union{Tuple,NamedTuple}=NamedTuple())
    return _add_vertex!(layer; metadata=metadata)
end

"""
    _add_vertex!( layer::L; metadata::Union{Tuple, NamedTuple}= NamedTuple()) where {T, U, G, L <: Layer{T,U,G}}

Add vertex with metadata to layer `layer`.
"""
function _add_vertex!(
    layer::L; metadata::Union{Tuple,NamedTuple}=NamedTuple()
) where {T,U,G,L<:Layer{T,U,G}}
    return __add_vertex!(layer.graph; metadata=metadata)
end

"""
    rem_vertex!(layer::Layer, mv::MultilayerVertex) 

Remove vertex `mv` from layer `layer`.
"""
function Graphs.rem_vertex!(layer::Layer, mv::MultilayerVertex)
    (isnothing(mv.layer) || mv.layer == layer.name) || return false
    return rem_vertex!(layer, mv.node)
end

"""
    rem_vertex!(layer::Layer, n::Node)

Remove node `n` from `layer`. Modify `layer.v_N_associations` according to how `rem_vertex!` works in [Graph.jl](https://juliagraphs.org/Graphs.jl/dev/core_functions/simplegraphs/#Graphs.SimpleGraphs.rem_vertex!-Tuple{Graphs.SimpleGraphs.AbstractSimpleGraph,%20Integer}).
"""
function Graphs.rem_vertex!(layer::Layer, n::Node)
    !has_node(layer, n) && return false
    success = rem_vertex!(layer, layer.v_V_associations(MV(n, layer.name)))

    if success
        # Get the key of the node to be removed
        v = layer.v_V_associations(MV(n, layer.name))
        # Get the last node and its key
        last_v = maximum(domain(layer.v_V_associations))

        if v != last_v
            last_V = layer.v_V_associations[last_v]
            # Delete the node to be removed
            delete!(layer.v_V_associations, v)
            # Delete the last node
            delete!(layer.v_V_associations, last_v)
            # Re-add the last node with the key of the node to be removed
            layer.v_V_associations[v] = last_V
            true
        else
            # Delete the node to be removed
            delete!(layer.v_V_associations, v)
            true
        end
    else
        return false
    end
end

"""
    rem_vertex!(layer::L, v::T) where {T, L <: Layer{T}}    

Remove vertex `v` from layer `layer`.
"""
Graphs.rem_vertex!(layer::L, v::T) where {T,L<:Layer{T}} = rem_vertex!(layer.graph, v)

"""
    get_v(layer::Layer, V::MultilayerVertex{nothing})

Return `v` associated with `V`, when `V` does not specify its layer. Only for `Layer`s. 
"""
function get_v(layer::Layer, V::MultilayerVertex{nothing})
    # Convert V to a bare vertex
    bare_V = get_bare_mv(V)
    # Check if subgraph has this vertex
    has_vertex(layer, bare_V) || return nothing
    # Get the list of edges
    return layer.v_V_associations(MV(node(V), name(layer)))
end

"""
    has_edge( layer::Layer, s::MultilayerVertex{nothing}, d::MultilayerVertex{nothing})

Return `true` if there is an edge between `s` and `d`, `false` otherwise.
"""
function Graphs.has_edge(
    layer::Layer, s::MultilayerVertex{nothing}, d::MultilayerVertex{nothing}
)
    return has_edge(
        layer,
        get_v(layer, MV(node(s), name(layer))),
        get_v(layer, MV(node(d), name(layer))),
    )
end

"""
    add_edge!(layer::L, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...) where {L <: Layer} 

Add edge from vertex `src` to vertex `dst` to layer `layer`. Returns true if succeeds. This method supports the uniform and transparent interfaces. See the [Edges](@ref edges_tut_subg) section of the Tutorial.
"""
function Graphs.add_edge!(
    layer::Layer, src::MultilayerVertex, dst::MultilayerVertex, args...; kwargs...
)
    # Check if the vertices exist
    !has_vertex(layer, src) &&
        throw(ErrorException("Vertex $(src) does not belong to the layer."))
    !has_vertex(layer, dst) &&
        throw(ErrorException("Vertex $(dst) does not belong to the layer."))

    # Check if the edge already exists
    if !has_edge(layer, src, dst)
        # If the edge doesn't exist, add it
        # If the user does not specify any arguments, we use the default values
        if isempty(args) &&
            length(kwargs) == 2 &&
            issetequal(Set([:weight, :metadata]), Set(keys(kwargs)))
            success = add_edge_standard!(
                layer,
                src,
                dst;
                weight=values(kwargs).weight,
                metadata=values(kwargs).metadata,
            )
            # If the user only specifies the weight, we use the default metadata
        elseif isempty(args) &&
            length(kwargs) == 1 &&
            issetequal(Set([:weight]), Set(keys(kwargs)))
            success = add_edge_standard!(
                layer,
                src,
                dst;
                weight=values(kwargs).weight,
                metadata=layer.default_edge_metadata(src, dst),
            )
            # If the user only specifies the metadata, we use the default weight
        elseif isempty(args) &&
            length(kwargs) == 1 &&
            issetequal(Set([:metadata]), Set(keys(kwargs)))
            success = add_edge_standard!(
                layer,
                src,
                dst;
                weight=layer.default_edge_weight(src, dst),
                metadata=values(kwargs).metadata,
            )
            # If the user does not specify any arguments, we use the default values
        elseif length(args) == length(kwargs) == 0
            success = add_edge_standard!(
                layer,
                src,
                dst;
                weight=layer.default_edge_weight(src, dst),
                metadata=layer.default_edge_metadata(src, dst),
            )
        else
            # If the user specifies arguments, we use those instead of the defaults
            success = add_edge!(
                layer.graph, get_v(layer, src), get_v(layer, dst), args...; kwargs...
            )
        end
        return success
    else
        return false
    end
end

# Base overloads
"""
    Base.(==)(x::Layer, y::Layer)

Overload equality for `Layer`s.
"""
function Base.:(==)(x::Layer, y::Layer)
    for field in fieldnames(Layer)
        if @eval $x.$field != $y.$field
            return false
        end
    end
    return true
end

"""
    Base.getproperty(layer::L, f::Symbol) where {L<:Layer}

Overload `getproperty` for `Layer`s.
"""
function Base.getproperty(layer::L, f::Symbol) where {L<:Layer}
    if f ∈ (:descriptor, :graph, :v_V_associations)
        Base.getfield(layer, f)
    elseif f ∈ (
        :name,
        :null_graph,
        :default_vertex_metadata,
        :default_edge_weight,
        :default_edge_metadata,
    )
        Base.getfield(layer.descriptor, f)
    end
end

# Console print utilities
function to_string(x::Layer)
    parameters = typeof(x).parameters
    return """
           Layer\t$(name(x))
           underlying_graph: $(typeof(graph(x)))
           vertex type: $(parameters[1])
           weight type: $(parameters[2]) 
           nv = $(nv(x))
           ne = $(ne(x))
           """
end
Base.show(io::IO, x::Layer) = print(io, to_string(x))
