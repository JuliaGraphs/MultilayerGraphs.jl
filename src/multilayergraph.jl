"""
    MultilayerGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultilayerGraph{T,U}

A concrete type that can represent a general multilayer graph. Its internal fields aren't meant to be modified by the user. Please prefer the provided API.

# FIELDS

- `layers::Vector{LayerDescriptor{T,U}}`: vector containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
- `interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}} `: the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
- `v_V_associations::Bijection{ T, <: MultilayerVertex}`: A Bijection from Bijections.jl that associates numeric vertices to `MultilayerVertex`s.
- `idx_N_associations::Bijection{ Int64 , Node}`: A Bijection from Bijections.jl that associates Int64 to `Node`s.
- `fadjlist::Vector{Vector{HalfEdge{ <: MultilayerVertex, <: Union{Nothing, U}}}}`: the forward adjacency list of the MultilayerGraph. It is a vector of vectors of [`HalfEdge`](@ref)s. Its i-th element are the[`HalfEdge`](@ref)s that originate from `v_V_associations[i]`.
- `v_metadata_dict::Dict{T, <: Union{ <: Tuple, <: NamedTuple}}`: A Dictionary that associates numeric vertices to their metadata.
"""
mutable struct MultilayerGraph{T,U} <: AbstractMultilayerUGraph{T,U}
    layers::Vector{LayerDescriptor{T,U}}
    interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}} 
    v_V_associations::Bijection{ T, <: MultilayerVertex}
    idx_N_associations::Bijection{ Int64 , Node}
    fadjlist::Vector{Vector{HalfEdge{ <: MultilayerVertex, <: Union{Nothing, U}}}}
    v_metadata_dict::Dict{T, <: Union{ <: Tuple, <: NamedTuple}}
end

# Traits
@traitimpl IsWeighted{MultilayerGraph}

# Constructors
"""
    MultilayerGraph(layers::Vector{ <: Layer{T,U }}, specified_interlayers::Vector{ <: Interlayer{T,U}};  default_interlayer::String  = "multiplex") where {T, U}      

Construct a MultilayerGraph with layers given by `layers`. The interlayers will be constructed by default according to `default_interlayer` where only `"multiplex"` and `"empty"` are allowed, except for those specified in `specified_interlayers`. `default_interlayer = "multiplex"` will imply that unspecified interlayers will have only diagonal couplings, while  `default_interlayer = "multiplex"` will produced interlayers that have no couplings.

# ARGUMENTS

- `layers::Vector{<:Layer{T,U}}`: The (ordered) list of layers the multilayer graph will have;
- `specified_interlayers::Vector{<:Interlayer{T,U}}`: The list of interlayers specified by the user. Note that the user does not need to specify all interlayers, as the unspecified ones will be automatically constructed using the indications given by the `default_interlayers_null_graph` and `default_interlayers_structure` keywords.;
- `default_interlayers_null_graph::H = SimpleGraph{T}()`: Sets the underlying graph for the interlayers that are to be automatically specified. Defaults to `SimpleGraph{T}()`. See the `Layer` constructors for more information.;
- `default_interlayers_structure::String = "multiplex"`: Sets the structure of the interlayers that are to be automatically specified. May be "multiplex" for diagonally coupled interlayers, or "empty" for empty interlayers (no edges).  "multiplex". See the `Interlayer` constructors for more information.;
"""
function MultilayerGraph(
    layers::Vector{<:Layer{T,U}},
    specified_interlayers::Vector{<:Interlayer{T,U}};
    default_interlayers_null_graph::H = SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U, H <: AbstractGraph{T}}

    multilayergraph = MultilayerGraph(T, U)

    for layer in deepcopy(layers)
        add_layer!(multilayergraph, layer; default_interlayers_null_graph = default_interlayers_null_graph, default_interlayers_structure = default_interlayers_structure)
    end

    if !isnothing(specified_interlayers)
        for interlayer in deepcopy(specified_interlayers)
            specify_interlayer!(multilayergraph, interlayer) 
        end
    end

    return multilayergraph
end

# General MultilayerGraph Utilities
fadjlist(mg::MultilayerGraph) = mg.fadjlist

# Nodes

# Vertices

# Edges
    
"""
    add_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerUGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add MultilayerEdge `me` to the MultilayerGraph `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.add_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerUGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}
    _src = get_bare_mv(src(me))
    _dst = get_bare_mv(dst(me))
    has_vertex(mg, _src) || throw(ErrorException("Vertex $_src does not belong to the multilayer graph."))
    has_vertex(mg, _dst) || throw(ErrorException("Vertex $_dst does not belong to the multilayer graph."))

    # Add edge to `edge_dict`
    src_V_idx = get_v(mg, _src)
    dst_V_idx = get_v(mg, _dst)

    _weight = isnothing(weight(me)) ? one(U) : weight(me)
    _metadata = metadata(me)

    if !has_edge(mg, _src, _dst)
        if src_V_idx != dst_V_idx
            push!(mg.fadjlist[src_V_idx], HalfEdge(_dst, _weight, _metadata))
            push!(mg.fadjlist[dst_V_idx], HalfEdge(_src, _weight, _metadata))
        else
            push!(mg.fadjlist[src_V_idx], HalfEdge(_dst, _weight, _metadata))
        end
        return true
    else
#=         # Should we modify weight and metadata or should we return false? This may be something to decide ecosystem-wise
        set_weight!(mg, src, dst, _weight)
        set_metadata!(mg, src, dst, _metadata) =#
        return false
    end
end

"""
    rem_edge!(mg::M, src::T, dst::T) where {T,U, M <: MultilayerGraph{T,U}}

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.rem_edge!(mg::M, src::MultilayerVertex, dst::MultilayerVertex) where  {T,U,M <: AbstractMultilayerUGraph{T,U}}
    # Perform routine checks
    has_vertex(mg, src) || throw(ErrorException("Vertex $_src does not belong to the multilayer graph."))
    has_vertex(mg, dst) || throw(ErrorException("Vertex $_dst does not belong to the multilayer graph."))

    has_edge(mg, src, dst) || return false

    src_V_idx = get_v(mg, src)
    dst_V_idx = get_v(mg, dst)

    _src = get_bare_mv(src)
    _dst = get_bare_mv(dst)

    if get_bare_mv(src) != get_bare_mv(dst)
        src_idx_tbr = findfirst(halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx])
        deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)

        dst_idx_tbr = findfirst(halfedge -> halfedge.vertex == _src, mg.fadjlist[dst_V_idx])
        deleteat!(mg.fadjlist[dst_V_idx], dst_idx_tbr)
    else
        src_idx_tbr = findfirst(halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx])
        deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)
    end

    return true
end

# Layers and Interlayers

# Graphs.jl's extensions

# Multilayer-specific methods
# "empty graph" could be the correct way of calling a graph with no edges: https://math.stackexchange.com/questions/320859/what-is-the-term-for-a-graph-on-n-vertices-with-no-edges
"""
    MultilayerGraph(empty_layers::Vector{<:Layer{T,U}}, empty_interlayers::Vector{<:Interlayer{T,U}}, degree_distribution::UnivariateDistribution; allow_self_loops::Bool = false, default_interlayers_null_graph::H = SimpleGraph{T}() ) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

Return a random MultilayerGraph that has `empty_layers` as layers and `empty_interlayers` as specified interlayers. `empty_layers` and `empty_interlayers` must respectively be `Layer`s and `Interlayer`s with whatever number of vertices but no edges (if any edge is found, an error is thrown). The  degree distribution of the returned random `MultilayerGraph` is given by `degree_distribution`, which must have a support that only contains positive numbers for obvious reasons. `allow_self_loops = true` allows for self loops t be present in the final random MultilayerGraph. `default_interlayers_null_graph` controls the `null_graph` argument passed to automatically-generated interlayers. 
"""
function MultilayerGraph(
    empty_layers::Vector{<:Layer{T,U}},
    empty_interlayers::Vector{<:Interlayer{T,U}},
    degree_distribution::UnivariateDistribution;
    allow_self_loops::Bool = false,
    default_interlayers_null_graph::H = SimpleGraph{T}(),
) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

    !allow_self_loops || throw(ErrorException("`allow_self_loops` must currently be set to `false`. The configuration model algorithm does not support self-loops yet."))

    minimum(support(degree_distribution)) >= 0 ||  throw(ErrorException("The `degree_distribution` must have positive support. Found $(support(degree_distribution))."))
    
    empty_multilayergraph = MultilayerGraph(empty_layers, empty_interlayers; default_interlayers_null_graph = default_interlayers_null_graph, default_interlayers_structure = "empty")

    n = nv(empty_multilayergraph)

    degree_sequence = nothing
    acceptable = false

    while !acceptable
        degree_sequence = round.(Ref(Int), rand(degree_distribution, nv(empty_multilayergraph)))

        acceptable = all(0 .<= degree_sequence .< n) 

        if acceptable
            acceptable = isgraphical(degree_sequence) 
        end

    end
    return MultilayerGraph(empty_multilayergraph, degree_sequence; allow_self_loops = false, perform_checks = false)
end

"""
    MultilayerGraph(empty_multilayergraph::MultilayerGraph{T,U}, degree_sequence::Vector{<:Integer}; allow_self_loops::Bool = false, perform_checks::Bool = false) where {T,U}

Return a random `MultilayerGraph` with degree sequence `degree_sequence`. `allow_self_loops` controls the presence of self-loops, while if `perform_checks` is true, the `degree_sequence` os checked to be graphical.
"""
function MultilayerGraph(empty_multilayergraph::MultilayerGraph{T,U}, degree_sequence::Vector{<:Integer}; allow_self_loops::Bool = false, perform_checks::Bool = true) where {T,U}

    (allow_self_loops && perform_checks) && @warn "Checks for graphicality and coherence with the provided `empty_multilayergraph` are currently performed without taking into account self-loops. Thus said checks may fail event though the provided `degree_sequence` may be graphical when one allows for self-loops within the multilayer graph to be present. If you are sure that the provided `degree_sequence` is indeed graphical under those circumstances, you may want to disable checks by setting `perform_checks = false`. We apologize for the inconvenient."
    
    _multilayergraph = deepcopy(empty_multilayergraph)

    ne(_multilayergraph) == 0 || throw(ErrorException("The `empty_multilayergraph` argument should be an empty MultilayerGraph. Found $(ne(_multilayergraph)) edges."))

    if perform_checks
        n = nv(_multilayergraph)
        n == length(degree_sequence) || throw(ErrorException("The number of vertices of the provided empty MultilayerGraph does not match the length of the degree sequence. Found $(nv(_multilayergraph)) and $(length(degree_sequence))."))

        isgraphical(degree_sequence) || throw(ArgumentError("degree_sequence must be graphical.")) 
    end

    edge_list = _random_undirected_configuration(_multilayergraph, degree_sequence, allow_self_loops)

    for edge in edge_list
        add_edge!(_multilayergraph, edge)
    end

    return _multilayergraph
end

"""
    MultilayerGraph(n_nodes::Int64, T::Type{ <: Number}, U::Type{ <: Number} )

Return a null MultilayerGraph with with vertex type `T` weighttype `U`. Use this constructor and then add Layers and Interlayers via the `add_layer!` and `specify_interlayer!` methods.
"""
MultilayerGraph(T::Type{<:Number}, U::Type{<:Number}) =  MultilayerGraph{T,U}(
                                                                                LayerDescriptor{T,U}[],
                                                                                OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}}(),
                                                                                Bijection{ T , MultilayerVertex}(),
                                                                                Bijection{ Int64 , Node}(),
                                                                                Vector{HalfEdge{MultilayerVertex, <: Union{Nothing, U}}}[],
                                                                                Dict{T, Union{Tuple,NamedTuple}}()
)

"""
    MultilayerGraph(layers::Vector{<:Layer{T,U}}; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String="multiplex") where {T,U, H <: AbstractGraph{T}}

Construct a MultilayerGraph with layers `layers` and all interlayers with structure `default_interlayers_structure` (only "multiplex" and "empty" are allowed) and type `default_interlayers_null_graph`.
"""
MultilayerGraph(layers::Vector{<:Layer{T,U}}; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String="multiplex") where {T,U, H <: AbstractGraph{T}} = MultilayerGraph(layers, Interlayer{T,U}[]; default_interlayers_null_graph = default_interlayers_null_graph, default_interlayers_structure = default_interlayers_structure)

# Base overloads
"""
    Base.getproperty(mg::M, f::Symbol) where { M <: MultilayerGraph }
"""
function Base.getproperty(mg::MultilayerGraph, f::Symbol)#  where {T,U,M<:AbstractMultilayerUGraph{T,U}}
    if f in (:v_V_associations, :fadjlist, :idx_N_associations, :layers, :interlayers, :v_metadata_dict) # :weight_tensor, :supra_weight_matrix, 
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