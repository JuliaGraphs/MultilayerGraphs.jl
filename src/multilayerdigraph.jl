"""
    MultilayerDiGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultilayerGraph{T,U}

A concrete type that can represent a general multilayer graph. Its internal fields aren't meant to be modified by the user. Please prefer the provided API.
"""
mutable struct MultilayerDiGraph{T,U} <: AbstractMultilayerDiGraph{T,U}
    layers::Vector{LayerDescriptor{T,U}} # vector containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
    interlayers::OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}} #  the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
    v_V_associations::Bijection{ T, <: MultilayerVertex} # A Bijection from Bijections.jl that associates numeric vertices to `MultilayerVertex`s.
    idx_N_associations::Bijection{ Int64 , Node} # A Bijection from Bijections.jl that associates Int64 to `Node`s.
    fadjlist::Vector{Vector{HalfEdge{ <: MultilayerVertex, <: Union{Nothing, U}}}} # the forward adjacency list of the MultilayerDiGraph. It is a vector of vectors of `HalfEdge`s. Its i-th element are the `HalfEdge`s that originate from `v_V_associations[i]`.
    badjlist::Vector{Vector{HalfEdge{ <: MultilayerVertex, <: Union{Nothing, U}}}} # the bacward adjacency list of the MultilayerDiGraph. It is a vector of vectors of `HalfEdge`s. Its i-th element are the `HalfEdge`s that insost on `v_V_associations[i]`.
    v_metadata_dict::Dict{T, <: Union{ <: Tuple, <: NamedTuple}} # A Dictionary that associates numeric vertices to their metadata
end

# Traits
@traitimpl IsWeighted{MultilayerDiGraph}
@traitimpl IsDirected{MultilayerDiGraph}

# Constructors
"""
    MultilayerDiGraph(
        layers::Vector{<:Layer{T,U}},
        specified_interlayers::Vector{<:Interlayer{T,U}};
        default_interlayers_null_graph::H = SimpleGraph{T}(),
        default_interlayers_structure::String="multiplex",
    ) where {T,U, H <: AbstractGraph{T}}    

Construct a MultilayerDiGraph with layers given by `layers`. The interlayers will be constructed by default according to `default_interlayer` where only `"multiplex"` and `"empty"` are allowed, except for those specified in `specified_interlayers`. `default_interlayer = "multiplex"` will imply that unspecified interlayers will have only diagonal couplings, while  `default_interlayer = "multiplex"` will produced interlayers that have no couplings.
"""
function MultilayerDiGraph(
    layers::Vector{<:Layer{T,U}},
    specified_interlayers::Vector{<:Interlayer{T,U}};
    default_interlayers_null_graph::H = SimpleGraph{T}(),
    default_interlayers_structure::String="multiplex",
) where {T,U, H <: AbstractGraph{T}}

    multilayerdigraph = MultilayerDiGraph(T, U)

    for layer in deepcopy(layers)
        add_layer!(multilayerdigraph, layer; default_interlayers_null_graph = default_interlayers_null_graph, default_interlayers_structure = default_interlayers_structure)
    end

    if !isnothing(specified_interlayers)
        for interlayer in deepcopy(specified_interlayers)
            specify_interlayer!(multilayerdigraph, interlayer) 
        end
    end

    return multilayerdigraph
end

# General MultilayerDiGraph Utilities
fadjlist(mg::MultilayerDiGraph) = mg.fadjlist
badjlist(mg::MultilayerDiGraph) = mg.badjlist

# Nodes

# Vertices

# Edges    
"""
    add_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerDiGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}

Add MultilayerEdge `me` to the MultilayerDiGraph `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.add_edge!(mg::M, me::E) where {T,U, M <: AbstractMultilayerDiGraph{T,U}, E <: MultilayerEdge{ <: Union{U,Nothing}}}
    
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
        push!(mg.fadjlist[src_V_idx], HalfEdge(_dst, _weight, _metadata))
        push!(mg.badjlist[dst_V_idx], HalfEdge(_src, _weight, _metadata))
    else

        return false
    end
#=     else
        push!(mg.fadjlist[src_V_idx], HalfEdge(_dst, _weight, _metadata))
    end =#

    return true
end

"""
    rem_edge!(mg::MultilayerDiGraph, src::MultilayerVertex, dst::MultilayerVertex)

Remove edge from `src` to `dst` from `mg`. Return true if succeeds, false otherwise.
"""
function Graphs.rem_edge!(mg::MultilayerDiGraph, src::MultilayerVertex, dst::MultilayerVertex)

    # Perform routine checks
    has_vertex(mg, src) || throw(ErrorException("Vertex $_src does not belong to the multilayer graph."))
    has_vertex(mg, dst) || throw(ErrorException("Vertex $_dst does not belong to the multilayer graph."))

    has_edge(mg, src, dst) || return false

    src_V_idx = get_v(mg, src)
    dst_V_idx = get_v(mg, dst)

    _src = get_bare_mv(src)
    _dst = get_bare_mv(dst)

    src_idx_tbr = findfirst(halfedge -> vertex(halfedge) == _dst, mg.fadjlist[src_V_idx])
    deleteat!(mg.fadjlist[src_V_idx], src_idx_tbr)

    dst_idx_tbr = findfirst(halfedge -> halfedge.vertex == _src, mg.badjlist[dst_V_idx])
    deleteat!(mg.badjlist[dst_V_idx], dst_idx_tbr)

    return true
end

# Layers and Interlayers

# Graphs.jl's extensions

# Multilayer-specific methods
# "empty graph" could be the correct way of calling a graph with no edges: https://math.stackexchange.com/questions/320859/what-is-the-term-for-a-graph-on-n-vertices-with-no-edges
"""
    MultilayerDiGraph(
        empty_layers::Vector{<:Layer{T,U}},
        empty_interlayers::Vector{<:Interlayer{T,U}},
        indegree_distribution::UnivariateDistribution,
        outdegree_distribution::UnivariateDistribution;
        allow_self_loops::Bool = false,
        default_interlayers_null_graph::H = SimpleGraph{T}(),
    ) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

Return a random MultilayerDiGraph that has `empty_layers` as layers and `empty_interlayers` as specified interlayers. `empty_layers` and `empty_interlayers` must respectively be `Layer`s and `Interlayer`s with whatever number of vertices but no edges (if any edge is found, an error is thrown). The  degree distribution of the returned random `MultilayerDiGraph` is given by `degree_distribution`, which must have a support that only contains positive numbers for obvious reasons. `allow_self_loops = true` allows for self loops t be present in the final random MultilayerDiGraph. `default_interlayers_null_graph` controls the `null_graph` argument passed to automatically-generated interlayers. 
"""
function MultilayerDiGraph(
    empty_layers::Vector{<:Layer{T,U}},
    empty_interlayers::Vector{<:Interlayer{T,U}},
    indegree_distribution::UnivariateDistribution,
    outdegree_distribution::UnivariateDistribution;
    allow_self_loops::Bool = false,
    default_interlayers_null_graph::H = SimpleGraph{T}(),
) where {T <: Integer, U <: Real, H <: AbstractGraph{T}}

    !allow_self_loops || throw(ErrorException("`allow_self_loops` must currently be set to `false`. The configuration model algorithm does not support self-loops yet."))

    minimum(support(indegree_distribution)) >= 0 ||  throw(ErrorException("Both the `indegree_distribution` and the `outdegree_distribution` must have positive support. Found $(support(indegree_distribution)) and $(support(outdegree_distribution))."))
    
    empty_multilayerdigraph = MultilayerDiGraph(empty_layers, empty_interlayers; default_interlayers_null_graph = default_interlayers_null_graph, default_interlayers_structure = "empty")

    n = nv(empty_multilayerdigraph)

    indegree_sequence  = Vector{Int64}(undef, n)
    outdegree_sequence = Vector{Int64}(undef, n)
    acceptable = false

    @info "Trying to sample a digraphical sequence from the two provided distributions..."
    while !acceptable
        indegree_sequence  .= round.(Ref(Int), rand(indegree_distribution, nv(empty_multilayerdigraph)))
        outdegree_sequence .= round.(Ref(Int), rand(outdegree_distribution, nv(empty_multilayerdigraph)))

        acceptable = all(0 .<= vcat(indegree_sequence,outdegree_sequence)  .< n) 

        if acceptable
            acceptable = isdigraphical(indegree_sequence, outdegree_sequence)  
        end

    end

    return MultilayerDiGraph(empty_multilayerdigraph, indegree_sequence, outdegree_sequence; allow_self_loops = false,perform_checks = false )
end

"""
    MultilayerDiGraph(
        empty_multilayerdigraph::MultilayerDiGraph{T,U}, 
        indegree_sequence::Vector{<:Integer},
        outdegree_sequence::Vector{<:Integer};
        allow_self_loops::Bool = false,
        perform_checks::Bool = false
    ) where {T,U}

Return a random `MultilayerDiGraph` with degree sequence `degree_sequence`. `allow_self_loops` controls the presence of self-loops, while if `perform_checks` is true, the `degree_sequence` os checked to be graphical.
"""
function MultilayerDiGraph(
    empty_multilayerdigraph::MultilayerDiGraph{T,U}, 
    indegree_sequence::Vector{<:Integer},
    outdegree_sequence::Vector{<:Integer};
    allow_self_loops::Bool = false,
     perform_checks::Bool = false
    ) where {T,U}

    (allow_self_loops && perform_checks) && @warn "Checks for graphicality and coherence with the provided `empty_multilayerdigraph` are currently performed without taking into account self-loops. Thus said checks may fail event though the provided `indegree_sequence` and `outdegree_sequence` may be graphical when one allows for self-loops within the directed multilayer graph to be present. If you are sure that the provided `indegree_sequence` and `outdegree_sequence` are indeed graphical under those circumstances, you may want to disable checks by setting `perform_checks = false`. We apologize for the inconvenient."
    
    _multilayerdigraph = deepcopy(empty_multilayerdigraph)

    ne(_multilayerdigraph) == 0 || throw(ErrorException("The `empty_multilayerdigraph` argument should be an empty MultilayerDiGraph. Found $(ne(_multilayerdigraph)) edges."))

    if perform_checks
        n = nv(_multilayerdigraph)
        n == length(indegree_sequence) == length(outdegree_sequence)  || throw(ErrorException("The number of vertices of the provided empty MultilayerDiGraph does not match the length of the `indegree_sequence` or the `outdegree_sequence`. Found $(nv(_multilayerdigraph)) , $(length(indegree_sequence)) and  $(length(outdegree_sequence))"))

        sum(indegree_sequence) == sum(outdegree_sequence) || throw(ErrorException("The sum of the `indegree_sequence` and the `outdegree_sequence` must match. Found $(sum(indegree_sequence)) and $(sum(outdegree_sequence))"))

        isdigraphical(degree_sequence) || throw(ArgumentError("`indegree_sequence` and `outdegree_sequence` must be digraphical"))
        
    end
    
    # edge_list = _random_directed_configuration(_multilayerdigraph, indegree_sequence, outdegree_sequence, allow_self_loops)
    equivalent_graph = kleitman_wang_graph_generator(indegree_sequence, outdegree_sequence)

    edge_list = [ME(_multilayerdigraph.v_V_associations[src(edge)], _multilayerdigraph.v_V_associations[dst(edge)]) for edge in edges(equivalent_graph) ]

    for edge in edge_list
        add_edge!(_multilayerdigraph, edge)
    end

    return _multilayerdigraph
end

"""
    MultilayerDiGraph(n_nodes::Int64, T::Type{ <: Number}, U::Type{ <: Number} )

Return a null MultilayerDiGraph with with vertex type `T` weighttype `U`. Use this constructor and then add Layers and Interlayers via the `add_layer!` and `specify_interlayer!` methods.
"""
MultilayerDiGraph(T::Type{<:Number}, U::Type{<:Number}) =  MultilayerDiGraph{T,U}(
                                                                                LayerDescriptor{T,U}[],
                                                                                OrderedDict{Set{Symbol},InterlayerDescriptor{T,U}}(),
                                                                                Bijection{ T , MultilayerVertex}(),
                                                                                Bijection{ Int64 , Node}(),
                                                                                Vector{HalfEdge{MultilayerVertex, <: Union{Nothing, U}}}[],
                                                                                Vector{HalfEdge{MultilayerVertex, <: Union{Nothing, U}}}[],
                                                                                Dict{T, Union{Tuple,NamedTuple}}()
)

"""
    MultilayerDiGraph(layers::Vector{<:Layer{T,U}}; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String="multiplex") where {T,U, H <: AbstractGraph{T}}

Construct a MultilayerDiGraph with layers `layers` and all interlayers with structure `default_interlayers_structure` (only "multiplex" and "empty" are allowed) and type `default_interlayers_null_graph`.
"""
MultilayerDiGraph(layers::Vector{<:Layer{T,U}}; default_interlayers_null_graph::H = SimpleGraph{T}(), default_interlayers_structure::String="multiplex") where {T,U, H <: AbstractGraph{T}} =  MultilayerDiGraph(layers, Interlayer{T,U}[]; default_interlayers_null_graph = default_interlayers_null_graph, default_interlayers_structure = default_interlayers_structure)

# Base overloads
"""
    Base.getproperty(mg::M, f::Symbol) where { M <: MultilayerDiGraph }
"""
function Base.getproperty(mg::MultilayerDiGraph, f::Symbol) # where {T,U,M<:AbstractMultilayerDiGraph{T,U}}
    if f in (:v_V_associations, :fadjlist, :badjlist, :idx_N_associations, :layers, :interlayers, :v_metadata_dict) # :weight_tensor, :supra_weight_matrix, 
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