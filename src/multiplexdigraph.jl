"""
    MultiplexDiGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultiplexUGraph{T,U}

A concrete type that can represent a general directed MultiplexDiGraph graph. That is, a multilayer graph whos einterlayers are all trivial in the sense that they only have links between vertices that represent the same node.

# FIELDS

- `adjacency_tensor::Array{U, 4}`: the 4-dimensional tensor that encodes all (weighted) connections within the graph. adjacency_tensor[1,2,3,4] encodes the strength of the (directed or undirected) link between node 1 in layer 3 and node 2 in layer 4.
- `layers::OrderedDict{ Tuple{Int64,Int64}, Layer{T,U,G}}`: the ordered dictionary containing all the layers of the MultiplexDiGraph graph. Their underlying graphs must be all undirected.
- `interlayers::OrderedDict{ Tuple{Int64,Int64}, Interlayer{T,U}}`: the ordered dictionary containing all the interlayers of the MultiplexDiGraph graph. 
"""
mutable struct MultiplexDiGraph{T,U} <: AbstractMultiplexDiGraph{T,U}
    adjacency_tensor::Array{U,4}
    layers::OrderedDict{Tuple{Int64,Int64},Layer{T,U}}
    interlayers::OrderedDict{Tuple{Int64,Int64},Interlayer{T,U}}

    function MultiplexDiGraph(adjacency_tensor::Array{U,4}, layers::OrderedDict{Tuple{Int64,Int64},Layer{T,U}}, interlayers::OrderedDict{Tuple{Int64,Int64},Interlayer{T,U}}) where {T,U}

        
        @assert all(is_multiplex_interlayer.(values(interlayers)))
        
        
        return new{T,U}(adjacency_tensor, layers, interlayers)

    end
    
end

@traitimpl IsDirected{MultiplexDiGraph}
@traitimpl IsWeighted{MultiplexDiGraph}

"""
    MultiplexDiGraph(num_layers::Int64, n_nodes::Int64, min_edges::Int64, max_edges::Int64, graph_types::Vector{DataType})

Return a random `MultiplexDiGraph` with `num_layers` layers, `n_nodes` nodes and each `Layer` and `Interlayer` has a random number of edges between `min_edges` and `max_edges`. `Layers` underlying graph type is randomly chosen among `graph_types`.
"""
function MultiplexDiGraph(
    num_layers::Int64,
    n_nodes::Int64,
    min_edges::Int64,
    max_edges::Int64,
    graph_types::Vector{DataType},
)

    #TODO: this method requires that all Graphs.jl extensions are such that:
    # The first type parameter of graphs types is the type of the nodes, like the T in SimpleGraph{T};
    # all graph types have a random constructor with signature (nv::Int64, ne::Int64)
    # all graph types have an empty constructor with signature ()
    all([istrait(IsDirected{graph_type}) for graph_type in graph_types]) || throw(ErrorException("Not all graph types in argument `graph_types` are directed."))
    T = get_common_type([graph_type.parameters[1] for graph_type in graph_types])
    U = get_common_type(
        eltype.([adjacency_matrix(graph_type()) for graph_type in graph_types])
    )

    layers = Layer{T,U}[]
    for i in 1:num_layers
        graph_type = rand(graph_types)
        push!(
            layers,
            Layer(Symbol("layer_$i"), graph_type(n_nodes, rand(min_edges:max_edges)); U=U),
        )
    end

    return MultiplexDiGraph(layers)::MultiplexDiGraph{T,U}
end

"""
    MultiplexDiGraph(n_nodes::Int64, T::Type{ <: Number}, U::Type{ <: Number} )

Return a MultiplexDiGraph with `n_nodes` of type `T` nodes and an adjacency tensor eltype `U`. Use this constructor and then add Layers and Interlayers via the `add_layer!` and `add_interlayer!` methods.
"""
function MultiplexDiGraph(n_nodes::Int64, T::Type{<:Number}, U::Type{<:Number})
    adjacency_tensor = zeros(U, n_nodes, n_nodes, 0, 0)
    return MultiplexDiGraph(
        adjacency_tensor,
        OrderedDict{Tuple{Int64,Int64},Layer{T,U}}(),
        OrderedDict{Tuple{Int64,Int64},Interlayer{T,U}}(),
    )
end

#= """
    MultiplexDiGraph(layers::Vector{ <: Layer{T,U}};  default_interlayer::String  = "multiplex") where {T,U} 

Construct a MultiplexDiGraph with layers `layers` and all interlayers of type `default_interlayer` (only "multiplex" is allowed).
"""
function MultiplexDiGraph(
    layers::Vector{<:Layer{T,U}}; default_interlayer::String="multiplex"
) where {T,U}
    return MultiplexDiGraph(
        layers, Interlayer{T,U}[]; default_interlayer=default_interlayer
    )
end
 =#
"""
    MultiplexDiGraph(layers::Vector{ <: Layer{T,U }}, specified_interlayers::Vector{ <: Interlayer{T,U}};  default_interlayer::String  = "multiplex") where {T, U} 

Construct a MultiplexDiGraph with layers given by `layers`. The interlayers will be constructed by default according to `default_interlayer` (only `"multiplex"` is allowed), except for those specified in `specified_interlayers`.
"""
function MultiplexDiGraph(
    layers::Vector{<:Layer{T,U}},
) where {T,U}
    # Check that all layers and specified interlayers are directed
    all([
        istrait(IsDirected{typeof(layer.graph)}) for layer in layers
    ]) || throw(
        ErrorException(
            "Not all the underlying layers' graphs are directed."
        ),
    )

    num_nodes = nv(layers[1])
    multilayerdigraph = MultiplexDiGraph(num_nodes, T, U)

    for layer in layers
        add_layer!(
            multilayerdigraph, layer
        )
    end

    return multilayerdigraph
end