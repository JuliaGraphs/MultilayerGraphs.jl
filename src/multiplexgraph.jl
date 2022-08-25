"""
    MultiplexGraph{T, U, G <: AbstractGraph{T}} <: AbstractMultiplexUGraph{T,U}

A concrete type that can represent a general multilayer graph.

# FIELDS

- `adjacency_tensor::Array{U, 4}`: the 4-dimensional tensor that encodes all (weighted) connections within the graph. adjacency_tensor[1,2,3,4] encodes the strength of the (directed or undirected) link between node 1 in layer 3 and node 2 in layer 4.
- `layers::OrderedDict{ Tuple{Int64,Int64}, Layer{T,U,G}}`: the ordered dictionary containing all the layers of the multilayer graph. Their underlying graphs must be all undirected.
- `interlayers::OrderedDict{ Tuple{Int64,Int64}, Interlayer{T,U}}`: the ordered dictionary containing all the interlayers of the multilayer graph. Their underlying graphs must be all undirected.
"""
mutable struct MultiplexGraph{T,U} <: AbstractMultiplexUGraph{T,U}
    adjacency_tensor::Array{U,4}
    layers::OrderedDict{Tuple{Int64,Int64},Layer{T}}
    interlayers::OrderedDict{Tuple{Int64,Int64},Interlayer{T}}
    
    function MultiplexGraph(adjacency_tensor::Array{U,4}, layers::OrderedDict{Tuple{Int64,Int64},Layer{T}}, interlayers::OrderedDict{Tuple{Int64,Int64},Interlayer{T}}) where {T,U}

        
        @assert all(is_multiplex_interlayer.(values(interlayers)))
        
        
        return new{T,U}(adjacency_tensor, layers, interlayers)

    end
end

"""
    MultiplexGraph(num_layers::Int64, n_nodes::Int64, min_edges::Int64, max_edges::Int64, graph_types::Vector{DataType})

Return a random `MultiplexGraph` with `num_layers` layers, `n_nodes` nodes and each `Layer` has a random number of edges between `min_edges` and `max_edges`. `Layers` have parametric type `graph_type`.
"""
function MultiplexGraph(
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
    all([!istrait(IsDirected{graph_type}) for graph_type in graph_types]) || throw(ErrorException("Not all graph types in argument `graph_types` are undirected."))
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

    return MultiplexGraph(layers)
end

@traitimpl IsWeighted{MultiplexGraph}

"""
    MultiplexGraph(n_nodes::Int64, T::Type{ <: Number}, U::Type{ <: Number} )

Return a MultiplexGraph with `n_nodes` of type `T` nodes and an adjacency tensor eltype `U`. Use this constructor and then add Layers via the `add_layer!` method.
"""
function MultiplexGraph(n_nodes::Int64, T::Type{<:Number}, U::Type{<:Number})
    adjacency_tensor = zeros(U, n_nodes, n_nodes, 0, 0)
    return MultiplexGraph(
        adjacency_tensor,
        OrderedDict{Tuple{Int64,Int64},Layer{T}}(),
        OrderedDict{Tuple{Int64,Int64},Interlayer{T}}(),
    )
end

"""
    MultiplexGraph(layers::Vector{ <: Layer{T,U }}) where {T, U}      

Construct a MultiplexGraph with layers given by `layers`.
"""
function MultiplexGraph(
    layers::Vector{<:Layer{T,U}},
) where {T,U}
    # Check that all layers are undirected
    # specified_subgraphs = layers
    all([
        !istrait(IsDirected{typeof(subgraph.graph)}) for subgraph in layers
    ]) || throw(
        ErrorException(
            "Not all the underlying layers' graphs are undirected."
        ),
    )

    num_nodes = nv(layers[1])
    multiplexgraph = MultiplexGraph(num_nodes, T, U)

    for layer in layers
        add_layer!(multiplexgraph, layer)
    end

    return multiplexgraph
end