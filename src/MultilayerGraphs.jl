module MultilayerGraphs

# Write your package code here.

export MultilayerGraph_ten, MultilayerGraph_m, MultilayerGraph_mcsc, MultilayerGraph_csc, Layer, InterLayer,MultiplexGraph

using LinearAlgebra,SparseArrays
using DataStructures
using Graphs, SimpleWeightedGraphs, MetaGraphs


include("layer.jl")
include("interlayer.jl")
include("utilities.jl")

"""
    abstract type MultilayerGraph{T <: Integer}

An abstract type for multilayer graphs, it must contain a fied `adjacency_tensor`. It is a subtype of AbstractGraph
"""
abstract type AbstractMultilayerGraph{T <: Integer} <: AbstractGraph{T} end

"""
AbstractMultiplexGraph{T}
"""
abstract type AbstractMultiplexGraph{T} <: AbstractMultilayerGraph{T} end


"""
"""
mutable struct MultilayerGraph_ten{T} <: AbstractMultilayerGraph{T}
    adjacency_tensor::Array{Float64, 4}
    layers::Vector{ <: Layer}
    interlayers::Vector{<: InterLayer} #Union{Nothing,OrderedDict{Symbol, Union{ <: AbstractGraph, <: InterLayer}} }
end

"""
    MultilayerGraph_ten(layers::Vector{ <: Layer }, specified_interlayers::Union{Vector{ <: InterLayer}};  default_interlayer::Symbol  = :multiplex) 
"""
function MultilayerGraph_ten(layers::Vector{ <: Layer }, specified_interlayers::Union{Vector{ <: InterLayer}};  default_interlayer::Symbol  = :multiplex) #where {G <: AbstractGraph, H <: InterLayer}#Symmetric{Pair{Symbol,H}, Matrix{Pair{Symbol,H}}}
    # Check that all the layers and specified_interlayers have the same (vertex) type  parameter
    @assert check_unique([typeof(graph).parameters[1] for graph in getproperty.(vcat(layers,specified_interlayers), Ref(:graph))] )
    vertex_type = typeof(layers[1].graph).parameters[1]
    # Check that the vertices are the same in all layers and specified_interlayers
    @assert check_unique(vertices.(layers)) == 1 # length(unique(vertices.(layers))) == 1
    n_nodes = nv(layers[1])
    # Check that all layers have unique names
    @assert check_unique(getproperty.(layers, Ref(:name)))
    @assert check_unique(getproperty.(specified_interlayers, Ref(:name))) #length(unique([interlayer.name for layer in specified_interlayers])) == length(specified_interlayers)
    # Check that the correct amount of specified_interlayers has been given
    num_layers = length(layers)
    @assert length(specified_interlayers) <= factorial(num_layers-1)

    # Initialize the adjacency tensor
    adjacency_tensor            = zeros(Float64, nv(layers[1]),nv(layers[1]), num_layers,num_layers ) #Array{Float64}(undef, nv(layers[1]),nv(layers[1]), num_layers,num_layers )
    # Array of tuples representing all the possible pair of values that the last two indexes of `adjacency_tensor` can take .
    all_layers_interlayers_idxs = vec(collect(Iterators.product(1:num_layers, 1:num_layers )))
    # Fill interlayer adjacency matrices
    for (i,layer) in enumerate(layers)
        adjacency_tensor[:, :, i, i] .= Matrix(adjacency_matrix(layer.graph))
        deleteat!(all_layers_interlayers_idxs, findfirst(x -> all(x .== (i,i)), all_layers_interlayers_idxs))
    end
    # Fill intralayer adjacency matrices (they are symmetric)
    layers_symbols = getproperty.(layers, Ref(:name)) #collect(keys(layers))
    for interlayer in specified_interlayers
        row_idx = findfirst(x -> x == interlayer.layer_1, layers_symbols)
        col_idx = findfirst(x -> x == interlayer.layer_2, layers_symbols)
        adjacency_tensor[:,:, row_idx, col_idx] .= adjacency_tensor[:,:, col_idx, row_idx] .= Matrix(adjacency_matrix(interlayer.graph))
        deleteat!(all_layers_interlayers_idxs, findfirst(x -> all(x .== (row_idx,col_idx)), all_layers_interlayers_idxs))
        deleteat!(all_layers_interlayers_idxs, findfirst(x -> all(x .== (col_idx,row_idx)), all_layers_interlayers_idxs))
    end
    

    # Intralayers that are not specified are filled with diagonal matrices (multiplex) or null matrices
    interlayers = deepcopy(specified_interlayers)
    if length(all_layers_interlayers_idxs) != 0
        for unspecified_interlayers_idxs in all_layers_interlayers_idxs
            if default_interlayer == :multiplex
                adjacency_matrix = Matrix{Float64}(I(n_nodes))
                adjacency_tensor[:,:, unspecified_interlayers_idxs[1], unspecified_interlayers_idxs[2]] .= adjacency_tensor[:,:, unspecified_interlayers_idxs[2], unspecified_interlayers_idxs[1]] .= adjacency_matrix
                push!(interlayers, InterLayer(Symbol("interlayer_$(unspecified_interlayers_idxs[1])_$(unspecified_interlayers_idxs[2])"), layers[unspecified_interlayers_idxs[1]].name, layers[unspecified_interlayers_idxs[2]].name, SimpleWeightedDiGraph(adjacency_matrix) ))
            elseif default_interlayer == :null
                adjacency_matrix = Matrix{Float64}(zeros(n_nodes, n_nodes))
                adjacency_tensor[:,:, unspecified_interlayers_idxs[1], unspecified_interlayers_idxs[2]] .= adjacency_tensor[:,:, unspecified_interlayers_idxs[2], unspecified_interlayers_idxs[1]] .= adjacency_matrix
                push!(interlayers, InterLayer(Symbol("interlayer_$(unspecified_interlayers_idxs[1])_$(unspecified_interlayers_idxs[2])"), layers[unspecified_interlayers_idxs[1]].name, layers[unspecified_interlayers_idxs[2]].name, SimpleDiGraph(adjacency_matrix) ))
            else
                error("`default_interlayer` may either be `:multiplex` or `null`: $(default_interlayer) is invalid")
            end
        end
    else
        println("All specified_interlayers have been specified, so none will be automatically assigned")
    end
    
    return MultilayerGraph_ten{vertex_type}(adjacency_tensor, layers, interlayers ) #convert(OrderedDict{Symbol, Union{ <: AbstractGraph, <: InterLayer}} , merge(layers, interlayers) )

end

"""
MultiplexGraph(layers::Vector{ <: Layer })
"""
MultiplexGraph(layers::Vector{ <: Layer }) = MultilayerGraph_ten(layers, InterLayer[])


"""
"""
function Base.getproperty(g::G, f::Symbol) where { G <: AbstractMultilayerGraph }
    if f == :adjacency_tensor
        Base.getfield(g, :adjacency_tensor)
    elseif f == :layers
        Base.getfield(g, :layers)
    elseif f == :interlayers
        Base.getfield(g, :interlayers)
    else
        try
            g.layers[findfirst(layer -> layer.name == f, g.layers)]
        catch
            g.interlayers[findfirst(interlayer -> interlayer.name == f, g.interlayers)]
        end
        # layer_idx      = findfirst(layer -> layer.name == f, g.layers)
        # interlayer_idx = findfirst(interlayer -> interlayer.name == f, g.interlayer)
        # Base.getfield(g, :layers_interlayers)[f]
    end
end






mutable struct MultilayerGraph_m{T} <: AbstractMultilayerGraph{T}
    adjecency_tensor::Matrix{Union{Float64,Int64}}
    layers_interlayers::Union{Nothing, Tuple{Vararg{<: AbstractGraph}}}
end

mutable struct MultilayerGraph_mcsc{T} <: AbstractMultilayerGraph{T}
    adjecency_tensor::Matrix{SparseMatrixCSC{Union{Float64,Int64}, Int64}}
    layers_interlayers::Union{Nothing, Tuple{Vararg{<: AbstractGraph}}}
end

mutable struct MultilayerGraph_csc{T} <: AbstractMultilayerGraph{T}
    adjecency_tensor::SparseMatrixCSC{Union{Float64,Int64}, Int64}
    layers_interlayers::Union{Nothing, Tuple{Vararg{<: AbstractGraph}}}
end




# TODO
# write performant loops and adjacency_tensor structures following https://docs.huihoo.com/julia/0.3/manual/performance-tips/index.html#access-arrays-in-memory-order-along-columns


# mutable struct MultiplexGraph{T} <: AbstractMultiplexGraph{T}
#     graphs::Tuple{Vararg{<:AbstractGraph}}
# end


end