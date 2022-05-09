# function MultilayerGraph_ten(layers::OrderedDict{Symbol, G}, interlayers::OrderedDict{Symbol, H}) where {G <: AbstractGraph, H <: InterLayer}#Symmetric{Pair{Symbol,H}, Matrix{Pair{Symbol,H}}}
#     # Check that the vertices are the same in all layers and interlayers
#     @assert length(unique(vertices.(values(layers)))) == 1
#     # Check that all layers have unique names
#     @assert length(unique(keys(layers))) == length(layers)
#     # Check that the correct amount of interlayers has been given
#     num_layers = length(layers)
#     @assert length(interlayers) <= factorial(num_layers-1)

#     # types = unique([typeof(layer).parameters[1] for layer in vcat(collect(values(layers)), collect(values(interlayers))) ])
#     # t = length(types) == 1 ? types[1] : get_common_type(types)


#     adjacency_tensor = zeros(Float64, nv(layers[1]),nv(layers[1]), num_layers,num_layers ) #Array{Float64}(undef, nv(layers[1]),nv(layers[1]), num_layers,num_layers )
#     # Fill interlayer adjacency matrices
#     for (i,layer) in enumerate(values(layers))
#         adjacency_tensor[:, :, i, i] .= Matrix(adjacency_matrix(layer))
#     end
#     # Fill intralayer adjacency matrices (they are symmetric)
#     layers_symbols = collect(keys(layers))
#     for interlayer in values(interlayers)
#         row_idx = findfirst(x -> x == interlayer.layer_1, layers_symbols)
#         col_idx = findfirst(x -> x == interlayer.layer_2, layers_symbols)
#         @assert row_idx != col_idx
#         adjacency_tensor[:,:, row_idx, col_idx] .= adjacency_tensor[:,:, col_idx, row_idx] .= Matrix(adjacency_matrix(interlayer.graph))
#     end


#     # for i in 1:(num_layers-1)
#     #     for j in (i+1):num_layers
#     #         adjacency_tensor[:, :, i, j] .= Matrix(adjacency_matrix(interlayers[i,j])) 
#     #         adjacency_tensor[:, :, i+1, i] .= adjacency_tensor[:, :, i, i+1]
#     #     end
#     # end
#     println(typeof(layers[1]).parameters[1])
#     #println(adjacency_tensor)
#     println(merge(layers,interlayers))
#     return MultilayerGraph_ten{typeof(layers[1]).parameters[1]}(adjacency_tensor, convert(OrderedDict{Symbol, Union{ <: AbstractGraph, <: InterLayer}}  , merge(layers,interlayers)))

# end