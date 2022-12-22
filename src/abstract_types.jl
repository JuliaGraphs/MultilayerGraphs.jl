"""
    abstract type AbstractSubGraph{T <: Integer,U <: Real,G <: AbstractGraph{T}}

An abstract type representing a subgraph (i.e. a layer or an interlayer).
"""
abstract type AbstractSubGraph{T<:Integer,U<:Real,G<:AbstractGraph{T}} end


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
    AbstractInterlayer{T,U,G}

An abstract type representing a generic Interlayer.

# PARAMETRIC TYPES

- `T`: the node type;
- `U`: the adjacency matrix/tensor eltype;
- `G`: the underlying graph type.
"""
abstract type AbstractInterlayer{T,U,G} <: AbstractSubGraph{T,U,G} end


"""
    AbstractMultilayerGraph{T <: Integer, U <: Real} <: AbstractGraph{T}

An abstract type for multilayer graphs. It is a subtype of AbstractGraph and its concrete subtypes may extend Graphs.jl.
"""
abstract type AbstractMultilayerGraph{T<:Integer,U<:Real} <: AbstractGraph{T} end

"""
    AbstractMultilayerUGraph{T,U} <: AbstractMultilayerGraph{T,U} 

Abstract type representing an undirected multilayer graph.
"""
abstract type AbstractMultilayerUGraph{T,U} <: AbstractMultilayerGraph{T,U} end

"""
    AbstractMultilayerDiGraph{T,U} <: AbstractMultilayerGraph{T,U} 

Abstract type representing an undirected multilayer graph.
"""
abstract type AbstractMultilayerDiGraph{T,U} <: AbstractMultilayerGraph{T,U} end



"""
    AbstractTensorRepresentation{U}

An abstract type encoding a generic tensorial representation of the links and metadata of a multilayer graph. 

Concrete subtypes must have an `array` field (a 4-dimensional tensor of eltype U, indexes as [source_node_idx, destination_node_idx, source_layer_idx, destination_layer_idx ]).

# PARAMETRIC TYPES

- `U`: the weight type of the multilayer graph.
"""
abstract type AbstractTensorRepresentation end


"""
    AbstractMatrixRepresentation{T,U}

An abstract type encoding a generic matrix representation of the links and metadata of a multilayer graph. 

Concrete subtypes must have an `array` field (a matrix of eltype U) and a `v_V_associations` (a `Bijection{T, Union{MissingVertex, MultilayerVertex}}`).

# PARAMETRIC TYPES

- `T`: the type of the internal representation of vertices;
- `U`: the weight type of the multilayer graph.
"""
abstract type AbstractMatrixRepresentation end
