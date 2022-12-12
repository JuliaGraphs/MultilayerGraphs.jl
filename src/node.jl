"""
    abstract type AbstractNode

An abstract type representing a node.
"""
abstract type AbstractNode end

"""
    struct Node <: AbstractNode

A custom concrete type representing a node of a multilayer graph.
"""
struct Node <: AbstractNode
    id::String
end

"""
    id(n::Node)

Return the id of `n`.
"""
id(n::Node) = n.id

# Base overrides
Base.:(==)(lhs::Node, rhs::Node) = (lhs.id == rhs.id)
Base.isequal(lhs::Node, rhs::Node) = (lhs.id == rhs.id)