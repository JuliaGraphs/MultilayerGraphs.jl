"""
    abstract type AbstractNode

An abstract type representing a node.
"""
abstract type AbstractNode end

"""
    struct Node <: AbstractNode

A custom concrete type representing a node of a multilayer graph.

# FIELDS

-`id::String`: the node's id i.e. a String stating its name.
"""
struct Node <: AbstractNode
    id::String
end

# Base overrides
Base.:(==)(lhs::Node, rhs::Node) = (lhs.id == rhs.id)
Base.isequal(lhs::Node, rhs::Node) = (lhs.id == rhs.id)