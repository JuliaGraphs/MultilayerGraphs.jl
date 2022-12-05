# FIXME: Agents.jl functions add_node! and rem_node! should be renames in add_vertex! and rem_vertex!. There is no reason for not following graphs.jl name conventions. Moreover, add_edge! is extended while rem_edge! isn't. We should open and issue after releasing it.

"""
    add_node!(model::ABM{<:GraphSpace}, node::Node)
Add a new `Node` to the model's graph.
"""
function Agents.add_node!(model::ABM{<:GraphSpace}, node::Node)
    add_node!(model.space.graph, node)
end

"""
    rem_node!(model::ABM{<:GraphSpace{ <: AbstractMultilayerGraph}}, n::Int)

Internal fallback of `rem_node!(model::ABM{<:GraphSpace{ <: AbstractMultilayerGraph}}, node::Node)`
"""
function Agents.rem_node!(model::ABM{<:GraphSpace{ <: AbstractMultilayerGraph}}, n::Int)
    n in domain(model.space.graph.idx_N_associations) || return false
    rem_node!(model, model.space.graph.idx_N_associations[n])
end
"""
    rem_node!(model::ABM{<: GraphSpace}, n::Int)
Remove `Node` `node` from the model's underlying multilayer graph. All agents in vertices representing `node` are killed.
"""
function Agents.rem_node!(model::ABM{<:GraphSpace{ <: AbstractMultilayerGraph}}, node::Node)
    has_node(model.space.graph, node) || return false

    for v in vertices(model.space.graph)
        mv = model.space.graph.v_V_associations[v]
        if mv.node == node
            rem_vertex!(model, v)
        end
    end

    idx_tbr = mg.idx_N_associations(n)
    delete!(model.space.graph.idx_N_associations, idx_tbr)

    return true
end