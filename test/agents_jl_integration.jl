################################################
####### TEST INTEGRATION WITH AGENTS.JL ########
################################################

# Create a multilayer graph 
graph = MultilayerGraph(all_layers_u, all_interlayers_u)
# Define the agent type
@agent EigenAgent GraphAgent begin
    previous_value::Float64
    old_value::Float64
    new_value::Float64
end
# Instantiate the agent
function EigenAgent(id, pos; initial_opinion::Float64)
    return EigenAgent(id, pos, initial_opinion, initial_opinion, -1.0)
end
# Define the agent-based model
EigenABM = ABM(EigenAgent, GraphSpace(graph))
# Add agents to the ABM so that agent i is located in vertex i
for (i, mv) in enumerate(mv_vertices(graph))
    initial_opinion = 1.0
    add_agent!(i, EigenABM; initial_opinion=initial_opinion)
end
# Define the individual-level dynamics (micro-dynamics)
function agent_step!(agent::EigenAgent, model)
    agent.previous_value = agent.old_value
    return agent.new_value = sum([
        outneighbor_agent.old_value for
        outneighbor_agent in nearby_agents(agent, model; neighbor_type=:all)
    ])
end
# Define the system-level dynamics (macro-dynamics)
function model_step!(model)
    tot_edges_weight = sqrt(sum(abs2, [agent.new_value for agent in allagents(model)]))
    for agent in allagents(model)
        agent.new_value = agent.new_value / tot_edges_weight
        agent.old_value = agent.new_value
    end
end
# Set the rule to stop the model simulation  
function terminate(model, s)
    # println(s, typeof(s))
    if any(
        !isapprox(a.previous_value, a.new_value; rtol=1e-18) for a in allagents(model)
    ) && !(s > 10000)
        return false
    else
        return true
    end
end
# Simulate the model 
agent_data, _ = run!(
    EigenABM, agent_step!, model_step!, terminate; adata=[:new_value], when=terminate
)
# Compute the eigenvector centrality of the surrounding multilayer graph 
eig_centr_swm, err_swm = eigenvector_centrality(
    EigenABM.space.graph; weighted=false, tol=1e-18, norm="null"
)

for (val_1, val_2) in zip(eig_centr_swm, agent_data.new_value)
    @test val_1 ≈ val_2
end
