using Revise
using MultilayerGraphs

# The users are represented as `Node`s, and everyuser is od course represented in every layer
const n_users = 100
const users = [Node("user_$i") for i in 1:n_users]

# The layers (follower, retweet, like and mention networks) will be instantiated as random layers. We define some useful variables for the instantiation of random layers hereafter:
const weighttype = Int64                   # The retweet, like and mention networks will be represented by weighted graphs with integer weights
const min_edges = 50                       # Minimum number of edges in each layer
const max_edges = 1000                     # Maximum number of edges in each layer


# Next we deine the layers
## Follower network
const follower_network = layer_simpledigraph( 
                                            :follower_network,          # Name of the layer
                                            users,                      # The `Node`s it will represent
                                            rand(min_edges:max_edges),   # The number of edges that will be randomly diestributed
                                            weighttype = weighttype
)

## Retweet network
const retweet_network = layer_simpleweighteddigraph(
                                                    :retweet_network, 
                                                    users, 
                                                    rand(min_edges:max_edges), 
                                                    default_edge_weight = (src,dst) -> rand(1:10),  # Function that assigns a weight to each edge upon creation
                                                    weighttype = weighttype                         # The weight type
)

## Like network
const like_network = layer_simpleweighteddigraph(
                                                :like_network, 
                                                users, 
                                                rand(min_edges:max_edges), 
                                                default_edge_weight = (src,dst) -> rand(1:10),
                                                weighttype = weighttype                       
)

## Retweet network
const mention_network = layer_simpleweighteddigraph(
                                                :mention_network, 
                                                users, 
                                                rand(min_edges:max_edges), 
                                                default_edge_weight = (src,dst) -> rand(1:10),
                                                weighttype = weighttype                       
)

# Instantiate "twitter" as a multiplex graph made up of the layers above
const twitter = MultilayerDiGraph(
                                [follower_network, retweet_network, like_network, mention_network],
                                default_interlayers_structure = "multiplex" # Interlayers will be automatically specified as containing diagonal couplings only
)


# Add a newly-subscribed user
## The user is represented by a new node
const new_user = Node("user_101")
## Add the user to the MultilayerGraph
add_node!(twitter, new_user)
## Add the corresponding vertex to each layer
for layer_name in twitter.layers_names
    add_vertex!(twitter, MV(new_user, layer_name))
end

# Evaluate some metrics

