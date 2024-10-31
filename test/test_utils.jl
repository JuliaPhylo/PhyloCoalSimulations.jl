net = readTopology("((C:0.9,(B:0.2)#H1:0.7::0.6):0.6,(#H1:0.6,A:1):0.5);");
PN.nameinternalnodes!(net, "i"); # "i" is a prefix to name internal nodes

rng = StableRNG(7)

tree = readTopology(writeTopology(simulatecoalescent(rng, net,1,1; nodemapping=true)[1]))

tree.node[1].name = "a name not in the species phylogeny" # give an incorrect name
message = "The gene and species phylogeny have different sets of node names"
@test_throws ErrorException(message) PCS.encode_edges!(tree,net)

tree.node[1].name = "B"

PCS.encode_edges!(tree,net)