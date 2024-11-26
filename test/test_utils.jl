@testset "utils" begin
net = readnewick("((C:0.9,(B:0.2)#H1:0.7::0.6)i1:0.6,(#H1:0.6::0.4,A:1.0)i2:0.5)i3;")
rng = StableRNG(7)
tree = readnewick(writenewick(simulatecoalescent(rng, net,1,1; nodemapping=true)[1]))

tree.node[1].name = "a name not in the species phylogeny" # give an incorrect name
message = "The gene and species phylogeny have different sets of node names"
@test_throws ErrorException(message) gene_edgemapping!(tree,net)

tree.node[1].name = "B"
gene_edgemapping!(tree,net)
@test (x-> x.inte1).(tree.edge) == [2,5,7,8,1,4,8,6,7,8,8]
end
