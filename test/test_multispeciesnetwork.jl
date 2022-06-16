@testset "multi-species network" begin

# on a tree
net = PN.readTopology("(A:1,B:1);")
Random.seed!(432)
genetree = simulatecoalescent(net, 1, 2)
@test length(genetree) == 1
@test sort(tipLabels(genetree[1])) == ["A_1","A_2","B_1","B_2"]
# differing number of individuals / species
genetree = simulatecoalescent(net, 1, Dict("A"=>2, "B"=>1))
@test sort(tipLabels(genetree[1])) == ["A_1","A_2","B"]

end
