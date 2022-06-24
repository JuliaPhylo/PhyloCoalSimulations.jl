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
  
#((t5:1.227715515,((#H17:0.7352338122)#H14:0.0002739851344,(t7:0.11848377,t2:0.11848377):1.0947791542999998):0.01445259019):0.383743467,(#H14:0.1395317581,(t8:0.477755127)#H17:0.8747655703):0.2589382843);

end
