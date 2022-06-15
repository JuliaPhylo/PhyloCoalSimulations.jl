@testset "one population" begin

@test isempty(PCS.simulatecoal_onepopulation(PN.Edge[], 2.0, 5))

e1 = PCS.initializetip("s","1",1,"",0.1)
f1 = PCS.initializetip(e1.node[1],1,3)
@test length(f1) == 1
@test f1[1].node[1].name == "s1"
f2 = PCS.initializetip(e1.node[1],2,4)
@test length(f2) == 2
@test [e.node[1].name for e in f2] == ["s1_1", "s1_2"]
@test [e.node[1].number for e in f2] == [4,5]

e2 = PCS.initializetip("s","2",2,"",0.2)
Random.seed!(432)
f1 = PCS.simulatecoal_onepopulation([e1,e2], 0.1, 3)
@test length(f1) == 2
@test all(f1 .=== (e1,e2))
@test e1.length ≈ 0.2
@test e2.length ≈ 0.3

f2 = PCS.simulatecoal_onepopulation(f1, Inf, 4)
@test length(f2) == 1
net = PCS.convert2tree!(f2[1].node[1])
@test PN.writeTopology(net, round=true, digits=1) == "(s2:1.1,s1:1.0);"

end
