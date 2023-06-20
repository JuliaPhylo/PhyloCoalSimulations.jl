@testset "one population" begin

f = PN.Edge[]
@test PCS.simulatecoal_onepopulation!(f, 2.0, 5) == 5
@test isempty(f)

e1 = PCS.initializetip("s","1",1,"")
f1 = PCS.initializetipforest(e1.node[1],1,3)
@test length(f1) == 1
@test f1[1].node[1].name == "s1"
f2 = PCS.initializetipforest(e1.node[1],2,4)
@test length(f2) == 2
@test [e.node[1].name for e in f2] == ["s1_1", "s1_2"]
@test [e.node[1].number for e in f2] == [4,5]

e2 = PCS.initializetip("s","2",2,"")
Random.seed!(432)
f1 = [e1,e2]
@test PCS.simulatecoal_onepopulation!(f1, 0.1, 3) == 3
@test length(f1) == 2
@test all(f1 .=== (e1,e2))
@test e1.length ≈ 0.1
@test e2.length ≈ 0.1

@test PCS.simulatecoal_onepopulation!(f1, Inf, 4) == 5
@test length(f1) == 1
net = PCS.convert2tree!(f1[1].node[1])
testgenetree = (VERSION < v"1.6"  ?   "(s1:0.5,s2:0.5);" : # for julia v1.5
    ( v"1.7" <= VERSION < v"1.10-" ? "(s2:0.9,s1:0.9);" :   # for julia v1.7-v1.9
    "random genetree depends on RNG")) # for later
@test PN.writeTopology(net, round=true, digits=1) == testgenetree

end
