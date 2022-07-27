@testset "multi-species network" begin

# on a tree
net = PN.readTopology("(A:1,B:1);")
Random.seed!(432)
genetree = PCS.simulatecoalescent(net, 1, 2)
@test length(genetree) == 1
@test sort(tipLabels(genetree[1])) == ["A_1","A_2","B_1","B_2"]
# differing number of individuals / species
genetree = PCS.simulatecoalescent(net, 1, Dict("A"=>2, "B"=>1))
@test sort(tipLabels(genetree[1])) == ["A_1","A_2","B"]

# on a network with hybrid ladder (and a 3_1-cycle). only t8 below reticulations
net = PN.readTopology("((t5:1.228,((#H17:0.735::0.4)#H14:0.0,(t7:0.118,t2:0.118):1.095):0.014):0.384,(#H14:0.14::0.3,(t8:0.478)#H17:0.875):0.259);")
Random.seed!(432)
genetree = PCS.simulatecoalescent(net, 2, 1)
@test length(genetree) == 2

# simple checks:
# degree-2 nodes have a name, degree-3 nodes don't have a name but have inCycle
# (john to cecile: can you remind me why we wanted it this way, and what the name/inCycle describe?)
# inCycle values should be a number from a species edge
Random.seed!(1624)
genetree = PCS.simulatecoalescent(net, 1, 1)[1]
for inode in 1:length(genetree.node)
    print(inode)
    node = genetree.node[inode]
    isroot = node == genetree.node[genetree.root]
    degree = length(node.edge)

    # inCycle values should be a number from a species edge
    speciesEdgeNumbers = Int64[]
    for speciesEdge in net.edge
        push!(speciesEdgeNumbers, speciesEdge.number)
    end
    @test node.inCycle == -1 || node.inCycle in speciesEdgeNumbers

    if degree == 1
        @test node.leaf && node.name != ""
    elseif degree == 2
        @test isroot || node.name != "" # degree-2 nodes have a name
    elseif degree == 3
        @test node.name == "" && node.inCycle != -1 # degree-3 nodes don't have a name but have inCycle
    else
        @test false # nodes should be degree 1, 2, or 3
    end 
end
# names should be species node names or string from species node number
# (john: ok, but i don't remember why this was useful.)

# to do: simulate with 1000 gene trees, then check a few things, e.g.:
# calculate qCFs, then check that they aren't too far from expectations
# get distances between t2-t7, or t2-t5, or t5-t7: those should be minimum value + exponential.

end
