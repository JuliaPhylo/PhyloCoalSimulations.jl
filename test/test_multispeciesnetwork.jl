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

#= basic check on node names and inCycle values
inCycle: numberID of net's edge that the node maps to for degree-3 nodes.
      -1 for nodes of degree 1 or 2 (which map to a node, not an edge)
      max(edge numberID)+1 for degree-3 nodes that map to the infinite root edge.
name: name of numberID (as string) of net's node that the node maps to,
      for gene nodes of degree 1 or 2. No name for degree-3 nodes or root.
=#
speciesedgenumbers = [e.number for e in net.edge]
push!(speciesedgenumbers, maximum(speciesedgenumbers)+1)
speciesleafnames = tipLabels(net)
speciesintnodenames = String[]
for n in net.node
    n.leaf && continue
    push!(speciesintnodenames, (n.name == "" ? replace(string(n.number), "-" => "minus") : n.name))
end
function checknodeattributes(genetree)
  nd2 = 0 # number of degree-2 nodes in gene tree
  for node in genetree.node
    isroot = node === genetree.node[genetree.root]
    degree = length(node.edge)
    @test degree in [1,2,3]
    if degree == 2 nd2 += 1; end
    if degree == 1
        @test node.leaf
        @test replace(node.name, r"_\d*$" => "") in speciesleafnames # rm "_individualnumber" from tip name
        @test node.inCycle == -1
    elseif degree == 2 && !isroot
        @test node.name in speciesintnodenames # degree-2 nodes have a name
        @test node.inCycle == -1
    else # root (of degree 2) or degree == 3
        @test node.name == ""
        @test node.inCycle in speciesedgenumbers
    end
  end
  return nd2
end
myseed = 1624
Random.seed!(myseed)
gt1 = simulatecoalescent(net, 1, 3)[1]
ndegree2 = checknodeattributes(gt1)
@test ndegree2 == 1 # root only, when nodemapping=false
Random.seed!(myseed)  # same as before
gt2 = simulatecoalescent(net, 1, 3; nodemapping=true)[1]
ndegree2 = checknodeattributes(gt2)
@test ndegree2 >= 9
# fuse degree-2 nodes in gt2, then check that gt2 == gt1
PN.removedegree2nodes!(gt2, true)
@test hardwiredClusterDistance(gt1, gt2, true)==0
# ideally: also check for equal edge lengths. Not done here

#= simulate 1000 gene trees, 1 indiv/species, then check:
- correct quartet CFs (approximately)
- correct distribution of pairwise distances: mixture of shifted exponential(s)
on same net as earlier (4-species net, but rotated)
=#
net = PN.readTopology("((t5:1.228,((t7:0.118,t2:0.118):1.095,(#H17:0.735::0.4)#H14:0.0::0.7):0.014):0.384,(#H14:0.14::0.3,(t8:0.478)#H17:0.875::0.6):0.259);")
nsim = 1000
Random.seed!(1602)
genetrees = PCS.simulatecoalescent(net, nsim, 1)
α = 0.05 # to test "approximately" correct results
#= expected CF:
plot(net, showedgenumber=true, showgamma=true);
t2,t7 sister after cut edge 6, t8 only descendant of either hybrid node, so:
expCFminor = exp(-net.edge[4].length) * ((0.6 + 0.4*0.3)*exp(-net.edge[7].length) + 0.4*0.7)/3
=#
expCFminor = 0.11039698101750112 # for t2t5|t7t8 and t2t8|t5t7
expCF = [expCFminor, 1-2*expCFminor, expCFminor]
obsCF, taxa = PN.countquartetsintrees(genetrees) # taxa alphabetically: t2,t5,t7,t8
obsCF = Int.(obsCF[1].data[1:3]*nsim) # change format to what HypothesisTests expects
qCFtest = HypothesisTests.ChisqTest(obsCF, expCF)
@test pvalue(qCFtest) >= α
# get distances between t2-t7, or t2-t5, or t5-t7: those should be minimum value + exponential.
# (note for later: there is repetition below that could probably be vectorized, or nested inside for-loops)
d_t2t7_min = 0.236 # sum(net.edge[i].length for i in [2,3])
d_t2t5_min = 2.455 # sum(net.edge[i].length for i in [1,3,4,7])
d_t5t7_min = 2.455 # sum(net.edge[i].length for i in [1,2,4,7])
#= d_t2t8 ~ mixture((gam1 * Dirac(min1) + gam2 * Dirac(min2)) + Exp(2))
   but not sure how to code this mixture of convolutions with Distributions.
d_t2t8_min1 = 2.426 # sum(net.edge[i].length for i in [3,4, 10,5,6])
d_t2t8_min2 = 3.223 # sum(net.edge[i].length for i in [3,4,7,8, 10,11,12])
d_t2t8_min2 = 3.223 # sum(net.edge[i].length for i in [3,4,7,8, 10,5,9,12])
d_t2t8_gam1 = 0.28  # 0.4*0.7: prod(net.edge[i].gamma for i in [5,6])
d_t2t8_gam2 = 0.72
d_t2t8 = Float64[]
=#
d_t2t7 = Float64[]; d_t2t5 = Float64[]; d_t5t7 = Float64[]
for tree in genetrees
    distances = PN.pairwiseTaxonDistanceMatrix(tree)
    o = sortperm(PN.tipLabels(tree)) # order to get taxa alphabetically: t2, t5, t7, t8
    push!(d_t2t7, distances[o[1], o[3]])
    push!(d_t2t5, distances[o[1], o[2]])
    push!(d_t5t7, distances[o[2], o[3]])
    # push!(d_t2t8, distances[o[2], o[4]])
end
expdist2 = Distributions.Exponential(2)
@test pvalue(HypothesisTests.OneSampleADTest(d_t2t7 .- d_t2t7_min, expdist2)) >= α
@test pvalue(HypothesisTests.OneSampleADTest(d_t2t5 .- d_t2t5_min, expdist2)) >= α
@test pvalue(HypothesisTests.OneSampleADTest(d_t5t7 .- d_t5t7_min, expdist2)) >= α

# on a tree with #generations + Ne dictionary
net = PN.readTopology("(A:1000,B:1000);")
Ne = Dict(1=>20, 2=>300, 3=>300)
Random.seed!(639)
genetree = PCS.simulatecoalescent(net, 1, 2, Ne)[1]
@test all(sort!([e.length for e in genetree.edge]) .> [1,1, 15,15, 800,800])
end
