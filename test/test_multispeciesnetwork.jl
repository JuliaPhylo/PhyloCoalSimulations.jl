@testset "multi-species network" begin

# on a tree, method taking the default RNG
net = PN.readnewick("(A:1,B:1);")
Random.seed!(432)
genetree = simulatecoalescent(net, 1, 2)
@test length(genetree) == 1
@test sort(tiplabels(genetree[1])) == ["A_1","A_2","B_1","B_2"]
@test isempty(collect(PCS.mappingnodes(genetree[1])))
# differing number of individuals / species
genetree = simulatecoalescent(net, 1, Dict("A"=>2, "B"=>1))
@test sort(tiplabels(genetree[1])) == ["A_1","A_2","B"]
# nodemapping and related utilities
genetree = simulatecoalescent(net, 1, 3; nodemapping=true)[1]
tmp = map(n -> n.name, PCS.mappingnodes(genetree))
@test !isempty(tmp) && all(tmp .== "minus2")
@test Set(unique(population_mappedto.(PN.getchildedge.(PCS.mappingnodes(genetree))))) == Set((1,2))
@test Set(unique(population_mappedto.(genetree.node))) == Set((1,2,3,nothing))
@test Set(unique(population_mappedto.(genetree.edge))) == Set((1,2,3))

# on a network with hybrid ladder (and a 3_1-cycle). only t8 below reticulations
net = PN.readnewick("((t5:1.228,((#H17:0.735::0.4)#H14:0.0,(t7:0.118,t2:0.118):1.095):0.014):0.384,(#H14:0.14::0.3,(t8:0.478)#H17:0.875):0.259);")
Random.seed!(432)
genetree = simulatecoalescent(net, 2, 1)
@test length(genetree) == 2

#= basic check on node names and intn1 values
intn1: numberID of net's edge that the node maps to for degree-3 nodes.
      -1 for nodes of degree 1 or 2 (which map to a node, not an edge)
      max(edge numberID)+1 for degree-3 nodes that map to the infinite root edge.
name: name of numberID (as string) of net's node that the node maps to,
      for gene nodes of degree 1 or 2. No name for degree-3 nodes or root.
=#
speciesedgenumbers = [e.number for e in net.edge]
push!(speciesedgenumbers, maximum(speciesedgenumbers)+1)
speciesleafnames = tiplabels(net)
speciesintnodenames = String[]
for n in net.node
    n.leaf && continue
    push!(speciesintnodenames, (n.name == "" ? replace(string(n.number), "-" => "minus") : n.name))
end
function checknodeattributes(genetree)
  nd2 = 0 # number of degree-2 nodes in gene tree
  for node in genetree.node
    isroot = node === genetree.node[genetree.rooti]
    degree = length(node.edge)
    @test degree in [1,2,3]
    if degree == 2 nd2 += 1; end
    if degree == 1
        @test node.leaf
        @test replace(node.name, r"_\d*$" => "") in speciesleafnames # rm "_individualnumber" from tip name
        @test node.intn1 == -1
    elseif degree == 2 && !isroot
        @test node.name in speciesintnodenames # degree-2 nodes have a name
        @test node.intn1 == -1
    else # root (of degree 2) or degree == 3
        @test node.name == ""
        @test node.intn1 in speciesedgenumbers
    end
  end
  return nd2
end
myseed = 1624
rng = StableRNG(myseed) # Random.seed!(myseed)
gt1 = simulatecoalescent(rng, net, 1, 3)[1]
ndegree2 = checknodeattributes(gt1)
@test ndegree2 == 1 # root only, when nodemapping=false
rng = StableRNG(myseed)  # same as before
gt2 = simulatecoalescent(rng, net, 1, 3; nodemapping=true)[1]
ndegree2 = checknodeattributes(gt2)
@test ndegree2 >= 9
# check edge lengths: okay bc stable RNG
@test writenewick(gt2, round=true, digits=2) == "((((((t2_2:0.12)minus6:1.09)minus4:0.01)minus3:0.38)minus2:0.45,(((((t8_1:0.35,(t8_2:0.07,t8_3:0.07):0.28):0.13)H17:0.74)H14:0.14)minus7:0.26)minus2:0.45):0.55,(((((((t7_3:0.12)minus6:0.06,(t7_1:0.12)minus6:0.06):0.57,((t2_1:0.12)minus6:0.32,((t7_2:0.12)minus6:0.05,(t2_3:0.12)minus6:0.05):0.27):0.31):0.46)minus4:0.01)minus3:0.18,((t5_1:0.2,(t5_3:0.04,t5_2:0.04):0.16):1.03)minus3:0.18):0.2)minus2:1.0);"
# fuse degree-2 nodes in gt2, then check that gt2 == gt1
PN.removedegree2nodes!(gt2, true)
@test hardwiredclusterdistance(gt1, gt2, true)==0

#= simulate 1000 gene trees, 1 indiv/species, then check:
- correct quartet CFs (approximately)
- correct distribution of pairwise distances: mixture of shifted exponential(s)
on same net as earlier (4-species net, but rotated)
=#
net = PN.readnewick("((t5:1.228,((t7:0.118,t2:0.118):1.095,(#H17:0.735::0.4)#H14:0.0::0.7):0.014):0.384,(#H14:0.14::0.3,(t8:0.478)#H17:0.875::0.6):0.259);")
nsim = 1000
rng = StableRNG(1602)
genetrees = simulatecoalescent(rng, net, nsim, 1)
α = 0.40 # to test "approximately" correct results
#= expected CF:
plot(net, showedgenumber=true, showgamma=true);
t2,t7 sister after cut edge 6, t8 only descendant of either hybrid node, so:
expCFminor = exp(-net.edge[4].length) * ((0.6 + 0.4*0.3)*exp(-net.edge[7].length) + 0.4*0.7)/3
=#
expCFminor = 0.11039698101750112 # for t2t5|t7t8 and t2t8|t5t7
expCF = [expCFminor, 1-2*expCFminor, expCFminor]
obsCF, taxa = countquartetsintrees(genetrees) # taxa alphabetically: t2,t5,t7,t8
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
    distances = PN.pairwisetaxondistancematrix(tree)
    o = sortperm(PN.tiplabels(tree)) # order to get taxa alphabetically: t2, t5, t7, t8
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
net = PN.readnewick("(A:1000,B:1000);")
Ne = Dict(1=>20, 2=>300, 3=>300)
rng = StableRNG(639)
genetree = simulatecoalescent(rng, net, 1, 2, Ne)[1]
@test all(sort!([e.length for e in genetree.edge]) .> [1,1, 15,15, 800,800])
# modify the tree to have same population size on all branches
Ne = 20
genetree = simulatecoalescent(rng, net, 1, 2, Ne)[1]
@test all(sort!([e.length for e in genetree.edge]) .> [1,1, 5,5, 800,800])

# test edge cases (like incorrect input)
# (because i am hungry for 100% code coverage)

nindividuals = "this is neither an integer nor a dict of integers"
message = "nindividuals should be an integer or dictionary of integers"
@test_throws ErrorException(message) simulatecoalescent(net, 1, nindividuals)

Ne = "this is neither a number nor a dictionary"
message = "populationsize should be a number or dictionary"
@test_throws ErrorException(message) simulatecoalescent(net, 1, 2, Ne)

end


@testset "correlated inheritance" begin

net = PN.readnewick("((((a:0.01)#H1:0.21::0.6,(#H1:0.1::0.4)#H3:0.11::0.6)i1:0.2)#H2:22.4::0.6,(#H2:11.1,#H3:11.41)i2:11.3)i3;")
# plot(net, showedgelength=true, shownodelabel=true, useedgelength=true);
rng = StableRNG(1234)
res = simulatecoalescent(rng, net, 10, 2; inheritancecorrelation=0.99)
@test all(all(e.length < 5 for e in t.edge) for t in res)
rng = StableRNG(544)
res = simulatecoalescent(rng, net, 10, 2; inheritancecorrelation=0.01)
@test sum(all(e.length < 10 for e in t.edge) for t in res) <= 5
# res = simulatecoalescent(net, 1, 5; inheritancecorrelation=0.99, nodemapping=true)[1]
# plot(res, shownodelabel=true); # see identical mapping from the root to all tips
net = PN.readnewick("((t:0.0)#H1:200::0.6,#H1:200)r;")
tips_alltogether(n) = all(e -> e.length < 100, n.edge)
rng = StableRNG(581)
@testset "correlated inheritance" for rho in [0.8,0.3,0.1]
  res = simulatecoalescent(rng, net, 100, 2; inheritancecorrelation=rho);
  @test isapprox(sum(tips_alltogether.(res))/100,  1 - (2*0.6*0.4)*(1-rho), atol=0.06)
end
# in Ne, and option to not round #generations in gene tree
Ne=0.1
gt = (@test_logs simulatecoalescent(rng, net, 1, 4, Ne; inheritancecorrelation=0.5, round_generationnumber=false))
@test all( 0<e.length<1 for e in gt[1].edge )

end
