@testset "Trait with ILS" begin
## Test covariance formula on a simple tree with BM
#= R code to produce the tree Matrix
devtools::install_github("pbastide/phylolm", ref = "ils") # branch with ILS transformation
library(phylolm)
tree <- read.tree(text = "((B:.5,A:.5):1.5,O:2)r;")
treeils <- transf.branch.lengths(tree, model = "ILS", parameters = list("lambda_ILS" = 0.1))
vcv(treeils$tree)
          B         A   O
B 2.1000000 0.8008171 0.0
A 0.8008171 2.1000000 0.0
O 0.0000000 0.0000000 2.1
=#
#= Simulation and test using Julia
net = PN.readnewick("((B:.5,A:.5):1.5,O:2)r;");
Random.seed!(1848)
m0 = 0.0; v0 = 0.1; sigmaBM = 1;
root_prior = NormalCanon(m0, 1 / v0)
transition_distribution(x,t) = Normal(x, sigmaBM*sqrt(t))
trait,lab,gt = simulate_polygenictrait(net, 100000, 1, root_prior, transition_distribution);
using Statistics
empcov = Statistics.cov(trait)
#3×3 Matrix{Float64}:
#  2.10979      0.799914    -0.00729777
#  0.799914     2.08758     -0.00950007
# -0.00729777  -0.00950007   2.10296
theocovbm = Matrix(PN.vcv(net));
theocovq = 1 .- exp.(-theocovbm);
theocovq[diagind(theocovq)] .= 0.0;
theocov = sigmaBM^2 .* (theocovbm - theocovq) + v0 .* (theocovq + Matrix(I, 3, 3))
#3×3 Matrix{Float64}:
# 2.1       0.800817  0.0
# 0.800817  2.1       0.0
# 0.0       0.0       2.1
@test empcov ≈ theocov rtol=1e-2
=#

v0 = 0.1; rootprior = NormalCanon(0, 1/v0)
transition(x,t) = Normal(x, sqrt(t))

# 1 replicate, 1 locus, default rng: basic checks
# 3 taxa, long single edge root --> MRCA, so P(1 ancestor at origin pop) ≈ 1
net = readnewick("((((B:.5,A:.5):1.5,O:2):100)r);");
_,lab,gt = simulate_polygenictrait(net,1,1,rootprior,transition;
    nindividuals=Dict("A"=>1, "B"=>2, "O"=>2, "extra"=>10));
@test lab == ["B_1","B_2","A","O_1","O_2"]
@test length(gt) == 1
@test all(n.booln5 for n in gt[1].node)
@test all(n.fvalue != -1 for n in gt[1].node)

nind = "this is neither an integer nor a dictionary"
message = "nindividuals should be an integer or dictionary of integers"
@test_throws ErrorException(message) simulate_polygenictrait(net,1,1,rootprior,transition; nindividuals=nind);

# 100 replicates, 1 locus, 3 individuals, 1 taxon + 2-cycle above it
ell = 2.0; # in coalescent units
net = readnewick("((t:0.0)#H1:$ell::0.6,#H1:$ell)r;")
v = 2.1 # ell + v0 # theoretical derivation: variance at the tip t
q = 0.8646647167633873 # 1-exp(-ell);
r = 0.5676676416183064 # 1-q/ell
c = 0.6353369125547349 # (0.6^2 +0.4^2)*(q*v0 + ell*r) # covariance at tip t
rng = StableRNG(771)
x,lab,_ = simulate_polygenictrait(rng, net,100,1,rootprior,transition; nindividuals=3);
@test lab == ["t_1","t_2","t_3"]
@test all(-.3 .≤ Statistics.mean(x) .≤ .3) # 0 ± v/sqrt(100) = ± 0.21
empcov = Statistics.cov(x) # compare to [v c c; c v c; c c v]
# quantile(Chisq(99), [.10, .90]) ./99 # 0.8227197247573989, 1.1859281129977668
@test all(0.81 .≤ diag(empcov) ./ v .≤ 1.21)
@test all(-.25 .≤ empcov[[2,3,6]] .- c .≤ .25)

#= longer simulations to check that the theoretical covariance
#  from PhyloTraits.gaussiancoalescent_covariancematrix(net,v0/1)
#  matches the covariance in simulated traits
# 5 taxa, level-2, subnet on A,B,C is a tree, d1 below 1 hybrid, d2 below 2
net = readnewick("(((C:2,#H1:.1):0.3,(((d1:1,#H2:.1):.8)#H1:.7::.6,(d2:.5)#H2:1::.7):.4):.3,(B:1,A:.5):2);");
# plot(net, showedgelength=true, showgamma=true);
# v0=0.1;
# cM, eV = PhyloTraits.gaussiancoalescent_covariancematrix(net,v0/1)
# tiplabels(eV) == tiplabels(cM) == ["C","d1","d2","B","A"]
cT = [2.7,3.02,2.396,3.1,2.6] # eV[:tips] + diag(cM[:tips])
v11=1.7668462203929007; v12=0.11761402816197718; v13=0.08199968747807541
v22=2.038420771136349; v23=0.443417290678788; v33=1.247326232249089
v44=2.144808361531078; v45=1.2218017549129516; v55=1.6738764987615091
cM = [v11 v12 v13 0 0; v12 v22 v23 0 0; v13 v23 v33 0 0;
    0 0 0 v44 v45; 0 0 0 v45 v55] # cM[:tips]
nrep=1_000_000
x,lab,_ = simulate_polygenictrait(net,nrep,1,rootprior,transition; nindividuals=2);
@test lab == ["C_1","C_2","d1_1","d1_2","d2_1","d2_2","B_1","B_2","A_1","A_2"]
mean(x) # should be 0 ± repeat(cT ./ sqrt(nrep), inner=2)
extrema(mean(x)) # (-0.0030091839177305334, 0.0034590473513944927): all good!
empcov = Statistics.cov(x)
#=
10×10 Matrix{Float64}:
 2.70302     1.76553      0.118256    0.119325      0.0843356     0.0835598     0.00710998    0.00399126   0.00457252    0.0037936
 1.76553     2.69819      0.117498    0.117525      0.0875738     0.0864092     0.00344701    0.00267235   0.00262211    0.00326495
 0.118256    0.117498     3.02681     2.0452        0.448045      0.450111     -0.00083566    0.00446083   0.00498976    0.00586739
 0.119325    0.117525     2.0452      3.0192        0.449673      0.450464      0.000317299   0.00486496   0.00129635    0.0028426
 0.0843356   0.0875738    0.448045    0.449673      2.39993       1.25083      -0.00213595   -0.00221455   0.000847208  -0.000626026
 0.0835598   0.0864092    0.450111    0.450464      1.25083       2.39465       0.000600053  -0.000335371  0.003998      0.0026118
 0.00710998  0.00344701  -0.00083566  0.000317299  -0.00213595    0.000600053   3.10008       2.14406      1.22063       1.22157
 0.00399126  0.00267235   0.00446083  0.00486496   -0.00221455   -0.000335371   2.14406       3.09624      1.22052       1.21953
 0.00457252  0.00262211   0.00498976  0.00129635    0.000847208   0.003998      1.22063       1.22052      2.59601       1.67167
 0.0037936   0.00326495   0.00586739  0.0028426    -0.000626026   0.0026118     1.22157       1.21953      1.67167       2.59521
=#
extrema(diag(empcov) ./ repeat(cT, inner=2) .- 1.0)
# (-0.001842170480209404, 0.002255885223152454) : good!
maximum(abs.(LinearAlgebra.tril(empcov - repeat(cM, inner=(2,2)), -1)))
# 0.0071099844807377575: good
=#

# 3 taxa, multivariate
P0 = [5.0  0.0  0.5; 0.0  10.0  3.2; 0.5  3.2  7.0]
rootpriormv = MvNormalCanon([0, 0, 0], P0)
transitionmv(x,t) = MvNormal(x, sqrt(t).*[1 0 0; 0 1 0; 0 0 1])
net = readnewick("((((B:.5,A:.5):1.5,O:2):1)r);");
X,_,_ = simulate_polygenictrait(net,5,10,rootpriormv,transitionmv;
    nindividuals=Dict("A"=>1, "B"=>2, "O"=>2, "extra"=>10));
@test length(X) == 5
@test size(X[2]) == (5,3)

end
