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

# 1 replicate, 1 locus: basic checks about labels etc.
# 3 taxa, long single edge root --> MRCA, so P(1 ancestor at origin pop) ≈ 1
net = readnewick("((((B:.5,A:.5):1.5,O:2):100)r);");
_,lab,gt = simulate_polygenictrait(net,1,1,rootprior,transition;
    nindividuals=Dict("A"=>1, "B"=>2, "O"=>2, "extra"=>10));
@test lab == ["B_1","B_2","A","O_1","O_2"]
@test length(gt) == 1
@test all(n.booln5 for n in gt[1].node)
@test all(n.fvalue != -1 for n in gt[1].node)

# 100 replicates, 1 locus, 3 individuals, 1 taxon + 2-cycle above it
ell = 2.0; # in coalescent units
net = readnewick("((t:0.0)#H1:$ell::0.6,#H1:$ell)r;")
v = 2.1 # ell + v0 # theoretical derivation: variance at the tip t
q = 0.8646647167633873 # 1-exp(-ell);
r = 0.5676676416183064 # 1-q/ell
c = 0.6353369125547349 # (0.6^2 +0.4^2)*(q*v0 + ell*r) # covariance at tip t
rng = StableRNG(771)
x,lab,_ = simulate_polygenictrait(net,100,1,rootprior,transition; nindividuals=3);
@test lab == ["t_1","t_2","t_3"]
@test all(-.22 .≤ Statistics.mean(x) .≤ .22) # 0 ± v/sqrt(100) = ± 0.21
empcov = Statistics.cov(x) # compare to [v c c; c v c; c c v]
# quantile(Chisq(99), [.10, .90]) ./99 # 0.8227197247573989, 1.1859281129977668
@test all(0.81 .≤ diag(empcov) ./ v .≤ 1.21)
@test all(-.25 .≤ empcov[[2,3,6]] .- c .≤ .25)
end
