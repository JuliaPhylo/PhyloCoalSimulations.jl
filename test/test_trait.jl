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
trait = simulate_polygenictrait(net, 100000, 1, root_prior, transition_distribution);
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
end