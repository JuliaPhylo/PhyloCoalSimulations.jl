using PhyloCoalSimulations

using Pkg
Pkg.add(PackageSpec(url="https://github.com/JuliaPhylo/SNaQ.jl", rev="dev"))
# fixit: remove line above after SNaQ is registered
# instead update dependencies in test/Project.toml:
# remove Pkg, but add SNaQ:
# SNaQ = "c2bf7a07-44d2-4d11-b5e2-76dc5aa199d6"

using Distributions
using HypothesisTests
using PhyloNetworks
using Random
using StableRNGs
using Test
using SNaQ

const PN = PhyloNetworks
const PCS = PhyloCoalSimulations

@testset "PhyloCoalSimulations.jl" begin
    include("test_onepopulation.jl")
    include("test_multispeciesnetwork.jl")
    include("test_utils.jl")
end
