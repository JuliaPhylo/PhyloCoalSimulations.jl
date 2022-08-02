using PhyloCoalSimulations

using Distributions
using HypothesisTests
using PhyloNetworks
using Random
using Test
const PN = PhyloNetworks
const PCS = PhyloCoalSimulations

@testset "PhyloCoalSimulations.jl" begin
    include("test_onepopulation.jl")
    include("test_multispeciesnetwork.jl")
end
