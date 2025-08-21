using PhyloCoalSimulations

using Distributions
using HypothesisTests
using LinearAlgebra: diag
using PhyloNetworks
using Random
using StableRNGs
using Statistics
using Test

const PN = PhyloNetworks
const PCS = PhyloCoalSimulations

@testset "PhyloCoalSimulations.jl" begin
    include("test_onepopulation.jl")
    include("test_multispeciesnetwork.jl")
    include("test_utils.jl")
end
