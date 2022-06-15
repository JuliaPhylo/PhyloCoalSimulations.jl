using PhyloCoalSimulations
import PhyloCoalSimulations as PCS
using Test
using Random
import PhyloNetworks as PN

@testset "PhyloCoalSimulations.jl" begin
    include("test_onepopulation.jl")
end
